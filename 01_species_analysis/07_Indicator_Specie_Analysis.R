# ===================== Indicator Species (with robust ID normalization) =====================
set.seed(123)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(indicspecies); library(permute); library(stringi)
})

# ----------------------------- user controls ----------------------------------
detect_limit <- 0        # TPM > detect_limit -> presence (0, 0.1, or 1 are common)
min_present  <- 2        # keep species present in >= this many samples
min_per_grp  <- 3        # min samples required per group; groups with < min_per_grp are dropped
topN_plot    <- 60       # how many significant indicators to plot
# ------------------------------------------------------------------------------

# --- Helper: normalize IDs aggressively but safely (no semantic changes needed) ---
normalize_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- stri_trans_nfc(x)                    # canonical Unicode form
  x <- gsub("[\u00A0\t\r\n ]+", "", x)      # remove non-breaking space & whitespace inside
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")  # strip non-ASCII if present
  x <- gsub("[\"'`]", "", x)                # strip stray quotes/backticks
  x <- gsub("[^A-Za-z0-9_.-]", "", x)       # keep only safe chars
  # R sometimes prefixes columns with 'X' (not here, but just in case someone re-saved):
  x <- sub("^X(?=(ERR|SRR|DRR)\\d+$)", "", x, perl = TRUE)
  x
}

# --- Helper: nearest neighbor suggestions when overlap is empty ---
nearest_match_table <- function(source_ids, target_ids, k = 1) {
  if (length(source_ids) == 0 || length(target_ids) == 0) return(tibble())
  m <- adist(source_ids, target_ids, partial = TRUE, ignore.case = FALSE)
  best <- apply(m, 1, function(row) {
    j <- which.min(row)
    tibble(Target = target_ids[j], edit_distance = row[j])
  })
  out <- bind_cols(tibble(Source = source_ids), bind_rows(best))
  arrange(out, edit_distance)
}

# ----------------------------- 1) Load -----------------------------------------

# 1️⃣ Load TPM and metadata
tpm <- read.csv("CRC_species_matrix_filtered.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("CRC_metadata_aligned.csv")
stopifnot(all(c("Sample_ID","Group") %in% names(meta)))

# Ensure orientation: species x samples
if (nrow(tpm) >= ncol(tpm)) tpm <- t(tpm)  # if rows look like samples, flip
tpm <- t(tpm)  # now rows = species, cols = samples (explicit)

# Replace NAs
tpm[is.na(tpm)] <- 0

# Optional: clean species labels a bit
rownames(tpm) <- gsub("_", " ", rownames(tpm))

# ----------------------------- 2) Normalize IDs --------------------------------
meta$Sample_ID_raw <- meta$Sample_ID
col_ids_raw        <- colnames(tpm)

meta$Sample_ID     <- normalize_ids(meta$Sample_ID_raw)
meta$Group         <- trimws(meta$Group)
colnames(tpm)      <- normalize_ids(col_ids_raw)

# Diagnostics about normalization
id_diag <- tibble(
  tpm_col_before  = col_ids_raw,
  tpm_col_after   = colnames(tpm)
)
meta_diag <- tibble(
  meta_id_before  = meta$Sample_ID_raw,
  meta_id_after   = meta$Sample_ID
)

write.csv(id_diag,  "indicator_IDs_TPM_normalization.csv",  row.names = FALSE)
write.csv(meta_diag,"indicator_IDs_META_normalization.csv", row.names = FALSE)

# ----------------------------- 3) Align samples --------------------------------
common_samples <- intersect(colnames(tpm), meta$Sample_ID)

if (length(common_samples) == 0) {
  # try a softer “core” match heuristic: drop trailing run parts like ".1", "_R1", etc.
  soften <- function(x) sub("^([A-Za-z]+\\d+).*", "\\1", x)
  cs_soft <- intersect(soften(colnames(tpm)), soften(meta$Sample_ID))
  if (length(cs_soft) > 0) {
    # remap both to soft core and de-duplicate
    core_tpm  <- soften(colnames(tpm))
    core_meta <- soften(meta$Sample_ID)
    keep_tpm  <- !duplicated(core_tpm)
    keep_meta <- !duplicated(core_meta)
    tpm  <- tpm[, keep_tpm, drop = FALSE]
    colnames(tpm) <- core_tpm[keep_tpm]
    meta <- meta[keep_meta, , drop = FALSE]
    meta$Sample_ID <- core_meta[keep_meta]
    common_samples <- intersect(colnames(tpm), meta$Sample_ID)
  }
}

if (length(common_samples) == 0) {
  # write nearest match suggestions and stop clearly
  sugg1 <- nearest_match_table(colnames(tpm), meta$Sample_ID, k = 1)
  write.csv(sugg1, "indicator_ID_nearest_matches.csv", row.names = FALSE)
  stop("❌ Still no overlapping Sample_IDs after normalization. See 'indicator_ID_nearest_matches.csv'.")
}

# Keep only overlapping samples in the same order
tpm  <- tpm[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$Sample_ID), , drop = FALSE]
stopifnot(all(colnames(tpm) == meta$Sample_ID))

# ----------------------------- 4) Presence/absence ------------------------------
species_pa <- (tpm > detect_limit) * 1
keep_sp    <- rowSums(species_pa) >= min_present
species_pa <- species_pa[keep_sp, , drop = FALSE]

# Drop groups with too few samples
grp_tbl      <- table(meta$Group)
drop_groups  <- names(grp_tbl[grp_tbl < min_per_grp])
keep_cols    <- !(meta$Group %in% drop_groups)

species_pa <- species_pa[, keep_cols, drop = FALSE]
meta_keep  <- meta[keep_cols, , drop = FALSE]
grp_tbl2   <- table(meta_keep$Group)

diag_summary <- tibble(
  n_species_in       = nrow(tpm),
  n_samples_in       = ncol(tpm),
  n_species_kept     = nrow(species_pa),
  n_samples_kept     = ncol(species_pa),
  groups_dropped     = paste(drop_groups, collapse = ", "),
  groups_kept        = paste(names(grp_tbl2), collapse = ", "),
  detect_limit_TPM   = detect_limit,
  min_present_filter = min_present,
  min_per_grp        = min_per_grp
)
write.csv(diag_summary, "indicator_DIAGNOSTICS.csv", row.names = FALSE)

if (length(grp_tbl2) < 2) {
  write.csv(data.frame(Group = names(grp_tbl2), n = as.integer(grp_tbl2)),
            "indicator_remaining_groups.csv", row.names = FALSE)
  stop("❌ Need at least 2 environments with sufficient samples after filtering.")
}

# ----------------------------- 5) IndVal.g -------------------------------------
comm <- as.data.frame(t(species_pa))  # samples x species
stopifnot(nrow(comm) == nrow(meta_keep))

ctrl <- how(nperm = 999)
indval <- multipatt(comm, meta_keep$Group, func = "IndVal.g", control = ctrl)

# ----------------------------- 6) Tidy results ---------------------------------
sign_tbl <- indval$sign %>% as.data.frame() %>% rownames_to_column("Species")
if (!"p.value" %in% names(sign_tbl) && !is.null(indval$p.value)) {
  sign_tbl$p.value <- indval$p.value
}
if (!"p.value" %in% names(sign_tbl)) stop("p-values not found in multipatt result.")

sign_tbl$p.adj <- p.adjust(sign_tbl$p.value, method = "BH")

sel_cols <- grep("^s\\.", names(sign_tbl), value = TRUE)
sign_tbl$Environment <- apply(sign_tbl[, sel_cols, drop = FALSE], 1, function(x) {
  labs <- sub("^s\\.", "", names(x)[x == 1])
  if (length(labs) == 0) NA_character_ else paste(labs, collapse = "+")
})

sign_tbl <- sign_tbl %>% arrange(p.adj, desc(stat))

write.csv(sign_tbl, "indicator_ALL_raw.csv", row.names = FALSE)
sig_tbl <- sign_tbl %>% filter(!is.na(p.adj) & p.adj < 0.05)
write.csv(sig_tbl, "indicator_SIGNIFICANT_BH.csv", row.names = FALSE)

# ----------------------------- 7) Plot (top N) ---------------------------------
if (nrow(sig_tbl) > 0) {
  top_tbl <- sig_tbl %>%
    slice_head(n = min(topN_plot, nrow(sig_tbl))) %>%
    mutate(Environment_simple = sub("\\+.*$", "", Environment))
  
  
  # Palette that safely maps whatever groups remain
  env_levels <- sort(unique(top_tbl$Environment_simple))
  
  base_cols <- c(
    "Healthy" = "#ffcc00",
    "Adenoma" = "#984ea3",
    "Cancer" = "#4daf4a"
  )
  
  

  # fallback colors for unexpected names
  fallback <- setNames(RColorBrewer::brewer.pal(max(3, length(env_levels)), "Set2"),
                       env_levels)
  
  # merge with preference for base_cols
  env_cols <- fallback
  env_cols[names(base_cols)] <- base_cols
  
  
  p <- ggplot(top_tbl, aes(x = stat, y = reorder(Species, stat), color = Environment_simple)) +
    geom_point(size = 2.5) +
    scale_color_manual(values = env_cols) +
    labs(x = "Indicator value (IndVal.g)", y = NULL, color = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y      = element_text(size = 8, color = "black"),
      axis.text.x      = element_text(size = 9, color = "black"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = .8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.text      = element_text(size = 8)
    )
  
  tiff("Top_Indicators.tiff", width = 6, height = 8, units = "in", res = 600, compression = "lzw")
  print(p); dev.off()
} else {
  message("No significant indicators at BH q<0.05. Check 'indicator_ALL_raw.csv' for near-significant taxa.")
}
