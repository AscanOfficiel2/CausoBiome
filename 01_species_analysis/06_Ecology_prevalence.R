######################################################################
# Genome-Resolved Metagenomic Analysis
# Alpha & Beta Diversity, Prevalence, WHO Pathogen Ecology, Networks
######################################################################

set.seed(42)

# ---- Packages ----
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(FSA)         # Dunn test
  library(tidyr)
  library(tibble)
  library(parallel)
  library(pheatmap)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(ggraph)
  library(patchwork)
})



# ============================================================
# Align TPM matrix with metadata after filtering
# ============================================================

# 1️⃣ Load TPM and metadata
tpm <- read.csv("CRC_species_matrix_filtered.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("CRC_Cleaned_MetaData.csv")

# 2️⃣ Transpose TPM if needed (samples must be columns)
###########if (nrow(tpm) > ncol(tpm)) tpm <- t(tpm)

# 3️⃣ Clean sample names (avoid hidden differences)
colnames(tpm) <- trimws(as.character(colnames(tpm)))
meta$Sample_ID <- trimws(as.character(meta$Sample_ID))

# 4️⃣ Keep only samples that exist in both datasets
shared_samples <- intersect(colnames(tpm), meta$Sample_ID)

# 5️⃣ Subset TPM and metadata
tpm <- tpm[, shared_samples]
meta_aligned <- meta[match(shared_samples, meta$Sample_ID), ]

# 6️⃣ Confirm alignment
stopifnot(all(meta_aligned$Sample_ID == colnames(tpm)))

# 7️⃣ Save aligned metadata for reuse
write.csv(meta_aligned, "CRC_metadata_aligned.csv", row.names = FALSE)

cat("✅ Metadata aligned successfully! Saved as 'CRC_metadata_aligned.csv'\n")
cat("   Total matched samples:", length(shared_samples), "\n")


######################################################################
# 1. ALPHA DIVERSITY (raw TPM)
######################################################################

tpm <- read.csv("CRC_species_matrix_filtered.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("CRC_metadata_aligned.csv")

if (nrow(tpm) > ncol(tpm)) tpm <- t(tpm)
tpm[is.na(tpm)] <- 0

stopifnot("Sample_ID" %in% colnames(meta))
rownames(meta) <- meta$Sample_ID
meta <- meta[rownames(tpm), , drop = FALSE]

alpha_df <- data.frame(
  Sample_ID = rownames(tpm),
  Shannon   = diversity(tpm, "shannon"),
  Simpson   = diversity(tpm, "simpson"),
  Richness  = specnumber(tpm)
)

thresholds <- c(0.1, 1, 10)
for (thr in thresholds) {
  alpha_df[[paste0("Richness_TPM≥", thr)]] <- rowSums(tpm >= thr)
}
alpha_df <- left_join(alpha_df, meta, by = "Sample_ID")

# plotting function
group_colors <- c(
  "Healthy" = "#ffcc00",            
  "Adenoma" = "#984ea3", 
  "Cancer" = "#4daf4a"
)

plot_alpha <- function(df, yvar, ylab, file) {
  p <- ggplot(df, aes(x = Group, y = .data[[yvar]], fill = Group)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.8) +
    #geom_jitter(width = 0.2, size = 0.1, alpha = 0.5) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    theme_minimal(base_size = 11) +
    labs(x = "", y = ylab) +
    theme(
      legend.position   = "none",
      axis.text.x       = element_text(angle = 45, hjust = 1, color = "black", size = 4),
      axis.text.y       = element_text(color = "black", size = 4),
      axis.title.y      = element_text(size = 4, face = "bold"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.ticks        = element_line(color = "black", linewidth = 0.2)
    )
  ggsave(file, p, width = 1, height = 2, dpi = 600)
}

plot_alpha(alpha_df, "Shannon", "Shannon Index", "Alpha_Shannon.tiff")
plot_alpha(alpha_df, "Simpson", "Simpson Index", "Alpha_Simpson.tiff")
plot_alpha(alpha_df, "Richness", "Richness", "Alpha_Richness.tiff")
for (thr in thresholds) {
  yvar <- paste0("Richness_TPM≥", thr)
  plot_alpha(alpha_df, yvar, paste0("Richness (TPM ≥ ", thr, ")"),
             paste0("Alpha_", yvar, ".tiff"))
}

write.csv(alpha_df, "AlphaDiversity_TPM_allmetrics.csv")

#####################################################################
#         Alpha Diversity Statistics
#####################################################################
# Stats: Kruskal–Wallis + Dunn ----
stat_list <- list()

# KW for Shannon, Simpson, Raw Richness
kw_shannon  <- kruskal.test(Shannon ~ Group, data = alpha_df)
kw_simpson  <- kruskal.test(Simpson ~ Group, data = alpha_df)
kw_richness <- kruskal.test(Richness ~ Group, data = alpha_df)

kw_df <- tibble(
  Test        = "Kruskal-Wallis",
  Metric      = c("Shannon", "Simpson", "Richness"),
  Comparison  = NA_character_,
  Z           = NA_real_,
  P.unadj     = NA_real_,
  P.adj       = NA_real_,
  Chi_squared = c(kw_shannon$statistic, kw_simpson$statistic, kw_richness$statistic),
  df          = c(kw_shannon$parameter, kw_simpson$parameter, kw_richness$parameter),
  p_value     = c(kw_shannon$p.value, kw_simpson$p.value, kw_richness$p.value)
)

stat_list[["KW"]] <- kw_df

# Dunn’s post-hoc for Shannon, Simpson, Richness
dunn_shannon  <- dunnTest(Shannon ~ Group, data = alpha_df, method = "bh")$res %>%
  mutate(Metric = "Shannon")
dunn_simpson  <- dunnTest(Simpson ~ Group, data = alpha_df, method = "bh")$res %>%
  mutate(Metric = "Simpson")
dunn_richness <- dunnTest(Richness ~ Group, data = alpha_df, method = "bh")$res %>%
  mutate(Metric = "Richness")

dunn_all <- bind_rows(dunn_shannon, dunn_simpson, dunn_richness) %>%
  select(Metric, Comparison, Z, P.unadj, P.adj) %>%
  mutate(Test = "Dunn")

# Threshold richness tests
for (thr in thresholds) {
  metric <- paste0("Richness_TPM≥", thr)
  kw <- kruskal.test(alpha_df[[metric]] ~ alpha_df$Group)
  dunn <- dunnTest(alpha_df[[metric]] ~ alpha_df$Group, method = "bh")$res %>%
    mutate(Metric = metric)
  
  # Append to KW
  kw_df_thr <- tibble(
    Test        = "Kruskal-Wallis",
    Metric      = metric,
    Comparison  = NA_character_,
    Z           = NA_real_,
    P.unadj     = NA_real_,
    P.adj       = NA_real_,
    Chi_squared = kw$statistic,
    df          = kw$parameter,
    p_value     = kw$p.value
  )
  stat_list[[paste0("KW_", metric)]] <- kw_df_thr
  
  # Append to Dunn
  dunn_all <- bind_rows(dunn_all,
                        dunn %>% select(Comparison, Z, P.unadj, P.adj) %>%
                          mutate(Test = "Dunn", Metric = metric))
}

# Combine everything
alpha_stats <- bind_rows(stat_list[["KW"]], dunn_all)
for (thr in thresholds) {
  alpha_stats <- bind_rows(alpha_stats, stat_list[[paste0("KW_", paste0("Richness_TPM≥", thr))]])
}

write.csv(alpha_stats, "AlphaDiversity_statistics.csv", row.names = FALSE)


######################################################################
# 2. BETA DIVERSITY (CLR-corrected TPM)
######################################################################

species_clr <- read.csv("TPM_clr_batch_corrected_v1.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("CRC_metadata_aligned.csv")

if (nrow(species_clr) > ncol(species_clr)) species_clr <- t(species_clr)
rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[rownames(species_clr), , drop = FALSE]

grouping <- metadata$Group

dist_matrix <- vegdist(species_clr, method = "euclidean")
anosim_res  <- anosim(dist_matrix, metadata$Group, permutations = 999, parallel = detectCores() - 1)
print(anosim_res)
# Set color vector based on ordering of the groups + "Between"
ordered_groups <- levels(as.factor(grouping))  # Ensure consistent order
anosim_colors <- c("black",                     # Between
                   group_colors["Healthy"],
                   group_colors["Adenoma"],
                   group_colors["Cancer"])

tiff("ANOSIM_CLR_Euclidean.tiff", width = 7, height = 4, units = "in", res = 600)
plot(anosim_res,
     col = anosim_colors,
     ylab = "Dissimilarity Rank Value", xlab = "",
     cex.lab = 0.8, cex.axis = 0.8,
     lwd = 1)

box(lwd = 5)  # Border thickness
dev.off()


NMDS <- metaMDS(species_clr, k = 2, distance = "euclidean", trymax = 50)
NMDS$stress  # check stress

nmds_df <- as.data.frame(NMDS$points) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(metadata, by="Sample_ID")

tiff("NMDS_CLR_Euclidean_AllGroups.tiff", width = 4, height = 4, units = "in", res = 600)
ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(size = 2.3, alpha = 0.6) +
  scale_color_manual(values = group_colors) +
  theme_classic(base_size = 12) +
  labs(x = "NMDS1", y = "NMDS2", title = "") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  )
dev.off()




# ----  Test homogeneity of dispersion ----
dispersion <- betadisper(dist_matrix, grouping, type = "centroid")
anova_disp <- anova(dispersion)
print(anova_disp)

write.csv(as.data.frame(anova_disp), "MultivariateDispersion_ANOVA.csv")

# ----  PERMANOVA ----
perm <- adonis2(dist_matrix ~ Group, data = metadata, permutations = 999)
print(perm)
write.csv(as.data.frame(perm), "PERMANOVA_CLR_Euclidean.csv")

# ----  Pairwise PERMANOVA ----
pairwise_adonis2 <- function(data, grouping, method = "euclidean", permutations = 999) {
  groups <- unique(grouping)
  comparisons <- combn(groups, 2, simplify = FALSE)
  results <- data.frame()
  
  for (comp in comparisons) {
    subset_rows <- grouping %in% comp
    dist_sub <- vegdist(data[subset_rows, ], method = method)
    meta_sub <- data.frame(Group = factor(grouping[subset_rows]))
    
    ad <- adonis2(dist_sub ~ Group, data = meta_sub, permutations = permutations)
    
    results <- rbind(results, data.frame(
      Comparison = paste(comp[1], "vs", comp[2]),
      R2 = round(ad$R2[1], 4),
      F_model = round(ad$F[1], 3),
      p_value = ad$`Pr(>F)`[1]
    ))
  }
  results$Significance <- cut(results$p_value,
                              breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                              labels = c("***", "**", "*", "ns"))
  return(results)
}

pairwise_results <- pairwise_adonis2(species_clr, metadata$Group)
write.csv(pairwise_results, "PERMANOVA_Pairwise_CLR_Euclidean.csv", row.names = FALSE)
print(pairwise_results)



# Extract centroid coordinates from betadisper
centroids_df <- data.frame(Group = rownames(scores(dispersion, display = "centroids")),
                           scores(dispersion, display = "centroids"))

# Rename columns to match NMDS data frame (MDS1, MDS2)
colnames(centroids_df)[2:3] <- c("MDS1", "MDS2")

# Save centroid coordinates
write.csv(centroids_df, "MultivariateDispersion_Centroids.csv", row.names = FALSE)


######################################################################
# 3. ENVIRONMENT-SPECIFIC PREVALENCE (all species)
######################################################################


if (nrow(tpm) > ncol(tpm)) tpm <- t(tpm)
tpm[is.na(tpm)] <- 0
rownames(meta) <- meta$Sample_ID
tpm <- tpm[rownames(meta), ]

pa_matrix <- (tpm > 0) * 1
pa_long <- as.data.frame(pa_matrix) %>%
  rownames_to_column("Sample_ID") %>%
  pivot_longer(-Sample_ID, names_to = "Species", values_to = "Present") %>%
  left_join(meta, by = "Sample_ID")

prevalence_env <- pa_long %>%
  group_by(Group, Species) %>%
  summarise(Prevalence = mean(Present) * 100, .groups="drop")

write.csv(prevalence_env, "Species_Prevalence_ByEnv_all_species.csv", row.names = FALSE)


### Prevalence at setting Threshold #######
prevalence_by_env <- pa_long %>%
  group_by(Group, Species) %>%
  summarise(Prevalence = mean(Present) * 100, .groups = "drop") %>%
  mutate(Category = case_when(
    Prevalence >= 80 ~ "Core",          # relaxed core
    Prevalence >= 50 ~ "Secondary",     # relaxed secondary
    TRUE ~ "Peripheral"
  ))


# Save full table
write.csv(prevalence_by_env, "Species_Prevalence_ByEnvironment_threshold.csv", row.names = FALSE)

# ----  Summary counts per environment with grouping ----
summary_counts <- prevalence_by_env %>%
  group_by(Group, Category) %>%
  summarise(n_species = n(), .groups = "drop")

write.csv(summary_counts, "Species_Category_Counts_ByEnvironment_threshold.csv", row.names = FALSE)

####----------------------------------------------


# ---- 3) Global prevalence ----
prevalence_global <- colSums(pa_matrix) / nrow(pa_matrix) * 100

species_categories <- data.frame(
  Species = names(prevalence_global),
  Prevalence = prevalence_global
) %>%
  mutate(Category = case_when(
    Prevalence >= 50 ~ "Core",
    Prevalence >= 30 ~ "Secondary",
    TRUE ~ "Peripheral"
  ))

# Make Category a factor
species_categories$Category <- factor(species_categories$Category,
                                      levels = c("Core", "Secondary", "Peripheral"))

# Summarise counts
summary_counts_global <- species_categories %>%
  group_by(Category) %>%
  summarise(n_species = n(), .groups = "drop")

# Add percentages
total_species <- sum(summary_counts_global$n_species)
summary_counts_global <- summary_counts_global %>%
  mutate(Percent = round(n_species / total_species * 100, 1))  # 1 decimal place




# Make Category a factor
prevalence_by_env$Category <- factor(prevalence_by_env$Category,
                                     levels = c("Core", "Secondary", "Peripheral"))

summary_counts_env <- prevalence_by_env %>%
  group_by(Group, Category) %>%
  summarise(n_species = n(), .groups = "drop")

# ---- 5) Barplots (Global vs Environment-specific) ----
category_colors <- c("Core" = "#e377c2", "Secondary" = "#969696", "Peripheral" = "#6baed6")


# ---- Global prevalence plot ----
p_global <- ggplot(summary_counts_global, aes(x = Category, y = n_species, fill = Category)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste0(Percent, "%")),
            vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = category_colors) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(x = "", y = "Number of Species", title = "")

# Save separately
tiff("Prevalence_Global_relaxed.tiff", width = 5, height = 5, units = "in", res = 600)
print(p_global)
dev.off()

# ---- Environment-specific prevalence plot ----
p_env <- ggplot(summary_counts_env, aes(x = Group, y = n_species, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = category_colors) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "Environment", y = "Number of Species", title = "")

# Save separately
tiff("Prevalence_Environment_relaxed.tiff", width = 6, height = 5, units = "in", res = 600)
print(p_env)
dev.off()

#####################################################

# ---- Heatmap of top core species per environment ----

# Pick top N species (by prevalence) within each environment’s "Core" category
top_core_env <- prevalence_by_env %>%
  filter(Category == "Core") %>%
  group_by(Group) %>%
  arrange(desc(Prevalence)) %>%
  slice_head(n = 50) %>%        # top 20 per environment
  ungroup()

# Build matrix: species (rows) x environment (columns)
top_core_matrix <- prevalence_by_env %>%
  filter(Species %in% top_core_env$Species) %>%
  select(Group, Species, Prevalence) %>%
  pivot_wider(names_from = Group, values_from = Prevalence, values_fill = 0) %>%
  column_to_rownames("Species")



# Ccolours
heat_colors <- colorRampPalette(c("white", "#6baed6", "#08306b"))(100)

tiff("Heatmap_EnvSpecific_Core_Species.tiff", 
     width = 7, height = 8, units = "in", res = 600, compression = "lzw")
pheatmap(as.matrix(top_core_matrix),
         cluster_rows = TRUE,
         cluster_cols = FALSE,        # keep environment order
         color = heat_colors,
         fontsize_row = 8,            # adjust for species labels
         fontsize_col = 12,
         angle_col = "45",            # slant environment names
         border_color = NA,
         main = "",
         cellwidth = 25,
         cellheight = 10)
dev.off()
#############################################################

# ---- Top core species per environment (relaxed thresholds) ----
top_core_env <- prevalence_by_env %>%
  filter(Category == "Core") %>%
  group_by(Group) %>%
  arrange(desc(Prevalence)) %>%
  slice_head(n = 20) %>%    # top 20 per environment
  ungroup()


tiff("Barplot_Core_Species_PerEnv.tiff", width = 10, height = 8, units = "in", res = 600, compression = "lzw")
ggplot(top_core_env, aes(x = reorder(Species, Prevalence), y = Prevalence, fill = Group)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  labs(x = "Species", y = "Prevalence (%)",
       title = "") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  )
dev.off()

###################################### 
library(UpSetR)

# ---- Build species × environment matrix ----
species_env_matrix <- prevalence_by_env %>%
  filter(Prevalence > 0) %>%         # keep only species that occur at least once
  distinct(Group, Species) %>%
  mutate(Present = 1) %>%
  pivot_wider(names_from = Group, values_from = Present, values_fill = 0)

# Set species as rownames
rownames(species_env_matrix) <- species_env_matrix$Species
species_env_matrix$Species <- NULL

species_env_matrix <- as.data.frame(species_env_matrix)  # convert tibble → data.frame

tiff("UpSet_Species_Shared_Environments.tiff", width = 10, height = 6, units = "in", res = 600, compression = "lzw")
UpSetR::upset(
  species_env_matrix,
  sets = colnames(species_env_matrix),
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = group_colors[colnames(species_env_matrix)], # colored set bars
  main.bar.color = "black",   # black bars for intersections
  text.scale = 1.5,
  point.size = 4,
  line.size = 1
)
dev.off()




#----- Environment Overlap Network to complement Upset----------------
library(igraph)


# ---- 1) Build environment overlap matrix ----
envs <- colnames(species_env_matrix)

edges <- data.frame(from = character(), to = character(), 
                    weight = numeric(), percent = numeric())

for(i in 1:(length(envs)-1)) {
  for(j in (i+1):length(envs)) {
    overlap <- sum(species_env_matrix[, envs[i]] == 1 & species_env_matrix[, envs[j]] == 1)
    min_size <- min(colSums(species_env_matrix)[c(envs[i], envs[j])])
    perc <- round(100 * overlap / min_size, 1)   # % overlap relative to smaller set
    edges <- rbind(edges, data.frame(from = envs[i], to = envs[j], 
                                     weight = overlap, percent = perc))
  }
}

# ---- 2) Node attributes ----
nodes <- data.frame(
  name = envs,
  richness = colSums(species_env_matrix)
)

# ---- 3) Create igraph object ----
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# ---- 4) Plot upgraded network ----
tiff("Network_Env_Species_Sharing.tiff", width = 7, height = 7, units = "in", res = 600, compression = "lzw")
ggraph(g, layout = "fr") +
  # Edges: thickness by shared species, color by % overlap
  geom_edge_link(aes(width = weight, color = percent), alpha = 0.7) +
  
  # Edge labels using edge midpoints
  geom_edge_link(aes(label = paste0(weight, " (", percent, "%)")),
                 angle_calc = "along", label_dodge = unit(2.5, "mm"),
                 label_size = 3, color = "black", show.legend = FALSE,
                 check_overlap = TRUE, repel = TRUE) +
  
  # Nodes
  geom_node_point(aes(size = richness, color = name)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
  
  scale_edge_width(range = c(0.5, 4)) +
  scale_edge_color_gradient(low = "#9ecae1", high = "#08306b") +  # light → dark blue
  scale_size(range = c(6, 15)) +
  
  scale_color_manual(values = c(
    "Ambulance" = "#1f77b4",
    "Hosp_env" = "#aec7e8",
    "Hosp_sewage" = "#ff7f0e",
    "Public_transp" = "#fdbf6f"
  )) +
  
  theme_void() +
  labs(title = "",
       edge_width = "Shared species (n)",
       edge_color = "% overlap (relative to smaller env)",
       size = "Species richness") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()
###################################

# ============================================================
# Identify unique and shared species among groups
# ============================================================

# Each column = Group (Healthy, Adenoma, Cancer)
# Each row = Species, value = 1 (present) or 0 (absent)

# 1️⃣ Species unique to each group
unique_species <- list()

for (group in colnames(species_env_matrix)) {
  others <- setdiff(colnames(species_env_matrix), group)
  # Species present in group but absent in others
  unique_species[[group]] <- rownames(species_env_matrix)[
    (species_env_matrix[, group] == 1) &
      (rowSums(species_env_matrix[, others, drop = FALSE]) == 0)
  ]
  cat(sprintf("✅ %d species unique to %s\n",
              length(unique_species[[group]]), group))
}

# Save each list to file
for (g in names(unique_species)) {
  write.csv(data.frame(Species = unique_species[[g]]),
            paste0("Unique_Species_", g, ".csv"), row.names = FALSE)
}


# 2️⃣ Species shared by all groups (core)
core_species <- rownames(species_env_matrix)[
  rowSums(species_env_matrix) == ncol(species_env_matrix)
]
cat(sprintf("\n✅ %d species shared by all groups (core microbiome)\n",
            length(core_species)))
write.csv(data.frame(Species = core_species),
          "Shared_Species_AllGroups_Core.csv", row.names = FALSE)


# 3️⃣ Species shared by exactly two groups
pairwise_shared <- list()

groups <- colnames(species_env_matrix)
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    pair <- paste(g1, g2, sep = "_")
    
    pairwise_shared[[pair]] <- rownames(species_env_matrix)[
      (species_env_matrix[, g1] == 1) &
        (species_env_matrix[, g2] == 1) &
        (rowSums(species_env_matrix[, setdiff(groups, c(g1, g2)), drop = FALSE]) == 0)
    ]
    
    cat(sprintf("✅ %d species shared by %s and %s (but not others)\n",
                length(pairwise_shared[[pair]]), g1, g2))
  }
}

# Save all pairwise shared species
for (p in names(pairwise_shared)) {
  write.csv(data.frame(Species = pairwise_shared[[p]]),
            paste0("Shared_Species_", p, ".csv"), row.names = FALSE)
}



