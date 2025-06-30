# plot_model_performance.R
set.seed(43)

# 1) Load libraries
library(tidyverse)
library(tidyr)      # for pivot_longer
library(ggplot2)

# 2) Read in your foldwise scores
df <- read.csv("foldwise_scores_microbial_only.csv", stringsAsFactors = FALSE)

# Make sure names are syntactically valid (F1 Macro → F1.Macro)
names(df) <- make.names(names(df))

# 3) If there is a Data column, filter to the “Microbial Only” set
if("Data" %in% colnames(df)) {
  df <- df %>% filter(Data == "Genes Only")
}

# 4) Pivot to long format for the four metrics
df_long <- df %>%
  pivot_longer(
    cols      = c(Accuracy, F1.Macro, Kappa, MCC),
    names_to  = "Metric",
    values_to = "Score"
  )

# 5) Rename metrics for prettier facet labels
df_long <- df_long %>%
  mutate(
    Metric = recode(
      Metric,
      "F1.Macro" = "F1 Macro",
      "Accuracy" = "Mean Balanced Accuracy"
    )
  )

# 6) Order the Model factor by descending mean Kappa
model_order <- df_long %>%
  filter(Metric == "Kappa") %>%
  group_by(Model) %>%
  summarise(mean_k = mean(Score, na.rm = TRUE)) %>%
  arrange(desc(mean_k)) %>%
  pull(Model)

df_long$Model <- factor(df_long$Model, levels = model_order)

# 7) Identify top-3 models by average Kappa and flag them
top3 <- df_long %>%
  filter(Metric == "Kappa") %>%
  group_by(Model) %>%
  summarise(mean_k = mean(Score, na.rm = TRUE)) %>%
  arrange(desc(mean_k)) %>%
  slice(1:3) %>%
  pull(Model)

df_long <- df_long %>%
  mutate(
    Highlight = if_else(Model %in% top3, "Top 3", "Other")
  )

# 8) Save plot as high-res TIFF
tiff(
  filename    = "model_comparison_plot_with_metadata.tiff",
  width       = 8,       # inches
  height      = 6,       # inches
  units       = "in",
  res         = 600,     # DPI
  compression = "lzw"
)

ggplot(df_long, aes(x = Score, y = Model)) +
  geom_boxplot(
    width         = 0.5,
    outlier.shape = NA,
    color         = "black",
    fill          = "transparent"
  ) +
  stat_summary(
    aes(fill = Highlight),
    fun   = mean,
    geom  = "point",
    shape = 22,
    size  = 4
  ) +
  scale_fill_manual(values = c("Top 3" = "deeppink", "Other" = "blue")) +
  facet_wrap(~ Metric, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(x = "Score", y = NULL) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 1.5),
    axis.line        = element_line(color = "black", size = 1),
    axis.text.x      = element_text(size = 8, color = "black"),
    axis.text.y      = element_text(size = 8, color = "black"),
    axis.title.x     = element_text(size = 8, color = "black", margin = margin(t = 10)),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#f0f0f0", color = "black", size = 1),
    legend.position  = "none"
  )

dev.off()
