# ============================================================
# DAGitty Causal Graph for ARG/VF → CRC Stage
# Author: AbdulAziz Ascandari
# Date: November 2025
# ============================================================

# ---- Load required library ----
if(!require(dagitty)) install.packages("dagitty", dependencies = TRUE)
library(dagitty)

# ---- Define the DAG ----
dag_str <- "
dag {
bb=\"-3.58,-3.211,3.965,3.913\"
ARG_VF_Abundance [exposure,pos=\"-0.667,2.486\"]
Age [pos=\"-2.085,-1.707\"]
BMI [pos=\"-0.832,-2.816\"]
CRC_Stage [outcome,pos=\"1.431,-1.529\"]
Diet [latent,pos=\"2.421,-0.484\"]
Medication [latent,pos=\"2.518,1.976\"]
Sex [pos=\"-3.036,-0.293\"]
ARG_VF_Abundance -> CRC_Stage
Age -> ARG_VF_Abundance
Age -> CRC_Stage
BMI -> ARG_VF_Abundance
BMI -> CRC_Stage
Diet -> ARG_VF_Abundance
Medication -> ARG_VF_Abundance
Sex -> ARG_VF_Abundance
Sex -> CRC_Stage
}
"

# ---- Load DAG into R ----
dag_model <- dagitty(dag_str)

# ---- Summarize DAG ----
print(dag_model)

# ---- Check DAG validity ----
if(is.dagitty(dag_model)) {
  cat("\n✅ DAG successfully loaded and parsed.\n")
} else {
  stop("❌ Error: DAG not recognized. Check syntax.")
}

# ---- Find Minimal Adjustment Set ----
adj_set <- adjustmentSets(dag_model, exposure = "ARG_VF_Abundance", outcome = "CRC_Stage")
cat("\n📊 Minimal Adjustment Set(s):\n")
print(adj_set)

# ---- Identify Testable Implications ----
implications <- impliedConditionalIndependencies(dag_model)
cat("\n🔍 Testable Implications:\n")
print(implications)

# ---- Plot the DAG (Publication Style) ----
pdf("Causal_DAG_ARG_VF_to_CRC.pdf", width = 7, height = 5)
plot(dag_model,
     coords = TRUE,
     main = "Causal DAG: ARG/VF Abundance → CRC Stage",
     cex = 1.2,
     col = "black")
dev.off()
cat("\n✅ DAG plot saved as 'Causal_DAG_ARG_VF_to_CRC.pdf'\n")

# ---- Export adjustment set to file ----
if (length(adj_set) > 0) {
  writeLines(paste("Minimal Adjustment Set(s):", 
                   paste(unlist(adj_set[[1]]), collapse = ", ")),
             con = "Minimal_Adjustment_Set.txt")
  cat("✅ Adjustment set saved to 'Minimal_Adjustment_Set.txt'\n")
}

# ---- (Optional) Visualize DAG interactively ----
# library(ggdag)
# ggdag(dag_model, layout = "nicely") +
#   theme_dag(base_size = 14) +
#   ggtitle("Causal Structure: ARG/VF → CRC Stage")

