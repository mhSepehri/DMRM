library(data.table)
library(dplyr)
library(ggplot2)

BC_Norm_Beta <- data.matrix(fread("results/BC_Norm_Beta.csv", data.table = FALSE, row.names = 1))
targets <- read.csv("data/phenotype_data.csv")

control_cols <- which(targets$Sample_Group == "Control")
ms_cols <- which(targets$Sample_Group == "MS")

compute_ratios <- function(beta_matrix, num_cpgs, den_cpgs) {
  num_mat <- beta_matrix[rownames(beta_matrix) %in% num_cpgs, , drop = FALSE]
  den_mat <- beta_matrix[rownames(beta_matrix) %in% den_cpgs, , drop = FALSE]
  
  cpg_pairs <- expand.grid(Num = rownames(num_mat), Den = rownames(den_mat), stringsAsFactors = FALSE)
  
  ratio_matrix <- matrix(nrow = nrow(cpg_pairs), ncol = ncol(beta_matrix))
  rownames(ratio_matrix) <- paste(cpg_pairs$Num, cpg_pairs$Den, sep = "/")
  colnames(ratio_matrix) <- colnames(beta_matrix)
  
  for (i in seq_len(nrow(cpg_pairs))) {
    ratio_matrix[i, ] <- num_mat[cpg_pairs$Num[i], ] / den_mat[cpg_pairs$Den[i], ]
  }
  return(ratio_matrix)
}

Hyper_Healthy_CD4 <- readRDS("results/Hyper_Healthy_CD4.rds")
Hypo_Healthy_CD4 <- readRDS("results/Hypo_Healthy_CD4.rds")
Hyper_Healthy_CD8 <- readRDS("results/Hyper_Healthy_CD8.rds")
Hypo_Healthy_CD8 <- readRDS("results/Hypo_Healthy_CD8.rds")

Hyper_MS_CD4 <- readRDS("results/Hyper_MS_CD4.rds")
Hypo_MS_CD4 <- readRDS("results/Hypo_MS_CD4.rds")
Hyper_MS_CD8 <- readRDS("results/Hyper_MS_CD8.rds")
Hypo_MS_CD8 <- readRDS("results/Hypo_MS_CD8.rds")

DMRM1_1 <- compute_ratios(BC_Norm_Beta, Hyper_Healthy_CD4, Hypo_Healthy_CD4)
DMRM1_2 <- compute_ratios(BC_Norm_Beta, Hyper_Healthy_CD8, Hypo_Healthy_CD8)
DMRM1_3 <- compute_ratios(BC_Norm_Beta, Hyper_Healthy_CD4, Hypo_Healthy_CD8)
DMRM1_4 <- compute_ratios(BC_Norm_Beta, Hypo_Healthy_CD4, Hyper_Healthy_CD8)
Scenario_1 <- rbind(DMRM1_1, DMRM1_2, DMRM1_3, DMRM1_4)

Scenario_2 <- compute_ratios(BC_Norm_Beta, Hyper_MS_CD4, Hypo_MS_CD4)
Scenario_3 <- compute_ratios(BC_Norm_Beta, Hyper_MS_CD8, Hypo_MS_CD8)
Scenario_4 <- rbind(Scenario_2, Scenario_3)

Scenario_5 <- compute_ratios(BC_Norm_Beta, Hyper_MS_CD4, Hyper_MS_CD8)
Scenario_6 <- compute_ratios(BC_Norm_Beta, Hypo_MS_CD4, Hypo_MS_CD8)
Scenario_7 <- rbind(Scenario_5, Scenario_6)

Scenario_8 <- compute_ratios(BC_Norm_Beta, Hyper_MS_CD4, Hypo_MS_CD8)

dir.create("results/ratio_scenarios", showWarnings = FALSE)
write.csv(Scenario_8, "results/ratio_scenarios/ratio_markers_final_Scenario8.csv", row.names = TRUE)

log_ratio_markers <- log2(Scenario_8)

results <- data.frame(
  probename = rownames(log_ratio_markers),
  log2FoldChange = rowMeans(log_ratio_markers[, ms_cols]) - rowMeans(log_ratio_markers[, control_cols]),
  pvalue = apply(log_ratio_markers, 1, function(x) t.test(x[ms_cols], x[control_cols])$p.value)
)

results$Significant <- "NO"
results$Significant[abs(results$log2FoldChange) > 0.5 & results$pvalue < 0.05] <- "YES"
results$delabel <- NA
results$delabel[results$Significant == "YES"] <- results$probename[results$Significant == "YES"]

write.csv(results, "results/ratio_scenarios/results_Scenario8.csv", row.names = FALSE)

volcano <- ggplot(data = results, aes(x = log2FoldChange, y = -log10(pvalue), col = Significant)) +
  ggtitle("FC volcano Plot of ratio markers Scenario 8") +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  scale_color_manual(values = c("NO" = "darkgrey", "YES" = "darkmagenta"))

ggsave("results/ratio_scenarios/volcano_plot_scenario8.pdf", plot = volcano, width = 8, height = 6)