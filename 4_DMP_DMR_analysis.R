library(DMRcate)
library(DMRcatedata)
library(ChAMP)
library(data.table)

# PBMC DMP and DMR Analysis
BC_Norm_Beta_PBMC <- data.matrix(fread("results/BC_Norm_Beta.csv", data.table = FALSE, row.names = 1))
targets_PBMC <- read.csv("data/phenotype_data.csv")
design_PBMC <- model.matrix(~factor(targets_PBMC$Sample_Group))

dmp_DMRcate <- cpg.annotate(datatype = "array", 
                            BC_Norm_Beta_PBMC, 
                            what = "Beta", 
                            arraytype = "450K", 
                            analysis.type = "differential",
                            adjPVal = 0.05, 
                            design = design_PBMC, 
                            coef = 2)

DMP_DMRcate <- as.data.frame(dmp_DMRcate@ranges@elementMetadata@listData)
DMP_DMRcate <- cbind(Probe = dmp_DMRcate@ranges@ranges@NAMES, DMP_DMRcate)
write.csv(DMP_DMRcate, "results/DMP_DMRcate_norm_bc.csv", row.names = FALSE)

DMR_DMRcate <- dmrcate(dmp_DMRcate, lambda = 1000, C = 2, pcutoff = "fdr")
write.csv(as.data.frame(extractRanges(DMR_DMRcate, genome = "hg19")), "results/DMR_DMRcate_norm_bc.csv", row.names = FALSE)

# Extracting CD4+ MS vs HC DMPs for DMRM Scenarios
beta_CD4 <- data.matrix(fread("data/GSE130029_RAW_CD4_data/BC_Norm_Beta_CD4.csv", data.table = FALSE, row.names = 1))
targets_CD4 <- read.csv("data/GSE130029_RAW_CD4/Sample_Sheet.csv")
myDMP_CD4 <- champ.DMP(beta = beta_CD4, pheno = targets_CD4$Sample_Group, adjPVal = 0.05)
dmp_CD4_df <- myDMP_CD4[[1]]

saveRDS(rownames(dmp_CD4_df[dmp_CD4_df$logFC > 0, ]), "results/Hyper_MS_CD4.rds")
saveRDS(rownames(dmp_CD4_df[dmp_CD4_df$logFC < 0, ]), "results/Hypo_MS_CD4.rds")

# Extracting CD8+ MS vs HC DMPs for DMRM Scenarios
beta_CD8 <- data.matrix(fread("data/GSE130030_RAW_CD8_data/BC_Norm_Beta_CD8.csv", data.table = FALSE, row.names = 1))
targets_CD8 <- read.csv("data/GSE130030_RAW_CD8/Sample_Sheet.csv")
myDMP_CD8 <- champ.DMP(beta = beta_CD8, pheno = targets_CD8$Sample_Group, adjPVal = 0.05)
dmp_CD8_df <- myDMP_CD8[[1]]

saveRDS(rownames(dmp_CD8_df[dmp_CD8_df$logFC > 0, ]), "results/Hyper_MS_CD8.rds")
saveRDS(rownames(dmp_CD8_df[dmp_CD8_df$logFC < 0, ]), "results/Hypo_MS_CD8.rds")

# Extracting Healthy CD4+ vs Healthy CD8+ DMPs for DMRM Scenarios
hc_CD4_cols <- targets_CD4$Sample_Group == "Control"
hc_CD8_cols <- targets_CD8$Sample_Group == "Control"
beta_HC <- cbind(beta_CD4[, hc_CD4_cols], beta_CD8[, hc_CD8_cols])
pheno_HC <- c(rep("CD4", sum(hc_CD4_cols)), rep("CD8", sum(hc_CD8_cols)))

myDMP_HC <- champ.DMP(beta = beta_HC, pheno = pheno_HC, adjPVal = 0.05)
dmp_HC_df <- myDMP_HC[[1]]

saveRDS(rownames(dmp_HC_df[dmp_HC_df$logFC > 0, ]), "results/Hyper_Healthy_CD4.rds")
saveRDS(rownames(dmp_HC_df[dmp_HC_df$logFC < 0, ]), "results/Hypo_Healthy_CD4.rds")
saveRDS(rownames(dmp_HC_df[dmp_HC_df$logFC < 0, ]), "results/Hyper_Healthy_CD8.rds")
saveRDS(rownames(dmp_HC_df[dmp_HC_df$logFC > 0, ]), "results/Hypo_Healthy_CD8.rds")