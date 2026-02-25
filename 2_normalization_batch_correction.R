
library(ChAMP)

# Assuming 'Beta_Final' and 'targets' are available from the preprocessing step
# If running sequentially, you can load them:
# Beta_Final <- data.matrix(read.csv("results/Beta_Final.csv", row.names = 1))
# targets <- read.csv("data/gse106648/Sample_Sheet.csv")

# BMIQ Normalization
Norm_Beta <- champ.norm(beta = Beta_Final,
                        resultsDir = "results/normalization",
                        method = "BMIQ",
                        plotBMIQ = TRUE,
                        arraytype = "450K",
                        cores = 16)

write.csv(Norm_Beta, "results/Norm_Beta.csv", row.names = TRUE)

# Singular Value Decomposition (SVD) before batch correction
champ.SVD(beta = Norm_Beta,
          rgSet = NULL,
          pd = targets,
          PDFplot = TRUE,
          Rplot = TRUE)

# ComBat Batch Effect Correction
BC_Norm_Beta <- champ.runCombat(beta = Norm_Beta,
                                pd = targets,
                                variablename = "Sample_Group", 
                                batchname = c("Sentrix_ID", "Sentrix_Position"))

write.csv(BC_Norm_Beta, "results/BC_Norm_Beta.csv", row.names = TRUE)
BC_Norm_Beta <- data.frame(BC_Norm_Beta)

# SVD after batch correction to verify effectiveness
champ.SVD(beta = BC_Norm_Beta,
          rgSet = NULL,
          pd = targets,  
          PDFplot = TRUE,
          Rplot = TRUE)