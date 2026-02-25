library(ChAMP)
library(MethylCIBERSORT)

# Load batch-corrected normalized beta values
Mat <- read.csv("results/BC_Norm_Beta.csv", as.is = TRUE, header = TRUE, sep = ",")
rownames(Mat) <- Mat$X
Mat[,1] <- NULL
BC_Norm_Beta <- data.matrix(Mat)

# Houseman Deconvolution
myRefBase <- champ.refbase(beta = BC_Norm_Beta, arraytype = "450K")
write.csv(myRefBase[["CellFraction"]], "results/deconv_gse106648.csv", row.names = TRUE)

# MethylCIBERSORT Deconvolution
data("StromalMatrix_V1")
table(Stromal_v1.pheno)
data("StromalMatrix_V2")
table(Stromal_v2.pheno)

# Intersect probes between dataset and stromal matrix
Int <- intersect(rownames(BC_Norm_Beta), rownames(Stromal_v2))
Mat_sub <- BC_Norm_Beta[match(Int, rownames(BC_Norm_Beta)), ]
Stromal_v2_sub <- Stromal_v2[match(Int, rownames(Stromal_v2)), ]

# Feature Selection
Signature <- FeatureSelect.V4(CellLines.matrix = NULL,
                              Heatmap = FALSE,
                              export = TRUE,
                              sigName = "Sig106648",
                              Stroma.matrix = Stromal_v2_sub,
                              deltaBeta = 0.2,
                              FDR = 0.01,
                              MaxDMRs = 100,
                              Phenotype.stroma = Stromal_v2.pheno)

# Prepare for CIBERSORT
Prep.CancerType(Beta = BC_Norm_Beta, 
                Probes = rownames(Signature$SignatureMatrix), 
                fname = "results/methylCiber106648")