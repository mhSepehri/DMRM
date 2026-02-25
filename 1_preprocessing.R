library(data.table)
library(minfi)
library(ChAMP)

# Read intensity data
Methylated <- fread("data/gse106648/GSE106648_Matrix_methylated_signal_intensities.txt.gz")
rownames(Methylated) <- Methylated$ID_REF
Methylated[,1] <- NULL
Methylated <- data.matrix(Methylated)

Unmethylated <- fread("data/gse106648/GSE106648_Matrix_unmethylated_signal_intensities.txt.gz")
rownames(Unmethylated) <- Unmethylated$ID_REF
Unmethylated[,1] <- NULL
Unmethylated <- data.matrix(Unmethylated)

# Convert intensities to Beta values
intensities.to.beta <- function(M, U, alpha = 100){
  B <- M/(M+U+alpha)
  return(B)
}

B <- intensities.to.beta(Methylated, Unmethylated)
Beta_intensities <- round(B, digits = 2)
write.csv(Beta_intensities, "data/gse106648/beta_GSE106648.csv", row.names = TRUE)

# Reading IDAT data
path <- "data/MS-PBMC-Russia"
targets <- read.csv(file.path(path, "Sample_Sheet.csv"), as.is = TRUE, header = TRUE, sep = ",")
targets$ID <- file.path(path, targets$ID)

RGSet <- read.metharray(targets$ID, force=TRUE)
MSet <- preprocessRaw(RGSet) 
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
Beta_idat <- getBeta(RSet)
write.csv(Beta_idat, "data/MS-PBMC-Russia/beta_PBMC-Russia.csv", row.names = TRUE)

# Custom filtering setup
P <- fread('data/gse106648/GSE106648_p.txt', data.table = FALSE)
rownames(P) <- P$ID_REF
P$ID_REF <- NULL
names(P) <- gsub(".detectionPval", "", names(P)) 

bad.detected.probes <- function(P){
  return(unique(rownames(which(P > 0.05, arr.ind = TRUE))))
}

negative.probes <- function(M, U){
  return(union(unique(rownames(which(M < 0, arr.ind = TRUE))), unique(rownames(which(U < 0, arr.ind = TRUE)))))
}

rs_probes <- function(){
  rs <- read.csv('data/450k_manifest/rs_af0.05.csv')
  return(unique(rs$X))
}

non_specific.probes <- function(){
  ps <- read.csv('data/450k_manifest/non_specific_sites.csv')
  return(unique(ps[,2]))
}

# Aggregate probes to delete
deleted.probes <- function(M, U, P){
  high.p.probes <- bad.detected.probes(P)
  neg.probes <- negative.probes(M, U)
  nspec.probes <- non_specific.probes()
  snp.probes <- rs_probes()
  deleted.probes <- unique(c(high.p.probes, neg.probes, nspec.probes, snp.probes))
  return(deleted.probes)
}

# Apply custom filters
filtering_beta <- function(M, U, P, Beta_matrix){
  bad.probes <- deleted.probes(M, U, P)
  Beta_filtered <- Beta_matrix[!(rownames(Beta_matrix) %in% bad.probes), ]
  return(Beta_filtered)
}

Beta_custom_filtered <- filtering_beta(Methylated, Unmethylated, P, Beta_intensities)

# Final ChAMP filtering
Beta_Final <- champ.filter(Beta_custom_filtered, targets, M=NULL)