# Composition-independence of DMRMs (discovery cohort GSE106648)
# Variance explained by cell composition, composition-vs-disease association,
# adjusted-model survival, and nested LRT of disease vs composition.

library(data.table)

PROJ  <- "E:/backup/D/Uni/Thesis/Data/ratio scenarios/12/New folder"
BETA  <- file.path(PROJ, "gse106648_Norm_ber_beta.csv")
S8    <- file.path(PROJ, "s8_results.csv")
PHENO <- file.path(PROJ, "results/stepA2/phenotype_GSE106648_FINAL.csv")
DECON <- file.path(PROJ, "results/deconv_gse106648.csv")   # 5 fractions per sample
OUT   <- file.path(PROJ, "results/step6")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

b <- fread(BETA, data.table = FALSE); rownames(b) <- b[[1]]; b[[1]] <- NULL
B <- as.matrix(b); storage.mode(B) <- "numeric"
ph <- read.csv(PHENO, stringsAsFactors = FALSE)
ph <- ph[match(colnames(B), ph$Beta_ID), ]
keep <- !is.na(ph$MS); B <- B[, keep]; ph <- ph[keep, ]
dc <- read.csv(DECON, stringsAsFactors = FALSE)
dc <- dc[match(colnames(B), dc[[1]]), ]
# adjust these to the champ.refbase column names in deconv_gse106648.csv
frac <- as.matrix(dc[, c("CD4T","CD8T","NK","Bcell","Mono")])

s8 <- read.csv(S8, row.names = 1, stringsAsFactors = FALSE)
sp <- do.call(rbind, strsplit(trimws(rownames(s8)), "[ /]+"))
s8$Num <- sp[,1]; s8$Den <- sp[,2]
s8 <- s8[s8$Num %in% rownames(B) & s8$Den %in% rownames(B), ]
sig <- if ("Significant" %in% names(s8)) s8[s8$Significant == "YES", ] else
       s8[s8$pvalue < 0.05 & abs(s8$log2FoldChange) > 0.5, ]

lr <- function(n, d) log2(B[n, ]) - log2(B[d, ])
M  <- t(mapply(lr, sig$Num, sig$Den))          # markers x samples

# 1. variance explained by composition (per-marker R2)
r2 <- apply(M, 1, function(y) summary(lm(y ~ frac))$r.squared)
cat("median R2 (composition):", round(median(r2), 3), "\n")

# 2. composition features vs disease
y <- ph$MS
p_cd4  <- summary(glm(y ~ frac[,"CD4T"], family = binomial))$coef[2,4]
p_cd8  <- summary(glm(y ~ frac[,"CD8T"], family = binomial))$coef[2,4]
p_rat  <- summary(glm(y ~ I(frac[,"CD4T"]/frac[,"CD8T"]), family = binomial))$coef[2,4]
cat("CD4 p:", signif(p_cd4,3), "| CD8 p:", signif(p_cd8,3), "| CD4/CD8 p:", signif(p_rat,3), "\n")

# 3. markers surviving adjustment
covs <- data.frame(ph[, c("Age","Sex","Smoking")], batch = ph$Batch)
adj_p <- function(extra) apply(M, 1, function(m) {
  d <- data.frame(m = m, MS = y, frac, extra)
  summary(lm(m ~ ., data = d))$coef["MS", 4]
})
q_comp <- p.adjust(adj_p(data.frame()), "BH")
q_full <- p.adjust(adj_p(covs), "BH")
cat("survive composition adj (q<0.05):", sum(q_comp < 0.05), "\n")
cat("survive full adj (q<0.05):", sum(q_full < 0.05), "\n")

# 4. nested LRT: disease vs composition variance contribution
dev_comp <- dev_dis <- numeric(nrow(M))
for (i in seq_len(nrow(M))) {
  base <- lm(M[i, ] ~ frac)
  ext  <- lm(M[i, ] ~ frac + y)
  dev_comp[i] <- summary(base)$r.squared
  dev_dis[i]  <- summary(ext)$r.squared - summary(base)$r.squared
}
cat("disease/composition variance ratio:", round(median(dev_dis)/median(dev_comp), 2), "\n")

write.csv(data.frame(marker = rownames(M), r2, q_comp, q_full),
          file.path(OUT, "composition_independence.csv"), row.names = FALSE)
