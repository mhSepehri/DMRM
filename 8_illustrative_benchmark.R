# Illustrative low-plex classification benchmark (Supporting Information).
# Leakage-free: 10-marker DMRM and DMP panels selected and fitted on discovery
# (GSE106648), applied unchanged to three external cohorts. Plain logistic
# regression; features present in all cohorts (no imputation); per-cohort
# label-free operating point at the median predicted score. Outputs Sup. Table 2
# and the Figure 5B AUC barplot.

library(data.table); library(pROC); library(ggplot2)

PROJ    <- "E:/backup/D/Uni/Thesis/Data/ratio scenarios/12/New folder"
BETA    <- file.path(PROJ, "gse106648_Norm_ber_beta.csv")
S8      <- file.path(PROJ, "s8_results.csv")
DMP_106 <- "E:/backup/D/Uni/Thesis/Data/gse106648/DMP/normal_beRemoved/DMP_champ_norm_beRemoved.csv"
A2      <- file.path(PROJ, "results/stepA2")
A3B     <- file.path(PROJ, "results/stepA3b")
OUT     <- file.path(PROJ, "results/step8")
FIG     <- "E:/backup/D/Uni/Thesis/Data/figures"
K       <- 10
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
set.seed(2026)

norm_id <- function(x) gsub("[^A-Za-z0-9]", "", as.character(x))
Mv <- function(v) log2(v/(1-v))
lr <- function(n, d, M) log2(M[n, ]) - log2(M[d, ])

d <- fread(BETA, data.table = FALSE); rownames(d) <- d[[1]]; d[[1]] <- NULL
B <- as.matrix(d); storage.mode(B) <- "numeric"
ph <- read.csv(file.path(A2, "phenotype_GSE106648_FINAL.csv"), stringsAsFactors = FALSE)
ph <- ph[match(colnames(B), ph$Beta_ID), ]; ok <- !is.na(ph$MS)
B <- B[, ok]; grp <- ph$MS[ok]

load_ext <- function(nm) {
  f <- file.path(A3B, paste0(nm, "_beta.csv")); if (!file.exists(f)) return(NULL)
  b <- fread(f, data.table = FALSE); rownames(b) <- b[[1]]; b[[1]] <- NULL
  b <- as.matrix(b); storage.mode(b) <- "numeric"
  p <- read.csv(file.path(A3B, paste0(nm, "_pheno.csv")), stringsAsFactors = FALSE)
  key <- NULL
  for (cc in intersect(c("Sample_Label","ID","Sample","Sample_Name"), names(p)))
    if (sum(norm_id(colnames(b)) %in% norm_id(p[[cc]])) > .8*ncol(b)) { key <- cc; break }
  if (is.null(key)) return(NULL)
  p <- p[match(norm_id(colnames(b)), norm_id(p[[key]])), ]
  k <- !is.na(p$MS); list(B = b[, k, drop = FALSE], g = p$MS[k])
}
EXT <- list()
for (nm in c("GSE88824","Kulakova","Twins")) { e <- load_ext(nm); if (!is.null(e)) EXT[[nm]] <- e }
shared <- rownames(B); for (nm in names(EXT)) shared <- intersect(shared, rownames(EXT[[nm]]$B))

# diversified DMRM panel: one per numerator, largest |log2FC| first
s8 <- read.csv(S8, row.names = 1, stringsAsFactors = FALSE)
sp <- do.call(rbind, strsplit(trimws(rownames(s8)), "[ /]+"))
s8$Num <- sp[,1]; s8$Den <- sp[,2]
s8 <- s8[s8$Num %in% shared & s8$Den %in% shared, ]
s8$sig <- if ("Significant" %in% names(s8)) s8$Significant == "YES" else
          (s8$pvalue < 0.05 & abs(s8$log2FoldChange) > 0.5)
s8 <- s8[order(!s8$sig, -abs(s8$log2FoldChange)), ]
selR <- s8[!duplicated(s8$Num), ][1:K, c("Num","Den")]

# DMP panel: significant, largest |logFC|
dmp <- read.csv(DMP_106, stringsAsFactors = FALSE, check.names = TRUE)
pref <- sub("\\.logFC$", "", grep("\\.logFC$", names(dmp), value = TRUE)[1])
ids <- as.character(dmp[[1]]); fc <- dmp[[paste0(pref,".logFC")]]
pv <- dmp[[paste0(pref,".P.Value")]]; pa <- dmp[[paste0(pref,".adj.P.Val")]]
keep <- ids %in% shared; ids <- ids[keep]; fc <- fc[keep]; pv <- pv[keep]; pa <- pa[keep]
sm <- if (sum(pa < 0.05, na.rm = TRUE) >= K) pa < 0.05 else pv < 0.05
ds <- data.frame(id = ids, fc = fc)[sm & !is.na(sm), ]
selC <- ds$id[order(-abs(ds$fc))][1:K]

Xr <- t(apply(selR, 1, function(r) lr(r["Num"], r["Den"], B)))
rownames(Xr) <- paste(selR$Num, selR$Den, sep = "/")
Xc <- Mv(B[selC, , drop = FALSE])

fit_glm <- function(Xtr, y) {
  mu <- rowMeans(Xtr); s <- apply(Xtr, 1, sd); s[s == 0] <- 1
  Z <- scale(t(Xtr), mu, s)
  df <- data.frame(y = y, Z); colnames(df) <- c("y", paste0("f", seq_len(ncol(Z))))
  list(model = glm(y ~ ., data = df, family = binomial()), mu = mu, sd = s)
}
mR <- fit_glm(Xr, grp); mC <- fit_glm(Xc, grp)
predp <- function(m, X) {
  Z <- scale(t(X), m$mu, m$sd); df <- as.data.frame(Z)
  colnames(df) <- paste0("f", seq_len(ncol(df)))
  as.numeric(predict(m$model, df, type = "response"))
}
build_X <- function(feats, type, beta) {
  if (type == "DMRM") { sp <- do.call(rbind, strsplit(feats, "/"))
    X <- matrix(NA_real_, nrow(sp), ncol(beta), dimnames = list(feats, colnames(beta)))
    for (i in seq_len(nrow(sp))) if (all(sp[i,] %in% rownames(beta)))
      X[i,] <- log2(beta[sp[i,1],]) - log2(beta[sp[i,2],])
  } else { X <- matrix(NA_real_, length(feats), ncol(beta), dimnames = list(feats, colnames(beta)))
    for (i in seq_along(feats)) if (feats[i] %in% rownames(beta)) X[i,] <- Mv(beta[feats[i],]) }
  miss <- rowSums(is.na(X)) > 0
  for (i in which(miss)) X[i,] <- if (type == "DMRM") mR$mu[i] else mC$mu[i]
  X
}
met <- function(y, pp) {
  ci <- as.numeric(ci.auc(roc(y, pp, quiet = TRUE, direction = "<")))
  cl <- as.integer(pp >= median(pp))
  c(auc = ci[2], lo = ci[1], hi = ci[3], sens = mean(cl[y==1]==1), spec = mean(cl[y==0]==0))
}
rows <- list()
for (nm in names(EXT)) {
  E <- EXT[[nm]]
  a <- met(E$g, predp(mR, build_X(rownames(Xr), "DMRM", E$B)))
  b <- met(E$g, predp(mC, build_X(selC, "DMP", E$B)))
  rows[[nm]] <- data.frame(cohort = nm,
    dmrm_auc = a["auc"], dmrm_lo = a["lo"], dmrm_hi = a["hi"], dmrm_sens = a["sens"], dmrm_spec = a["spec"],
    dmp_auc = b["auc"], dmp_lo = b["lo"], dmp_hi = b["hi"], dmp_sens = b["sens"], dmp_spec = b["spec"])
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
write.csv(tab, file.path(OUT, "sup_table_2_benchmark.csv"), row.names = FALSE)
print(tab)

# Figure 5B: AUC barplot
df <- data.frame(
  cohort = factor(rep(tab$cohort, each = 2), levels = tab$cohort),
  panel  = factor(rep(c("DMRM","DMP"), nrow(tab)), levels = c("DMRM","DMP")),
  auc = as.vector(t(tab[, c("dmrm_auc","dmp_auc")])),
  lo  = as.vector(t(tab[, c("dmrm_lo","dmp_lo")])),
  hi  = as.vector(t(tab[, c("dmrm_hi","dmp_hi")])))
g <- ggplot(df, aes(cohort, auc, fill = panel)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_col(position = position_dodge(0.72), width = 0.62, colour = "grey20", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.72), width = 0.18, linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", auc)), position = position_dodge(0.72), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c(DMRM = "#C0392B", DMP = "#2E5FA3"), name = NULL) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "AUC") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        axis.text = element_text(colour = "black"))
tiff(file.path(FIG, "Fig5B_AUC_barplot.tiff"), width = 4.6, height = 3.6, units = "in",
     res = 300, compression = "lzw"); print(g); dev.off()
ggsave(file.path(FIG, "Fig5B_AUC_barplot.pdf"), g, width = 4.6, height = 3.6)
