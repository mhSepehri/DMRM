# Cross-cohort overlap of significant features under matched nominal p<0.05
# (Figure 5B Venn source; DMPs and DMRMs, discovery vs GSE88824).

library(ggVennDiagram); library(ggplot2)

DMP_DISC   <- "E:/backup/D/Uni/Thesis/Data/gse106648/DMP/normal_beRemoved/DMP_champ_norm_beRemoved.csv"
DMP_VALID  <- "E:/backup/D/Uni/Thesis/Data/ratio scenarios/validation/DMP_gse88824_norm_beR.csv"
DMRM_DISC  <- "E:/backup/D/Uni/Thesis/Data/ratio scenarios/12/New folder/dmrm_disc.csv"
DMRM_VALID <- "E:/backup/D/Uni/Thesis/Data/ratio scenarios/12/New folder/dmrm_gse88824.csv"
OUT        <- "E:/backup/D/Uni/Thesis/Data/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

P <- 0.05

sig_ids <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = TRUE)
  pcol <- grep("p[._]?val|pvalue|P\\.Value", names(d), ignore.case = TRUE, value = TRUE)
  pcol <- pcol[!grepl("adj|fdr|q", pcol, ignore.case = TRUE)][1]
  idcol <- grep("^X$|cpg|probe|marker|name|id", names(d), ignore.case = TRUE, value = TRUE)[1]
  ids <- if (is.na(idcol)) rownames(d) else as.character(d[[idcol]])
  unique(ids[which(as.numeric(d[[pcol]]) < P)])
}

dmp_d  <- sig_ids(DMP_DISC);  dmp_v  <- sig_ids(DMP_VALID)
dmrm_d <- sig_ids(DMRM_DISC); dmrm_v <- sig_ids(DMRM_VALID)

cat("DMPs  shared:", length(intersect(dmp_d, dmp_v)), "\n")
cat("DMRMs shared:", length(intersect(dmrm_d, dmrm_v)), "\n")

mk <- function(a, b, title)
  ggVennDiagram(list(GSE106648 = a, GSE88824 = b), label = "count", label_alpha = 0) +
    scale_fill_gradient(low = "#EAF2F8", high = "#AED6F1", guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    labs(title = title) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")

for (x in list(list(dmp_d, dmp_v, "Single-CpG DMPs", "Fig5B_venn_DMP"),
               list(dmrm_d, dmrm_v, "Candidate DMRMs", "Fig5B_venn_DMRM"))) {
  g <- mk(x[[1]], x[[2]], x[[3]])
  tiff(file.path(OUT, paste0(x[[4]], ".tiff")), width = 4.4, height = 4, units = "in",
       res = 300, compression = "lzw"); print(g); dev.off()
}
