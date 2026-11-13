## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
dd <- DATA
si <- read.csv(file.path(dd, "sample_info.csv"), stringsAsFactors = FALSE)
fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE)
x <- read.csv(file.path(EXT, "DBS_nACTB_analysis_data_20250808.csv"), check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(x)[1] <- "Accession"
getrow <- function(acc) { r <- x[x$Accession == acc, ]; v <- as.numeric(unlist(r[, si$SampleID])); names(v) <- si$SampleID; v }
a13 <- getrow("Q9UJS0"); a12 <- getrow("O75746"); actb_present <- any(x$Accession == "P60709")
cat("ACTB row present in table:", actb_present, "\n")
cat("SLC25A13 peptides / unique:", x[x$Accession=="Q9UJS0", c("Peptides","Unique peptides")][1,1], "/", x[x$Accession=="Q9UJS0", c("Peptides","Unique peptides")][1,2], "\n")
cat("Representation of non-detection in CSV: zeros =", sum(a13 == 0), ", NA =", sum(is.na(a13)), "\n")
si$SLC25A13 <- a13[si$SampleID]; si$SLC25A12 <- a12[si$SampleID]
si$SLC25A13_detected <- si$SLC25A13 > 0
cat("\n=== SLC25A13 detection by group x batch (detected / n)\n")
print(with(si, tapply(SLC25A13_detected, list(group, batch), function(z) sprintf("%d/%d", sum(z), length(z)))))
cat("\n=== SLC25A12 non-detections:", sum(si$SLC25A12 == 0 | is.na(si$SLC25A12)), "\n")
cat("\n=== SLC25A13 / ACTB (x1e4) by group: median (min–max), ratio to DD median\n")
med_dd <- median(si$SLC25A13[si$group == "DD"])
for (g in c("DD","WW","WS","SM","SS")) { v <- si$SLC25A13[si$group == g]
  cat(sprintf("%-3s n=%2d  %.2f (%.2f–%.2f)  ratio_to_DD=%.2f\n", g, length(v), 1e4*median(v), 1e4*min(v), 1e4*max(v), median(v)/med_dd)) }
cat("\n=== Ratio to control (DD+WW pooled) with bootstrap 95% CI of ratio of medians (B=5000)\n")
set.seed(1); ctrl <- si$SLC25A13[si$group %in% c("DD","WW")]
for (g in c("WS","SM")) { v <- si$SLC25A13[si$group == g]
  r <- replicate(5000, median(sample(v, replace=TRUE))/median(sample(ctrl, replace=TRUE)))
  cat(sprintf("%s vs control: ratio of medians %.2f (95%% CI %.2f–%.2f); Wilcoxon p=%.2g\n", g, median(v)/median(ctrl), quantile(r,.025), quantile(r,.975), wilcox.test(v, ctrl)$p.value)) }
cat(sprintf("SM vs WS: Wilcoxon p=%.2g\n", wilcox.test(si$SLC25A13[si$group=="SM"], si$SLC25A13[si$group=="WS"])$p.value))
# raw (non-normalized) check if raw_data.xlsx readable
if (requireNamespace("openxlsx", quietly = TRUE)) {
  sh <- openxlsx::getSheetNames(file.path(dd, "raw_data.xlsx")); cat("\nraw_data.xlsx sheets:", paste(sh, collapse=", "), "\n")
  rw <- openxlsx::read.xlsx(file.path(dd, "raw_data.xlsx"), sheet = 1, check.names = FALSE)
  names(rw)[1] <- "Accession"; r <- rw[rw$Accession == "Q9UJS0", ]
  if (nrow(r)) { v <- unlist(r[, intersect(si$SampleID, names(rw))]); cat("raw SLC25A13: n cols", length(v), "; zeros", sum(as.numeric(v)==0, na.rm=TRUE), "; NA", sum(is.na(v)), "\n"); print(v[grep("^SS", names(v))]) }
}
m <- merge(si, fm, by = "SampleID")
write.csv(m[, c("SampleID","group","batch","SLC25A13","SLC25A12","SLC25A13_detected")],
  file.path(REPO, "03_sensitivity", "two_stage", "results", "slc25a13_per_specimen_60.csv"), row.names = FALSE)
