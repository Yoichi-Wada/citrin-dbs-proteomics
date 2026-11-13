## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
# Supplementary Table: discriminative performance of SLC25A13 alone (ACTB-normalized), five comparisons,
# with family-cluster bootstrap 95% CIs; plus the same comparisons under the four alternative
# normalization strategies as sensitivity rows (point estimates). Descriptive, not screening performance.
suppressPackageStartupMessages({library(openxlsx)})
setwd(REPO)
D <- DATA
d  <- read.xlsx(file.path(D, "raw_data.xlsx"), sheet = "data", sep.names = " ")
si <- read.csv(file.path(D, "sample_info.csv")); fm <- read.csv(file.path(EXT, "family_map.csv"))
sc <- grep("^(SS|SM|WS|WW|DD)_[0-9]+$", names(d), value = TRUE)
acc <- d[["Master accession number"]]; up <- d[["Unique peptides"]]
M <- as.matrix(d[, sc]); storage.mode(M) <- "numeric"; M[is.na(M)] <- 0; rownames(M) <- acc
x <- M["Q9UJS0", ]; actb <- M["P60709", ]
g <- setNames(si$group, si$SampleID)[sc]; f <- setNames(fm$family_id, fm$SampleID)[sc]
# normalization variants (log2 scale), as in actb_normalization_sensitivity: complete-case proteins with >=3 unique peptides
cc <- up >= 3 & rowSums(M > 0) == 60
L <- log2(M); L[!is.finite(L)] <- NA
norm <- list(
  `ACTB (primary)`        = log2(x) - log2(actb),
  `Sample median`         = log2(x) - apply(L[cc, ], 2, median),
  `Total intensity`       = log2(x) - log2(colSums(M[cc, ])),
  `None (raw intensity)`  = log2(x))
pg <- read.delim(file.path(D, "DBS-Tohoku.pg_matrix.tsv"), check.names = FALSE)
pgs <- names(pg)[sapply(pg, is.numeric)]; v <- as.numeric(pg[pg$Protein.Group == "Q9UJS0", pgs]); names(v) <- sub(".*[\\\\/]", "", sub("\\.raw$|\\.d$", "", pgs))
v <- v[sc]; v[is.na(v)] <- 0; norm[["DIA-NN MaxLFQ"]] <- log2(v)
for (nm in names(norm)) { z <- norm[[nm]]; z[!is.finite(z)] <- -Inf; norm[[nm]] <- z }  # not quantified -> lowest score

auc <- function(s, y) { rk <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(rk[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
boot_auc <- function(score, keep, pos, B = 2000, seed = 20261106) {
  s <- -score[keep]; y <- as.integer(pos[keep]); ff <- f[keep]; obs <- auc(s, y)
  set.seed(seed); fams <- unique(ff)
  b <- replicate(B, { fs <- sample(fams, replace = TRUE); idx <- unlist(lapply(fs, function(k) which(ff == k))); yy <- y[idx]
    if (sum(yy) > 0 && sum(yy == 0) > 0) auc(s[idx], yy) else NA })
  c(obs = obs, lo = unname(quantile(b, 0.025, na.rm = TRUE)), hi = unname(quantile(b, 0.975, na.rm = TRUE)), valid = sum(!is.na(b)))
}
comps <- list(
  list(label = "LoF/LoF vs all other participants", keep = rep(TRUE, 60), pos = g == "SS", n1 = 16, n0 = 44),
  list(label = "Reduced (carriers + LoF/Mis) vs normal (disease controls + healthy controls)", keep = g != "SS", pos = g %in% c("WS", "SM"), n1 = 17, n0 = 27),
  list(label = "Any citrin deficiency (LoF/LoF + LoF/Mis) vs all other participants", keep = rep(TRUE, 60), pos = g %in% c("SS", "SM"), n1 = 23, n0 = 37),
  list(label = "LoF/Mis vs all non-LoF/LoF participants", keep = g != "SS", pos = g == "SM", n1 = 7, n0 = 37),
  list(label = "LoF/Mis vs carriers", keep = g %in% c("SM", "WS"), pos = g == "SM", n1 = 7, n0 = 10))
rows <- lapply(comps, function(cm) {
  pr <- boot_auc(norm[["ACTB (primary)"]], cm$keep, cm$pos)
  alt <- sapply(names(norm)[-1], function(nm) auc(-norm[[nm]][cm$keep], as.integer(cm$pos[cm$keep])))
  data.frame(Comparison = cm$label, `Positive group (n)` = cm$n1, `Negative group (n)` = cm$n0,
             `AUC, ACTB-normalized (95% CI)` = sprintf("%.3f (%.3f–%.3f)", pr["obs"], pr["lo"], pr["hi"]),
             `AUC, sample-median normalization` = sprintf("%.3f", alt["Sample median"]),
             `AUC, total-intensity normalization` = sprintf("%.3f", alt["Total intensity"]),
             `AUC, DIA-NN MaxLFQ` = sprintf("%.3f", alt["DIA-NN MaxLFQ"]),
             `AUC, unnormalized` = sprintf("%.3f", alt["None (raw intensity)"]), check.names = FALSE)
})
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)

title <- "Table S[x]. Discriminative performance of SLC25A13 abundance alone."
legend <- c(
  "Area under the ROC curve (AUC) for SLC25A13 abundance as a single marker, with lower abundance indicating the positive (disease or reduced) group. The primary column uses the ACTB-normalized abundance (log2 of the SLC25A13 intensity divided by the ACTB intensity in the same specimen); specimens in which SLC25A13 was not quantified (all 16 LoF/LoF specimens) were assigned the lowest rank. 95% confidence intervals are from a family-cluster bootstrap (families resampled with replacement, 2,000 replicates). The remaining columns give point estimates under the alternative normalization strategies described in the Methods (sample-median centering and total-intensity scaling computed on the 3,325 protein groups with at least three unique peptides quantified in all 60 specimens; DIA-NN MaxLFQ protein quantities; unnormalized raw intensity).",
  "The AUC of 1.000 for LoF/LoF versus all other participants, and for reduced versus normal groups, reflects complete separation by construction (SLC25A13 was quantified in 0 of 16 LoF/LoF specimens and in 44 of 44 other specimens) and is equivalent to the detection frequencies reported in the Results; the corresponding confidence intervals are degenerate. These values describe a clinically ascertained cohort in which the control groups differ from the affected groups in age and clinical state; they are not estimates of screening performance. The comparison of LoF/Mis individuals with carriers shows that SLC25A13 abundance alone does not separate these groups.")
dir.create(file.path(REPO, "04_SLC25A13_measurement", "results"), showWarnings = FALSE)
wb <- createWorkbook(); addWorksheet(wb, "Table S SLC25A13 AUC")
hs <- createStyle(textDecoration = "bold", border = "TopBottom", wrapText = TRUE, valign = "center")
writeData(wb, 1, title, startRow = 1); addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
writeData(wb, 1, tab, startRow = 3, headerStyle = hs)
for (i in seq_along(legend)) writeData(wb, 1, legend[i], startRow = 3 + nrow(tab) + 1 + i)
setColWidths(wb, 1, cols = 1:ncol(tab), widths = c(52, 12, 12, 26, 16, 16, 14, 14))
addStyle(wb, 1, createStyle(halign = "center"), rows = 4:(3 + nrow(tab)), cols = 2:ncol(tab), gridExpand = TRUE, stack = TRUE)
saveWorkbook(wb, file.path(REPO, "04_SLC25A13_measurement", "results", "TableS_SLC25A13_alone_AUC.xlsx"), overwrite = TRUE)
mdrow <- function(r) paste0("| ", paste(r, collapse = " | "), " |")
writeLines(c(paste0("# ", title), "", mdrow(names(tab)), paste0("|", paste(rep("---", ncol(tab)), collapse = "|"), "|"), apply(tab, 1, mdrow), "", legend),
           file.path(REPO, "04_SLC25A13_measurement", "results", "TableS_SLC25A13_alone_AUC.md"))
write.csv(tab, file.path(REPO, "04_SLC25A13_measurement", "results", "TableS_SLC25A13_alone_AUC_values.csv"), row.names = FALSE)
cat("written\n")
