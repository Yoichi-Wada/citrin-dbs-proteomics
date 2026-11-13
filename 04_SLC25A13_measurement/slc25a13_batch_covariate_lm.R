## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# SLC25A13/ACTB log2 ratio ~ group + batch (Run 2 / Run 3), 44 non-LoF/LoF specimens with SLC25A13 quantified.
# Reproduces the figures quoted in Response R2 Major 4 "Findings with batch modelled as a covariate"
# (carrier -1.55, fold 0.34; LoF/Mis -1.84, fold 0.28; Run 3 -0.005, P = 0.98). Reference group = disease + healthy controls.
# Written 2026-11-09 to document provenance; original values were computed inline on 2026-11-05 and recomputed after the WW fix on 2026-11-07 (no script was saved then).
suppressPackageStartupMessages(library(openxlsx))
dd <- DATA
si <- read.csv(file.path(dd, "sample_info.csv"), stringsAsFactors = FALSE)
rw <- read.xlsx(file.path(dd, "raw_data.xlsx"), sheet = "data", check.names = FALSE, sep.names = " ")
M <- as.matrix(rw[, si$SampleID]); storage.mode(M) <- "numeric"
acc <- rw[["Master accession number"]]; sym <- rw[["Master Gene Symbol"]]
slc <- which(sym == "SLC25A13"); actb <- which(sym == "ACTB")
cat("SLC25A13 idx:", slc, acc[slc], " ACTB idx:", actb, acc[actb], "\n")
r <- log2(M[slc[1], ] / M[actb[1], ])
d <- data.frame(SampleID = si$SampleID, group = si$group, batch = si$batch, lr = r)
d44 <- d[!is.na(d$lr) & is.finite(d$lr) & d$group != "SS", ]
cat("n quantified non-SS:", nrow(d44), "\n"); print(table(d44$group, d44$batch))
d44$grp <- factor(ifelse(d44$group %in% c("DD","WW"), "control", d44$group), levels = c("control","WS","SM"))
d44$batch <- factor(d44$batch)
cat("\n--- lm(log2 ratio ~ group + batch) ---\n")
m1 <- lm(lr ~ grp + batch, d44); print(round(summary(m1)$coefficients, 4)); cat("fold:", round(2^coef(m1)[2:3], 3), "\n")
cat("\n--- unadjusted ---\n")
m0 <- lm(lr ~ grp, d44); print(round(summary(m0)$coefficients, 4)); cat("fold:", round(2^coef(m0)[2:3], 3), "\n")
cat("\n--- 4-level group + batch ---\n")
d44$g4 <- factor(d44$group, levels=c("DD","WW","WS","SM")); print(round(summary(lm(lr ~ g4 + batch, d44))$coefficients,4))
cat("\n--- run-stratified medians x1e4 ---\n"); print(aggregate(2^lr*1e4 ~ group + batch, d44, function(x) round(median(x),2)))
