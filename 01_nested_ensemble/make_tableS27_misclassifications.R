## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
# Supplementary Table S27 (computed columns): false positives of the primary consensus-30 ensemble at the complete-inner-
# sensitivity threshold and the three true positives with the smallest margins, with the ACTB-normalized SLC25A13 intensity
# as a percentage of the median of the 27 disease and healthy controls. The clinical columns of the printed table (diagnosis,
# genotype, age, sex) were added from the clinical tables and are not generated here.
suppressPackageStartupMessages(library(openxlsx))
p <- read.csv(file.path(REPO, "01_nested_ensemble", "results", "outer_predictions.csv"), stringsAsFactors = FALSE)
p <- p[p$strategy == "consensus30" & p$policy == "full_inner_recall", ]; stopifnot(nrow(p) == 44, sum(p$truth) == 7)
si <- read.csv(file.path(DATA, "sample_info.csv"), stringsAsFactors = FALSE)
rw <- read.xlsx(file.path(DATA, "raw_data.xlsx"), sheet = "data", check.names = FALSE, sep.names = " ")
M <- as.matrix(rw[, si$SampleID]); storage.mode(M) <- "numeric"; sym <- rw[["Master Gene Symbol"]]
ratio <- M[which(sym == "SLC25A13")[1], ] / M[which(sym == "ACTB")[1], ]; names(ratio) <- si$SampleID
ctrl_med <- median(ratio[si$group %in% c("DD", "WW")])
p$margin <- p$probability - p$threshold
fp <- p[p$truth == 0 & p$predicted_class == 1, ]; fp <- fp[order(-fp$margin), ]
tp <- p[p$truth == 1 & p$predicted_class == 1, ]; tp <- tp[order(tp$margin), ][1:3, ]; tp <- tp[order(tp$SampleID), ]
out <- rbind(fp, tp)
tab <- data.frame(Specimen = out$SampleID, Group = out$group, Analytical_run = out$batch,
  SLC25A13_ACTB_percent_of_control_median = round(100 * ratio[out$SampleID] / ctrl_med), Out_of_fold_probability = round(out$probability, 3),
  Fold_specific_threshold = round(out$threshold, 3), Classification = ifelse(out$truth == 1, "true positive", "false positive"),
  Probability_minus_threshold = round(out$margin, 3), stringsAsFactors = FALSE)
print(tab, row.names = FALSE); cat("false positives:", nrow(fp), "; true positives:", sum(p$truth == 1 & p$predicted_class == 1), "of 7; false negatives:", sum(p$truth == 1 & p$predicted_class == 0), "\n")
write.csv(tab, file.path(REPO, "01_nested_ensemble", "results", "TableS27_misclassifications_values.csv"), row.names = FALSE)
