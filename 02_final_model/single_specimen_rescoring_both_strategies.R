## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
# Single-specimen rescoring of the 16 LoF/LoF specimens with the frozen final ensemble,
# for BOTH strategies (consensus30, diverse30). Each specimen is scored alone so that its
# MinProb draws do not depend on other specimens or on input order (same rule as the
# reproducibility package). No refitting; uses the functions of the user's scoring script.
# Run from the R project root (directory containing 20261029_rebuttal.Rproj).
suppressPackageStartupMessages({library(glmnet); library(openxlsx)})
ROOT <- REPO
OUT  <- file.path(REPO, "02_final_model", "results")
setwd(ROOT)
src <- file.path(ROOT, "02_final_model", "citrin_SLC25A13_LoFLoF_frozen_ensemble_scoring_20261104.R")
lines <- readLines(src, warn = FALSE, encoding = "UTF-8")
stopifnot(identical(trimws(lines[length(lines)]), "main()"))
eval(parse(text = lines[-length(lines)]))   # define functions and constants without running main()

root <- ROOT; data_dir <- DATA; model_file <- file.path(ROOT, "02_final_model", "FinalModel_reproducibility", "final_training_model_FROZEN_NOT_VALIDATED.rds")
model <- readRDS(model_file); validate_model(model)
input <- load_ss_features(data_dir)
cat("[MODEL] ", model_file, " md5 ", unname(tools::md5sum(model_file)), "\n", sep = "")

# joint (batch) scoring as in the original run
joint <- score_frozen_model(model, input)$scores
# single-specimen scoring
alone <- do.call(rbind, lapply(seq_len(nrow(input$ss)), function(i) {
  one <- input
  one$ss <- input$ss[i, , drop = FALSE]
  one$X  <- input$X[i, , drop = FALSE]
  for (nm in c("slc_ratio", "slc_detected", "n_features")) if (!is.null(input[[nm]]) && length(input[[nm]]) == nrow(input$ss)) one[[nm]] <- input[[nm]][i]
  score_frozen_model(model, one)$scores
}))
key <- c("SampleID", "strategy")
m <- merge(joint[, c(key, "batch", "probability", "threshold_full_inner_recall", "threshold_inner_spec90")],
           alone[, c(key, "probability")], by = key, suffixes = c("_joint", "_alone"))
m$diff <- m$probability_alone - m$probability_joint
m$call_full_joint <- as.integer(m$probability_joint >= m$threshold_full_inner_recall)
m$call_full_alone <- as.integer(m$probability_alone >= m$threshold_full_inner_recall)
m$call_spec90_joint <- as.integer(m$probability_joint >= m$threshold_inner_spec90)
m$call_spec90_alone <- as.integer(m$probability_alone >= m$threshold_inner_spec90)
m$call_changed <- (m$call_full_joint != m$call_full_alone) | (m$call_spec90_joint != m$call_spec90_alone)
m <- m[order(m$strategy, as.integer(sub("SS_", "", m$SampleID))), ]
write.csv(m, file.path(OUT, "CHECK_single_specimen_rescoring_both_strategies.csv"), row.names = FALSE)
# compare with the earlier consensus30 scored_alone column
chk <- read.csv(file.path(OUT, "CHECK_seed_sensitivity_consensus30.csv"))
cc <- m[m$strategy == "consensus30", ]
cat("consensus30: max |alone - earlier scored_alone| =", max(abs(round(cc$probability_alone, 3) - chk$scored_alone[match(cc$SampleID, chk$SampleID)])), "\n")
for (st in c("consensus30", "diverse30")) {
  x <- m[m$strategy == st, ]
  cat(sprintf("%s: max|diff| %.4f; calls changed: %d; positives full joint/alone %d/%d; spec90 joint/alone %d/%d; median joint %.3f alone %.3f\n",
      st, max(abs(x$diff)), sum(x$call_changed), sum(x$call_full_joint), sum(x$call_full_alone), sum(x$call_spec90_joint), sum(x$call_spec90_alone),
      median(x$probability_joint), median(x$probability_alone)))
  print(x[x$call_changed, c("SampleID", "probability_joint", "probability_alone")], row.names = FALSE)
}
# per-run medians (alone)
print(aggregate(probability_alone ~ strategy + batch, m, function(v) round(median(v), 3)))
