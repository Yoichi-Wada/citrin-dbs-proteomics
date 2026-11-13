## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
# Two-stage evaluation on all 60 specimens.
# Step 1: SLC25A13 (Q9UJS0) not quantified (ACTB-normalised intensity == 0) -> call LoF/LoF.
# Step 2: specimens with SLC25A13 quantified -> LoF/Mis classifier out-of-fold call.
# Classifiers: (A) fixed30 primary, per seed (10 seeds, stage1 no-ComBat); (B) ensemble consensus30 primary; (C) ensemble diverse30 (secondary).
res <- file.path(REPO, "03_sensitivity", "two_stage", "results")
cl <- Sys.getenv("CITRIN_MULTISEED_RUNS", "")   # folder holding the 40 single-model runs (not distributed); fixed30 rows are skipped if empty
s60 <- read.csv(file.path(res, "slc25a13_per_specimen_60.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(s60) == 60)
s60$step1_call <- ifelse(s60$SLC25A13_detected, "proceed_to_classifier", "LoF/LoF")
cat("Step 1 (SLC25A13 not quantified -> LoF/LoF):\n"); print(table(truth = s60$group, step1 = s60$step1_call))

two_stage <- function(pred, label) {
  # pred: SampleID, predicted_class (1 = LoF/Mis) for the 44 quantified specimens
  d <- merge(s60, pred[, c("SampleID", "predicted_class")], by = "SampleID", all.x = TRUE)
  d$final_call <- ifelse(d$step1_call == "LoF/LoF", "LoF/LoF",
                   ifelse(d$predicted_class == 1, "LoF/Mis", "not_citrin_deficient"))
  d$truth3 <- ifelse(d$group == "SS", "LoF/LoF", ifelse(d$group == "SM", "LoF/Mis", "not_citrin_deficient"))
  cm <- table(truth = factor(d$truth3, levels = c("LoF/LoF", "LoF/Mis", "not_citrin_deficient")),
              call  = factor(d$final_call, levels = c("LoF/LoF", "LoF/Mis", "not_citrin_deficient")))
  ss_sens <- cm["LoF/LoF", "LoF/LoF"] / 16; sm_sens <- cm["LoF/Mis", "LoF/Mis"] / 7
  any_cd_sens <- (cm["LoF/LoF", "LoF/LoF"] + cm["LoF/LoF", "LoF/Mis"] + cm["LoF/Mis", "LoF/LoF"] + cm["LoF/Mis", "LoF/Mis"]) / 23
  spec <- cm["not_citrin_deficient", "not_citrin_deficient"] / 37
  acc <- sum(diag(cm)) / 60
  data.frame(classifier = label, SS_detected = cm["LoF/LoF", "LoF/LoF"], SM_detected = cm["LoF/Mis", "LoF/Mis"],
             controls_FP = 37 - cm["not_citrin_deficient", "not_citrin_deficient"],
             sens_LoFLoF = ss_sens, sens_LoFMis = sm_sens, sens_any_citrin_deficiency = round(any_cd_sens, 3),
             specificity_60 = round(spec, 3), accuracy_3class = round(acc, 3), stringsAsFactors = FALSE)
}
rows <- list(); cms <- list()
files <- if (nzchar(cl)) sort(Sys.glob(file.path(cl, "runs/44_fixed30_sample_seed*/result/two_run_nested_cv/*/stage1_no_combat/outer_predictions.csv"))) else character(0)
cat("\nfixed30 prediction files:", length(files), "\n")
for (f in files) { seed <- sub(".*seed(\\d+)/.*", "\\1", f); p <- read.csv(f, stringsAsFactors = FALSE)
  rows[[length(rows) + 1]] <- two_stage(p, paste0("fixed30_seed", seed)) }
ens <- read.csv(file.path(REPO, "01_nested_ensemble", "results", "outer_predictions.csv"), stringsAsFactors = FALSE)
ens <- ens[ens$policy == "full_inner_recall", ]
for (s in c("consensus30", "diverse30", "adaptive")) rows[[length(rows) + 1]] <- two_stage(ens[ens$strategy == s, ], paste0("ensemble_", s))
tab <- do.call(rbind, rows); rownames(tab) <- NULL
print(tab, row.names = FALSE)
fx <- tab[grepl("fixed30", tab$classifier), ]
if (nrow(fx)) cat(sprintf("\nfixed30 across 10 seeds: SM detected median %.0f (range %d–%d); controls FP median %.0f (range %d–%d); spec median %.3f; any-CD sens median %.3f (range %.3f–%.3f); 3-class accuracy median %.3f\n",
  median(fx$SM_detected), min(fx$SM_detected), max(fx$SM_detected), median(fx$controls_FP), min(fx$controls_FP), max(fx$controls_FP),
  median(fx$specificity_60), median(fx$sens_any_citrin_deficiency), min(fx$sens_any_citrin_deficiency), max(fx$sens_any_citrin_deficiency), median(fx$accuracy_3class)))
write.csv(tab, file.path(res, "two_stage_60_results.csv"), row.names = FALSE)
# full 3x3 matrices for the two primaries (median fixed30 seed by SM_detected, and ensemble consensus30)
p <- if (length(files)) read.csv(files[1], stringsAsFactors = FALSE) else NULL
d <- merge(s60, ens[ens$strategy == "consensus30", c("SampleID", "predicted_class")], by = "SampleID", all.x = TRUE)
d$final_call <- ifelse(d$step1_call == "LoF/LoF", "LoF/LoF", ifelse(d$predicted_class == 1, "LoF/Mis", "not_citrin_deficient"))
d$truth3 <- ifelse(d$group == "SS", "LoF/LoF", ifelse(d$group == "SM", "LoF/Mis", "not_citrin_deficient"))
cat("\n3x3 confusion, ensemble consensus30:\n"); print(table(truth = d$truth3, call = d$final_call))
write.csv(d[, intersect(c("SampleID", "group", "batch", "SLC25A13", "step1_call", "predicted_class", "final_call"), names(d))], file.path(res, "two_stage_60_per_specimen_consensus30.csv"), row.names = FALSE)
