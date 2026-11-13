## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Regenerate (WW-corrected model, 2026-11-07) the LoF/LoF reference files that the 2026-11-04 ad hoc run produced:
# CHECK_single_specimen_rescoring_both_strategies.csv, CHECK_seed_sensitivity_consensus30.csv,
# REFERENCE_OOF_vs_LoFLoF_comparison.csv, TableSx_LoFLoF_frozen_ensemble_scores.csv. No refitting.
suppressPackageStartupMessages({library(glmnet); library(openxlsx)})
ROOT <- REPO
OUT  <- file.path(REPO, "02_final_model", "results")
RES  <- file.path(OUT, "LoFLoF_frozen_ensemble")
PROD <- file.path(REPO, "01_nested_ensemble", "results")
setwd(ROOT)
src <- file.path(ROOT, "02_final_model", "citrin_SLC25A13_LoFLoF_frozen_ensemble_scoring_20261104.R")
lines <- readLines(src, warn = FALSE, encoding = "UTF-8"); stopifnot(identical(trimws(lines[length(lines)]), "main()"))
eval(parse(text = lines[-length(lines)]))
root <- ROOT; data_dir <- DATA; model_file <- file.path(ROOT, "02_final_model", "FinalModel_reproducibility", "final_training_model_FROZEN_NOT_VALIDATED.rds")
model <- readRDS(model_file); validate_model(model); input <- load_ss_features(data_dir)
cat("[MODEL] ", model_file, " md5 ", unname(tools::md5sum(model_file)), "\n", sep = "")
score_alone <- function(input) do.call(rbind, lapply(seq_len(nrow(input$ss)), function(i) {
  one <- input; one$ss <- input$ss[i, , drop = FALSE]; one$X <- input$X[i, , drop = FALSE]
  for (nm in c("slc_ratio", "slc_detected", "n_features")) if (!is.null(input[[nm]]) && length(input[[nm]]) == nrow(input$ss)) one[[nm]] <- input[[nm]][i]
  score_frozen_model(model, one)$scores }))
joint <- score_frozen_model(model, input)$scores
alone <- score_alone(input)
key <- c("SampleID", "strategy")
m <- merge(joint[, c(key, "batch", "probability", "threshold_full_inner_recall", "threshold_inner_spec90")], alone[, c(key, "probability")], by = key, suffixes = c("_joint", "_alone"))
m$diff <- m$probability_alone - m$probability_joint
m$call_full_joint <- as.integer(m$probability_joint >= m$threshold_full_inner_recall); m$call_full_alone <- as.integer(m$probability_alone >= m$threshold_full_inner_recall)
m$call_spec90_joint <- as.integer(m$probability_joint >= m$threshold_inner_spec90); m$call_spec90_alone <- as.integer(m$probability_alone >= m$threshold_inner_spec90)
m$call_changed <- (m$call_full_joint != m$call_full_alone) | (m$call_spec90_joint != m$call_spec90_alone)
m <- m[order(m$strategy, as.integer(sub("SS_", "", m$SampleID))), ]
# seed sensitivity (consensus30, scored alone, imputation seed offsets)
orig_seed_fn <- model_scoring_seed
offs <- c(0, 1000, 2000, 3000, 4000); seedmat <- NULL
for (o in offs) {
  model_scoring_seed <- function(member_seed) orig_seed_fn(member_seed) + o
  a <- score_alone(input); a <- a[a$strategy == "consensus30", ]
  seedmat <- cbind(seedmat, setNames(a$probability[match(paste0("SS_", 1:16), a$SampleID)], paste0("seed_off_", o)))
}
model_scoring_seed <- orig_seed_fn
cs <- data.frame(SampleID = paste0("SS_", 1:16), batch = input$ss$batch[match(paste0("SS_", 1:16), input$ss$SampleID)], round(seedmat, 3))
cs$sd <- round(apply(seedmat, 1, sd), 3); cs$scored_alone <- round(m$probability_alone[m$strategy == "consensus30"][match(cs$SampleID, m$SampleID[m$strategy == "consensus30"])], 3)
# reference comparison with WW-corrected out-of-fold predictions
oof <- read.csv(file.path(PROD, "outer_predictions.csv"))
auc <- function(pos, neg) mean(outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b)))
ref <- do.call(rbind, lapply(c("consensus30", "diverse30"), function(st) {
  o1 <- oof[oof$strategy == st & oof$policy == "full_inner_recall", ]; o2 <- oof[oof$strategy == st & oof$policy == "inner_spec90", ]
  thr <- m[m$strategy == st, ][1, c("threshold_full_inner_recall", "threshold_inner_spec90")]
  ctrl <- o1[o1$truth == 0, ]; pos <- o1[o1$truth == 1, ]
  a <- m[m$strategy == st, ]
  row <- function(group, x, npf, nps) data.frame(strategy = st, group = group, n = length(x), median = round(median(x), 3), min = round(min(x), 3), max = round(max(x), 3),
                                                n_positive_full_recall = npf, n_positive_spec90 = nps, threshold_full = round(thr[[1]], 4), threshold_spec90 = round(thr[[2]], 4))
  rbind(row("OOF_controls_nestedCV", ctrl$probability, sum(ctrl$predicted_class), sum(o2$predicted_class[o2$truth == 0])),
        row("OOF_LoF-Mis_positives_nestedCV", pos$probability, sum(pos$predicted_class), sum(o2$predicted_class[o2$truth == 1])),
        row("LoF-LoF_frozen_model_all", a$probability_alone, sum(a$call_full_alone), sum(a$call_spec90_alone)),
        do.call(rbind, lapply(c("Run1", "Run2", "Run3"), function(b) { x <- a[a$batch == b, ]; row(paste0("LoF-LoF_", b), x$probability_alone, sum(x$call_full_alone), sum(x$call_spec90_alone)) })),
        data.frame(strategy = st, group = c("AUC_LoFLoF_vs_OOFcontrols", "AUC_OOFLoFMis_vs_LoFLoF", "AUC_OOFLoFMis_vs_OOFcontrols"), n = NA,
                   median = round(c(auc(a$probability_alone, ctrl$probability), auc(pos$probability, a$probability_alone), auc(pos$probability, ctrl$probability)), 3),
                   min = NA, max = NA, n_positive_full_recall = NA, n_positive_spec90 = NA, threshold_full = NA, threshold_spec90 = NA))
}))
# Table S19 base
sc <- alone[alone$strategy == "consensus30", ]
tab <- data.frame(Specimen = paste0("SS_", 1:16), Genotype = "LoF/LoF")
i <- match(tab$Specimen, sc$SampleID)
tab$Acquisition_run <- sc$batch[i]; tab$Run_present_in_training <- ifelse(sc$batch[i] %in% c("Run2", "Run3"), "yes", "no")
tab$SLC25A13_detected <- ifelse(sc$SLC25A13_raw_quantified[i], "yes", "no"); tab$Quantified_proteins <- sc$number_of_quantified_eligible_features[i]
for (st in c("consensus30", "diverse30")) { x <- m[m$strategy == st, ]; j <- match(tab$Specimen, x$SampleID)
  tab[[paste0(st, "_probability")]] <- round(x$probability_alone[j], 3)
  tab[[paste0(st, "_full_recall_call")]] <- ifelse(x$call_full_alone[j] == 1, "positive", "negative")
  tab[[paste0(st, "_spec90_call")]] <- ifelse(x$call_spec90_alone[j] == 1, "positive", "negative") }
for (d in c(OUT, RES)) {
  write.csv(m, file.path(d, "CHECK_single_specimen_rescoring_both_strategies.csv"), row.names = FALSE)
  write.csv(cs, file.path(d, "CHECK_seed_sensitivity_consensus30.csv"), row.names = FALSE)
  write.csv(ref, file.path(d, "REFERENCE_OOF_vs_LoFLoF_comparison.csv"), row.names = FALSE)
  write.csv(tab, file.path(d, "TableSx_LoFLoF_frozen_ensemble_scores.csv"), row.names = FALSE)
}
for (f in c("SS16_frozen_ensemble_scores.csv", "SS16_frozen_scores_summary.csv", "SS16_frozen_member_scores.csv", "SS16_scoring_methods_and_caveats.txt")) file.copy(file.path(RES, f), file.path(OUT, f), overwrite = TRUE)
writeLines(c(paste("Regenerated", Sys.time()), paste("model", model_file, unname(tools::md5sum(model_file))), paste("OOF", file.path(PROD, "outer_predictions.csv"))), file.path(OUT, "REGENERATED_reference_files_log.txt"))
print(ref, row.names = FALSE)
for (st in c("consensus30", "diverse30")) { x <- m[m$strategy == st, ]
  cat(sprintf("%s: max|alone-joint| %.4f; calls changed %d (%s); alone positives full %d spec90 %d; median alone %.3f\n", st, max(abs(x$diff)), sum(x$call_changed), paste(x$SampleID[x$call_changed], collapse = ","), sum(x$call_full_alone), sum(x$call_spec90_alone), median(x$probability_alone))) }
cat("seed SD max:", max(cs$sd), "\n"); print(cs[, c("SampleID", "batch", "sd", "scored_alone")], row.names = FALSE)
