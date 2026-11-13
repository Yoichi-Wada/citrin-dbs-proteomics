## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(openxlsx))
setwd(file.path(REPO, "02_final_model", "results"))
d <- read.csv("TableSx_LoFLoF_frozen_ensemble_scores.csv", check.names = FALSE)
ck <- read.csv("CHECK_single_specimen_rescoring_both_strategies.csv")
for (st in c("consensus30", "diverse30")) {
  x <- ck[ck$strategy == st, ]; i <- match(d$Specimen, x$SampleID)
  d[[paste0(st, "_probability")]] <- x$probability_alone[i]
  d[[paste0(st, "_full_recall_call")]] <- ifelse(x$call_full_alone[i] == 1, "positive", "negative")
  d[[paste0(st, "_spec90_call")]] <- ifelse(x$call_spec90_alone[i] == 1, "positive", "negative")
}
ref <- read.csv("REFERENCE_OOF_vs_LoFLoF_comparison.csv")
sm <- read.csv("SS16_frozen_scores_summary.csv")
d$Specimen <- factor(d$Specimen, levels = paste0("SS_", 1:16)); d <- d[order(d$Specimen), ]

tab <- data.frame(
  Specimen = as.character(d$Specimen),
  Genotype = d$Genotype,
  `Acquisition run` = d$Acquisition_run,
  `Run contained training specimens` = ifelse(d$Run_present_in_training == "yes", "Yes", "No (extrapolation)"),
  `SLC25A13 detected` = ifelse(d$SLC25A13_detected == "yes", "Yes", "No"),
  `Quantified proteins (of 5,199)` = d$Quantified_proteins,
  `Consensus-30 probability` = sprintf("%.3f", d$consensus30_probability),
  `Consensus-30 call, full inner-OOF recall threshold (0.034)` = d$consensus30_full_recall_call,
  `Consensus-30 call, inner-OOF specificity threshold (0.210)` = d$consensus30_spec90_call,
  `Diverse-30 probability` = sprintf("%.3f", d$diverse30_probability),
  `Diverse-30 call, full inner-OOF recall threshold (0.090)` = d$diverse30_full_recall_call,
  `Diverse-30 call, inner-OOF specificity threshold (0.517)` = d$diverse30_spec90_call,
  check.names = FALSE)

grp_label <- function(g) switch(g,
  OOF_controls_nestedCV = "Reference: 37 controls, out-of-fold (nested CV)",
  `OOF_LoF-Mis_positives_nestedCV` = "Reference: 7 LoF/Mis, out-of-fold (nested CV)",
  `LoF-LoF_frozen_model_all` = "LoF/LoF, all runs (n = 16), frozen model", g)
summ <- do.call(rbind, lapply(c("consensus30", "diverse30"), function(st) {
  r <- ref[ref$strategy == st & ref$group %in% c("OOF_controls_nestedCV", "OOF_LoF-Mis_positives_nestedCV"), ]
  s <- sm[sm$strategy == st & sm$batch != "ALL", ]
  rbind(
    data.frame(Strategy = st, Group = sapply(r$group, grp_label),
               `Median probability (range)` = sprintf("%.3f (%.3f–%.3f)", r$median, r$min, r$max),
               `Positive at full-recall threshold` = "—",
               `Positive at specificity threshold` = "—", check.names = FALSE),
    do.call(rbind, lapply(c("ALL", "Run1", "Run2", "Run3"), function(b) {
      x <- ck[ck$strategy == st & (b == "ALL" | ck$batch == b), ]
      data.frame(Strategy = st, Group = if (b == "ALL") "LoF/LoF, all runs (n = 16), frozen model" else sprintf("LoF/LoF, %s (n = %d), frozen model", b, nrow(x)),
                 `Median probability (range)` = sprintf("%.3f (%.3f–%.3f)", median(x$probability_alone), min(x$probability_alone), max(x$probability_alone)),
                 `Positive at full-recall threshold` = sprintf("%d/%d", sum(x$call_full_alone), nrow(x)),
                 `Positive at specificity threshold` = sprintf("%d/%d", sum(x$call_spec90_alone), nrow(x)), check.names = FALSE) })))
}))
rownames(summ) <- NULL
summ$Strategy <- ifelse(summ$Strategy == "consensus30", "Consensus-30", "Diverse-30")

title <- "Table S19. Application of the frozen final classifier to the 16 LoF/LoF specimens."
legend <- c(
  "The final model (consensus 30-protein panel with ten elastic-net members; secondary diverse-panel strategy) was trained on the 44 non-LoF/LoF specimens and applied without modification. Columns give the acquisition run, whether that run contained training specimens (Run 1 did not; scores for Run 1 are extrapolations), SLC25A13 detection (undetectable in all 16), the number of quantified proteins among the 5,199 eligible proteins, and for each strategy the ensemble-mean predicted probability and the calls at the two thresholds derived from the inner family-leave-one-family-out validation of the 44-specimen training cohort (full inner-OOF recall: 0.034 consensus / 0.090 diverse; inner-OOF specificity ≥ 0.90: 0.210 / 0.517).",
  "For reference, the summary block gives the out-of-fold probabilities of the 37 controls and 7 LoF/Mis specimens in the nested cross-validation (medians 0.030 and 0.131 for consensus; 0.048 and 0.152 for diverse). The thresholds were estimated on the training cohort by inner validation, do not describe the final model's sensitivity or specificity on that cohort, and differ from the fold-specific thresholds used in the cross-validation; threshold counts are therefore not given for the reference groups. Missing values were imputed from each specimen's own 1% quantile with the training-derived scale and a fixed seed; probabilities were insensitive to the imputation seed (within-specimen SD ≤ 0.020 over five seeds). The reference scores were generated by fold-specific models in the nested cross-validation, whereas the LoF/LoF scores were generated by the final model trained on all 44 development specimens; no AUC, test or performance estimate is derived from comparisons between them. Each specimen was scored alone, so that its imputation draws do not depend on other specimens or on input order; scoring the 16 specimens together changed the probabilities by at most 0.031 and altered the full-recall call of one specimen under the diverse strategy (SS_7) and none under the consensus strategy.")

wb <- createWorkbook(); addWorksheet(wb, "Table S LoF-LoF")
hs <- createStyle(textDecoration = "bold", border = "TopBottom", wrapText = TRUE, valign = "center")
writeData(wb, 1, title, startRow = 1); addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
writeData(wb, 1, tab, startRow = 3, headerStyle = hs)
r2 <- 3 + nrow(tab) + 2
writeData(wb, 1, "Summary", startRow = r2); addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = r2, cols = 1)
writeData(wb, 1, summ, startRow = r2 + 1, headerStyle = hs)
r3 <- r2 + 1 + nrow(summ) + 2
for (i in seq_along(legend)) writeData(wb, 1, legend[i], startRow = r3 + i - 1)
setColWidths(wb, 1, cols = 1:ncol(tab), widths = c(10, 9, 11, 18, 10, 12, 12, 18, 18, 12, 18, 18))
addStyle(wb, 1, createStyle(halign = "center"), rows = 4:(3 + nrow(tab)), cols = 3:ncol(tab), gridExpand = TRUE, stack = TRUE)
saveWorkbook(wb, "TableS_LoFLoF_frozen_model_scores.xlsx", overwrite = TRUE)

mdrow <- function(x) paste0("| ", paste(x, collapse = " | "), " |")
md <- c(paste0("# ", title), "",
        mdrow(names(tab)), paste0("|", paste(rep("---", ncol(tab)), collapse = "|"), "|"), apply(tab, 1, mdrow), "",
        "**Summary**", "", mdrow(names(summ)), paste0("|", paste(rep("---", ncol(summ)), collapse = "|"), "|"), apply(summ, 1, mdrow), "", legend)
writeLines(md, "TableS_LoFLoF_frozen_model_scores.md")
cat("written\n"); print(summ, row.names = FALSE)
