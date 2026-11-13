## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Supplementary Table: protein-exclusion sensitivity analysis (ablation) of the nested consensus-30 ensemble.
# Scenarios: baseline (production), all4, no_PRL, no_MTTP, no_ARHGEF10, no_SEC31A. Complete nested
# re-development per scenario (exclusion before any selection step). Primary = consensus30, complete-inner-
# sensitivity threshold. CIs: family-cluster bootstrap (AUC, BA; 2,000 replicates), exact binomial (sens, spec).
suppressPackageStartupMessages({library(openxlsx)})
setwd(REPO)
R <- REPO
A <- file.path(REPO, "03_sensitivity", "ablation", "results")
scen <- c(baseline = file.path(R, "01_nested_ensemble", "results"),
          no_PRL = file.path(A, "no_PRL"), no_MTTP = file.path(A, "no_MTTP"),
          no_ARHGEF10 = file.path(A, "no_ARHGEF10"), no_SEC31A = file.path(A, "no_SEC31A"),
          all4 = file.path(A, "all4"))
label <- c(baseline = "None (primary analysis)", no_PRL = "PRL", no_MTTP = "MTTP", no_ARHGEF10 = "ARHGEF10", no_SEC31A = "SEC31A",
           all4 = "PRL, MTTP, ARHGEF10 and SEC31A")
auc <- function(s, y) { r <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cp <- function(x, n) binom.test(x, n)$conf.int
rows <- list(); fps <- list(); panels <- list(); preds <- list()
for (k in names(scen)) {
  d <- add_family(read.csv(file.path(scen[k], "outer_primary_predictions.csv"))); stopifnot(nrow(d) == 44, sum(d$truth) == 7)
  pf <- read.csv(file.path(scen[k], "outer_performance_all_policies.csv")); pf <- pf[pf$strategy == "consensus30" & pf$policy == "full_inner_recall", ]
  set.seed(20261106); fams <- unique(d$family_id)
  bs <- replicate(2000, { f <- sample(fams, replace = TRUE); x <- do.call(rbind, lapply(f, function(j) d[d$family_id == j, ]))
    if (sum(x$truth) > 0 && sum(x$truth == 0) > 0) c(auc(x$probability, x$truth), (mean(x$predicted_class[x$truth == 1] == 1) + mean(x$predicted_class[x$truth == 0] == 0)) / 2) else c(NA, NA) })
  qa <- quantile(bs[1, ], c(.025, .975), na.rm = TRUE); qb <- quantile(bs[2, ], c(.025, .975), na.rm = TRUE)
  if (k == "baseline") { qa <- c(0.602, 0.971); qb <- c(0.848, 0.974) }  # WWfix production CIs (Table 2, 2026-11-07; same bootstrap method, seed of the production summary)
  cs <- cp(pf$TP, 7); cq <- cp(pf$TN, 37)
  feat <- read.csv(file.path(scen[k], "outer_consensus30_features.csv")); nd <- length(unique(feat$FeatureID))
  rows[[k]] <- data.frame(`Excluded protein(s)` = label[k], `AUC (95% CI)` = sprintf("%.3f (%.3f–%.3f)", pf$AUC, qa[1], qa[2]),
    `Sensitivity (95% CI)` = sprintf("%d/7 (%.2f–%.2f)", pf$TP, cs[1], cs[2]), `Specificity (95% CI)` = sprintf("%d/37 (%.2f–%.2f)", pf$TN, cq[1], cq[2]),
    `Balanced accuracy (95% CI)` = sprintf("%.3f (%.3f–%.3f)", (pf$TP / 7 + pf$TN / 37) / 2, qb[1], qb[2]), `LoF/Mis families detected` = sprintf("%d/5", pf$positive_families_detected),
    `Brier score` = sprintf("%.3f", pf$Brier), `Distinct proteins in the 35 consensus panels` = nd, check.names = FALSE)
  fps[[k]] <- d$SampleID[d$truth == 0 & d$predicted_class == 1]; preds[[k]] <- d
  panels[[k]] <- table(feat$FeatureID)
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
print(tab, row.names = FALSE)
# per-specimen: which controls are FP in each scenario; which LoF/Mis missed
ids <- preds$baseline$SampleID[order(preds$baseline$group, preds$baseline$SampleID)]
call_mat <- sapply(names(scen), function(k) { d <- preds[[k]]; setNames(d$predicted_class, d$SampleID)[ids] })
prob_mat <- sapply(names(scen), function(k) { d <- preds[[k]]; round(setNames(d$probability, d$SampleID)[ids], 3) })
per <- data.frame(SampleID = ids, group = setNames(preds$baseline$group, preds$baseline$SampleID)[ids], truth = setNames(preds$baseline$truth, preds$baseline$SampleID)[ids],
                  call_mat, check.names = FALSE); colnames(per)[4:9] <- paste0("call_", names(scen))
per2 <- data.frame(per, prob_mat, check.names = FALSE); colnames(per2)[10:15] <- paste0("prob_", names(scen))
fp_summary <- data.frame(Scenario = label[names(scen)], `False positives (controls called positive)` = sapply(names(scen), function(k) paste(sort(fps[[k]]), collapse = ", ")),
                         `Missed LoF/Mis` = sapply(names(scen), function(k) { d <- preds[[k]]; m <- d$SampleID[d$truth == 1 & d$predicted_class == 0]; if (length(m)) paste(m, collapse = ", ") else "none" }), check.names = FALSE)
rownames(fp_summary) <- NULL; print(fp_summary, row.names = FALSE)
# core proteins presence across scenarios (in all 35 panels)
core <- c(SEC31A = "O94979", ARHGEF10 = "O15013", NAPB = "Q9H115", TESC = "Q96BS2", RAB2B = "Q8WUD1", BLTP3B = "A0JNW5", DOCK11 = "Q5JSL3", PRL = "P01236", MTTP = "P55157", CPNE1 = "Q99829", PCYOX1 = "Q9UHG3", PPAT = "Q06203")
core_tab <- data.frame(Protein = names(core), sapply(names(scen), function(k) { x <- panels[[k]]; v <- as.integer(x[core]); v[is.na(v)] <- 0; v }), check.names = FALSE)
print(core_tab, row.names = FALSE)

title <- "Table S15. Protein-exclusion sensitivity analysis of the exploratory consensus-30 ensemble."
legend <- c(
  "For each scenario the complete nested model-development procedure (training-only imputation, stability selection, consensus-panel formation, member fitting and inner-loop threshold derivation; family-grouped leave-one-family-out, 35 outer folds, ten members) was re-run after removing the indicated protein(s) from the 5,199-protein candidate matrix before any selection step. Values are family-held-out out-of-fold results of the primary strategy (consensus-30, complete-inner-sensitivity threshold) on the same 44 specimens (7 LoF/Mis from 5 families; 37 controls). AUC and balanced accuracy: family-cluster bootstrap 95% CI (2,000 replicates, conditional on the fitted models; the baseline intervals are those reported in Table 2, obtained with the same method); sensitivity and specificity: exact binomial 95% CI (specimen-level, not accounting for within-family correlation). The last column gives the number of distinct proteins that entered the consensus candidate panel in at least one of the 35 outer folds.",
  "The scenarios (joint exclusion of the four proteins highlighted by the reviewer, and each protein individually) were fixed before the runs and are all reported. These are post hoc sensitivity analyses within the development cohort: they indicate whether a protein is individually indispensable for obtaining discrimination in this cohort, but do not establish its contribution in the primary model, do not show that the remaining signal is free of confounding, and are not independent validation. AUC, Brier score and threshold-based sensitivity and specificity are distinct aspects of performance: because the inner-loop threshold is re-derived in every scenario, threshold-based counts can change substantially with little change in AUC. Selection counts of core proteins that remain at 35/35 do not imply that the rest of the panel or the fitted coefficients were unchanged, since a removed protein's place in each fold's panel is taken by another protein and the stochastic imputation differs once a feature is removed. Sheet 'Per-specimen calls' gives each specimen's out-of-fold call and probability under every scenario; sheet 'Core proteins' gives, for selected proteins, the number of outer folds (of 35) in which they entered the consensus panel under each scenario.")
dir.create(A, showWarnings = FALSE)
wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold", border = "TopBottom", wrapText = TRUE, valign = "center")
addWorksheet(wb, "Table S ablation"); writeData(wb, 1, title, startRow = 1); addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
writeData(wb, 1, tab, startRow = 3, headerStyle = hs); r <- 3 + nrow(tab) + 2
writeData(wb, 1, "Misclassified specimens by scenario", startRow = r); addStyle(wb, 1, createStyle(textDecoration = "bold"), rows = r, cols = 1)
writeData(wb, 1, fp_summary, startRow = r + 1, headerStyle = hs); r2 <- r + 1 + nrow(fp_summary) + 2
for (i in seq_along(legend)) writeData(wb, 1, legend[i], startRow = r2 + i - 1)
setColWidths(wb, 1, cols = 1:8, widths = c(30, 22, 20, 20, 22, 14, 10, 18))
addWorksheet(wb, "Per-specimen calls"); writeData(wb, 2, per2, headerStyle = hs); freezePane(wb, 2, firstRow = TRUE, firstCol = TRUE)
addWorksheet(wb, "Core proteins"); writeData(wb, 3, core_tab, headerStyle = hs)
saveWorkbook(wb, file.path(A, "TableS_ablation_consensus30.xlsx"), overwrite = TRUE)
write.csv(tab, file.path(A, "TableS_ablation_summary.csv"), row.names = FALSE); write.csv(per2, file.path(A, "ablation_per_specimen_calls.csv"), row.names = FALSE); write.csv(core_tab, file.path(A, "ablation_core_protein_panel_counts.csv"), row.names = FALSE)
mdrow <- function(x) paste0("| ", paste(x, collapse = " | "), " |")
writeLines(c(paste0("# ", title), "", mdrow(names(tab)), paste0("|", paste(rep("---", ncol(tab)), collapse = "|"), "|"), apply(tab, 1, mdrow), "",
             "**Misclassified specimens by scenario**", "", mdrow(names(fp_summary)), paste0("|", paste(rep("---", 3), collapse = "|"), "|"), apply(fp_summary, 1, mdrow), "",
             "**Consensus-panel membership of selected proteins (outer folds of 35)**", "", mdrow(names(core_tab)), paste0("|", paste(rep("---", ncol(core_tab)), collapse = "|"), "|"), apply(core_tab, 1, mdrow), "", legend),
           file.path(A, "TableS_ablation_consensus30.md"))
cat("written\n")
