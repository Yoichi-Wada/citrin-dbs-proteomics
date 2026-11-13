## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
suppressPackageStartupMessages({library(openxlsx)})
P <- file.path(REPO, "01_nested_ensemble", "results")
pr <- add_family(read.csv(file.path(P,"outer_predictions.csv"))); cat("outer_predictions columns:", paste(names(pr), collapse=", "), "\n")
perf <- read.csv(file.path(P,"outer_performance_all_policies.csv"))
ci <- read.csv(file.path(P, "ensemble_v1_1_performance_with_CI.csv")); ci <- ci[ci$comparison=="SM vs all controls",]  # canonical AUC CIs (used in all Response drafts)
auc <- function(s,y){ r <- rank(s); n1 <- sum(y==1); n0 <- sum(y==0); (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
cp <- function(x,n) binom.test(x,n)$conf.int
rows <- list(c("consensus30","full_inner_recall"), c("diverse30","full_inner_recall"), c("adaptive","full_inner_recall"), c("consensus30","inner_spec90"), c("diverse30","inner_spec90"), c("adaptive","inner_spec90"))
set.seed(20261105); B <- 2000
out <- do.call(rbind, lapply(rows, function(rw){
  st <- rw[1]; po <- rw[2]
  d <- pr[pr$strategy==st & pr$policy==po,]; if(!nrow(d)) d <- pr[pr$strategy==st,]
  stopifnot(nrow(d)==44, sum(d$truth==1)==7)
  pf <- perf[perf$strategy==st & perf$policy==po,]
  fams <- unique(d$family_id); ba <- function(x){ (mean(x$predicted_class[x$truth==1]==1)+mean(x$predicted_class[x$truth==0]==0))/2 }
  bs <- replicate(B, { f <- sample(fams, replace=TRUE); x <- do.call(rbind, lapply(f, function(k) d[d$family_id==k,])); if(sum(x$truth==1)>0 && sum(x$truth==0)>0) c(auc(x$probability,x$truth), ba(x)) else c(NA,NA) })
  qa <- quantile(bs[1,], c(.025,.975), na.rm=TRUE); qb <- quantile(bs[2,], c(.025,.975), na.rm=TRUE)
  TP <- pf$TP; FN <- pf$FN; TN <- pf$TN; FP <- pf$FP
  cs <- cp(TP,7); cq <- cp(TN,37)
  data.frame(Strategy=st, Threshold_rule=po, AUC=sprintf("%.3f (%s)", pf$AUC, ci$AUC_CI95_family_bootstrap[ci$strategy==st]), Family_AUC=sprintf("%.2f", pf$family_auc_descriptive),
    Sensitivity=sprintf("%d/7 = %.3f (%.3f–%.3f)", TP, TP/7, cs[1], cs[2]), Specificity=sprintf("%d/37 = %.3f (%.3f–%.3f)", TN, TN/37, cq[1], cq[2]),
    Balanced_accuracy=sprintf("%.3f (%.3f–%.3f)", (TP/7+TN/37)/2, qb[1], qb[2]), Positive_families_detected=sprintf("%d/5", pf$positive_families_detected),
    Brier=sprintf("%.3f", pf$Brier), boot_valid=sum(!is.na(bs[1,])), stringsAsFactors=FALSE)
}))
print(out[, 1:9], row.names=FALSE)
write.csv(out, file.path(P, "Table2_nested_performance_values.csv"), row.names=FALSE)
## formatted table
lab_st <- c(consensus30="Consensus-30 ensemble", diverse30="Diverse-30 ensemble", adaptive="Adaptive ensemble")
lab_po <- c(full_inner_recall="Complete inner sensitivity (primary rule)", inner_spec90="Inner specificity ≥ 0.90")
tab <- data.frame(`Strategy`=unname(lab_st[out$Strategy]), `Threshold rule`=unname(lab_po[out$Threshold_rule]), `AUC (95% CI)`=out$AUC, `Family-level AUC`=out$Family_AUC,
  `Sensitivity (95% CI)`=out$Sensitivity, `Specificity (95% CI)`=out$Specificity, `Balanced accuracy (95% CI)`=out$Balanced_accuracy, `LoF/Mis families detected`=out$Positive_families_detected, `Brier score`=out$Brier, check.names=FALSE)
tab$Strategy[1] <- paste0(tab$Strategy[1], " (primary)")
fn <- c("Table 2. Out-of-fold performance of the exploratory multiprotein ensemble in fully nested, family-grouped leave-one-family-out cross-validation.",
 "Analysis set: 44 specimens from 35 families (7 LoF/Mis individuals from 5 families; 37 controls from 30 families). All members of a family were held out together in the outer loop (35 folds); imputation, stability selection, consensus-panel formation, standardization, λ selection, member fitting and threshold derivation were repeated within each outer training set, and thresholds were derived in an inner family-leave-one-family-out loop (34 folds). The panel size was fixed at 30. No batch correction was applied.",
 "Strategies: Consensus-30, one common 30-protein candidate panel formed from ten members (random seeds 2025–2034) by ranking proteins on the number of members in whose top-30 lists they appeared (ties resolved by mean stability-selection frequency), each member an elastic-net model, probabilities averaged; Diverse-30, each member retains its own 30 proteins; Adaptive, the strategy with fewer inner false positives at complete inner sensitivity chosen per fold. The consensus strategy with the complete-inner-sensitivity threshold was designated primary before the final run; the other rows are alternatives specified for the revised ensemble evaluation and are reported irrespective of outcome. All strategies remain exploratory because their development was informed by the same cohort; the ensemble specification was fixed after exploratory single-seed analyses of these specimens.",
 "AUC and balanced accuracy: 95% confidence intervals by family-cluster bootstrap (families resampled with replacement from the out-of-fold predictions, 2,000 replicates; conditional on the fitted models). Sensitivity and specificity: exact binomial (Clopper–Pearson) 95% confidence intervals; counts at the fold-specific thresholds. These are descriptive specimen-level intervals and do not account for within-family correlation. Family-level AUC: one score per family (maximum member score), descriptive. A LoF/Mis family was counted as detected if at least one member was classified positive. Brier score: mean squared difference between out-of-fold probability and outcome. AUC and Brier score do not depend on the threshold rule and are identical between the two rows of each strategy.",
 "These estimates are internal; the ensemble was developed on this cohort and has not been validated in an independent cohort. Single-model seed-sensitivity analyses are given in Table S[x].")
wb <- createWorkbook(); addWorksheet(wb, "Table 2")
writeData(wb, 1, tab, startRow=2, headerStyle=createStyle(textDecoration="bold", border="TopBottom", wrapText=TRUE, valign="center"))
writeData(wb, 1, fn[1], startRow=1); addStyle(wb, 1, createStyle(textDecoration="bold"), rows=1, cols=1)
for(i in seq_along(fn)[-1]) writeData(wb, 1, fn[i], startRow=2+nrow(tab)+i)
setColWidths(wb, 1, cols=1:ncol(tab), widths=c(26,30,22,12,26,26,22,14,10))
addStyle(wb, 1, createStyle(border="bottom"), rows=2+nrow(tab), cols=1:ncol(tab), gridExpand=TRUE, stack=TRUE)
saveWorkbook(wb, file.path(P, "Table2_nested_performance.xlsx"), overwrite=TRUE)
md <- c("# Table 2 (main text)", "", fn[1], "", paste0("| ", paste(names(tab), collapse=" | "), " |"), paste0("|", paste(rep("---", ncol(tab)), collapse="|"), "|"), apply(tab, 1, function(r) paste0("| ", paste(r, collapse=" | "), " |")), "", fn[-1])
writeLines(md, file.path(P, "Table2_nested_performance.md")); cat("written\n")
