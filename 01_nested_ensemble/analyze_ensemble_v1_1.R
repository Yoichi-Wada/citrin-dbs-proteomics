## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(pROC))
out <- file.path(REPO, "01_nested_ensemble", "results")
pr <- read.csv(file.path(out, "outer_predictions.csv"), stringsAsFactors = FALSE)
pr <- pr[pr$policy == "full_inner_recall", ]
cat("\n=== Subgroup AUC (outer OOF probabilities, full_inner_recall rows)\n")
for (s in c("consensus30", "diverse30", "adaptive")) {
  d <- pr[pr$strategy == s, ]
  a_all <- as.numeric(auc(d$truth, d$probability, quiet = TRUE))
  res <- c(all = a_all)
  for (g in c("DD", "WS", "WW")) {
    dd <- d[d$group %in% c("SM", g), ]
    res[g] <- as.numeric(auc(dd$truth, dd$probability, quiet = TRUE))
  }
  cat(sprintf("%-12s AUC all=%.3f  vsDD=%.3f  vsWS=%.3f  vsWW=%.3f\n", s, res["all"], res["DD"], res["WS"], res["WW"]))
}
cat("\n=== diverse30 full_inner_recall: positives and false positives\n")
d <- pr[pr$strategy == "diverse30", ]
print(d[d$truth == 1, c("SampleID", "probability", "threshold", "predicted_class", "inner_false_positives_full")], row.names = FALSE)
print(d[d$truth == 0 & d$predicted_class == 1, c("SampleID", "group", "probability", "threshold")], row.names = FALSE)
cat("\n=== adaptive full_inner_recall: FP list\n")
d <- pr[pr$strategy == "adaptive", ]
print(d[d$truth == 0 & d$predicted_class == 1, c("SampleID", "group", "fitted_strategy", "probability", "threshold")], row.names = FALSE)
cat("\n=== Post hoc frontier for diverse30 (NOT validated): FP needed for 7/7 and 6/7\n")
d <- pr[pr$strategy == "diverse30", ]
pos <- sort(d$probability[d$truth == 1]); neg <- d$probability[d$truth == 0]
cat(sprintf("7/7 -> FP=%d ; 6/7 -> FP=%d ; 5/7 -> FP=%d\n", sum(neg >= pos[1]), sum(neg >= pos[2]), sum(neg >= pos[3])))
d <- pr[pr$strategy == "consensus30", ]
pos <- sort(d$probability[d$truth == 1]); neg <- d$probability[d$truth == 0]
cat(sprintf("consensus30: 7/7 -> FP=%d ; 6/7 -> FP=%d ; 5/7 -> FP=%d\n", sum(neg >= pos[1]), sum(neg >= pos[2]), sum(neg >= pos[3])))
cat("\n=== Control ranking (consensus30 / diverse30) top 8\n")
c30 <- pr[pr$strategy == "consensus30", c("SampleID", "group", "probability")]
d30 <- pr[pr$strategy == "diverse30", c("SampleID", "probability")]
m <- merge(c30, d30, by = "SampleID", suffixes = c("_c30", "_d30"))
m <- m[order(-m$probability_c30), ]
print(head(m[!grepl("^SM", m$SampleID), ], 8), row.names = FALSE)

cat("\n=== Consensus30 panel stability across 35 outer folds\n")
f <- read.csv(file.path(out, "outer_consensus30_features.csv"), stringsAsFactors = FALSE)
tab <- sort(table(f$FeatureID), decreasing = TRUE)
cat("distinct proteins ever in an outer consensus30 panel:", length(tab), "\n")
cat("in all 35 folds:", sum(tab == 35), " in >=30:", sum(tab >= 30), " in >=18:", sum(tab >= 18), "\n")
fin <- read.csv(file.path(out, "final_consensus30_proteins_NOT_INDEPENDENTLY_VALIDATED.csv"), stringsAsFactors = FALSE)
fin$outer_folds_in_panel <- as.integer(tab[fin$FeatureID]); fin$outer_folds_in_panel[is.na(fin$outer_folds_in_panel)] <- 0L
print(fin[, c("rank", "FeatureID", "members_top30", "outer_folds_in_panel", "in_original15_exploratory")], row.names = FALSE)
write.csv(fin, file.path(out, "final_consensus30_with_outer_recurrence.csv"), row.names = FALSE)
cat("original15 members in final consensus30:", sum(fin$in_original15_exploratory), "\n")

cat("\n=== Member audit: fallbacks / failures\n")
a <- read.csv(file.path(out, "outer_member_audit.csv"), stringsAsFactors = FALSE)
cat("members:", nrow(a), " stability_success min:", min(a$stability_success), " fallbacks total:", sum(a$stability_lambda_fallbacks), " failures total:", sum(a$stability_failures), "\n")
print(table(a$consensus_lambda_method)); print(table(a$diverse_lambda_method))

cat("\n=== Inner AUC per held-out family (consensus30)\n")
t <- read.csv(file.path(out, "outer_threshold_audit.csv"), stringsAsFactors = FALSE)
tc <- t[t$strategy == "consensus30", ]
cat(sprintf("median inner AUC %.3f (range %.3f-%.3f); median full-recall inner FP %d of ~36 (range %d-%d)\n",
  median(tc$inner_auc), min(tc$inner_auc), max(tc$inner_auc), median(tc$full_inner_FP), min(tc$full_inner_FP), max(tc$full_inner_FP)))
print(tc[tc$held_family %in% c("13","14","17","21","22"), c("held_family","inner_auc","full_threshold","full_inner_FP","full_inner_specificity")], row.names = FALSE)
cat("adaptive chose diverse30 in", sum(tc$chosen_adaptive == "diverse30"), "of 35 outer folds\n")
