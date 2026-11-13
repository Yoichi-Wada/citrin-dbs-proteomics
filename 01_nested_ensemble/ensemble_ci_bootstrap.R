## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
suppressPackageStartupMessages(library(pROC))
set.seed(20261103)
out <- file.path(REPO, "01_nested_ensemble", "results")
res_dir <- out
pr <- add_family(read.csv(file.path(out, "outer_predictions.csv"), stringsAsFactors = FALSE))
pr <- pr[pr$policy == "full_inner_recall", ]
B <- 2000
fam_boot_auc <- function(d, B) {
  fams <- split(seq_len(nrow(d)), d$family_id); nf <- length(fams); v <- numeric(0)
  for (b in seq_len(B)) {
    idx <- unlist(fams[sample.int(nf, nf, replace = TRUE)], use.names = FALSE)
    y <- d$truth[idx]; if (length(unique(y)) < 2) next
    v <- c(v, as.numeric(auc(y, d$probability[idx], quiet = TRUE)))
  }
  c(lower = unname(quantile(v, .025)), upper = unname(quantile(v, .975)), n_valid = length(v))
}
cp <- function(x, n) { ci <- binom.test(x, n)$conf.int; sprintf("%d/%d = %.3f (%.3f–%.3f)", x, n, x/n, ci[1], ci[2]) }
rows <- list()
for (s in c("consensus30", "diverse30", "adaptive")) {
  d <- pr[pr$strategy == s, ]
  for (g in c("all", "DD", "WS", "WW")) {
    dd <- if (g == "all") d else d[d$group %in% c("SM", g), ]
    a <- as.numeric(auc(dd$truth, dd$probability, quiet = TRUE)); ci <- fam_boot_auc(dd, B)
    tp <- sum(dd$truth == 1 & dd$predicted_class == 1); fn <- sum(dd$truth == 1 & dd$predicted_class == 0)
    tn <- sum(dd$truth == 0 & dd$predicted_class == 0); fp <- sum(dd$truth == 0 & dd$predicted_class == 1)
    pos_fam <- unique(dd$family_id[dd$truth == 1]); det_fam <- unique(dd$family_id[dd$truth == 1 & dd$predicted_class == 1])
    neg_fam <- unique(dd$family_id[dd$truth == 0]); fp_fam <- unique(dd$family_id[dd$truth == 0 & dd$predicted_class == 1])
    rows[[length(rows) + 1]] <- data.frame(strategy = s, comparison = if (g == "all") "SM vs all controls" else paste("SM vs", g),
      n = nrow(dd), n_families = length(unique(dd$family_id)),
      AUC = round(a, 3), AUC_CI95_family_bootstrap = sprintf("%.3f–%.3f", ci["lower"], ci["upper"]), boot_valid = ci["n_valid"],
      sensitivity = cp(tp, tp + fn), specificity = cp(tn, tn + fp),
      positive_families_detected = sprintf("%d/%d", length(det_fam), length(pos_fam)),
      negative_families_with_FP = sprintf("%d/%d", length(fp_fam), length(neg_fam)),
      Brier = round(mean((dd$probability - dd$truth)^2), 3), stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
print(tab, row.names = FALSE, right = FALSE)
write.csv(tab, file.path(res_dir, "ensemble_v1_1_performance_with_CI.csv"), row.names = FALSE)
cat("\nNotes: conditional CIs computed on fixed outer OOF probabilities and calls (pipeline not re-run); AUC CI = percentile family-cluster bootstrap, B=2000, seed 20261103; sens/spec CI = Clopper-Pearson exact.\n")
