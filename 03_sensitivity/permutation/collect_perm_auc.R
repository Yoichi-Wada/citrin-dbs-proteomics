## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Collect the consensus30 (primary) and diverse30 (auxiliary) out-of-fold AUCs of the 100 label permutations.
W <- Sys.getenv("CITRIN_PERMUTATION_RUNS", file.path(REPO, "03_sensitivity", "permutation", "runs"))   # folder holding runs/perm_001 ... perm_100 (regenerable; not distributed)
RES <- file.path(REPO, "03_sensitivity", "permutation", "results")
obs <- read.csv(file.path(file.path(REPO, "01_nested_ensemble", "results"), "outer_performance_all_policies.csv"))
obs_auc <- setNames(obs$AUC[obs$policy == "full_inner_recall"], obs$strategy[obs$policy == "full_inner_recall"])
runs <- sort(list.files(W, pattern = "^perm_[0-9]{3}$"))
rows <- lapply(runs, function(p) {
  f <- file.path(W, p, "outer_performance_all_policies.csv"); if (!file.exists(f)) return(NULL)
  d <- read.csv(f); d <- d[d$policy == "full_inner_recall", ]
  data.frame(perm = p, consensus30_AUC = d$AUC[d$strategy == "consensus30"], diverse30_AUC = d$AUC[d$strategy == "diverse30"],
             consensus30_family_AUC = d$family_auc_descriptive[d$strategy == "consensus30"])
})
null <- do.call(rbind, rows); write.csv(null, file.path(RES, "permutation_null_AUC.csv"), row.names = FALSE)
summ <- function(x, o) c(n = length(x), observed = o, null_median = median(x), null_q025 = quantile(x, .025), null_q975 = quantile(x, .975),
                         null_max = max(x), n_ge_observed = sum(x >= o), P_plus_one = (sum(x >= o) + 1) / (length(x) + 1))
out <- rbind(consensus30 = summ(null$consensus30_AUC, obs_auc["consensus30"]), diverse30 = summ(null$diverse30_AUC, obs_auc["diverse30"]))
write.csv(out, file.path(RES, "permutation_summary.csv")); print(round(out, 4))
