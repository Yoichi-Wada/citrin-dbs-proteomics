## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Constrained family-level label permutation plan for the select pipeline null.
# Constraints (fixed before any permutation run):
#  - positives only in Run3 (as observed); Run2 specimens are never positive
#  - exactly 5 positive families with sizes {1,1,1,2,2} = 7 positive specimens (as observed)
#  - family structure, batch, SS exclusion unchanged; only the SM label moves
#  - observed configuration excluded; 200 distinct draws; generator seed 20261103
set.seed(20261103)
dd <- DATA
out <- file.path(Sys.getenv("CITRIN_PERMUTATION_RUNS", file.path(REPO, "03_sensitivity", "permutation", "runs")), "plan")   # per-permutation input folders (not distributed)
si <- read.csv(file.path(dd, "sample_info.csv"), stringsAsFactors = FALSE)
fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE)
m <- merge(si, fm[, c("SampleID", "family_id")], by = "SampleID")
m44 <- m[m$group != "SS", ]
stopifnot(nrow(m44) == 44, sum(m44$group == "SM") == 7)
r3 <- m44[m44$batch == "Run3", ]
fs <- table(r3$family_id)
single <- names(fs)[fs == 1]; pair <- names(fs)[fs == 2]
cat("Run3 families:", length(fs), "| size1:", length(single), "| size2:", length(pair), "| size3:", sum(fs == 3), "\n")
obs <- sort(unique(m44$family_id[m44$group == "SM"]))
key <- function(f) paste(sort(f), collapse = "_")
total <- choose(length(single), 3) * choose(pair |> length(), 2)
cat("configurations:", total, " (observed excluded -> ", total - 1, ")\n")
seen <- key(obs); plans <- list(); i <- 0
while (length(plans) < 200) {
  f <- c(sample(single, 3), sample(pair, 2)); k <- key(f)
  if (k %in% seen) next
  seen <- c(seen, k); i <- i + 1
  pos_ids <- m44$SampleID[m44$family_id %in% as.integer(f)]
  stopifnot(length(pos_ids) == 7)
  si_p <- si
  si_p$group[si_p$group == "SM"] <- "DD"          # placeholder negative label for the true SM
  si_p$group[si_p$SampleID %in% pos_ids] <- "SM"
  stopifnot(sum(si_p$group == "SM") == 7, all(si_p$batch[si_p$group == "SM"] == "Run3"), sum(si_p$group == "SS") == 16)
  pd <- file.path(out, sprintf("perm_%03d", i)); dir.create(pd, showWarnings = FALSE)
  write.csv(si_p, file.path(pd, "sample_info.csv"), row.names = FALSE, quote = FALSE)
  if (!file.exists(file.path(pd, "raw_data.xlsx"))) file.symlink(file.path(dd, "raw_data.xlsx"), file.path(pd, "raw_data.xlsx"))
  plans[[i]] <- data.frame(perm_id = i, positive_families = k, positive_specimens = paste(sort(pos_ids), collapse = ";"),
                           n_true_SM_relabelled_positive = sum(pos_ids %in% m44$SampleID[m44$group == "SM"]),
                           groups_of_pseudo_positives = paste(sort(table(m44$group[m44$SampleID %in% pos_ids])), collapse = ";"))
}
plan <- do.call(rbind, plans)
plan$groups_of_pseudo_positives <- sapply(plan$positive_specimens, function(s) { ids <- strsplit(s, ";")[[1]]; g <- m44$group[match(ids, m44$SampleID)]; paste(names(table(g)), table(g), sep = "=", collapse = ";") })
write.csv(plan, file.path(out, "permutation_plan.csv"), row.names = FALSE)
write.csv(plan, file.path(file.path(REPO, "03_sensitivity", "permutation", "results"), "permutation_plan.csv"), row.names = FALSE)
cat("distinct plans:", length(unique(plan$positive_families)), "\n")
print(table(plan$n_true_SM_relabelled_positive))
print(head(plan, 5))
cat("observed key:", key(obs), "\n")
