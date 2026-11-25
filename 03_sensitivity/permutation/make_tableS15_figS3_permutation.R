#!/usr/bin/env Rscript
## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Table S13 and Fig. S4 for the label-permutation null of the primary consensus-30 ensemble AUC (100 permutations, completed 2026-11-08 14:09).
suppressPackageStartupMessages({library(ggplot2); library(openxlsx)})
W <- file.path(REPO, "03_sensitivity", "permutation", "results"); O <- W
null <- read.csv(file.path(W, "permutation_null_AUC.csv")); summ <- read.csv(file.path(W, "permutation_summary.csv"), row.names = 1)
plan <- read.csv(file.path(W, "permutation_plan.csv"))
null$perm_id <- as.integer(sub("perm_", "", null$perm))
null <- merge(null, plan[, c("perm_id", "n_true_SM_relabelled_positive")], by = "perm_id", all.x = TRUE)
stopifnot(nrow(null) == 100)
obs <- c(consensus30 = summ["consensus30", "observed.consensus30"], diverse30 = summ["diverse30", "observed.consensus30"])
q <- function(x) c(median = median(x), q025 = unname(quantile(x, .025)), q25 = unname(quantile(x, .25)), q75 = unname(quantile(x, .75)), q975 = unname(quantile(x, .975)), max = max(x))
sc <- q(null$consensus30_AUC); sd_ <- q(null$diverse30_AUC); sf <- q(null$consensus30_family_AUC)
b_c <- sum(null$consensus30_AUC >= obs["consensus30"]); b_d <- sum(null$diverse30_AUC >= obs["diverse30"]); b_f <- sum(null$consensus30_family_AUC >= 0.773)
# ---- Table S13 workbook ----
wb <- createWorkbook(); addWorksheet(wb, "Table S15")
title <- "Table S15. Label-permutation reference distribution for the out-of-fold AUC of the primary consensus-30 ensemble (100 family-structured label permutations)."
writeData(wb, "Table S15", title, startRow = 1, startCol = 1)
hdr <- data.frame(Permutation = sprintf("perm_%03d", null$perm_id),
  `Consensus-30 out-of-fold AUC (primary statistic)` = round(null$consensus30_AUC, 4),
  `Diverse-30 out-of-fold AUC (alternative)` = round(null$diverse30_AUC, 4),
  `Consensus-30 family-level AUC (descriptive)` = round(null$consensus30_family_AUC, 4),
  `Permuted positive specimens coinciding with true LoF/Mis specimens (n of 7)` = null$n_true_SM_relabelled_positive, check.names = FALSE)
writeData(wb, "Table S15", hdr, startRow = 3, startCol = 1)
r <- 3 + nrow(hdr) + 2
notes <- c(
  sprintf("Observed values (original labels): consensus-30 AUC %.3f; diverse-30 AUC %.3f; consensus-30 family-level AUC 0.773.", obs["consensus30"], obs["diverse30"]),
  sprintf("Consensus-30 (primary, pre-specified): null median %.3f, interquartile range %.3f–%.3f, 2.5th–97.5th percentiles %.3f–%.3f, maximum %.3f; %d of 100 permutations ≥ observed; one-sided P = (%d + 1)/(100 + 1) = %.3f (Monte Carlo resolution 0.01).", sc["median"], sc["q25"], sc["q75"], sc["q025"], sc["q975"], sc["max"], b_c, b_c, (b_c + 1) / 101),
  sprintf("Diverse-30 (alternative; not the statistic specified in advance): null median %.3f, 2.5th–97.5th percentiles %.3f–%.3f, maximum %.3f; %d of 100 ≥ observed; P = %.3f.", sd_["median"], sd_["q025"], sd_["q975"], sd_["max"], b_d, (b_d + 1) / 101),
  sprintf("Consensus-30 family-level AUC (descriptive, not pre-specified): null median %.3f, 2.5th–97.5th percentiles %.3f–%.3f; %d of 100 ≥ observed 0.773; P = %.3f.", sf["median"], sf["q025"], sf["q975"], b_f, (b_f + 1) / 101),
  "Design (fixed on 2026-11-07 before any permutation result was inspected): in each permutation the seven positive labels were assigned to five randomly chosen Run 3 families of sizes 1, 1, 1, 2 and 2 (measurements, family map and run composition unchanged), and the complete outer family-grouped procedure (training-partition imputation, stability selection with ten members, consensus-30 panel, standardization, λ tuning, elastic-net fitting, averaging of member probabilities) was repeated in every outer fold. The inner leave-one-family-out calibration, used only for decision thresholds and the adaptive choice, was omitted; on the original labels this omission reproduced the 44 out-of-fold probabilities of the primary analysis exactly (maximum absolute difference 0). The number of permutations (100) was fixed in advance; P = (b + 1)/(N + 1).",
  "The last column gives, for each permutation, how many of the seven permuted positive specimens coincide with true LoF/Mis specimens (0–7); permutations with high overlap are expected to yield higher AUC and are retained because they are part of the pre-specified family-structured null.",
  "Source: 03_sensitivity/permutation (citrin_nested_ensemble_familyLOFO_v1_1_permAUC.R, engine v2_2 md5 fbff0fa82f0b2da1dafb553bc5564d27; collect_perm_auc.R).")
for (i in seq_along(notes)) writeData(wb, "Table S15", notes[i], startRow = r + i - 1, startCol = 1)
setColWidths(wb, "Table S15", cols = 1:5, widths = c(14, 22, 22, 22, 26))
saveWorkbook(wb, file.path(O, "TableS15_permutation_ensemble_AUC.xlsx"), overwrite = TRUE)
write.csv(hdr, file.path(O, "TableS15_values.csv"), row.names = FALSE)
# ---- Fig. S4 ----
d <- rbind(data.frame(strategy = "Consensus-30 (primary)", AUC = null$consensus30_AUC, obs = obs["consensus30"], b = b_c),
           data.frame(strategy = "Diverse-30 (alternative)", AUC = null$diverse30_AUC, obs = obs["diverse30"], b = b_d))
lab <- data.frame(strategy = c("Consensus-30 (primary)", "Diverse-30 (alternative)"),
  text = c(sprintf("observed AUC %.3f\n%d/100 permutations ≥ observed\nP = %.3f", obs["consensus30"], b_c, (b_c + 1) / 101),
           sprintf("observed AUC %.3f\n%d/100 permutations ≥ observed\nP = %.3f", obs["diverse30"], b_d, (b_d + 1) / 101)))
g <- ggplot(d, aes(x = AUC)) +
  geom_histogram(breaks = seq(0, 1.0, 0.05), fill = "#B6C7D6", colour = "white", linewidth = .3) +
  geom_vline(aes(xintercept = obs), colour = "#BD6059", linewidth = .8) +
  geom_text(data = lab, aes(x = 0.02, y = 19, label = text), hjust = 0, vjust = 1, size = 3, lineheight = .95, colour = "#17212B") +
  facet_wrap(~ strategy, ncol = 1) +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 20), expand = expansion(0)) +
  labs(x = "Out-of-fold AUC under 100 family-structured label permutations", y = "Number of permutations") +
  theme_classic(base_size = 10) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0))
ggsave(file.path(O, "FigureS3_permutation_null_ensemble_AUC.pdf"), g, width = 6.5, height = 6, device = cairo_pdf)
ggsave(file.path(O, "FigureS3_permutation_null_ensemble_AUC.png"), g, width = 6.5, height = 6, dpi = 200, bg = "white")
cat(sprintf("consensus30: med %.3f IQR %.3f-%.3f 95%% %.3f-%.3f max %.3f b=%d P=%.3f | diverse30: med %.3f 95%% %.3f-%.3f max %.3f b=%d P=%.3f | family: med %.3f 95%% %.3f-%.3f b=%d P=%.3f\n",
  sc["median"], sc["q25"], sc["q75"], sc["q025"], sc["q975"], sc["max"], b_c, (b_c+1)/101, sd_["median"], sd_["q025"], sd_["q975"], sd_["max"], b_d, (b_d+1)/101, sf["median"], sf["q025"], sf["q975"], b_f, (b_f+1)/101))
print(table(null$n_true_SM_relabelled_positive))
