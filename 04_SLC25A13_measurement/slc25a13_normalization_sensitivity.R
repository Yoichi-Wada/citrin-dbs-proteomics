## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
# ============================================================================
# SLC25A13 conclusions under alternative normalization strategies (R2 Major 6, R4 Major 4)
# and stability of the ACTB normaliser across groups and batches.
# Inputs (R/data): raw_data.xlsx (raw DIA-NN protein intensities, 0 = not quantified),
#   DBS-Tohoku.pg_matrix.tsv (DIA-NN MaxLFQ, DIA-NN RT-dependent normalisation),
#   sample_info.csv, Batch_20250809.xlsx (not needed: batch in sample_info).
# Normalisations compared (log2 scale): ACTB (manuscript), sample median, total intensity,
#   DIA-NN MaxLFQ, none (raw). Median/total factors use the complete-case protein set
#   (quantified in all 60 specimens) so that missingness does not bias the factor.
# Outputs: this folder (results CSVs, figure PDF).
# ============================================================================
suppressPackageStartupMessages({library(openxlsx); library(ggplot2)})
set.seed(20261104); B <- 5000
dd <- DATA
out <- file.path(REPO, "04_SLC25A13_measurement", "results")
si <- read.csv(file.path(dd, "sample_info.csv"), stringsAsFactors = FALSE)
si$group <- factor(si$group, levels = c("DD","WW","WS","SM","SS"))
rw <- read.xlsx(file.path(dd, "raw_data.xlsx"), sheet = "data", check.names = FALSE, sep.names = " ")
stopifnot(all(si$SampleID %in% names(rw)))
acc <- rw[["Master accession number"]]; upep <- as.numeric(rw[["Unique peptides"]])
M <- as.matrix(rw[, si$SampleID]); storage.mode(M) <- "numeric"; rownames(M) <- acc
M[!is.finite(M) | M <= 0] <- NA
L <- log2(M)
cat("proteins:", nrow(L), "; with >=3 unique peptides:", sum(upep >= 3, na.rm = TRUE), "\n")
complete <- rowSums(is.na(L)) == 0 & upep >= 3
cat("complete-case proteins (all 60, >=3 unique peptides) used for median/total factors:", sum(complete), "\n")
actb <- L["P60709", ]; stopifnot(all(is.finite(actb)))
f_median <- apply(L[complete, ], 2, median)
f_total  <- log2(colSums(M[complete, ]))
pg <- read.delim(file.path(dd, "DBS-Tohoku.pg_matrix.tsv"), check.names = FALSE)
lfq <- as.numeric(pg[pg$Protein.Group == "Q9UJS0", si$SampleID]); names(lfq) <- si$SampleID
lfq[!is.finite(lfq) | lfq <= 0] <- NA
a13 <- L["Q9UJS0", ]
norm <- list(
  "ACTB (manuscript)"      = a13 - actb,
  "Sample median"          = a13 - f_median,
  "Total intensity"        = a13 - f_total,
  "DIA-NN MaxLFQ"          = log2(lfq),
  "None (raw log2)"        = a13)
stopifnot(all(is.na(a13[si$group == "SS"])), all(is.na(lfq[si$group == "SS"])))

ratio_ci <- function(v, ctrl) {   # ratio of medians on the linear scale, bootstrap percentile CI
  v <- v[is.finite(v)]; ctrl <- ctrl[is.finite(ctrl)]
  obs <- 2^(median(v) - median(ctrl))
  bs <- replicate(B, 2^(median(sample(v, replace = TRUE)) - median(sample(ctrl, replace = TRUE))))
  c(ratio = obs, lo = unname(quantile(bs, .025)), hi = unname(quantile(bs, .975)))
}
rows <- list()
for (nm in names(norm)) {
  v <- norm[[nm]]; g <- si$group
  ctrl <- v[g %in% c("DD","WW")]
  ws <- ratio_ci(v[g == "WS"], ctrl); sm <- ratio_ci(v[g == "SM"], ctrl)
  smws <- ratio_ci(v[g == "SM"], v[g == "WS"])
  rows[[nm]] <- data.frame(normalisation = nm,
    median_DD = 2^median(v[g=="DD"]), median_WW = 2^median(v[g=="WW"]), median_WS = 2^median(v[g=="WS"]), median_SM = 2^median(v[g=="SM"]),
    n_quantified_SS = sum(is.finite(v[g=="SS"])),
    WS_vs_control = ws["ratio"], WS_lo = ws["lo"], WS_hi = ws["hi"], p_WS_vs_control = wilcox.test(v[g=="WS"], ctrl)$p.value,
    SM_vs_control = sm["ratio"], SM_lo = sm["lo"], SM_hi = sm["hi"], p_SM_vs_control = wilcox.test(v[g=="SM"], ctrl)$p.value,
    SM_vs_WS = smws["ratio"], SMWS_lo = smws["lo"], SMWS_hi = smws["hi"], p_SM_vs_WS = wilcox.test(v[g=="SM"], v[g=="WS"])$p.value,
    p_DD_vs_WW = wilcox.test(v[g=="DD"], v[g=="WW"])$p.value,
    auc_low_vs_normal = { lo <- v[g %in% c("WS","SM")]; hi <- v[g %in% c("DD","WW")]; mean(outer(lo, hi, "<")) + 0.5*mean(outer(lo, hi, "==")) },
    row.names = NULL)
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(out, "slc25a13_normalisation_sensitivity.csv"), row.names = FALSE)
cat("\n=== SLC25A13: ratios of medians (linear scale) with bootstrap 95% CI, by normalisation\n")
print(data.frame(normalisation = res$normalisation,
  WS_vs_ctrl = sprintf("%.2f (%.2f-%.2f) p=%.1e", res$WS_vs_control, res$WS_lo, res$WS_hi, res$p_WS_vs_control),
  SM_vs_ctrl = sprintf("%.2f (%.2f-%.2f) p=%.1e", res$SM_vs_control, res$SM_lo, res$SM_hi, res$p_SM_vs_control),
  SM_vs_WS   = sprintf("%.2f (%.2f-%.2f) p=%.2f", res$SM_vs_WS, res$SMWS_lo, res$SMWS_hi, res$p_SM_vs_WS),
  AUC_low_vs_normal = sprintf("%.3f", res$auc_low_vs_normal), p_DD_vs_WW = sprintf("%.2f", res$p_DD_vs_WW)), row.names = FALSE)
# rank concordance between normalisations
V <- sapply(norm, identity); V <- V[si$group != "SS", ]
cat("\n=== Spearman correlation of SLC25A13 values across normalisations (44 specimens with quantification)\n")
print(round(cor(V, method = "spearman", use = "pairwise"), 3))

# ---- ACTB stability -------------------------------------------------------
stab <- data.frame(SampleID = si$SampleID, group = si$group, batch = si$batch,
                   ACTB_log2 = actb, ACTB_minus_median = actb - f_median, ACTB_minus_total = actb - f_total,
                   median_factor = f_median, total_log2 = f_total)
write.csv(stab, file.path(out, "actb_per_specimen.csv"), row.names = FALSE)
summ <- function(x) sprintf("%.2f (%.2f-%.2f)", median(x), min(x), max(x))
cat("\n=== ACTB log2 raw intensity by group x batch: median (range); CV on linear scale\n")
tab <- aggregate(cbind(ACTB_log2) ~ group + batch, stab, function(x) c(n = length(x), med = median(x), min = min(x), max = max(x), cv = sd(2^x)/mean(2^x)))
tab <- do.call(data.frame, tab); names(tab) <- c("group","batch","n","median_log2","min_log2","max_log2","CV_linear"); print(tab, digits = 3)
write.csv(tab, file.path(out, "actb_stability_by_group_batch.csv"), row.names = FALSE)
cat("\nKruskal-Wallis ACTB log2 across groups: p =", signif(kruskal.test(ACTB_log2 ~ group, stab)$p.value, 2),
    "; across batches: p =", signif(kruskal.test(ACTB_log2 ~ batch, stab)$p.value, 2), "\n")
cat("ACTB relative to sample median (log2 difference) across groups: p =", signif(kruskal.test(ACTB_minus_median ~ group, stab)$p.value, 2),
    "; range of log2(ACTB) - median =", sprintf("%.2f to %.2f", min(stab$ACTB_minus_median), max(stab$ACTB_minus_median)), "\n")
cat("Spearman rho(ACTB log2, sample median factor) =", round(cor(stab$ACTB_log2, stab$median_factor, method = "spearman"), 3),
    "; rho(ACTB, total) =", round(cor(stab$ACTB_log2, stab$total_log2, method = "spearman"), 3), "\n")
cat("Overall ACTB linear CV across 60 specimens:", round(sd(2^actb)/mean(2^actb), 3), "; log2 range", sprintf("%.2f", diff(range(actb))), "\n")

# ---- Figure ---------------------------------------------------------------
glab <- c(DD = "Disease control", WW = "Healthy", WS = "Carrier", SM = "LoF/Mis", SS = "LoF/LoF")
stab$Group <- factor(glab[as.character(stab$group)], levels = glab)
pal <- c(Run1 = "#56B4E9", Run2 = "#E69F00", Run3 = "#0072B2")
p1 <- ggplot(stab, aes(Group, ACTB_log2)) + geom_boxplot(outlier.shape = NA, width = .6, colour = "grey40") +
  geom_jitter(aes(colour = batch), width = .15, size = 1.6) + scale_colour_manual(values = pal, name = "Batch") +
  labs(y = expression(log[2]~"ACTB raw intensity"), x = NULL, title = "A  ACTB raw intensity by group and batch") +
  theme_classic(base_size = 10) + theme(legend.position = "top")
p2 <- ggplot(stab, aes(Group, ACTB_minus_median)) +
  geom_boxplot(outlier.shape = NA, width = .6, colour = "grey40") + geom_jitter(aes(colour = batch), width = .15, size = 1.6) +
  scale_colour_manual(values = pal, guide = "none") +
  labs(y = expression(log[2]~"ACTB" - "sample median"), x = NULL, title = "B  ACTB relative to the sample median") + theme_classic(base_size = 10)
long <- do.call(rbind, lapply(names(norm), function(nm) data.frame(normalisation = nm, SampleID = si$SampleID, group = si$group, value = norm[[nm]])))
long <- long[long$group != "SS", ]; long$Group <- factor(glab[as.character(long$group)], levels = glab)
long$normalisation <- factor(long$normalisation, levels = names(norm))
ctrl_med <- aggregate(value ~ normalisation, long[long$group %in% c("DD","WW"), ], median)
long <- merge(long, setNames(ctrl_med, c("normalisation","ctrl")), by = "normalisation"); long$rel <- 2^(long$value - long$ctrl)
p3 <- ggplot(long, aes(Group, rel)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_boxplot(outlier.shape = NA, width = .6, colour = "grey40") + geom_jitter(width = .15, size = 1.2, alpha = .8) +
  facet_wrap(~ normalisation, nrow = 1) + scale_y_log10(breaks = c(0.1, 0.2, 0.3, 0.5, 1, 2), labels = c("0.1", "0.2", "0.3", "0.5", "1", "2")) +
  labs(y = "SLC25A13 relative to pooled control median", x = NULL, title = "C  SLC25A13 under five normalisation strategies (LoF/LoF: not quantified under any)") +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ---- Panel D (added 2026-11-06): quantified protein groups per specimen, by group and batch ----
# Counts over the 7,474 protein groups reported by DIA-NN (as cited in the Response); counts within the
# 5,199-protein analysis matrix (>=3 unique peptides) are written to the CSV alongside.
qc <- data.frame(SampleID = si$SampleID, group = si$group, batch = si$batch,
                 n_quantified_7474 = colSums(!is.na(M)), n_quantified_5199 = colSums(!is.na(M[upep >= 3 & !is.na(upep), ])))
qc$Group <- factor(glab[as.character(qc$group)], levels = glab)
write.csv(qc, file.path(out, "quantified_protein_groups_per_specimen.csv"), row.names = FALSE)
p4 <- ggplot(qc, aes(Group, n_quantified_7474)) + geom_boxplot(outlier.shape = NA, width = .6, colour = "grey40") +
  geom_jitter(aes(colour = batch), width = .15, size = 1.6) + scale_colour_manual(values = pal, guide = "none") +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(y = "Quantified protein groups per specimen", x = NULL,
       title = "D  Quantified protein groups per specimen by group and batch (of 7,474 reported)") +
  theme_classic(base_size = 10)
pdf(file.path(out, "FigureS1_ACTB_stability_and_SLC25A13_normalisation.pdf"), width = 11, height = 12, useDingbats = FALSE)
grid::grid.newpage(); lay <- grid::grid.layout(3, 2, heights = grid::unit(c(1, 1.15, 1), "null"))
grid::pushViewport(grid::viewport(layout = lay))
print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)); print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(p4, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1)); dev.off()
png(file.path(out, "FigureS1_ACTB_stability_and_SLC25A13_normalisation.png"), width = 11, height = 12, units = "in", res = 200)
grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = lay))
print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)); print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(p4, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1)); dev.off()
print(aggregate(cbind(n_quantified_7474, n_quantified_5199) ~ Group + batch, qc, function(v) c(median = median(v), min = min(v), max = max(v))))
cat("[OUTPUT]", file.path(out, "FigureS1_ACTB_stability_and_SLC25A13_normalisation.pdf"), "\n")
## summary statistics of the ACTB anchor quoted in the Methods and in the footnote of Supplementary Table S9
write.csv(data.frame(statistic = c("complete_case_proteins", "ACTB_linear_CV", "ACTB_log2_range", "ACTB_fold_range", "KruskalWallis_P_groups", "KruskalWallis_P_runs"),
  value = c(sum(complete), round(sd(2^actb)/mean(2^actb), 3), round(diff(range(actb)), 2), round(2^diff(range(actb)), 1),
            signif(kruskal.test(ACTB_log2 ~ group, stab)$p.value, 2), signif(kruskal.test(ACTB_log2 ~ batch, stab)$p.value, 2))),
  file.path(out, "actb_summary_statistics.csv"), row.names = FALSE)
