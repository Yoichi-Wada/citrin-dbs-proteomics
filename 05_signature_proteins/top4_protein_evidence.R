#!/usr/bin/env Rscript
## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
# Measurement evidence for the four most influential consensus30 proteins
# (PRL P01236, MTTP P55157, ARHGEF10 O15013, SEC31A O94979):
#  1. non-imputed values (ACTB ratio, log2) per specimen and group, zeros shown as "not quantified"
#  2. missing (not quantified) rate per group; Fisher test SM vs 37 controls of the 44 cohort
#  3. peptide evidence: unique peptides, precursors per specimen (DIA-NN pr_matrix), q-values (DIA-NN report summary)
#  4. specimen / family dependency of the single-protein discrimination (SM 7 vs 37 controls):
#     AUC with zeros ranked lowest (what MinProb imputation effectively does), leave-one-specimen-out and
#     leave-one-family-out AUC, detection-only Wilcoxon.
suppressPackageStartupMessages({library(openxlsx); library(ggplot2)})
D  <- DATA
OUT <- file.path(REPO, "05_signature_proteins", "results", "top4"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ACC <- c(P01236 = "PRL", P55157 = "MTTP", O15013 = "ARHGEF10", O94979 = "SEC31A")

d  <- read.xlsx(file.path(D, "raw_data.xlsx"), sheet = "data", sep.names = " ")
si <- read.csv(file.path(D, "sample_info.csv")); fm <- read.csv(file.path(EXT, "family_map.csv"))
sc <- grep("^(SS|SM|WS|WW|DD)_[0-9]+$", names(d), value = TRUE)
stopifnot(length(sc) == 60L, all(sc %in% si$SampleID), all(sc %in% fm$SampleID))
grp <- setNames(si$group, si$SampleID)[sc]; fam <- setNames(fm$family_id, fm$SampleID)[sc]
actb <- as.numeric(d[d[["Master accession number"]] == "P60709", sc]); stopifnot(all(actb > 0))
rows <- d[match(names(ACC), d[["Master accession number"]]), ]
stopifnot(!anyNA(rows$Number))
raw <- as.matrix(rows[, sc]); storage.mode(raw) <- "numeric"; rownames(raw) <- names(ACC)

# ---- 1/2. long table of non-imputed values ------------------------------------
long <- do.call(rbind, lapply(names(ACC), function(a) data.frame(
  accession = a, gene = ACC[[a]], SampleID = sc, group = grp, family = fam,
  raw_intensity = raw[a, ], quantified = raw[a, ] > 0,
  log2_ACTB_ratio = ifelse(raw[a, ] > 0, log2(raw[a, ] / actb), NA_real_))))

gorder <- c("SM", "DD", "WS", "WW", "SS")
miss <- do.call(rbind, lapply(names(ACC), function(a) {
  q <- raw[a, ] > 0
  cbind(data.frame(accession = a, gene = ACC[[a]], unique_peptides = rows[["Unique peptides"]][match(a, rows[["Master accession number"]])]),
        as.data.frame(t(sapply(gorder, function(g) sprintf("%d/%d", sum(q[grp == g]), sum(grp == g))))))
}))
names(miss)[4:8] <- paste0("quantified_", gorder)
# 44-cohort: SM vs controls (DD+WS+WW)
ctrl <- grp %in% c("DD", "WS", "WW"); pos <- grp == "SM"
miss$fisher_p_SM_vs_ctrl <- sapply(names(ACC), function(a) {
  q <- factor(raw[a, ] > 0, levels = c(FALSE, TRUE)); g <- factor(pos, levels = c(FALSE, TRUE))
  fisher.test(table(q[pos | ctrl], g[pos | ctrl]))$p.value })
miss$median_log2ratio_SM_quantified <- sapply(names(ACC), function(a) median(long$log2_ACTB_ratio[long$accession == a & long$group == "SM"], na.rm = TRUE))
miss$median_log2ratio_ctrl_quantified <- sapply(names(ACC), function(a) median(long$log2_ACTB_ratio[long$accession == a & long$group %in% c("DD","WS","WW")], na.rm = TRUE))
miss$wilcoxon_p_quantified_only <- sapply(names(ACC), function(a) {
  x <- long[long$accession == a & long$quantified, ]
  s <- x$log2_ACTB_ratio[x$group == "SM"]; c0 <- x$log2_ACTB_ratio[x$group %in% c("DD","WS","WW")]
  if (length(s) >= 2 && length(c0) >= 2) wilcox.test(s, c0)$p.value else NA_real_ })
write.csv(miss, file.path(OUT, "top4_missingness_and_group_difference.csv"), row.names = FALSE)

# ---- 3. precursor evidence per specimen from DIA-NN pr_matrix --------------------
pr <- read.delim(file.path(EXT, "DBS-Tohoku.pr_matrix.tsv"), check.names = FALSE, quote = "", na.strings = c("", "NA"))
prec <- do.call(rbind, lapply(names(ACC), function(a) {
  x <- pr[pr$Protein.Group == a & pr$Protein.Ids == a, ]
  m <- as.matrix(x[, sc]); storage.mode(m) <- "numeric"; det <- is.finite(m) & m > 0
  n_per <- colSums(det)
  data.frame(accession = a, gene = ACC[[a]], n_precursors_matrix = nrow(x),
             n_distinct_peptides_matrix = length(unique(x$Stripped.Sequence)),
             t(sapply(gorder, function(g) sprintf("%d/%d", sum(n_per[grp == g] > 0), sum(grp == g)))),
             median_precursors_per_quantified_specimen = median(n_per[n_per > 0]),
             specimens_with_1_precursor_only = sum(n_per == 1L),
             specimens_with_ge2_precursors = sum(n_per >= 2L), check.names = FALSE)
}))
names(prec)[5:9] <- paste0("specimens_with_precursor_", gorder)
write.csv(prec, file.path(OUT, "top4_precursor_evidence_pr_matrix.csv"), row.names = FALSE)
# per-specimen precursor counts (for the figure annotation)
prec_per <- do.call(rbind, lapply(names(ACC), function(a) {
  x <- pr[pr$Protein.Group == a & pr$Protein.Ids == a, ]
  m <- as.matrix(x[, sc]); storage.mode(m) <- "numeric"
  data.frame(accession = a, SampleID = sc, n_precursors = colSums(is.finite(m) & m > 0)) }))
long <- merge(long, prec_per, by = c("accession", "SampleID"))
long <- long[order(match(long$accession, names(ACC)), match(long$SampleID, sc)), ]
write.csv(long, file.path(OUT, "top4_nonimputed_values_60specimens.csv"), row.names = FALSE)

# ---- 4. specimen / family dependency of single-protein discrimination ----------
auc <- function(score, y) { # zeros ranked lowest; ties handled by rank
  r <- rank(score); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
keep <- pos | ctrl; y44 <- pos[keep]; f44 <- fam[keep]; s44 <- sc[keep]
dep <- do.call(rbind, lapply(names(ACC), function(a) {
  x <- ifelse(raw[a, keep] > 0, log2(raw[a, keep] / actb[keep]), -Inf)
  full <- auc(x, y44)
  loso <- sapply(seq_along(x), function(i) auc(x[-i], y44[-i]))
  lofo <- sapply(unique(f44), function(f) { k <- f44 != f; if (sum(y44[k]) < 2) NA else auc(x[k], y44[k]) })
  # direction-aware: if full < 0.5 the protein is lower in SM; report AUC in the informative direction
  dirn <- if (full >= 0.5) "higher in SM" else "lower in SM"
  orient <- function(v) if (full >= 0.5) v else 1 - v
  data.frame(accession = a, gene = ACC[[a]], direction = dirn,
             AUC_single_protein = round(orient(full), 3),
             LOSO_min = round(min(orient(loso)), 3), LOSO_max = round(max(orient(loso)), 3),
             LOSO_min_specimen = s44[which.min(orient(loso))],
             LOFO_min = round(min(orient(lofo), na.rm = TRUE), 3), LOFO_max = round(max(orient(lofo), na.rm = TRUE), 3),
             LOFO_min_family = unique(f44)[which.min(orient(lofo))],
             n_SM_quantified = sum(raw[a, keep][y44] > 0), n_ctrl_quantified = sum(raw[a, keep][!y44] > 0))
}))
write.csv(dep, file.path(OUT, "top4_specimen_family_dependency.csv"), row.names = FALSE)

# ---- figure: non-imputed values by group, not-quantified shown as open marks at the floor
long$group <- factor(long$group, levels = gorder, labels = c("LoF/Mis (SM)", "Disease ctrl (DD)", "Carrier (WS)", "Healthy (WW)", "LoF/LoF (SS)"))
long$gene <- factor(long$gene, levels = ACC)
floor_y <- tapply(long$log2_ACTB_ratio, long$gene, min, na.rm = TRUE) - 1.2
long$y <- ifelse(long$quantified, long$log2_ACTB_ratio, floor_y[as.character(long$gene)])
long$status <- factor(ifelse(long$quantified, "Quantified (non-imputed value)", "Not quantified (would be MinProb-imputed)"),
                      levels = c("Quantified (non-imputed value)", "Not quantified (would be MinProb-imputed)"))
set.seed(1)
p <- ggplot(long, aes(x = group, y = y)) +
  geom_jitter(aes(shape = status, colour = status), width = 0.18, height = 0, size = 1.9, stroke = 0.7) +
  stat_summary(data = long[long$quantified, ], fun = median, geom = "crossbar", width = 0.5, linewidth = 0.35, colour = "#17212B") +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  scale_colour_manual(values = c("#4269D0", "#9AA5B1"), name = NULL) +
  labs(x = NULL, y = "log2(protein / ACTB), non-imputed",
       caption = "Open circles: specimens in which the protein was not quantified (0 in raw_data.xlsx); plotted at an arbitrary floor. Bars: median of quantified values.") +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5), legend.position = "top",
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold"), plot.caption = element_text(size = 7, hjust = 0, colour = "#5B6B7C"))
ggsave(file.path(OUT, "FigureS_top4_proteins_nonimputed_values.pdf"), p, width = 9.5, height = 3.6, useDingbats = FALSE)
ggsave(file.path(OUT, "FigureS_top4_proteins_nonimputed_values.png"), p, width = 9.5, height = 3.6, dpi = 300, bg = "white")

cat("\n[missingness / group difference]\n"); print(miss, row.names = FALSE)
cat("\n[precursor evidence]\n"); print(prec, row.names = FALSE)
cat("\n[dependency]\n"); print(dep, row.names = FALSE)
