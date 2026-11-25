#!/usr/bin/env Rscript
## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Supplementary QC figure/table for the 15 proteins with non-zero coefficients in the submitted full-data fit.
#  A  fraction of specimens in which the protein was quantified (raw intensity > 0, before MinProb imputation), by group
#  B  non-imputed log2(protein/ACTB) per specimen by group; not-quantified specimens shown as open circles at a floor
#  C  selection stability of each protein across the 35 family-held-out outer folds of the primary consensus30 ensemble
#     (tier A/B/C as in main Fig. 3A; source ensemble_v1_1_results/tableS8_extended_ensemble_tiers.csv)
# Protein order (rows / facets) = tier table order = Fig. 3A row order.
# Extends protein_evidence_top4/top4_protein_evidence.R (same data sources, same definitions) to all 15 proteins.
suppressPackageStartupMessages({library(openxlsx); library(ggplot2); library(patchwork)})
D   <- DATA
P   <- REPO
OUT <- file.path(REPO, "05_signature_proteins", "results")

tier <- read.csv(file.path(OUT, "final30_ensemble_tiers.csv"), stringsAsFactors = FALSE)
ACC  <- setNames(tier$gene, tier$accession)            # order = tier order (A -> C)

d  <- read.xlsx(file.path(D, "raw_data.xlsx"), sheet = "data", sep.names = " ")
si <- read.csv(file.path(D, "sample_info.csv")); fm <- read.csv(file.path(EXT, "family_map.csv"))
sc <- grep("^(SS|SM|WS|WW|DD)_[0-9]+$", names(d), value = TRUE)
stopifnot(length(sc) == 60L, all(sc %in% si$SampleID), all(sc %in% fm$SampleID))
grp <- setNames(si$group, si$SampleID)[sc]; fam <- setNames(fm$family_id, fm$SampleID)[sc]
actb <- as.numeric(d[d[["Master accession number"]] == "P60709", sc]); stopifnot(all(actb > 0))
rows <- d[match(names(ACC), d[["Master accession number"]]), ]
stopifnot(!anyNA(rows$Number), identical(rows[["Master accession number"]], names(ACC)))
stopifnot(length(ACC) == 30L)
raw <- as.matrix(rows[, sc]); storage.mode(raw) <- "numeric"; rownames(raw) <- names(ACC)

gorder <- c("SM", "DD", "WS", "WW", "SS")
glab   <- c(SM = "LoF/Mis", DD = "Disease control", WS = "Carrier", WW = "Healthy", SS = "LoF/LoF")

# ---- B table: non-imputed values per specimen ------------------------------------------------
long <- do.call(rbind, lapply(names(ACC), function(a) data.frame(
  accession = a, gene = ACC[[a]], SampleID = sc, group = grp, family = fam,
  raw_intensity = raw[a, ], quantified = raw[a, ] > 0,
  log2_ACTB_ratio = ifelse(raw[a, ] > 0, log2(raw[a, ] / actb), NA_real_), stringsAsFactors = FALSE)))

# ---- A table: quantified fraction by group, Fisher SM vs 37 controls, medians of quantified values ----
ctrl <- grp %in% c("DD", "WS", "WW"); pos <- grp == "SM"
miss <- do.call(rbind, lapply(names(ACC), function(a) {
  q <- raw[a, ] > 0
  out <- data.frame(accession = a, gene = ACC[[a]],
                    unique_peptides = rows[["Unique peptides"]][match(a, rows[["Master accession number"]])],
                    stringsAsFactors = FALSE)
  for (g in gorder) out[[paste0("quantified_", g)]] <- sprintf("%d/%d", sum(q[grp == g]), sum(grp == g))
  for (g in gorder) out[[paste0("n_quantified_", g)]] <- sum(q[grp == g])
  out$quantified_44 <- sum(q[pos | ctrl])
  qf <- factor(q, levels = c(FALSE, TRUE)); gf <- factor(pos, levels = c(FALSE, TRUE))
  out$fisher_p_quantified_SM_vs_ctrl <- fisher.test(table(qf[pos | ctrl], gf[pos | ctrl]))$p.value
  x <- long[long$accession == a & long$quantified, ]
  s <- x$log2_ACTB_ratio[x$group == "SM"]; c0 <- x$log2_ACTB_ratio[x$group %in% c("DD", "WS", "WW")]
  out$median_log2ratio_SM_quantified   <- if (length(s))  median(s)  else NA_real_
  out$median_log2ratio_ctrl_quantified <- if (length(c0)) median(c0) else NA_real_
  out$wilcoxon_p_quantified_only <- if (length(s) >= 2 && length(c0) >= 2) wilcox.test(s, c0)$p.value else NA_real_
  out }))

# ---- peptide / precursor evidence from DIA-NN pr_matrix ------------------------------------
pr <- read.delim(file.path(EXT, "DBS-Tohoku.pr_matrix.tsv"), check.names = FALSE, quote = "", na.strings = c("", "NA"))
prec <- do.call(rbind, lapply(names(ACC), function(a) {
  x <- pr[pr$Protein.Group == a & pr$Protein.Ids == a, ]
  m <- as.matrix(x[, sc]); storage.mode(m) <- "numeric"; det <- is.finite(m) & m > 0
  n_per <- colSums(det)
  data.frame(accession = a, n_precursors_matrix = nrow(x),
             n_distinct_peptides_matrix = length(unique(x$Stripped.Sequence)),
             median_precursors_per_quantified_specimen = if (any(n_per > 0)) median(n_per[n_per > 0]) else NA_real_,
             specimens_with_1_precursor_only = sum(n_per == 1L), stringsAsFactors = FALSE) }))
prec_per <- do.call(rbind, lapply(names(ACC), function(a) {
  x <- pr[pr$Protein.Group == a & pr$Protein.Ids == a, ]
  m <- as.matrix(x[, sc]); storage.mode(m) <- "numeric"
  data.frame(accession = a, SampleID = sc, n_precursors = colSums(is.finite(m) & m > 0)) }))
long <- merge(long, prec_per, by = c("accession", "SampleID"))
long <- long[order(match(long$accession, names(ACC)), match(long$SampleID, sc)), ]
write.csv(long, file.path(OUT, "final30_nonimputed_values_60specimens.csv"), row.names = FALSE)

# ---- single-specimen concentration of the raw signal in the 44-specimen modelling set ---------
# (same definition as Top_share_44 in data/citrin_original15_QC_summary.csv: share of the summed raw intensity
#  over the 44 non-LoF/LoF specimens contributed by the single largest specimen)
keep44 <- pos | ctrl
conc <- do.call(rbind, lapply(names(ACC), function(a) {
  v <- raw[a, keep44]; tot <- sum(v)
  data.frame(accession = a, top_specimen_44 = names(v)[which.max(v)],
             top_share_44 = if (tot > 0) max(v) / tot else NA_real_, stringsAsFactors = FALSE) }))
miss <- merge(miss, conc, by = "accession", sort = FALSE)
# ---- rule-based QC notes (replaces the QC_flag of the 2026-11-03 summary; criteria stated explicitly) ----
pr44 <- sapply(names(ACC), function(a) sum(long$n_precursors[long$accession == a & long$SampleID %in% sc[keep44]] > 0))
miss$precursor_detected_44 <- pr44[miss$accession]
miss$precursor_minus_protein_44 <- miss$precursor_detected_44 - miss$quantified_44
miss$qc_note <- sapply(seq_len(nrow(miss)), function(i) {
  n <- character(0)
  if (miss$quantified_44[i] / 44 < 0.75) n <- c(n, "quantified in <75% of the 44 modelling specimens")
  if (!is.na(miss$top_share_44[i]) && miss$top_share_44[i] >= 0.30) n <- c(n, sprintf("single specimen (%s) carries %.0f%% of the summed raw intensity", miss$top_specimen_44[i], 100 * miss$top_share_44[i]))
  if (miss$precursor_minus_protein_44[i] != 0) n <- c(n, sprintf("precursor-level vs protein-level detection differ by %+d specimens", miss$precursor_minus_protein_44[i]))
  if (length(n)) paste(n, collapse = "; ") else "" })
# ---- summary table (A + peptides + C) --------------------------------------------------------
summ <- merge(miss, prec, by = "accession", sort = FALSE)
summ <- merge(summ, tier[, c("accession", "final_members_nonzero_of_10", "final_mean_coefficient", "in_candidate_panel_of_35",
                             "folds_nonzero_any_member_of_35", "median_members_nonzero_of_10",
                             "absent_when_family_held_out", "tier")], by = "accession", sort = FALSE)
summ <- summ[match(names(ACC), summ$accession), ]
first <- c("tier", "gene", "accession", "final_members_nonzero_of_10", "final_mean_coefficient", "unique_peptides", "n_precursors_matrix",
           "n_distinct_peptides_matrix", paste0("quantified_", gorder), "quantified_44",
           "fisher_p_quantified_SM_vs_ctrl", "median_log2ratio_SM_quantified", "median_log2ratio_ctrl_quantified",
           "wilcoxon_p_quantified_only", "in_candidate_panel_of_35", "folds_nonzero_any_member_of_35",
           "median_members_nonzero_of_10", "absent_when_family_held_out",
           "median_precursors_per_quantified_specimen", "specimens_with_1_precursor_only", "top_specimen_44", "top_share_44", "precursor_detected_44", "precursor_minus_protein_44", "qc_note")
summ <- summ[, c(first, setdiff(names(summ), first))]
write.csv(summ, file.path(OUT, "TableS_final30_QC_summary.csv"), row.names = FALSE)

# (cross-checks against the 15-protein summaries omitted for the 30-protein panel)

# ---- Figure --------------------------------------------------------------------------------------
tier_col <- c(A = "#17212B", B = "#4269D0", C = "#9AA5B1")
lev_gene <- rev(ACC)                                   # top row = first protein
# A: quantified fraction tiles
A <- do.call(rbind, lapply(names(ACC), function(a) data.frame(
  gene = ACC[[a]], group = gorder, n = sapply(gorder, function(g) sum(raw[a, grp == g] > 0)),
  N = sapply(gorder, function(g) sum(grp == g)))))
A$frac <- A$n / A$N; A$gene <- factor(A$gene, levels = lev_gene)
A$group <- factor(A$group, levels = gorder, labels = sprintf("%s\n(n = %d)", glab[gorder], table(grp)[gorder]))
pA <- ggplot(A, aes(x = group, y = gene, fill = frac)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%d/%d", n, N), colour = frac > 0.55), size = 2.7, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#17212B")) +
  scale_fill_gradient(low = "#F4F6F8", high = "#4269D0", limits = c(0, 1), name = "Quantified\nfraction") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "A  Specimens with a quantified (non-imputed) value, by group") +
  theme_minimal(base_size = 9.5) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(face = "bold", colour = tier_col[rev(tier$tier)]),
        axis.text.x = element_text(size = 7.5), legend.position = "right", legend.key.height = unit(10, "pt"),
        plot.title = element_text(face = "bold", size = 10))

# C: outer-fold selection stability
Cd <- tier; Cd$gene <- factor(Cd$gene, levels = lev_gene)
Cd$lab <- ifelse(is.na(Cd$absent_when_family_held_out) | Cd$absent_when_family_held_out == "", "",
                 paste0("absent when fam ", gsub(";", ", ", Cd$absent_when_family_held_out), " held out"))
Cd$lab <- ifelse(Cd$in_candidate_panel_of_35 == 0, "never selected", Cd$lab)
pC <- ggplot(Cd, aes(y = gene)) +
  geom_col(aes(x = in_candidate_panel_of_35, fill = tier), width = 0.7) +
  geom_text(aes(x = pmax(in_candidate_panel_of_35, 0.5) + 0.6, label = lab), hjust = 0, size = 2.5, colour = "#5B6B7C") +
  scale_fill_manual(values = tier_col, name = "Stability tier", labels = c("A: all 35 folds", "B: 30–34 folds", "C: < 30 folds")) +
  scale_x_continuous(limits = c(0, 62), breaks = seq(0, 35, 5), expand = c(0, 0)) +
  labs(x = "Outer folds (of 35) in which the protein entered the consensus-30 panel",
       y = NULL, title = "C  Selection stability across the 35 family-held-out outer folds") +
  theme_minimal(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold", colour = tier_col[rev(tier$tier)]),
        legend.position = "bottom", legend.direction = "vertical", legend.key.size = unit(9, "pt"),
        plot.title = element_text(face = "bold", size = 10))

# B: non-imputed values, 15 facets in tier order
long$group_f <- factor(long$group, levels = gorder, labels = c("LoF/Mis", "Disease ctrl", "Carrier", "Healthy", "LoF/LoF"))
long$gene_f  <- factor(paste0(long$gene, "  (tier ", tier$tier[match(long$accession, tier$accession)], ")"),
                       levels = paste0(ACC, "  (tier ", tier$tier, ")"))
floor_y <- tapply(long$log2_ACTB_ratio, long$gene_f, function(v) if (all(is.na(v))) 0 else min(v, na.rm = TRUE) - 1.2)
long$y <- ifelse(long$quantified, long$log2_ACTB_ratio, floor_y[as.character(long$gene_f)])
long$status <- factor(ifelse(long$quantified, "Quantified (non-imputed value)", "Not quantified (MinProb-imputed in the model)"),
                      levels = c("Quantified (non-imputed value)", "Not quantified (MinProb-imputed in the model)"))
set.seed(1)
pB <- ggplot(long, aes(x = group_f, y = y)) +
  geom_jitter(aes(shape = status, colour = status), width = 0.2, height = 0, size = 1.4, stroke = 0.6) +
  stat_summary(data = long[long$quantified, ], fun = median, geom = "crossbar", width = 0.55, linewidth = 0.3, colour = "#17212B") +
  facet_wrap(~ gene_f, scales = "free_y", nrow = 5) +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  scale_colour_manual(values = c("#4269D0", "#9AA5B1"), name = NULL) +
  labs(x = NULL, y = "log2(protein / ACTB), non-imputed",
       title = "B  Non-imputed values by group (open circles: not quantified, drawn at an arbitrary floor; bars: median of quantified values)") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.8), legend.position = "bottom",
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", size = 8), plot.title = element_text(face = "bold", size = 10))

fig <- (pA | pC) / pB + plot_layout(heights = c(1, 1.45), widths = c(1, 1.15))
ggsave(file.path(OUT, "FigureS8_final30_QC.pdf"), fig, width = 13.5, height = 17, device = cairo_pdf)
ggsave(file.path(OUT, "FigureS8_final30_QC.png"), fig, width = 13.5, height = 17, dpi = 200, bg = "white")

# ---- xlsx + md ------------------------------------------------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "Summary"); writeData(wb, "Summary", summ)
addWorksheet(wb, "Non-imputed values"); writeData(wb, "Non-imputed values", long[, c("accession","gene","SampleID","group","family","raw_intensity","quantified","log2_ACTB_ratio","n_precursors")])
addWorksheet(wb, "README"); writeData(wb, "README", data.frame(note = c(
  "Summary: one row per protein of the final consensus-30 panel (Table S18), ordered by stability tier (A: all 35 outer folds; B: 30-34; C: < 30) and fold count.",
  "quantified_<group>: specimens with raw intensity > 0 in raw_data.xlsx (before MinProb imputation) / group size. SM LoF/Mis, DD disease control, WS carrier, WW healthy adult, SS LoF/LoF.",
  "fisher_p_quantified_SM_vs_ctrl: Fisher exact test of quantified vs not, 7 LoF/Mis vs 37 controls (44-specimen modelling set).",
  "median_log2ratio_*_quantified and wilcoxon_p_quantified_only use quantified specimens only (no imputed values).",
  "unique_peptides: DIA-NN protein-level unique peptide count (raw_data.xlsx). n_precursors_matrix / n_distinct_peptides_matrix: rows of DBS-Tohoku.pr_matrix.tsv assigned to the protein group.",
  "in_candidate_panel_of_35, folds_nonzero_any_member_of_35, median_members_nonzero_of_10, absent_when_family_held_out, tier: from the primary consensus30 ensemble (ensemble_v1_1_results/tableS8_extended_ensemble_tiers.csv); exploratory feature-selection stability in this cohort, not a biological classification.",
  "top_specimen_44 / top_share_44: specimen contributing the largest share of the summed raw intensity over the 44 non-LoF/LoF specimens, and that share (single-specimen concentration of the signal).",
  "Precursor rows are restricted to protein groups whose Protein.Group and Protein.Ids both equal the accession. For NAPB (Q9H115) three of eight peptides (AIAHYEQSADYYK, AIAHYEQSADYYKGEESNSSANK, LDQWLTTMLLR) are shared with NAPA (P54920) and assigned to the P54920 group; NAPB-specific precursors were detected in 40/44 modelling specimens whereas protein-level NAPB was quantified in 41/44.",
  "qc_note: rule-based flags — quantified in <75% of the 44 modelling specimens; a single specimen carrying >=30% of the summed raw intensity; precursor-level detection (>=1 precursor in the DIA-NN precursor matrix) differing from protein-level quantification.",
  "Non-imputed values: per specimen x protein; n_precursors = precursors with a finite positive value in the DIA-NN precursor matrix.")))
saveWorkbook(wb, file.path(OUT, "TableS_final30_QC_summary.xlsx"), overwrite = TRUE)

cat("\n[summary]\n"); print(summ[, c("tier","gene","unique_peptides", paste0("quantified_", gorder), "in_candidate_panel_of_35","absent_when_family_held_out")], row.names = FALSE)
cat("\nsessionInfo R:", R.version.string, "| ggplot2", as.character(packageVersion("ggplot2")), "| patchwork", as.character(packageVersion("patchwork")), "\n")
