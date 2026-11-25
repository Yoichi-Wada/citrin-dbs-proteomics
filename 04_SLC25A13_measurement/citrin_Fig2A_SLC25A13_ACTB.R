#!/usr/bin/env Rscript
## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# ============================================================================
# Revised Fig. 2A: SLC25A13/ACTB per specimen by group (replaces the "QC_SLC25A13_over_ACTB_ratio_strip"
# plot of Analysis_full_20260313.R, section "1) SLC25A13/ACTB plot").
#
# Reviewer requests addressed:
#  - R1 #3  : state how undetected values are handled and whether an axis scale factor is used
#             -> zeros (not quantified by DIA-NN) are shown explicitly as "not detected" and excluded from
#                ratio statistics; y axis is the raw ratio x 10^4 (stated on the axis; no other scaling).
#  - R2 Mj6 : "approximately half" vs ratio scale -> medians per group and ratios of medians with
#             bootstrap 95% CI (carriers / LoF/Mis vs pooled controls; LoF/Mis vs carriers) drawn on the figure.
#  - R2 Mn5 : Shapiro-Wilk / Welch t against an all-zero group removed -> Fisher's exact test on detection
#             (LoF/LoF 0/16 vs 44/44) and Wilcoxon rank-sum tests on quantified values; P = 0.27 for LoF/Mis vs
#             carriers is reported with the ratio and CI, not as evidence of similarity.
#  - R4 Mj4 : "not detected by the current DBS DIA workflow" wording (no claim of biological absence).
#  - R4 Mn e: standardised group labels; internal codes (RefRef_CT, LoFMS...) removed.
#  - R2 Mj4 / R4 Mj2: analytical run can be encoded as point shape (SHOW_BATCH_SHAPE); off by default to keep the
#             submitted point design (batch x group is given in the Supplement).
#
# Statistics reproduce the canonical values used in the Response (actb_normalization_sensitivity, 2026-11-04):
#   same data (raw_data.xlsx, ACTB P60709, SLC25A13 Q9UJS0), same bootstrap (seed 20261104, B = 5000, ratio of
#   medians on the linear scale, percentile CI, pooled controls = disease controls + healthy adults).
# Inputs : R/data/raw_data.xlsx (sheet "data"), R/data/sample_info.csv (group, batch)
# Outputs: <project>/Figures/Figure2A_SLC25A13_ACTB_rebuttal/
#           Figure2A_SLC25A13_ACTB_by_group.pdf/.png, per_specimen_values.csv, group_statistics.csv, legend_draft.txt
# ============================================================================
suppressPackageStartupMessages({library(openxlsx); library(ggplot2)})

PATH_RAW_XLSX   <- file.path(DATA, "raw_data.xlsx")
PATH_SAMPLEINFO <- file.path(DATA, "sample_info.csv")
OUT <- file.path(REPO, "04_SLC25A13_measurement", "results", "Figure2A")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ACTB_ACC <- "P60709"; SLC25A13_ACC <- "Q9UJS0"
BOOT_SEED <- 20261104; B <- 5000
SCALE <- 1e4   # y axis = ratio x 10^4 (display only; stated on the axis)

# ---- data ----------------------------------------------------------------------------------------
si <- read.csv(PATH_SAMPLEINFO, stringsAsFactors = FALSE)
rw <- read.xlsx(PATH_RAW_XLSX, sheet = "data", check.names = FALSE, sep.names = " ")
stopifnot(all(si$SampleID %in% names(rw)), nrow(si) == 60)
acc <- rw[["Master accession number"]]
getrow <- function(a) { i <- which(acc == a); stopifnot(length(i) == 1); v <- as.numeric(rw[i, si$SampleID]); names(v) <- si$SampleID; v }
actb <- getrow(ACTB_ACC); slc <- getrow(SLC25A13_ACC)
stopifnot(all(is.finite(actb) & actb > 0))
slc[!is.finite(slc)] <- 0

GROUP_LAB <- c(DD = "Disease controls", WW = "Healthy controls", WS = "Carriers", SM = "LoF/Mis", SS = "LoF/LoF")
GROUP_ORDER <- c("WW", "WS", "SS", "SM", "DD")   # x-axis order as in the submitted Fig. 2A (RefRef_CT, RefLoF, LoFLoF, LoFMS, RefRef_DD)
df <- data.frame(SampleID = si$SampleID, group = si$group, batch = si$batch,
                 SLC25A13_raw = slc, ACTB_raw = actb,
                 quantified = slc > 0,
                 ratio = ifelse(slc > 0, slc / actb, 0),
                 log2_ratio = ifelse(slc > 0, log2(slc / actb), NA_real_), stringsAsFactors = FALSE)
stopifnot(all(df$group %in% GROUP_ORDER))
n_by <- table(factor(df$group, levels = GROUP_ORDER))
df$group_f <- factor(df$group, levels = GROUP_ORDER, labels = sprintf("%s\nn = %d", GROUP_LAB[GROUP_ORDER], n_by[GROUP_ORDER]))
write.csv(df[, c("SampleID","group","batch","SLC25A13_raw","ACTB_raw","quantified","ratio","log2_ratio")],
          file.path(OUT, "per_specimen_values.csv"), row.names = FALSE)

# ---- statistics (identical procedure to the 2026-11-04 sensitivity analysis) --------------------------
g <- df$group; v <- df$log2_ratio
ctrl <- v[g %in% c("DD", "WW")]
set.seed(BOOT_SEED)
ratio_ci <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  obs <- 2^(median(a) - median(b))
  bs <- replicate(B, 2^(median(sample(a, replace = TRUE)) - median(sample(b, replace = TRUE))))
  c(ratio = obs, lo = unname(quantile(bs, .025)), hi = unname(quantile(bs, .975)))
}
ws   <- ratio_ci(v[g == "WS"], ctrl)          # order of calls matters for seed reproduction
sm   <- ratio_ci(v[g == "SM"], ctrl)
smws <- ratio_ci(v[g == "SM"], v[g == "WS"])
p_ws  <- wilcox.test(v[g == "WS"], ctrl)$p.value
p_sm  <- wilcox.test(v[g == "SM"], ctrl)$p.value
p_smws <- wilcox.test(v[g == "SM"], v[g == "WS"])$p.value
p_ddww <- wilcox.test(v[g == "DD"], v[g == "WW"])$p.value
det <- table(LoFLoF = g == "SS", quantified = df$quantified)
p_fisher <- fisher.test(det)$p.value
med <- tapply(df$ratio[df$quantified], g[df$quantified], median)
iqr <- tapply(df$ratio[df$quantified], g[df$quantified], quantile, probs = c(.25, .75))

stats <- data.frame(
  comparison = c("Carrier vs pooled controls (DD+WW)", "LoF/Mis vs pooled controls (DD+WW)", "LoF/Mis vs carrier",
                 "Disease control vs healthy adult", "LoF/LoF detection vs all other groups"),
  statistic = c("ratio of medians", "ratio of medians", "ratio of medians", "Wilcoxon only", "Fisher exact on quantified / not"),
  estimate = c(ws["ratio"], sm["ratio"], smws["ratio"], NA, NA),
  ci_low = c(ws["lo"], sm["lo"], smws["lo"], NA, NA), ci_high = c(ws["hi"], sm["hi"], smws["hi"], NA, NA),
  P = c(p_ws, p_sm, p_smws, p_ddww, p_fisher),
  test = c("Wilcoxon rank-sum", "Wilcoxon rank-sum", "Wilcoxon rank-sum", "Wilcoxon rank-sum", "Fisher exact (0/16 vs 44/44)"),
  row.names = NULL)
medtab <- data.frame(group = GROUP_ORDER, n = as.integer(n_by[GROUP_ORDER]),
                     n_quantified = as.integer(tapply(df$quantified, factor(g, levels = GROUP_ORDER), sum)),
                     median_ratio = as.numeric(med[GROUP_ORDER]),
                     q25 = sapply(GROUP_ORDER, function(k) if (is.null(iqr[[k]])) NA else iqr[[k]][1]),
                     q75 = sapply(GROUP_ORDER, function(k) if (is.null(iqr[[k]])) NA else iqr[[k]][2]), row.names = NULL)
write.csv(stats, file.path(OUT, "group_statistics.csv"), row.names = FALSE)
write.csv(medtab, file.path(OUT, "group_medians.csv"), row.names = FALSE)

# ---- figure ----------------------------------------------------------------------------------------
# Point design: filled circles with a dark outline, alpha 0.7, jitter width 0.18 (seed 1); group colours as in Fig. S7.
# Analytical run is NOT encoded in this panel (batch x group is tabulated in the Supplement); set SHOW_BATCH_SHAPE <- TRUE
# to encode it as point shape instead.
SHOW_BATCH_SHAPE <- FALSE
pal <- c(WW = "#FF725C", WS = "#EFB118", SS = "#6CC5B0", SM = "#97BBF5", DD = "#FF8AB7")
fmtP <- function(p) if (p < 1e-4) sprintf("P = %.1e", p) else sprintf("P = %.2f", p)
fmtR <- function(x) sprintf("%.2f (%.2f–%.2f)", x["ratio"], x["lo"], x["hi"])
ymax <- max(df$ratio) * SCALE
lab_df <- data.frame(group_f = levels(df$group_f)[match(c("WS", "SM", "SS"), GROUP_ORDER)],
                     y = ymax * c(1.11, 1.11, 0.40),
                     label = c(paste0("ratio of medians vs controls\n", fmtR(ws), "\n", fmtP(p_ws)),
                               paste0("ratio of medians vs controls\n", fmtR(sm), "\n", fmtP(p_sm),
                                      "\nvs carriers ", fmtR(smws), ", ", fmtP(p_smws)),
                               paste0("not detected by the DBS DIA workflow\n0/16 (all other groups 44/44)\nFisher ", fmtP(p_fisher))))
qdf <- df[df$quantified, ]
set.seed(1)
p <- ggplot(df, aes(x = group_f, y = ratio * SCALE)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
  # IQR as a thin bar, median as a thick bar (quantified specimens only)
  stat_summary(data = qdf, fun.min = function(z) quantile(z, .25), fun.max = function(z) quantile(z, .75), fun = median,
               geom = "linerange", linewidth = 0.6, colour = "#17212B") +
  stat_summary(data = qdf, fun = median, geom = "crossbar", width = 0.5, linewidth = 0.5, colour = "#17212B") +
  { if (SHOW_BATCH_SHAPE)
      geom_point(aes(fill = group, shape = batch), position = position_jitter(width = 0.18, height = 0, seed = 1),
                 size = 3, stroke = 0.6, colour = "black", alpha = 0.7)
    else
      geom_point(aes(fill = group), shape = 21, position = position_jitter(width = 0.18, height = 0, seed = 1),
                 size = 3, stroke = 0.6, colour = "black", alpha = 0.7) } +
  geom_text(data = lab_df, aes(x = group_f, y = y, label = label), size = 2.6, colour = "#17212B", lineheight = 0.95) +
  # explicit axis limits: the first trained layer (medians of quantified specimens) has no LoF/LoF, and the discrete scale
  # would otherwise drop that level and re-append it at the end of the axis
  scale_x_discrete(limits = levels(df$group_f)) +
  scale_fill_manual(values = pal, guide = "none") +
  { if (SHOW_BATCH_SHAPE) scale_shape_manual(values = c(Run1 = 24, Run2 = 22, Run3 = 21), name = "Analytical run") } +
  scale_y_continuous(limits = c(-0.1, ymax * 1.22), breaks = seq(0, 6, 1), expand = c(0, 0)) +
  labs(x = NULL, y = expression("SLC25A13/ACTB intensity ratio (" %*% 10^-4 * ")"),
       title = "SLC25A13 protein relative to ACTB in dried blood spots",
       subtitle = "Points: specimens (raw DIA-NN intensities, no imputation). Bars: median and interquartile range of quantified specimens.") +
  guides(shape = if (SHOW_BATCH_SHAPE) guide_legend(override.aes = list(fill = "grey60")) else "none") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = if (SHOW_BATCH_SHAPE) "bottom" else "none", plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8, colour = "#5B6B7C"),
        axis.text.x = element_text(size = 8.5))
ggsave(file.path(OUT, "Figure2A_SLC25A13_ACTB_by_group.pdf"), p, width = 8.2, height = 5.6, device = cairo_pdf)
ggsave(file.path(OUT, "Figure2A_SLC25A13_ACTB_by_group.png"), p, width = 8.2, height = 5.6, dpi = 300, bg = "white")

# ---- portrait version (width:height = 350:400) = main-figure panel (Fig. 2A): group names only, LoF/LoF and LoF/Mis
#      bracketed as citrin deficiency, no title and no statistics in the plot (they are in group_statistics.csv and the legend) ----
LAB_PANEL <- c(WW = "Healthy\ncontrols", WS = "Carriers", SS = "LoF/LoF", SM = "LoF/Mis", DD = "Disease\ncontrols")
df$group_p <- factor(df$group, levels = GROUP_ORDER, labels = LAB_PANEL[GROUP_ORDER])
qdf <- df[df$quantified, ]
xb <- match(c("SS", "SM"), GROUP_ORDER)            # bracket under the two citrin-deficiency groups
set.seed(1)
p2 <- ggplot(df, aes(x = group_p, y = ratio * SCALE)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
  stat_summary(data = qdf, fun.min = function(z) quantile(z, .25), fun.max = function(z) quantile(z, .75), fun = median,
               geom = "linerange", linewidth = 0.6, colour = "black") +
  stat_summary(data = qdf, fun = median, geom = "crossbar", width = 0.5, linewidth = 0.5, colour = "black") +
  geom_point(aes(fill = group), shape = 21, position = position_jitter(width = 0.18, height = 0, seed = 1),
             size = 3, stroke = 0.6, colour = "black", alpha = 0.7) +
  annotate("segment", x = xb[1] - 0.4, xend = xb[2] + 0.4, y = -0.46, yend = -0.46, linewidth = 0.4, colour = "black") +
  annotate("text", x = mean(xb), y = -0.52, label = "Citrin deficiency", size = 7.5 / .pt, vjust = 1, colour = "black") +
  scale_x_discrete(limits = levels(df$group_p)) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(breaks = 0:6, labels = function(b) ifelse(b %% 2 == 0, b, ""), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.1, ymax * 1.20), clip = "off") +   # limits on the coordinate system: the bracket is drawn below the panel
  labs(x = NULL, y = expression("SLC25A13/ACTB intensity ratio (" %*% 10^-4 * ")")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "none",
        axis.text = element_text(colour = "black"), axis.text.x = element_text(size = 7.5, lineheight = 0.95),
        axis.title.y = element_text(size = 9), plot.margin = margin(8, 8, 22, 8))
ggsave(file.path(OUT, "Figure2A_SLC25A13_ACTB_by_group_portrait.pdf"), p2, width = 5.25, height = 6.0, device = cairo_pdf)
ggsave(file.path(OUT, "Figure2A_SLC25A13_ACTB_by_group_portrait.png"), p2, width = 5.25, height = 6.0, dpi = 300, bg = "white")

# ---- legend draft ------------------------------------------------------------------------------------
leg <- sprintf(paste(
  "Fig. 2A. SLC25A13 protein relative to ACTB in dried blood spots. Each point is one specimen (raw DIA-NN intensity of SLC25A13 divided by that of ACTB in the same specimen, x 10^-4; no imputation, no other scaling);",
  if (SHOW_BATCH_SHAPE) "shapes indicate the analytical run." else "", "Bars show the median and interquartile range of specimens in which SLC25A13 was quantified.",
  "SLC25A13 was not detected by the DBS DIA workflow in any of the 16 LoF/LoF specimens and was quantified in all 44 other specimens (Fisher's exact test, P = %.1e);",
  "these specimens are shown at zero and were excluded from the ratio statistics. Relative to the pooled controls (disease controls and healthy adults, n = 27),",
  "the median SLC25A13/ACTB ratio was %.2f (bootstrap 95%% CI %.2f–%.2f; Wilcoxon rank-sum P = %.1e) in carriers and %.2f (%.2f–%.2f; P = %.1e) in LoF/Mis individuals;",
  "the LoF/Mis to carrier ratio was %.2f (%.2f–%.2f; P = %.2f). Disease controls and healthy adults did not differ (P = %.2f). Group medians (x 10^-4): disease controls %.2f, healthy adults %.2f, carriers %.2f, LoF/Mis %.2f.",
  "LoF, loss-of-function variant; Mis, missense variant. n: disease controls 21, healthy adults 6, carriers 10, LoF/Mis 7, LoF/LoF 16."),
  p_fisher, ws["ratio"], ws["lo"], ws["hi"], p_ws, sm["ratio"], sm["lo"], sm["hi"], p_sm, smws["ratio"], smws["lo"], smws["hi"], p_smws, p_ddww,
  med["DD"] * SCALE, med["WW"] * SCALE, med["WS"] * SCALE, med["SM"] * SCALE)
writeLines(leg, file.path(OUT, "legend_draft.txt"))

cat("\n[group medians]\n"); print(medtab, row.names = FALSE)
cat("\n[statistics]\n"); print(stats, row.names = FALSE)
cat("\n[legend]\n"); cat(leg, "\n")
cat("\nOutputs in", OUT, "\n")
