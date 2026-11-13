#!/usr/bin/env Rscript
# ============================================================================
# Citrin deficiency: precursor-level detection of SLC25A13-specific, shared
# SLC25A12/SLC25A13 and SLC25A12-specific peptides across all 60 DBS specimens
# (extended version of citrin_SLC25A13_precursor_detection_heatmap.R, 2026-11-04)
#
# Input:  R/data/DBS-Tohoku.pr_matrix.tsv (DIA-NN precursor matrix, 60 runs)
# Script: R/20261029_rebuttal/script/
# Output: R/20261029_rebuttal/Figures/FigureS_SLC25A13_SLC25A12_precursor_detection.{pdf,png}
#         + <stem>_summary.csv
#
# In RStudio, open 20261029_rebuttal.Rproj and run:
# source("script/citrin_SLC25A13_SLC25A12_precursor_detection_heatmap.R", encoding="UTF-8")
#
# BINARY DETECTION HEATMAP (not abundance). Positive finite precursor intensity
# in the DIA-NN matrix (1% precursor and protein-group FDR) = quantified;
# NA/blank/0 = not quantified. No imputation, ComBat or normalization.
#
# Three row blocks:
#   1. SLC25A13-specific: Protein.Group == Protein.Ids == "Q9UJS0", proteotypic
#      (8 precursors; GLLPQLLGVAPEK flagged as I/L-isobaric with SLC25A12).
#   2. Shared SLC25A12/SLC25A13: Protein.Ids contains both O75746 and Q9UJS0
#      (3 precursors: KIYSTLAGTR, LQVAGEITTGPR, YEGFFGLYR).
#   3. SLC25A12-specific positive controls: Protein.Group == Protein.Ids ==
#      "O75746", proteotypic, quantified in >= A12_MIN_FRACTION of specimens
#      (objective rule, stated in the caption; I/L-isobaric GLIPQLIGVAPEK excluded).
# Sample order within groups is as in the original script. Group order (WW, WS, SS, SM, DD)
# and group labels follow Fig. 2A (healthy adult, carrier, LoF/LoF, LoF/Mis, disease control; 2026-11-07).
# ============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop('Package ggplot2 is required. Run install.packages("ggplot2").')
}
suppressPackageStartupMessages(library(ggplot2))

A12_MIN_FRACTION <- 0.90   # SLC25A12-specific control rows: quantified in >= 90% of the 60 specimens

# --- Stable, project-relative paths -----------------------------------------
PROJECT_FILE <- "20261029_rebuttal.Rproj"
find_project_root <- function() {
  candidates <- character(0)
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    p <- tryCatch(rstudioapi::getActiveProject(), error = function(e) NULL)
    if (length(p) == 1L && !is.na(p) && nzchar(p)) candidates <- c(candidates, p)
  }
  p <- normalizePath(getwd(), mustWork = TRUE)
  for (i in seq_len(6L)) {
    candidates <- c(candidates, p)
    if (dirname(p) == p) break
    p <- dirname(p)
  }
  candidates <- unique(candidates)
  matches <- candidates[file.exists(file.path(candidates, PROJECT_FILE))]
  if (!length(matches)) stop("Open ", PROJECT_FILE, " in RStudio, then source this script.")
  normalizePath(matches[1L], mustWork = TRUE)
}

ROOT <- normalizePath(Sys.getenv("CITRIN_REPO", "."))   # repository root (public release)
INPUT <- file.path(Sys.getenv("CITRIN_EXTERNAL_DIR", file.path(ROOT, "data")), "DBS-Tohoku.pr_matrix.tsv")   # not distributed (jPOST deposit)
OUTPUT_DIR <- file.path(ROOT, "04_SLC25A13_measurement", "results")
OUTPUT_STEM <- "Figure2B_SLC25A13_SLC25A12_precursor_detection"
if (!file.exists(INPUT)) stop("Input TSV not found: ", INPUT)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Read DIA-NN precursor matrix -------------------------------------------
pr <- utils::read.delim(INPUT, sep = "\t", header = TRUE, check.names = FALSE,
                        stringsAsFactors = FALSE, quote = "", comment.char = "",
                        na.strings = c("", "NA"), fileEncoding = "UTF-8-BOM")
need <- c("Protein.Group", "Protein.Ids", "Genes", "Proteotypic",
          "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "Precursor.Id")
if (!all(need %in% names(pr))) {
  stop("Required DIA-NN column(s) absent: ", paste(setdiff(need, names(pr)), collapse = ", "))
}

sample_cols <- grep("^(SS|SM|WS|WW|DD)_[0-9]+$", names(pr), value = TRUE)
group_code <- sub("_.*$", "", sample_cols)
group_order <- c("WW", "WS", "SS", "SM", "DD")
expected_n <- c(SS = 16L, SM = 7L, WS = 10L, WW = 6L, DD = 21L)
observed_n <- table(factor(group_code, levels = group_order))
if (length(sample_cols) != 60L || !identical(as.integer(observed_n), as.integer(expected_n[group_order]))) {
  stop("Expected SS16/SM7/WS10/WW6/DD21 sample columns; got: ",
       paste(names(observed_n), observed_n, collapse = ", "))
}
sample_cols <- unlist(lapply(group_order, function(gr) {
  s <- sample_cols[group_code == gr]
  s[order(as.integer(sub("^[^_]+_", "", s)))]
}), use.names = FALSE)
ss_cols <- startsWith(sample_cols, "SS_")

presence_matrix <- function(df) {
  m <- as.matrix(df[, sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(m) <- "numeric")
  out <- is.finite(m) & m > 0
  rownames(out) <- df$Stripped.Sequence; colnames(out) <- sample_cols
  out
}
has_id <- function(ids, acc) grepl(paste0("(^|;)", acc, "(;|$)"), ids)
IL_signature <- function(s) gsub("[IL]", "J", toupper(s))

# --- Block 1: SLC25A13-specific (identical selection to the original script) -
slc13 <- pr[pr$Protein.Group == "Q9UJS0" & pr$Protein.Ids == "Q9UJS0" &
              as.numeric(pr$Proteotypic) == 1, , drop = FALSE]
if (nrow(slc13) != 8L || length(unique(slc13$Precursor.Id)) != 8L) {
  stop("Expected exactly 8 SLC25A13-only precursors; found ", nrow(slc13))
}
a12_all_seq <- unique(pr$Stripped.Sequence[has_id(pr$Protein.Ids, "O75746") & !is.na(pr$Stripped.Sequence)])
il_ambiguous13 <- IL_signature(slc13$Stripped.Sequence) %in% IL_signature(a12_all_seq)
if (!identical(sort(slc13$Stripped.Sequence[il_ambiguous13]), "GLLPQLLGVAPEK")) {
  warning("Unexpected SLC25A12 I/L-isobaric peptide set in block 1; inspect the annotation.")
}
present13 <- presence_matrix(slc13)

# --- Block 2: shared SLC25A12/SLC25A13 precursors ----------------------------
shared <- pr[has_id(pr$Protein.Ids, "O75746") & has_id(pr$Protein.Ids, "Q9UJS0"), , drop = FALSE]
if (nrow(shared) < 1L) stop("No shared SLC25A12/SLC25A13 precursors found.")
present_sh <- presence_matrix(shared)

# --- Block 3: SLC25A12-specific positive controls ----------------------------
slc12 <- pr[pr$Protein.Group == "O75746" & pr$Protein.Ids == "O75746" &
              as.numeric(pr$Proteotypic) == 1, , drop = FALSE]
slc13_all_seq <- unique(pr$Stripped.Sequence[has_id(pr$Protein.Ids, "Q9UJS0") & !is.na(pr$Stripped.Sequence)])
il_ambiguous12 <- IL_signature(slc12$Stripped.Sequence) %in% IL_signature(slc13_all_seq)
slc12 <- slc12[!il_ambiguous12, , drop = FALSE]                    # drop GLIPQLIGVAPEK etc.
present12_all <- presence_matrix(slc12)
keep12 <- rowSums(present12_all) >= ceiling(A12_MIN_FRACTION * length(sample_cols))
if (sum(keep12) < 1L) stop("No SLC25A12-specific precursor passes the ", A12_MIN_FRACTION, " rule.")
slc12 <- slc12[keep12, , drop = FALSE]
present12 <- present12_all[keep12, , drop = FALSE]

# --- Safeguards --------------------------------------------------------------
if (any(present13[, ss_cols, drop = FALSE])) stop("SLC25A13-specific detection found in SS: verify filtering.")
if (any(colSums(present13[, !ss_cols, drop = FALSE]) < 1L)) warning("A non-SS specimen lacks SLC25A13-specific precursors.")
if (any(colSums(rbind(present_sh, present12)[, ss_cols, drop = FALSE]) < 1L)) {
  warning("An SS specimen lacks shared/SLC25A12 control precursors: check specimen quality.")
}

# --- Long format ---------------------------------------------------------------
block_labels <- c(SLC13 = sprintf("SLC25A13-specific\n(n = %d)", nrow(present13)),
                  SH    = sprintf("Shared\nA12 / A13 (n = %d)", nrow(present_sh)),
                  SLC12 = sprintf("SLC25A12-specific\ncontrols (n = %d)", nrow(present12)))
make_long <- function(present, block, il_flag = rep(FALSE, nrow(present))) {
  peptides <- rownames(present)
  n_det <- rowSums(present)
  ord <- peptides[order(-n_det, peptides)]
  d <- expand.grid(peptide = peptides, SampleID = sample_cols, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  d$detected <- as.vector(present)
  d$block <- block
  d$label <- paste0(d$peptide, ifelse(il_flag[match(d$peptide, peptides)], "*", ""))
  d$row_order <- match(d$peptide, rev(ord))
  d
}
hmap <- rbind(make_long(present13, "SLC13", il_ambiguous13),
              make_long(present_sh, "SH"),
              make_long(present12, "SLC12"))
hmap$Detection <- factor(ifelse(hmap$detected, "Quantified", "Not quantified"),
                         levels = c("Quantified", "Not quantified"))
hmap$GroupCode <- sub("_.*$", "", hmap$SampleID)
group_labels <- c(WW = "Healthy adult\n(Ref/Ref), n = 6", WS = "Carrier\n(Ref/LoF), n = 10", SS = "LoF/LoF\nn = 16",
                  SM = "LoF/Mis\nn = 7", DD = "Disease control\n(Ref/Ref), n = 21")
hmap$Group <- factor(unname(group_labels[hmap$GroupCode]), levels = unname(group_labels[group_order]))
hmap$SampleID <- factor(hmap$SampleID, levels = sample_cols)
hmap$Block <- factor(unname(block_labels[hmap$block]), levels = unname(block_labels))
# a single y factor across blocks; facet_grid with free y keeps only the rows of each block
y_levels <- unlist(lapply(c("SLC13", "SH", "SLC12"), function(b) {
  d <- hmap[hmap$block == b, ]; u <- unique(d[, c("label", "row_order")]); u$label[order(u$row_order)]
}), use.names = FALSE)
hmap$label <- factor(hmap$label, levels = y_levels)

# --- Figure ------------------------------------------------------------------
p <- ggplot(hmap, aes(x = SampleID, y = label, fill = Detection)) +
  geom_tile(width = 0.94, height = 0.92, color = "#FFFFFF", linewidth = 0.22) +
  facet_grid(Block ~ Group, scales = "free", space = "free", switch = "y") +
  scale_fill_manual(values = c("Quantified" = "#4269D0", "Not quantified" = "#E8EDF3"), drop = FALSE, name = NULL) +
  scale_x_discrete(drop = TRUE, expand = expansion(add = 0.48)) +
  scale_y_discrete(expand = expansion(add = 0.30)) +
  labs(title = "Precursor-level detection of SLC25A13, shared and SLC25A12 peptides across 60 DBS specimens",
       subtitle = "DIA-NN precursors at 1% precursor and protein-group FDR; pre-imputation quantification",
       x = NULL, y = NULL,
       caption = paste0(
         "Detection: positive, finite precursor intensity in DIA-NN pr_matrix.tsv; blank/0: not quantified. ",
         "SLC25A12-specific control rows: proteotypic SLC25A12 precursors quantified in >= ",
         round(100 * A12_MIN_FRACTION), "% of specimens.\n",
         "* GLLPQLLGVAPEK is isobaric with SLC25A12 GLIPQLIGVAPEK (I/L ambiguity); DIA-NN assigns the two by retention time, ",
         "so this row is not used as evidence of sequence specificity.")) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12, color = "#17212B"),
        plot.subtitle = element_text(size = 9.3, color = "#445467"),
        plot.caption = element_text(size = 8, color = "#445467", hjust = 0, lineheight = 1.15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.50, size = 6.8, color = "#233147"),
        axis.text.y = element_text(size = 8.9, color = "#17212B"),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = grid::unit(0.20, "cm"),
        panel.spacing.y = grid::unit(0.35, "cm"),
        strip.background = element_rect(fill = "#EEF2F7", color = NA),
        strip.text.x = element_text(size = 8.4, face = "bold", color = "#17212B", lineheight = 0.95),
        strip.text.y.left = element_text(size = 8.2, lineheight = 0.95, face = "bold", color = "#17212B", angle = 90),
        strip.placement = "outside",
        legend.position = "top", legend.justification = "right",
        legend.text = element_text(size = 9),
        legend.key.width = grid::unit(0.48, "cm"), legend.key.height = grid::unit(0.35, "cm"),
        plot.margin = margin(8, 10, 8, 8))

n_rows <- nrow(present13) + nrow(present_sh) + nrow(present12)
fig_h <- 2.6 + 0.33 * n_rows
out_pdf <- file.path(OUTPUT_DIR, paste0(OUTPUT_STEM, ".pdf"))
out_png <- file.path(OUTPUT_DIR, paste0(OUTPUT_STEM, ".png"))
ggsave(out_pdf, plot = p, width = 13, height = fig_h, units = "in", device = "pdf", useDingbats = FALSE, bg = "white")
ggsave(out_png, plot = p, width = 13, height = fig_h, units = "in", dpi = 300, bg = "white")

# --- Summary table -------------------------------------------------------------
summ <- function(present, block, il_flag = rep(FALSE, nrow(present))) {
  data.frame(Block = block, Peptide = rownames(present), IL_ambiguous = il_flag,
             SS = rowSums(present[, startsWith(sample_cols, "SS_"), drop = FALSE]),
             SM = rowSums(present[, startsWith(sample_cols, "SM_"), drop = FALSE]),
             WS = rowSums(present[, startsWith(sample_cols, "WS_"), drop = FALSE]),
             WW = rowSums(present[, startsWith(sample_cols, "WW_"), drop = FALSE]),
             DD = rowSums(present[, startsWith(sample_cols, "DD_"), drop = FALSE]),
             Total = rowSums(present), check.names = FALSE, row.names = NULL)
}
summary_table <- rbind(summ(present13, "SLC25A13-specific", il_ambiguous13),
                       summ(present_sh, "Shared"), summ(present12, "SLC25A12-specific"))
cat("\n[Detection counts by group, per precursor]\n"); print(summary_table, row.names = FALSE)
utils::write.csv(summary_table, file.path(OUTPUT_DIR, paste0(OUTPUT_STEM, "_summary.csv")), row.names = FALSE)
cat("\nSS specimens: SLC25A13-specific precursors per specimen = ",
    paste(range(colSums(present13[, ss_cols, drop = FALSE])), collapse = "-"),
    "; shared + SLC25A12 control precursors per specimen = ",
    paste(range(colSums(rbind(present_sh, present12)[, ss_cols, drop = FALSE])), collapse = "-"), "\n", sep = "")
cat("[OUTPUT] ", out_pdf, "\n[OUTPUT] ", out_png, "\n", sep = "")
