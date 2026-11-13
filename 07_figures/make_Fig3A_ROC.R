#!/usr/bin/env Rscript
## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Panel B (ROC) at 1.78 x 1.54 in for the composite: sources the v1.8 figure script (recomputes everything from the production
# out-of-fold predictions; same values as Table 2), then re-saves pC without title, with compact legend labels and small fonts.
R <- REPO
OUT <- file.path(REPO, "07_figures", "results", "Figure3"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PROJECT_ROOT_OVERRIDE <- ""; setwd(R)
sink(tempfile()); source(file.path(R, "07_figures", "citrin_classifier_figure_v1_8_final30.R"), encoding = "UTF-8"); sink()
suppressPackageStartupMessages(library(ggplot2))
short_labels <- setNames(vapply(STRATEGIES, function(st) {
  zz <- metric_check[metric_check$strategy == st & metric_check$policy == POLICIES[1], ]
  sprintf("%s  %.3f (%s)", strategy_title[st], zz$AUC, AUC_CI95[[st]]) }, character(1)), STRATEGIES)
pB <- ggplot(roc_dat, aes(x = fpr, y = tpr, color = strategy, linetype = strategy, group = strategy)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dashed", color = "#B9BDC5") +
  geom_step(direction = "hv", linewidth = 0.6, lineend = "butt") +
  scale_color_manual(values = strategy_cols, labels = short_labels, breaks = STRATEGIES, drop = FALSE) +
  scale_linetype_manual(values = strategy_linetypes, labels = short_labels, breaks = STRATEGIES, drop = FALSE) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, .5, 1)) + scale_y_continuous(breaks = c(0, .5, 1)) +
  labs(x = "1 − Specificity", y = "Sensitivity") +
  guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(linewidth = 0.7)), linetype = guide_legend(ncol = 1)) +
  theme_classic(base_size = 5.5) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 4.6), legend.key.width = grid::unit(0.45, "cm"),
        legend.key.height = grid::unit(0.22, "cm"), legend.spacing.y = grid::unit(0, "cm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(1, "pt"),
        axis.line = element_line(linewidth = .25), axis.ticks = element_line(linewidth = .25), axis.text = element_text(size = 4.6), axis.title = element_text(size = 5),
        plot.margin = margin(2, 2, 1, 2))
ggsave(file.path(OUT, "Fig3A_ROC_1.78x1.54in.pdf"), pB, width = 1.78, height = 1.54, device = cairo_pdf)
ggsave(file.path(OUT, "Fig3A_ROC_1.78x1.54in.png"), pB, width = 1.78, height = 1.54, dpi = 600, bg = "white")
print(short_labels); cat("ROC done\n")
