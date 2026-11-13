#!/usr/bin/env Rscript
## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Panels at the physical sizes used in the user's composite newFigure3_20261108_yw.svg (7 x 9.8 in page), so that they can be
# placed at 100 % without scaling: A heatmap 3.28 x 3.16 in; C recurrence 1.04 x 1.59 in; D margins 2.37 x 1.27 in.
# (B ROC is re-exported by make_Fig3_composite_sized_ROC.R.)
suppressPackageStartupMessages({library(openxlsx); library(pheatmap); library(RColorBrewer); library(ggplot2); library(grid)})
P <- REPO
OUT <- file.path(REPO, "07_figures", "results", "Figure3"); SRC <- OUT
# ---------- A: heatmap (re-use the data preparation of make_Fig3B_E_final30_pheatmap_style.R) ----------
src <- readLines(file.path(P, "07_figures", "make_Fig3C_heatmap_Fig3D_selection_frequency.R"))
eval(parse(text = src[1:(grep("^ph <- function", src) - 1)]))
OUT <- file.path(REPO, "07_figures", "results", "Figure3")   # (the sourced prefix resets OUT)
hm <- function(file, w, h) pheatmap(zc, color = pal, breaks = brk, na_col = "#D9D9D9", border_color = "grey60",
  cluster_rows = hr, cluster_cols = hc, annotation_col = ann, annotation_colors = ann_cols, annotation_names_col = FALSE,
  show_colnames = FALSE, fontsize = 5.5, fontsize_row = 5.5, treeheight_row = 12, treeheight_col = 16, legend = TRUE,
  cellwidth = NA, main = "", filename = file, width = w, height = h)
hm(file.path(OUT, "Fig3C_heatmap_3.28x3.16in.pdf"), 3.28, 3.16)
hm(file.path(OUT, "Fig3C_heatmap_3.28x3.16in.png"), 3.28, 3.16)
# ---------- C: recurrence, narrow (1.04 x 1.59 in) ----------
rec <- read.csv(file.path(REPO, "07_figures", "results", "Figure3_classifier_panels", "PanelE_final30_feature_selection_recurrence.csv"), stringsAsFactors = FALSE)
rec <- rec[order(rec$selection_rate, decreasing = TRUE), ]; rec$Gene <- factor(rec$Gene, levels = rev(rec$Gene)); rec$lab <- sprintf("%.0f%%", 100 * rec$selection_rate)
pE <- ggplot(rec, aes(y = Gene)) +
  geom_vline(xintercept = c(.25, .5, .75, 1), colour = "#E7EBF0", linewidth = .2) +
  geom_col(aes(x = selection_rate), fill = "#FFCCAA", width = .7) +
  geom_point(aes(x = nonzero_rate), shape = 16, size = 1.1, colour = "white") +
  geom_point(aes(x = nonzero_rate), shape = 16, size = 0.8, colour = "#4269D0") +
  geom_text(aes(x = 1.04, label = lab), hjust = 0, size = 1.3, colour = "#1F2D3D") +
  scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep, breaks = c(0, .5, 1), labels = c("0%", "50%", "100%"), expand = expansion(0)) +
  coord_cartesian(xlim = c(0, 1), clip = "off", expand = FALSE) +
  labs(x = "Member fits", y = NULL) +
  theme_classic(base_size = 5) +
  theme(axis.text.y = element_text(size = 4.3, colour = "#17212B"), axis.text.x = element_text(size = 4.3), axis.title.x = element_text(size = 4.8),
        axis.line = element_line(linewidth = .25), axis.ticks = element_line(linewidth = .25), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(2, 13, 2, 2))
ggsave(file.path(OUT, "Fig3D_selection_frequency_1.04x1.59in.pdf"), pE, width = 1.04, height = 1.59, device = cairo_pdf)
ggsave(file.path(OUT, "Fig3D_selection_frequency_1.04x1.59in.png"), pE, width = 1.04, height = 1.59, dpi = 600, bg = "white")
# ---------- D: margins (2.37 x 1.27 in) ----------
E <- file.path(REPO, "01_nested_ensemble", "results")
p <- read.csv(file.path(E, "outer_predictions.csv"), stringsAsFactors = FALSE); p <- p[p$strategy == "consensus30", ]
p$margin <- p$probability - p$threshold
p$grp <- factor(c(SM = "LoF/Mis", DD = "Disease control", WS = "Carrier", WW = "Healthy adult")[p$group], levels = c("LoF/Mis", "Healthy adult", "Carrier", "Disease control"))
p$policy <- factor(p$policy, levels = c("full_inner_recall", "inner_spec90"), labels = c("Primary threshold", "Specificity ≥ 0.90 threshold"))
p <- p[order(p$policy, -p$margin), ]; p$pos <- ave(seq_len(nrow(p)), p$policy, FUN = seq_along); p$called <- p$margin >= 0
p$lab <- ifelse(p$group == "SM", p$SampleID, "")
p$stag <- 0; for (pl in levels(p$policy)) { ix <- which(p$policy == pl & p$lab != ""); if (length(ix) > 1) { st <- 0; for (k in seq_along(ix)) { if (k > 1 && p$pos[ix[k]] - p$pos[ix[k - 1]] <= 2) st <- 1 - st else st <- 0; p$stag[ix[k]] <- st } } }
cols <- c("LoF/Mis" = "#00DAE0", "Healthy adult" = "#F8766D", "Carrier" = "#A3A500", "Disease control" = "#E76BF3")
summ <- aggregate(cbind(tp = called & group == "SM", fp = called & group != "SM") ~ policy, data = p, FUN = sum)
fpn <- sapply(levels(p$policy), function(pl) paste(p$SampleID[p$policy == pl & p$called & p$group != "SM"], collapse = ", "))
wrapfp <- function(x) { v <- strsplit(x, ", ")[[1]]; if (length(v) > 4) paste(paste(v[1:3], collapse = ", "), paste(v[-(1:3)], collapse = ", "), sep = ",\n") else x }
summ$lab <- sprintf("LoF/Mis %d/7; FP %d/37\nFP: %s", summ$tp, summ$fp, sapply(fpn[as.character(summ$policy)], wrapfp))
pD <- ggplot(p, aes(x = pos, y = margin)) +
  geom_hline(yintercept = 0, colour = "#17212B", linewidth = .25) +
  geom_col(aes(fill = grp), width = .8, colour = "#17212B", linewidth = .08, alpha = 0.7) +
  geom_text(aes(label = lab, y = ifelse(margin >= 0, margin + 0.015 + 0.13 * stag, margin - 0.015 - 0.13 * stag), hjust = ifelse(margin >= 0, 0, 1)), size = 1.15, angle = 90, vjust = 0.5, colour = "#17212B") +
  geom_text(data = summ, aes(x = 44.5, y = 1.25, label = lab), hjust = 1, vjust = 1, size = 1.3, lineheight = 0.9, colour = "#17212B") +
  facet_wrap(~ policy, nrow = 1) +
  scale_fill_manual(values = cols, name = NULL) + guides(fill = guide_legend(override.aes = list(alpha = 0.7), nrow = 1)) +
  scale_y_continuous(limits = c(-0.8, 1.25), breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.9), labels = c("-0.6", "-0.3", "0.0", "0.3", "0.6", "0.9")) +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 1)) +
  labs(x = "44 specimens, ordered by margin", y = "OOF probability − threshold") +
  theme_minimal(base_size = 5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(linewidth = .15),
        legend.position = "bottom", legend.key.size = unit(4, "pt"), legend.text = element_text(size = 4.5), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "pt"),
        strip.text = element_text(face = "bold", size = 4.8), axis.title = element_text(size = 4.8), axis.text = element_text(size = 4.3),
        panel.spacing = unit(4, "pt"), plot.margin = margin(2, 3, 1, 2))
ggsave(file.path(OUT, "Fig3B_margins_2.37x1.27in.pdf"), pD, width = 2.37, height = 1.27, device = cairo_pdf)
ggsave(file.path(OUT, "Fig3B_margins_2.37x1.27in.png"), pD, width = 2.37, height = 1.27, dpi = 600, bg = "white")
cat("done\n")
