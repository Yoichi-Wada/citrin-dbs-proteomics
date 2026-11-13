#!/usr/bin/env Rscript
## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Panel D draft: out-of-fold score minus the fold-specific training-derived threshold for ALL 44 specimens,
# three strategies x two threshold rules; bars above zero are called positive. Source: WW-fixed production ensemble.
suppressPackageStartupMessages({library(ggplot2)})
E <- file.path(REPO, "01_nested_ensemble", "results"); O <- file.path(REPO, "07_figures", "results", "Figure3"); dir.create(O, recursive = TRUE, showWarnings = FALSE)
p <- read.csv(file.path(E, "outer_predictions.csv"), stringsAsFactors = FALSE)
p <- p[p$strategy == "consensus30", ]
p$margin <- p$probability - p$threshold
p$grp <- factor(c(SM = "LoF/Mis", DD = "Disease control", WS = "Carrier", WW = "Healthy adult")[p$group], levels = c("LoF/Mis", "Healthy adult", "Carrier", "Disease control"))
p$strategy <- factor(p$strategy, levels = c("consensus30", "diverse30", "adaptive"), labels = c("Consensus-30 (primary)", "Diverse-30", "Adaptive"))
p$policy <- factor(p$policy, levels = c("full_inner_recall", "inner_spec90"), labels = c("Primary threshold (complete inner sensitivity)", "Specificity-constrained threshold (inner specificity ≥ 0.90)"))
# order specimens within each facet by margin (descending)
p <- p[order(p$strategy, p$policy, -p$margin), ]
p$pos <- ave(seq_len(nrow(p)), p$strategy, p$policy, FUN = seq_along)
p$called <- p$margin >= 0
p$lab <- ifelse(p$group == "SM" | (p$called & p$group != "SM"), p$SampleID, "")
cols <- c("LoF/Mis" = "#00DAE0", "Healthy adult" = "#F8766D", "Carrier" = "#A3A500", "Disease control" = "#E76BF3")   # rFig2A colours; LoF/Mis cyan as in the heatmap
summ <- aggregate(cbind(tp = called & group == "SM", fp = called & group != "SM") ~ strategy + policy, data = p, FUN = sum)
summ$lab <- sprintf("LoF/Mis detected %d/7; false positives %d/37", summ$tp, summ$fp)
g <- ggplot(p, aes(x = pos, y = margin)) +
  geom_hline(yintercept = 0, colour = "#17212B", linewidth = .4) +
  geom_col(aes(fill = grp), width = .8, colour = "#17212B", linewidth = .15, alpha = 0.7) +   # lighter fill with dark outline, as in rFig2A
  geom_text(aes(label = lab, y = ifelse(margin >= 0, margin + 0.015, margin - 0.015), hjust = ifelse(margin >= 0, 0, 1)), size = 2.3, angle = 90, vjust = 0.5, colour = "#17212B") +
  geom_text(data = summ, aes(x = 44, y = 0.95, label = lab), hjust = 1, vjust = 1, size = 3, colour = "#17212B") +
  facet_wrap(~ policy, nrow = 1) +
  scale_fill_manual(values = cols, name = NULL) + guides(fill = guide_legend(override.aes = list(alpha = 0.7))) +
  scale_y_continuous(limits = c(-0.7, 1.0), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0), labels = c("-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 1)) +
  labs(x = "44 specimens, ordered by margin (consensus-30 ensemble)", y = "Out-of-fold probability minus fold-specific threshold") +
  theme_minimal(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 9), axis.title = element_text(size = 9), legend.text = element_text(size = 9))
ggsave(file.path(O, "Fig3B_margins_consensus30_all44.pdf"), g, width = 10, height = 4.8, device = cairo_pdf)
ggsave(file.path(O, "Fig3B_margins_consensus30_all44.png"), g, width = 10, height = 4.8, dpi = 200, bg = "white")
print(summ[, c("strategy", "policy", "tp", "fp")])
