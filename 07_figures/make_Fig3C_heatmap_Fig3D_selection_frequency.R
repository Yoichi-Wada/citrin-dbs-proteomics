#!/usr/bin/env Rscript
## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Fig. 3B (heatmap of the 30 final-panel proteins) in the style of the user's Fig3A_Heatmap_m30_nonzero_LoFMS_vs_Others_exclLoFLoF.pdf
# (pheatmap defaults: RdYlBu palette, grey60 cell borders, row/column dendrograms, annotation bars hm_group + group), and
# Fig. 3E (recurrence of the 30 proteins) with the composite-draft colours (bars #FFCCAA, non-zero dots #4269D0). 2026-11-08.
suppressPackageStartupMessages({library(openxlsx); library(pheatmap); library(RColorBrewer); library(ggplot2); library(grid)})
D   <- DATA
P   <- REPO
OUT <- file.path(REPO, "07_figures", "results", "Figure3"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
panel <- read.csv(file.path(P, "05_signature_proteins", "results", "final30_panel.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(panel) == 30L)
d  <- read.xlsx(file.path(D, "raw_data.xlsx"), sheet = "data", sep.names = " ")
si <- read.csv(file.path(D, "sample_info.csv"), stringsAsFactors = FALSE)
sc <- grep("^(SS|SM|WS|WW|DD)_[0-9]+$", names(d), value = TRUE); stopifnot(length(sc) == 60L)
grp <- setNames(si$group, si$SampleID)[sc]
keep <- sc[grp != "SS"]; stopifnot(length(keep) == 44L)
actb <- as.numeric(d[d[["Master accession number"]] == "P60709", keep]); stopifnot(all(actb > 0))
rows <- d[match(panel$FeatureID, d[["Master accession number"]]), ]
stopifnot(identical(rows[["Master accession number"]], panel$FeatureID))
raw <- as.matrix(rows[, keep]); storage.mode(raw) <- "numeric"; rownames(raw) <- panel$Gene; colnames(raw) <- keep
x <- raw; x[!is.finite(x) | x <= 0] <- NA
x <- sweep(log2(x), 2, log2(actb), "-")                       # non-imputed log2(protein/ACTB)
z <- t(apply(x, 1, function(v) { ok <- is.finite(v); if (sum(ok) < 3 || sd(v[ok]) == 0) return(rep(NA_real_, length(v))); (v - mean(v[ok])) / sd(v[ok]) }))
dimnames(z) <- dimnames(x)
z0 <- z; z0[!is.finite(z0)] <- 0                               # clustering only (NA -> 0, as in the v1.x figure script)
dr <- dist(z0); dc <- dist(t(z0)); hr <- hclust(dr); hc <- hclust(dc)   # complete linkage, euclidean (pheatmap defaults)
# column dendrogram re-ordered so that LoF/Mis specimens are placed on the LEFT (leaf weights: LoF/Mis = 0, others = 1; topology unchanged)
wts <- ifelse(grp[keep] == "SM", 0, 1)
hc <- as.hclust(reorder(as.dendrogram(hc), wts, agglo.FUN = mean))
zc <- z; zc[is.finite(zc) & zc > 4] <- 4; zc[is.finite(zc) & zc < -4] <- -4   # display range as in the user's figure (-4 .. 4); keeps dims
g44 <- grp[keep]
ann <- data.frame(group = factor(c(WW = "Healthy adult", WS = "Carrier", SM = "LoF/Mis", DD = "Disease control")[g44], levels = c("LoF/Mis", "Healthy adult", "Carrier", "Disease control")), row.names = keep)
ann_cols <- list(group = c("LoF/Mis" = "#00DAE0", "Healthy adult" = "#F8766D", "Carrier" = "#A3A500", "Disease control" = "#E76BF3"))   # rFig2A colours; LoF/Mis in the cyan formerly used for hm_group
pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
brk <- seq(-4, 4, length.out = 101)
ph <- function(file, w, h, title) {
  pheatmap(zc, color = pal, breaks = brk, na_col = "#D9D9D9", border_color = "grey60",
           cluster_rows = hr, cluster_cols = hc, annotation_col = ann, annotation_colors = ann_cols,
           show_colnames = FALSE, fontsize = 8, fontsize_row = 7.5, treeheight_row = 28, treeheight_col = 36,
           main = title, filename = file, width = w, height = h)
}
ttl <- "Final consensus-30 panel (row z-score; grey = not quantified)"
ph(file.path(OUT, "Fig3C_heatmap_final30.pdf"), 8.5, 8.2, ttl)
ph(file.path(OUT, "Fig3C_heatmap_final30.png"), 8.5, 8.2, ttl)
write.csv(data.frame(row_order = rownames(z)[hr$order]), file.path(OUT, "Fig3C_row_order.csv"), row.names = FALSE)
write.csv(data.frame(column_order = colnames(z)[hc$order], group = g44[hc$order]), file.path(OUT, "Fig3C_column_order.csv"), row.names = FALSE)
write.csv(round(z, 3), file.path(OUT, "Fig3C_row_zscores_nonimputed.csv"))

# ---- Fig. 3E: recurrence, composite-draft colours -----------------------------------------------
rec <- read.csv(file.path(REPO, "07_figures", "results", "Figure3_classifier_panels", "PanelE_final30_feature_selection_recurrence.csv"), stringsAsFactors = FALSE)
rec <- rec[order(rec$selection_rate, decreasing = TRUE), ]
rec$Gene <- factor(rec$Gene, levels = rev(rec$Gene)); rec$lab <- sprintf("%.0f%%", 100 * rec$selection_rate)
pE <- ggplot(rec, aes(y = Gene)) +
  geom_vline(xintercept = c(.25, .5, .75, 1), colour = "#E7EBF0", linewidth = .3) +
  geom_col(aes(x = selection_rate), fill = "#FFCCAA", width = .64) +
  geom_point(aes(x = nonzero_rate), shape = 16, size = 3.4, colour = "white") +
  geom_point(aes(x = nonzero_rate), shape = 16, size = 2.6, colour = "#4269D0") +
  geom_text(aes(x = 1.045, label = lab), hjust = 0, size = 3, colour = "#1F2D3D") +
  scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep, breaks = c(0, .25, .5, .75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), expand = expansion(0)) +
  coord_cartesian(xlim = c(0, 1), clip = "off", expand = FALSE) +
  labs(x = "Proportion of ensemble member fits", y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 8.6, colour = "#17212B"), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(6, 58, 6, 6))
ggsave(file.path(OUT, "Fig3D_selection_frequency_final30.pdf"), pE, width = 7.5, height = 7.6, device = cairo_pdf)
ggsave(file.path(OUT, "Fig3D_selection_frequency_final30.png"), pE, width = 7.5, height = 7.6, dpi = 200, bg = "white")
cat("row order:", paste(rownames(z)[hr$order], collapse = ", "), "\n")
cat("done:", OUT, "\n")
