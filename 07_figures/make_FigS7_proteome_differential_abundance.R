## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Regeneration of the submitted Supplementary Fig. S1 (now Supplementary Fig. S7):
# descriptive differential-abundance overview across SLC25A13 genotype groups.
# Data: WW-fixed raw_data.xlsx (DIA-NN protein groups, 60 DBS). No ComBat.
# Preprocessing as in the manuscript: log2 (0 -> NA), ACTB centering (P60709), Unique peptides >= 3,
# MinProb imputation (imputeLCMD::impute.MinProb, q = 0.01, seed 2025) on the full 60-specimen matrix.
# Model: limma, ~ 0 + group + run (analytical run as covariate); contrasts LoF/LoF vs Ref/Ref, LoF/Mis vs Ref/Ref, LoF/Mis vs LoF/LoF.
# Global test: F-test of the four group contrasts (ANOVA-like, batch-adjusted), BH; top 100 by adjusted P shown as z-score heatmap.
suppressPackageStartupMessages({library(openxlsx); library(limma); library(ggplot2); library(ggrepel); library(pheatmap); library(cowplot); library(grid)})
DD <- DATA
OUT <- file.path(REPO, "07_figures", "results", "FigureS7"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED <- 2025; MINPROB_Q <- 0.01; UNIQUE_PEPTIDE_MIN <- 3; FDR <- 0.05; LFC <- 1; TOPN <- 100
si <- read.csv(file.path(DD, "sample_info.csv"), stringsAsFactors = FALSE)
rw <- read.xlsx(file.path(DD, "raw_data.xlsx"), sheet = "data", check.names = FALSE, sep.names = " ")
quant_cols <- si$SampleID; stopifnot(all(quant_cols %in% names(rw)))
M <- as.matrix(rw[, quant_cols]); storage.mode(M) <- "numeric"
M[M == 0] <- NA; M <- log2(M)
acc <- toupper(trimws(as.character(rw[["Master accession number"]]))); gene <- as.character(rw[["Master Gene Symbol"]])
actb <- which(acc == "P60709"); stopifnot(length(actb) == 1)
M <- sweep(M, 2, M[actb, ], "-")                       # log2(X) - log2(ACTB)
keep <- as.numeric(rw[["Unique peptides"]]) >= UNIQUE_PEPTIDE_MIN
M <- M[keep, ]; acc <- acc[keep]; gene <- gene[keep]
dup <- duplicated(acc); M <- M[!dup, ]; acc <- acc[!dup]; gene <- gene[!dup]
rownames(M) <- acc
cat("features after filter:", nrow(M), " specimens:", ncol(M), "\n")
set.seed(SEED); Mi <- imputeLCMD::impute.MinProb(M, q = MINPROB_Q)
Mi <- Mi[rownames(Mi) != "P60709", ]                     # ACTB itself is constant after centering
gene <- gene[acc != "P60709"]; acc <- acc[acc != "P60709"]
grp_map <- c(WW = "RefRef_CT", WS = "RefLoF", SS = "LoFLoF", SM = "LoFMS", DD = "RefRef_DD")
si$grp <- factor(unname(grp_map[substr(si$SampleID, 1, 2)]), levels = names(sort(table(grp_map)))[0])  # placeholder
si$grp <- factor(unname(grp_map[substr(si$SampleID, 1, 2)]), levels = c("RefRef_CT","RefLoF","LoFLoF","LoFMS","RefRef_DD"))
si$run <- factor(si$batch)
print(table(si$grp, si$run))
feat <- data.frame(FeatureID = acc, Gene = ifelse(is.na(gene) | gene == "", acc, gene), stringsAsFactors = FALSE)
# ---- limma, group + run ----
design <- model.matrix(~ 0 + grp + run, data = si); colnames(design) <- sub("^grp", "", colnames(design))
fit <- lmFit(Mi, design)
contr <- makeContrasts(LoFLoF_vs_RefRef = LoFLoF - RefRef_CT, LoFMS_vs_RefRef = LoFMS - RefRef_CT, LoFMS_vs_LoFLoF = LoFMS - LoFLoF, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))
tabs <- list()
for (cn in colnames(contr)) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "P"); tt$FeatureID <- rownames(tt)
  tt <- merge(tt, feat, by = "FeatureID"); tt <- tt[order(tt$P.Value), ]
  tt$Significance <- ifelse(tt$adj.P.Val < FDR & tt$logFC >= LFC, "up", ifelse(tt$adj.P.Val < FDR & tt$logFC <= -LFC, "down", "ns"))
  write.csv(tt, file.path(OUT, paste0("DE_limma_groupRun_", cn, ".csv")), row.names = FALSE); tabs[[cn]] <- tt
  cat(sprintf("%s: adjP<%.2f & |logFC|>=%g: up %d, down %d (adjP<0.05 any: %d)\n", cn, FDR, LFC, sum(tt$Significance=="up"), sum(tt$Significance=="down"), sum(tt$adj.P.Val < FDR)))
}
# no-covariate variant (for the report only)
design0 <- model.matrix(~ 0 + grp, data = si); colnames(design0) <- sub("^grp", "", colnames(design0))
fit0 <- eBayes(contrasts.fit(lmFit(Mi, design0), makeContrasts(LoFLoF - RefRef_CT, LoFMS - RefRef_CT, LoFMS - LoFLoF, levels = design0)))
for (k in 1:3) { tt0 <- topTable(fit0, coef = k, number = Inf); cat(sprintf("  [no run covariate] contrast %d: up %d down %d\n", k, sum(tt0$adj.P.Val<FDR & tt0$logFC>=LFC), sum(tt0$adj.P.Val<FDR & tt0$logFC<= -LFC))) }
# ---- global F-test across groups (batch-adjusted): intercept model, F over the 4 group coefficients ----
designF <- model.matrix(~ grp + run, data = si)
fitF <- eBayes(lmFit(Mi, designF))
ttF <- topTable(fitF, coef = grep("^grp", colnames(designF)), number = Inf, sort.by = "F"); ttF$FeatureID <- rownames(ttF)
ttF <- merge(ttF, feat, by = "FeatureID"); ttF <- ttF[order(ttF$P.Value), ]
write.csv(ttF, file.path(OUT, "DE_limma_globalF_groupRun_all.csv"), row.names = FALSE)
depF <- ttF[ttF$adj.P.Val < FDR, ]; cat("global F adjP<0.05:", nrow(depF), "\n")
top_ids <- head(depF$FeatureID, TOPN)
# ---- plots ----
FAM <- "Arial"
lab_grp <- c(RefRef_CT = "Ref/Ref controls", RefLoF = "Ref/LoF", LoFLoF = "LoF/LoF", LoFMS = "LoF/Mis", RefRef_DD = "Ref/Ref disease controls")
grp_cols <- c(RefRef_CT = "#FF725C", RefLoF = "#EFB118", LoFLoF = "#6CC5B0", LoFMS = "#97BBF5", RefRef_DD = "#FF8AB7")
base_theme <- theme_classic(base_size = 8, base_family = FAM) + theme(plot.title = element_text(size = 8, face = "plain", hjust = 0), legend.position = "none", axis.text = element_text(size = 7, colour = "black"))
titles <- c(LoFLoF_vs_RefRef = "LoF/LoF vs Ref/Ref controls", LoFMS_vs_RefRef = "LoF/Mis vs Ref/Ref controls", LoFMS_vs_LoFLoF = "LoF/Mis vs LoF/LoF")
ups <- c(LoFLoF_vs_RefRef = "LoF/LoF", LoFMS_vs_RefRef = "LoF/Mis", LoFMS_vs_LoFLoF = "LoF/Mis")
sig_cols <- c(up = "#D7301F", down = "#2166AC", ns = "grey70")
ma_list <- list(); vol_list <- list()
for (cn in names(titles)) {
  tt <- tabs[[cn]]; tt$Significance <- factor(tt$Significance, levels = c("ns","down","up")); tt <- tt[order(tt$Significance), ]
  lab_df <- tt[tt$Significance != "ns", ]; lab_df <- head(lab_df[order(lab_df$adj.P.Val), ], 20)
  ma_list[[cn]] <- ggplot(tt, aes(AveExpr, logFC, colour = Significance)) + geom_hline(yintercept = c(-LFC, LFC), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_point(size = 0.7, alpha = 0.7) + scale_colour_manual(values = sig_cols) +
    geom_text_repel(data = lab_df, aes(label = Gene), size = 2, family = FAM, max.overlaps = 30, segment.size = 0.2, colour = "black") +
    labs(title = paste0("MA plot: ", titles[cn]), x = "Average log2 abundance (ACTB-normalized)", y = paste0("log2 fold change (", ups[cn], " vs comparator)")) + base_theme
  vol_list[[cn]] <- ggplot(tt, aes(logFC, -log10(adj.P.Val), colour = Significance)) + geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(FDR), linetype = "dashed", colour = "grey50", linewidth = 0.3) + geom_point(size = 0.7, alpha = 0.7) + scale_colour_manual(values = sig_cols) +
    geom_text_repel(data = lab_df, aes(label = Gene), size = 2, family = FAM, max.overlaps = 30, segment.size = 0.2, colour = "black") +
    labs(title = paste0("Volcano plot: ", titles[cn]), x = paste0("log2 fold change (", ups[cn], " vs comparator)"), y = "-log10 adjusted P") + base_theme
}
# heatmap
mat <- Mi[top_ids, , drop = FALSE]; matz <- t(scale(t(mat))); matz[!is.finite(matz)] <- 0
rownames(matz) <- feat$Gene[match(rownames(matz), feat$FeatureID)]
ann <- data.frame(Group = factor(lab_grp[as.character(si$grp)], levels = lab_grp), row.names = si$SampleID)
ann_cols <- list(Group = setNames(grp_cols, lab_grp))
hm <- pheatmap(matz, annotation_col = ann, annotation_colors = ann_cols, cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = FALSE, show_rownames = TRUE,
               fontsize = 7, fontsize_row = 7, fontfamily = FAM, cellwidth = 5.2, cellheight = 26, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), breaks = seq(-3, 3, length.out = 102),
               border_color = NA, treeheight_row = 15, treeheight_col = 15, silent = TRUE, annotation_names_col = FALSE,
               main = sprintf("Proteins differing among the five groups\n(F-test, adjusted P < 0.05; n = %d)", nrow(depF)))
left <- plot_grid(ma_list[[1]], vol_list[[1]], ma_list[[2]], vol_list[[2]], ma_list[[3]], vol_list[[3]], ncol = 2, labels = c("A","B","C","D","E","F"), label_size = 10, label_fontfamily = FAM)
right <- plot_grid(hm$gtable, labels = "G", label_size = 10, label_fontfamily = FAM)
fig <- plot_grid(left, right, ncol = 2, rel_widths = c(1, 1))
cairo_pdf(file.path(OUT, "FigureS7_proteome_differential_abundance.pdf"), width = 14, height = 9, family = FAM); print(fig); dev.off()
png(file.path(OUT, "FigureS7_proteome_differential_abundance.png"), width = 14, height = 9, units = "in", res = 200, type = "cairo", family = FAM); print(fig); dev.off()
writeLines(c(paste("features tested:", nrow(Mi)), capture.output(table(si$grp, si$run)), sapply(names(tabs), function(cn) sprintf("%s: up %d down %d", cn, sum(tabs[[cn]]$Significance=="up"), sum(tabs[[cn]]$Significance=="down"))), paste("globalF adjP<0.05:", nrow(depF)), capture.output(sessionInfo())), file.path(OUT, "run_log.txt"))
cat("done\n")
