## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
# ============================================================
# RNA-seq analysis (RESUBMISSION)
#   Derived from reference/Analysis_7_enriched_exclude_SS234578_20260317.R
#   (same structure; changes listed below)
#
#   - featureCounts -> DESeq2
#   - Sample set = manuscript contrast: LoF/LoF (n = 7) vs Ref/Ref (n = 6)
#     (carriers and the LoF/Mis specimen are NOT included; counts file holds
#      only these 13 samples)
#   - Designs:
#       PRIMARY     : ~ batch + sex + status   (Reviewer 2: sex as covariate)
#       SUBMITTED   : ~ batch + status         (design used at submission; sensitivity)
#     Age cannot be adjusted: all patients <= 14 y, all controls >= 32 y
#     (complete confounding) -> reported in Limitations, not modelled.
#   - Modes:
#       ENRICHED : capture-panel genes only; ranked list = captured genes
#                  -> gene sets evaluated AGAINST THE CAPTURED-GENE BACKGROUND
#                  (Reviewer 2: capture panel as background). PRIMARY MODE.
#       ALL      : whole transcriptome (off-target genes included); kept only
#                  as a reference run showing the capture-design artefact.
#   - VST / PCA (raw VST + batch-regressed residual PCA)
#   - Volcano / MA with sex-chromosome (chrX/chrY) and extreme-LFC labels
#   - GSEA (GO BP/MF/CC + Reactome) with NES and FDR in every table, fixed seed
#   - STRICT GSEA (padj < 0.05)
#   - PANEL-AWARE ORA for ENRICHED mode (universe = captured genes)
#   - SLC25A13 focus table (Fig. 3A/B), capture statistics (Methods numbers),
#     cross-run summary and comparison of previously reported terms
#     (ALL vs ENRICHED; Table S12 replacement material)
#
# CHANGES vs reference script
#   - sample metadata (genotype, age, sex, batch) read from data/info_20261103.xlsx
#   - two designs per mode (PRIMARY with sex, SUBMITTED without sex)
#   - set.seed + seed = TRUE in gseGO/gsePathway (reproducible NES/FDR)
#   - figures exported as PDF only (no PNG)
#   - new summary/ directory with Methods numbers and cross-run tables
#   - removed: renv snapshot, zip bundle of PNGs (bundle keeps script+sessionInfo+hashes)
#
# Outputs are organized under:
#   OUTDIR/
#     master_logs/
#     summary/                       (cross-run tables; Methods numbers)
#     <MODE>/<DESIGN>/
#       logs/ figures/ tables/ rdata/ submission_bundle/
#       tables/{metadata, DESeq2, matrices, PCA, Volcano_MA, GSEA, SLC25A13}
# ============================================================

# ============================================================
# 0) GLOBAL SETTINGS
# ============================================================
PROJECT_ROOT <- file.path(REPO, "06_RNAseq")   # expects data/ (count matrix, sample info, panel list; available on request) under this folder
OUTDIR <- file.path(PROJECT_ROOT, "results_resubmission")

MASTER_LOGDIR <- file.path(OUTDIR, "master_logs")
SUMMARY_DIR   <- file.path(OUTDIR, "summary")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MASTER_LOGDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUMMARY_DIR, showWarnings = FALSE, recursive = TRUE)

setwd(OUTDIR)
message("[OUTDIR] ", OUTDIR)

# -----------------------------
# Inputs
# -----------------------------
count_file <- file.path(PROJECT_ROOT, "data", "counts_featureCounts_20261103.txt")
info_xlsx  <- file.path(PROJECT_ROOT, "data", "info_20261103.xlsx")
panel_csv  <- file.path(PROJECT_ROOT, "data", "panel_symbols_from_xlsx.csv")

script_self_path <- file.path(PROJECT_ROOT, "script", "Analysis_resubmission_RNAseq_20261103.R")

# -----------------------------
# Analysis settings
# -----------------------------
RUN_MODES   <- c("ENRICHED", "ALL")             # ENRICHED = primary
DESIGNS     <- list(
  PRIMARY   = ~ batch + sex + status,
  SUBMITTED = ~ batch + status
)
PADJ_CUT    <- 0.05
LFC_CUT     <- 1
GSEA_SEED   <- 20261103
GSEA_MIN    <- 10
GSEA_MAX    <- 500
RUN_PANEL_ORA <- TRUE
SHOW_N      <- 15
EXTREME_LFC <- 5     # |log2FC| above which genes are labelled on Volcano/MA (Reviewer 2)
# Figure titles: the primary run (ENRICHED_PRIMARY; Supplementary Figs. S9-S11) is untagged; the other runs keep their run tag.
fig_title <- function(run_tag, txt) if (identical(run_tag, "ENRICHED_PRIMARY")) txt else paste0("[", run_tag, "] ", txt)

# GO/Reactome terms reported in the submitted manuscript (Fig. 3C/D, Table S12,
# Results/Discussion). They are re-evaluated in every run.
PREVIOUSLY_REPORTED_TERMS <- c(
  "GO:0002181", # cytoplasmic translation (Fig. 3C)
  "GO:0006119", # oxidative phosphorylation (Fig. 3C)
  "GO:0008380", # RNA splicing (Fig. 3C)
  "GO:0060271", # cilium assembly (Fig. 3C)
  "GO:0098803", # respiratory chain complex (mitochondrial function)
  "GO:0005509", # calcium ion binding
  "GO:0005516", # calmodulin binding
  "GO:0022832", # voltage-gated channel activity
  "GO:0034220", # monoatomic ion transmembrane transport
  "GO:0015075"  # monoatomic ion transmembrane transporter activity
)

# ============================================================
# 0A) MASTER LOGGING
# ============================================================
master_log_file <- file.path(
  MASTER_LOGDIR,
  paste0("master_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)
cat("[MASTER LOG START] ", format(Sys.time()), "\n", file = master_log_file, append = FALSE)

master_log <- function(...) {
  txt <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste(..., collapse = ""))
  message(txt)
  cat(txt, "\n", file = master_log_file, append = TRUE)
}

# ============================================================
# 0B) Libraries
# ============================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(matrixStats)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggrepel)
  library(patchwork)
  library(grid)

  library(clusterProfiler)
  library(enrichplot)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(AnnotationDbi)

  library(ggrastr)
  library(scales)
  library(openxlsx)
})

select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange
rename  <- dplyr::rename

# ============================================================
# 0C) Generic helpers
# ============================================================
write_csv_safe <- function(df, path, row.names = FALSE){
  utils::write.csv(df, file = path, row.names = row.names, quote = TRUE)
}

save_rds_safe <- function(obj, path){
  saveRDS(obj, file = path)
}

make_inputs_manifest <- function(files){
  files <- unique(files)
  info <- file.info(files)
  out <- data.frame(
    file = files,
    exists = file.exists(files),
    size = ifelse(file.exists(files), info$size, NA),
    mtime = ifelse(file.exists(files), as.character(info$mtime), NA),
    stringsAsFactors = FALSE
  )
  if (requireNamespace("digest", quietly = TRUE)) {
    out$sha256 <- ifelse(
      out$exists,
      vapply(out$file, function(f) digest::digest(f, algo = "sha256", file = TRUE), FUN.VALUE = character(1)),
      NA
    )
  } else {
    out$sha256 <- NA
  }
  out
}

make_run_dirs <- function(base_outdir, mode_tag, design_tag) {
  run_root <- file.path(base_outdir, mode_tag, design_tag)
  dirs <- list(
    ROOT      = run_root,
    LOGS      = file.path(run_root, "logs"),
    FIG       = file.path(run_root, "figures"),
    TABLES    = file.path(run_root, "tables"),
    RDATA     = file.path(run_root, "rdata"),
    BUNDLE    = file.path(run_root, "submission_bundle"),

    META      = file.path(run_root, "tables", "metadata"),
    DESEQ     = file.path(run_root, "tables", "DESeq2"),
    MAT       = file.path(run_root, "tables", "matrices"),
    PCA       = file.path(run_root, "tables", "PCA"),
    VOLMA     = file.path(run_root, "tables", "Volcano_MA"),
    GSEA      = file.path(run_root, "tables", "GSEA"),
    SLC       = file.path(run_root, "tables", "SLC25A13")
  )
  invisible(lapply(dirs, function(d) dir.create(d, showWarnings = FALSE, recursive = TRUE)))
  dirs
}

mode_log_factory <- function(log_file) {
  force(log_file)
  function(...) {
    txt <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste(..., collapse = ""))
    message(txt)
    cat(txt, "\n", file = log_file, append = TRUE)
  }
}

safe_mapIds <- function(keys, column, keytype, multiVals = "first") {
  out <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = keys,
    column = column,
    keytype = keytype,
    multiVals = multiVals
  )
  unname(out)
}

read_enriched_panel <- function(panel_csv) {
  stopifnot(file.exists(panel_csv))
  x <- utils::read.csv(panel_csv, stringsAsFactors = FALSE, check.names = FALSE)
  gene_col <- intersect(c("panel_symbol", "Gene", "gene", "SYMBOL", "Symbol", "symbol"), colnames(x))
  if (length(gene_col) == 0) stop("Could not find a gene-symbol column in ", panel_csv)
  gene_col <- gene_col[1]
  panel_symbols <- unique(trimws(as.character(x[[gene_col]])))
  panel_symbols <- panel_symbols[!is.na(panel_symbols) & panel_symbols != ""]

  panel_ensembl <- safe_mapIds(panel_symbols, "ENSEMBL", "SYMBOL", "first")
  panel_ensembl <- unique(panel_ensembl[!is.na(panel_ensembl) & panel_ensembl != ""])
  panel_entrez <- safe_mapIds(panel_symbols, "ENTREZID", "SYMBOL", "first")
  panel_entrez <- unique(panel_entrez[!is.na(panel_entrez) & panel_entrez != ""])

  list(panel_symbols = panel_symbols, panel_ensembl = panel_ensembl,
       panel_entrez = panel_entrez, raw_table = x, gene_col = gene_col)
}

# Sample info workbook: rows = Genotype group / Age / Sex / batch, columns = samples
read_sample_info <- function(info_xlsx) {
  stopifnot(file.exists(info_xlsx))
  x <- openxlsx::read.xlsx(info_xlsx, sheet = 1, colNames = TRUE)
  key <- trimws(as.character(x[[1]]))
  samples <- colnames(x)[-1]
  get_row <- function(pattern) {
    i <- which(grepl(pattern, key, ignore.case = TRUE))
    if (length(i) == 0) stop("Row not found in info xlsx: ", pattern)
    trimws(as.character(unlist(x[i[1], -1])))
  }
  genotype_group <- get_row("^Genotype")
  age_txt <- get_row("^Age")
  sex_txt <- tolower(get_row("^Sex"))
  batch   <- get_row("^batch")

  parse_age_years <- function(s) {
    s <- tolower(s)
    y <- str_match(s, "([0-9]+)\\s*y")[, 2]
    m <- str_match(s, "([0-9]+)\\s*m")[, 2]
    y <- ifelse(is.na(y), 0, as.numeric(y))
    m <- ifelse(is.na(m), 0, as.numeric(m))
    round(y + m / 12, 2)
  }
  sex <- ifelse(sex_txt %in% c("woman", "girl", "female", "f"), "F",
         ifelse(sex_txt %in% c("man", "boy", "male", "m"), "M", NA))
  if (any(is.na(sex))) stop("Unrecognised sex label(s): ", paste(sex_txt[is.na(sex)], collapse = ", "))

  status <- ifelse(genotype_group == "Ref/Ref", "control",
            ifelse(genotype_group == "LoF/LoF", "patient", NA))
  if (any(is.na(status))) stop("Only Ref/Ref and LoF/LoF are expected in the resubmission contrast; found: ",
                               paste(unique(genotype_group), collapse = ", "))

  data.frame(
    sample = samples,
    genotype_group = genotype_group,
    status = factor(status, levels = c("control", "patient")),
    age_label = age_txt,
    age_years = parse_age_years(age_txt),
    sex = factor(sex, levels = c("F", "M")),
    batch = factor(batch),
    row.names = samples,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 1) Read featureCounts, sample info, panel
# ============================================================
stopifnot(file.exists(count_file))
stopifnot(file.exists(info_xlsx))
stopifnot(file.exists(panel_csv))

cts_raw <- read.table(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)
cts_raw <- subset(cts_raw, Geneid != "")
cts_raw$Geneid_clean <- sub("\\.\\d+$", "", cts_raw$Geneid)
cts_raw$Chr_first <- sub(";.*$", "", cts_raw$Chr)

symbol_map <- safe_mapIds(cts_raw$Geneid_clean, "SYMBOL", "ENSEMBL", "first")

coldata_full <- read_sample_info(info_xlsx)
count_cols <- coldata_full$sample
missing_cols <- setdiff(count_cols, colnames(cts_raw))
if (length(missing_cols) > 0) {
  stop("These count columns are missing in count file: ", paste(missing_cols, collapse = ", "))
}
extra_cols <- setdiff(colnames(cts_raw), c("Geneid", "Chr", "Start", "End", "Strand", "Length",
                                           "Geneid_clean", "Chr_first", count_cols))
if (length(extra_cols) > 0) {
  master_log("[WARN] count file has columns not in sample info (ignored): ", paste(extra_cols, collapse = ", "))
}

cts_all_full <- data.frame(
  Geneid = cts_raw$Geneid,
  ENSEMBL = cts_raw$Geneid_clean,
  SYMBOL = symbol_map,
  Chr = cts_raw$Chr_first,
  Length = cts_raw$Length,
  cts_raw[, count_cols, drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(cts_all_full) <- cts_all_full$ENSEMBL

panel_info <- read_enriched_panel(panel_csv)

cts_all_full$IN_PANEL <- (
  (!is.na(cts_all_full$SYMBOL)  & cts_all_full$SYMBOL %in% panel_info$panel_symbols) |
    (cts_all_full$ENSEMBL %in% panel_info$panel_ensembl)
)
cts_all_panel <- cts_all_full[cts_all_full$IN_PANEL, , drop = FALSE]

cts_all_mode <- as.matrix(cts_all_full[, count_cols, drop = FALSE]); mode(cts_all_mode) <- "integer"
cts_enriched_mode <- as.matrix(cts_all_panel[, count_cols, drop = FALSE]); mode(cts_enriched_mode) <- "integer"

# -----------------------------
# Capture statistics (Methods numbers; Reviewer 2 "state how many genes were quantified")
# -----------------------------
on_target_frac <- colSums(cts_enriched_mode) / colSums(cts_all_mode)
capture_stats <- data.frame(
  item = c("panel_symbols_in_list",
           "panel_symbols_mapped_to_ensembl",
           "panel_symbols_mapped_to_entrez",
           "featureCounts_rows_total",
           "featureCounts_rows_matched_to_panel",
           "panel_rows_with_any_count",
           "panel_rows_with_total_count_ge10",
           "all_rows_with_any_count",
           "all_rows_with_total_count_ge10",
           "on_target_read_fraction_median",
           "on_target_read_fraction_min",
           "on_target_read_fraction_max",
           "n_samples", "n_patient_LoF_LoF", "n_control_Ref_Ref"),
  value = c(length(panel_info$panel_symbols),
            length(panel_info$panel_ensembl),
            length(panel_info$panel_entrez),
            nrow(cts_all_full),
            sum(cts_all_full$IN_PANEL),
            sum(rowSums(cts_enriched_mode) > 0),
            sum(rowSums(cts_enriched_mode) >= 10),
            sum(rowSums(cts_all_mode) > 0),
            sum(rowSums(cts_all_mode) >= 10),
            round(median(on_target_frac), 4),
            round(min(on_target_frac), 4),
            round(max(on_target_frac), 4),
            nrow(coldata_full),
            sum(coldata_full$status == "patient"),
            sum(coldata_full$status == "control")),
  stringsAsFactors = FALSE
)
write_csv_safe(capture_stats, file.path(SUMMARY_DIR, "capture_statistics.csv"))
write_csv_safe(data.frame(sample = names(on_target_frac), on_target_read_fraction = round(on_target_frac, 4)),
               file.path(SUMMARY_DIR, "on_target_fraction_per_sample.csv"))
write_csv_safe(data.frame(panel_symbol = panel_info$panel_symbols,
                          matched_in_counts = panel_info$panel_symbols %in% cts_all_panel$SYMBOL),
               file.path(SUMMARY_DIR, "panel_symbols_matching_audit.csv"))

# Participant table (Reviewer 2: ages, sexes; same IDs as the proteomic cohort)
participants <- coldata_full %>%
  mutate(age_group = ifelse(status == "patient", "paediatric", "adult")) %>%
  select(sample, genotype_group, status, age_label, age_years, age_group, sex, batch)
write_csv_safe(participants, file.path(SUMMARY_DIR, "participants_RNAseq.csv"))

master_log("[PANEL] n(panel symbols) = ", length(panel_info$panel_symbols),
           "; matched rows = ", sum(cts_all_full$IN_PANEL),
           "; on-target read fraction (median) = ", round(median(on_target_frac), 3))
master_log("[SAMPLES] ", nrow(coldata_full), " samples: ",
           paste(sprintf("%s(%s,%s,%s,b%s)", coldata_full$sample, coldata_full$status,
                         coldata_full$age_label, coldata_full$sex, coldata_full$batch), collapse = " "))
master_log("[DESIGN NOTE] age is completely confounded with status (patients <= 14 y, controls >= 32 y) and is not modelled")

analysis_inputs <- list(
  ENRICHED = list(cts = cts_enriched_mode, cts_all = cts_all_panel, coldata = coldata_full, panel_info = panel_info),
  ALL      = list(cts = cts_all_mode,      cts_all = cts_all_full,  coldata = coldata_full, panel_info = NULL)
)

# ============================================================
# 2) Main pipeline runner
# ============================================================
run_one_pipeline <- function(mode_tag, design_tag, design_formula,
                             cts_use, cts_all_use, coldata_use, dirs,
                             panel_info_use = NULL) {

  log_file <- file.path(dirs$LOGS, paste0("run_", mode_tag, "_", design_tag, "_",
                                          format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  mode_log <- mode_log_factory(log_file)
  run_tag <- paste0(mode_tag, "_", design_tag)

  mode_log("============================================================")
  mode_log("[RUN] mode = ", mode_tag, " | design = ", design_tag, " : ", deparse(design_formula))
  mode_log("============================================================")

  DIR_FIG    <- dirs$FIG
  DIR_RDATA  <- dirs$RDATA
  DIR_BUNDLE <- dirs$BUNDLE
  DIR_META   <- dirs$META
  DIR_DESEQ  <- dirs$DESEQ
  DIR_MAT    <- dirs$MAT
  DIR_PCA    <- dirs$PCA
  DIR_VOLMA  <- dirs$VOLMA
  DIR_GSEA   <- dirs$GSEA
  DIR_SLC    <- dirs$SLC

  # PDF only
  save_pdf <- function(p, prefix, w = 6, h = 6, units = "in") {
    out_pdf <- file.path(DIR_FIG, paste0(basename(prefix), ".pdf"))
    ggplot2::ggsave(out_pdf, p, width = w, height = h, units = units)
  }

  export_deseq2_results <- function(res, res_shr = NULL, model_tag, contrast_tag,
                                    padj_cut = PADJ_CUT, lfc_cut = LFC_CUT){
    df <- as.data.frame(res)
    df$ENSEMBL <- rownames(df)
    df$SYMBOL <- cts_all_use[df$ENSEMBL, "SYMBOL"]
    df$Chr <- cts_all_use[df$ENSEMBL, "Chr"]
    df <- df[, c("ENSEMBL", "SYMBOL", "Chr", setdiff(colnames(df), c("ENSEMBL", "SYMBOL", "Chr")))]
    df <- df[order(df$pvalue), , drop = FALSE]

    write_csv_safe(df, file.path(DIR_DESEQ, paste0("DESeq2_results_", model_tag, "_", contrast_tag, "_full.csv")))

    df_sig <- df[is.finite(df$padj) & df$padj < padj_cut &
                   is.finite(df$log2FoldChange) & abs(df$log2FoldChange) >= lfc_cut, , drop = FALSE]
    write_csv_safe(df_sig, file.path(DIR_DESEQ, paste0("DESeq2_results_", model_tag, "_", contrast_tag,
                                                      "_sig_padj", padj_cut, "_lfc", lfc_cut, ".csv")))

    if (!is.null(res_shr)) {
      df2 <- as.data.frame(res_shr)
      df2$ENSEMBL <- rownames(df2)
      df2$SYMBOL <- cts_all_use[df2$ENSEMBL, "SYMBOL"]
      df2$Chr <- cts_all_use[df2$ENSEMBL, "Chr"]
      df2 <- df2[, c("ENSEMBL", "SYMBOL", "Chr", setdiff(colnames(df2), c("ENSEMBL", "SYMBOL", "Chr")))]
      write_csv_safe(df2, file.path(DIR_DESEQ, paste0("DESeq2_results_", model_tag, "_", contrast_tag, "_shrunk.csv")))
    }

    capture.output(summary(res), file = file.path(DIR_DESEQ, paste0("DESeq2_summary_", model_tag, "_", contrast_tag, ".txt")))
    invisible(list(full = df, sig = df_sig))
  }

  is_gsea_nonempty <- function(gse_obj) {
    if (is.null(gse_obj)) return(FALSE)
    df <- tryCatch(as.data.frame(gse_obj), error = function(e) NULL)
    if (is.null(df)) return(FALSE)
    nrow(df) > 0
  }

  get_plot_xrange <- function(p) {
    b <- tryCatch(ggplot_build(p), error = function(e) NULL)
    if (is.null(b)) return(NULL)
    xs <- unlist(lapply(b$data, function(d) {
      out <- NULL
      if (!is.null(d$x))    out <- c(out, d$x)
      if (!is.null(d$xmin)) out <- c(out, d$xmin)
      if (!is.null(d$xmax)) out <- c(out, d$xmax)
      out
    }), use.names = FALSE)
    xs <- xs[is.finite(xs)]
    if (length(xs) < 2) return(NULL)
    r <- range(xs); w <- diff(r)
    if (!is.finite(w) || w <= 0) return(NULL)
    r
  }
  tighten_plot_x <- function(p, pad = 0.05) {
    r <- get_plot_xrange(p); if (is.null(r)) return(p)
    w <- diff(r); p + coord_cartesian(xlim = c(r[1] - pad * w, r[2] + pad * w), clip = "off")
  }
  apply_compact_margins <- function(p, left_margin_pts = 14, right_margin_pts = 6) {
    p + theme(plot.margin = margin(t = 6, r = right_margin_pts, b = 6, l = left_margin_pts))
  }
  wrap_y_labels <- function(p, wrap_width = 32) {
    p + scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = wrap_width))
  }

  gsea_table_df <- function(gse_obj, ont_tag) {
    if (!is_gsea_nonempty(gse_obj)) return(NULL)
    df <- as.data.frame(gse_obj)
    df$Ontology <- ont_tag
    df$mode <- mode_tag
    df$design <- design_tag
    df[, c("Ontology", "ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue",
           "rank", "leading_edge", "core_enrichment", "mode", "design")]
  }
  save_gsea_table <- function(gse_obj, tag_prefix, ont_tag) {
    df <- gsea_table_df(gse_obj, ont_tag)
    if (is.null(df)) return(invisible(NULL))
    df <- df[order(df$p.adjust, df$pvalue), , drop = FALSE]
    write_csv_safe(df, file.path(DIR_GSEA, paste0(tag_prefix, "_results.csv")))
    invisible(df)
  }

  save_gsea_dot_ridge <- function(gse_obj, tag_prefix, showN = SHOW_N, wrap_width = 32, pad = 0.05,
                                  left_margin_pts = 14, right_margin_pts = 6, w = 9, h = 6, subtitle = NULL) {
    if (!is_gsea_nonempty(gse_obj)) {
      mode_log("[GSEA] ", tag_prefix, " is NULL/empty -> skipped")
      return(invisible(NULL))
    }
    p_dot <- tryCatch({
      p <- enrichplot::dotplot(gse_obj, showCategory = showN) +
        ggtitle(paste0(tag_prefix, " dotplot"), subtitle = subtitle) +
        theme_bw() + theme(plot.title = element_text(face = "bold")) + coord_cartesian(clip = "off")
      p <- wrap_y_labels(p, wrap_width); p <- tighten_plot_x(p, pad)
      apply_compact_margins(p, left_margin_pts, right_margin_pts)
    }, error = function(e) { mode_log("[GSEA] dotplot failed (", tag_prefix, "): ", conditionMessage(e)); NULL })
    if (!is.null(p_dot)) save_pdf(p_dot, paste0(tag_prefix, "_dotplot"), w = w, h = h)

    p_rid <- tryCatch({
      p <- enrichplot::ridgeplot(gse_obj, showCategory = showN) +
        ggtitle(paste0(tag_prefix, " ridgeplot"), subtitle = subtitle) +
        theme_bw() + theme(plot.title = element_text(face = "bold")) + coord_cartesian(clip = "off")
      p <- wrap_y_labels(p, wrap_width); p <- tighten_plot_x(p, pad)
      apply_compact_margins(p, left_margin_pts, right_margin_pts)
    }, error = function(e) { mode_log("[GSEA] ridgeplot failed (", tag_prefix, "): ", conditionMessage(e)); NULL })
    if (!is.null(p_rid)) save_pdf(p_rid, paste0(tag_prefix, "_ridgeplot"), w = w, h = h)
    invisible(list(dot = p_dot, ridge = p_rid))
  }

  # NES bar chart of the top terms by FDR (Fig. 3C/D replacement; NES and FDR shown)
  save_nes_barplot <- function(df_all, tag_prefix, topN = SHOW_N, w = 8, h = 6) {
    if (is.null(df_all) || nrow(df_all) == 0) return(invisible(NULL))
    d <- df_all[order(df_all$p.adjust, df_all$pvalue), , drop = FALSE]
    d <- head(d, topN)
    d$label <- sprintf("%s (%s, n=%d)", d$Description, d$Ontology, d$setSize)
    d$label <- factor(d$label, levels = rev(d$label))
    d$sig <- ifelse(d$p.adjust < PADJ_CUT, "FDR < 0.05", "FDR >= 0.05")
    p <- ggplot(d, aes(x = NES, y = label, fill = sig)) +
      geom_col(width = 0.7, colour = "grey30") +
      geom_text(aes(label = sprintf("FDR %.2g", p.adjust),
                    x = ifelse(NES >= 0, NES + 0.05, NES - 0.05),
                    hjust = ifelse(NES >= 0, 0, 1)), size = 2.8) +
      geom_vline(xintercept = 0, colour = "grey40") +
      scale_fill_manual(values = c("FDR < 0.05" = "#4269D0", "FDR >= 0.05" = "white"), name = NULL, drop = FALSE) +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, 45)) +
      scale_x_continuous(expand = expansion(mult = 0.25)) +
      labs(title = paste0("[", run_tag, "] top ", topN, " gene sets by FDR"),
           subtitle = ifelse(mode_tag == "ENRICHED",
                             "Background = captured (panel) genes; patient vs control",
                             "Whole-transcriptome ranking (reference run, capture-design artefact)"),
           x = "Normalized enrichment score (NES)", y = NULL) +
      theme_bw(base_size = 10) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
    save_pdf(p, tag_prefix, w = w, h = h)
    invisible(p)
  }

  make_geneList_from_deseq <- function(res) {
    df <- as.data.frame(res)
    df$ENSEMBL <- rownames(df)
    df <- df %>% dplyr::filter(is.finite(stat))
    entrez <- safe_mapIds(df$ENSEMBL, "ENTREZID", "ENSEMBL", "first")
    rank_df <- data.frame(ENTREZID = entrez, stat = df$stat, stringsAsFactors = FALSE)
    rank_df <- rank_df[!is.na(rank_df$ENTREZID), , drop = FALSE]
    rank_df <- rank_df[order(abs(rank_df$stat), decreasing = TRUE), , drop = FALSE]
    rank_df <- rank_df[!duplicated(rank_df$ENTREZID), , drop = FALSE]
    geneList <- rank_df$stat
    names(geneList) <- rank_df$ENTREZID
    sort(geneList, decreasing = TRUE)
  }

  make_sig_entrez_from_deseq <- function(res, padj_cut = PADJ_CUT, lfc_cut = LFC_CUT) {
    df <- as.data.frame(res)
    df$ENSEMBL <- rownames(df)
    df <- df[is.finite(df$padj) & df$padj < padj_cut &
               is.finite(df$log2FoldChange) & abs(df$log2FoldChange) >= lfc_cut, , drop = FALSE]
    if (nrow(df) == 0) return(character(0))
    entrez <- safe_mapIds(df$ENSEMBL, "ENTREZID", "ENSEMBL", "first")
    sig_entrez <- unique(entrez)
    sig_entrez[!is.na(sig_entrez) & sig_entrez != ""]
  }

  safe_gseGO <- function(geneList, ont, tag, pvalueCutoff = 1) {
    if (length(geneList) < 10) {
      mode_log("[GSEA] geneList too short for ", tag, " / ", ont, " -> skipped"); return(NULL)
    }
    set.seed(GSEA_SEED)
    tryCatch({
      clusterProfiler::gseGO(geneList, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
                             pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff,
                             minGSSize = GSEA_MIN, maxGSSize = GSEA_MAX, verbose = FALSE, seed = TRUE)
    }, error = function(e) { mode_log("[GSEA] gseGO failed (", tag, " / ", ont, "): ", conditionMessage(e)); NULL })
  }

  safe_gsePathway <- function(geneList, tag, pvalueCutoff = 1) {
    if (length(geneList) < 10) {
      mode_log("[GSEA] geneList too short for Reactome ", tag, " -> skipped"); return(NULL)
    }
    set.seed(GSEA_SEED)
    tryCatch({
      ReactomePA::gsePathway(geneList, organism = "human", pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff,
                             minGSSize = GSEA_MIN, maxGSSize = GSEA_MAX, verbose = FALSE, seed = TRUE)
    }, error = function(e) { mode_log("[GSEA] gsePathway failed (", tag, "): ", conditionMessage(e)); NULL })
  }

  # STRICT = same objects filtered to padj < 0.05 (identical statistics; no re-run)
  strict_subset <- function(gse_obj, padj_cut = PADJ_CUT) {
    if (!is_gsea_nonempty(gse_obj)) return(NULL)
    df0 <- as.data.frame(gse_obj)
    keep_id <- df0$ID[is.finite(df0$p.adjust) & df0$p.adjust < padj_cut]
    if (length(keep_id) == 0) return(NULL)
    g2 <- gse_obj
    g2@result <- g2@result[g2@result$ID %in% keep_id, , drop = FALSE]
    g2
  }

  safe_enrichGO_panelORA <- function(sig_entrez, universe_entrez, ont, tag) {
    if (length(sig_entrez) == 0 || length(universe_entrez) == 0) return(NULL)
    tryCatch({
      enrichGO(gene = sig_entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
               ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
               minGSSize = GSEA_MIN, maxGSSize = GSEA_MAX, readable = TRUE)
    }, error = function(e) { mode_log("[ORA] enrichGO failed (", tag, " / ", ont, "): ", conditionMessage(e)); NULL })
  }
  safe_enrichPathway_panelORA <- function(sig_entrez, universe_entrez, tag) {
    if (length(sig_entrez) == 0 || length(universe_entrez) == 0) return(NULL)
    tryCatch({
      ReactomePA::enrichPathway(gene = sig_entrez, universe = universe_entrez, organism = "human",
                                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                                minGSSize = GSEA_MIN, maxGSSize = GSEA_MAX, readable = TRUE)
    }, error = function(e) { mode_log("[ORA] enrichPathway failed (", tag, "): ", conditionMessage(e)); NULL })
  }
  save_ora_table <- function(ora_obj, tag_prefix, ont_tag) {
    if (!is_gsea_nonempty(ora_obj)) return(invisible(NULL))
    df <- as.data.frame(ora_obj); df$Ontology <- ont_tag; df$mode <- mode_tag; df$design <- design_tag
    df <- df[order(df$p.adjust, df$pvalue), , drop = FALSE]
    write_csv_safe(df, file.path(DIR_GSEA, paste0(tag_prefix, "_results.csv")))
    invisible(df)
  }

  # -----------------------------
  # Metadata export
  # -----------------------------
  coldata_out <- as.data.frame(coldata_use)
  write_csv_safe(coldata_out, file.path(DIR_META, "coldata.csv"))
  write_csv_safe(data.frame(mode = mode_tag, design = design_tag, formula = deparse(design_formula),
                            n_samples = nrow(coldata_use),
                            n_patient = sum(coldata_use$status == "patient"),
                            n_control = sum(coldata_use$status == "control"),
                            n_genes_in_count_matrix = nrow(cts_use),
                            samples = paste(coldata_use$sample, collapse = ";"),
                            note_age = "age completely confounded with status; not modelled",
                            stringsAsFactors = FALSE),
                 file.path(DIR_META, "run_manifest.csv"))
  write_csv_safe(as.data.frame(cts_use), file.path(DIR_MAT, "counts_raw_gene_x_sample.csv"), row.names = TRUE)
  write_csv_safe(as.data.frame(cts_all_use[, c("Geneid", "ENSEMBL", "SYMBOL", "Chr", "Length", "IN_PANEL")]),
                 file.path(DIR_META, "gene_annotation.csv"))
  if (!is.null(panel_info_use)) {
    write_csv_safe(data.frame(panel_symbol = panel_info_use$panel_symbols), file.path(DIR_META, "panel_symbols.csv"))
  }
  write_csv_safe(make_inputs_manifest(c(count_file, info_xlsx, panel_csv)), file.path(DIR_META, "inputs_manifest.csv"))

  # -----------------------------
  # DESeq2
  # -----------------------------
  mode_log("[STEP] DESeq2 start")
  dds <- DESeqDataSetFromMatrix(countData = cts_use, colData = coldata_use, design = design_formula)
  dds$status <- relevel(dds$status, ref = "control")
  dds <- DESeq(dds, quiet = TRUE)

  res_pvc <- results(dds, contrast = c("status", "patient", "control"))
  # Shrunken LFC (DESeq2 normal estimator, as stated in the Methods and used for the deposited results).
  # Fixed to type = "normal" so that the result does not depend on whether apeglm happens to be installed.
  res_pvc_shr <- tryCatch({
    lfcShrink(dds, contrast = c("status", "patient", "control"), type = "normal", quiet = TRUE)
  }, error = function(e) { mode_log("[DESeq2] lfcShrink failed, using unshrunken LFC: ", conditionMessage(e)); res_pvc })
  shr_type <- "normal"
  mode_log("[DESeq2] LFC shrinkage type for plots: ", shr_type)

  de <- export_deseq2_results(res_pvc, res_shr = res_pvc_shr, model_tag = design_tag, contrast_tag = "patient_vs_control")
  save_rds_safe(dds, file.path(DIR_RDATA, "dds.rds"))

  # Extreme UNSHRUNKEN log2FC among DEGs (Reviewer 2: "log2 fold changes of -20 to -30 ... suggest
  # sex-chromosome genes"). Such values arise when one group has zero counts; per-group mean
  # normalized counts and chromosome are listed so the reader can see the origin of each value.
  nc <- counts(dds, normalized = TRUE)
  grp_mean <- function(ids, grp) rowMeans(nc[ids, coldata_use$status == grp, drop = FALSE])
  ext <- de$sig[abs(de$sig$log2FoldChange) >= EXTREME_LFC, , drop = FALSE]
  if (nrow(ext) > 0) {
    ext$mean_norm_count_control <- round(grp_mean(ext$ENSEMBL, "control"), 2)
    ext$mean_norm_count_patient <- round(grp_mean(ext$ENSEMBL, "patient"), 2)
    ext$n_control_nonzero <- rowSums(cts_use[ext$ENSEMBL, coldata_use$status == "control", drop = FALSE] > 0)
    ext$n_patient_nonzero <- rowSums(cts_use[ext$ENSEMBL, coldata_use$status == "patient", drop = FALSE] > 0)
    ext$sex_chromosome <- ext$Chr %in% c("chrX", "chrY")
    ext$shrunken_log2FC <- round(as.data.frame(res_pvc_shr)[ext$ENSEMBL, "log2FoldChange"], 3)
    ext <- ext[order(-abs(ext$log2FoldChange)), , drop = FALSE]
  }
  write_csv_safe(ext, file.path(DIR_DESEQ, paste0("DEG_extreme_unshrunken_LFC_ge", EXTREME_LFC, ".csv")))
  # DEG list with the |LFC| threshold applied to the SHRUNKEN estimate (robust to sparse low-count genes)
  shr_df <- as.data.frame(res_pvc_shr); shr_df$ENSEMBL <- rownames(shr_df)
  de_shr <- de$full
  de_shr$shrunken_log2FC <- shr_df$log2FoldChange[match(de_shr$ENSEMBL, shr_df$ENSEMBL)]
  de_shr <- de_shr[is.finite(de_shr$padj) & de_shr$padj < PADJ_CUT &
                     is.finite(de_shr$shrunken_log2FC) & abs(de_shr$shrunken_log2FC) >= LFC_CUT, , drop = FALSE]
  write_csv_safe(de_shr, file.path(DIR_DESEQ, paste0("DESeq2_results_", design_tag, "_patient_vs_control_sig_padj",
                                                     PADJ_CUT, "_shrunkenLFC", LFC_CUT, ".csv")))
  mode_log("[DESeq2] DEGs with padj < ", PADJ_CUT, " and SHRUNKEN |LFC| >= ", LFC_CUT, ": ", nrow(de_shr))
  mode_log("[DESeq2] DEGs with unshrunken |LFC| >= ", EXTREME_LFC, ": ", nrow(ext),
           " (on chrX/chrY: ", sum(ext$sex_chromosome), "; zero in one group: ",
           sum(ext$n_control_nonzero == 0 | ext$n_patient_nonzero == 0), ")")
  mode_log("[DESeq2] DEGs (padj < ", PADJ_CUT, ", |LFC| >= ", LFC_CUT, "): ", nrow(de$sig))

  # SLC25A13 focus (Fig. 3A/B). Default results apply Cook's outlier filtering (padj may be NA
  # when no cell of the design has >= 7 replicates); the unfiltered row is given alongside.
  res_pvc_nocook <- results(dds, contrast = c("status", "patient", "control"), cooksCutoff = FALSE)
  slc_id <- rownames(cts_all_use)[which(cts_all_use$SYMBOL == "SLC25A13")]
  if (length(slc_id) > 0) {
    slc_tab <- rbind(
      data.frame(filter = "default_cooks", as.data.frame(res_pvc[slc_id, , drop = FALSE])),
      data.frame(filter = "cooksCutoff_FALSE", as.data.frame(res_pvc_nocook[slc_id, , drop = FALSE]))
    )
    slc_tab$ENSEMBL <- slc_id; slc_tab$mode <- mode_tag; slc_tab$design <- design_tag
    write_csv_safe(slc_tab, file.path(DIR_SLC, "SLC25A13_DESeq2.csv"))
    norm_counts <- counts(dds, normalized = TRUE)[slc_id, ]
    slc_counts <- data.frame(sample = names(norm_counts), normalized_count = round(norm_counts, 2),
                             raw_count = cts_use[slc_id, names(norm_counts)],
                             coldata_use[names(norm_counts), c("status", "sex", "age_years", "batch")])
    write_csv_safe(slc_counts, file.path(DIR_SLC, "SLC25A13_counts_per_sample.csv"))
    p_slc <- ggplot(slc_counts, aes(x = status, y = normalized_count, colour = sex, shape = batch)) +
      geom_jitter(width = 0.12, size = 3, alpha = 0.9) +
      labs(title = fig_title(run_tag, "SLC25A13 normalized counts"),
           subtitle = sprintf("DESeq2 patient vs control: log2FC = %.2f, P = %.2g, FDR = %.2g%s",
                              res_pvc_nocook[slc_id, "log2FoldChange"], res_pvc_nocook[slc_id, "pvalue"],
                              res_pvc_nocook[slc_id, "padj"],
                              ifelse(is.na(res_pvc[slc_id, "padj"]),
                                     "\n(statistics computed without Cook's outlier filtering)", "")),
           x = NULL, y = "DESeq2 normalized count") +
      theme_bw() + theme(plot.subtitle = element_text(size = 8))
    save_pdf(p_slc, "SLC25A13_normalized_counts", w = 6, h = 4.5)
  } else {
    mode_log("[WARN] SLC25A13 not found in this count matrix")
  }
  mode_log("[STEP] DESeq2 done")

  # -----------------------------
  # VST / PCA
  # -----------------------------
  vsd <- if (identical(mode_tag, "ENRICHED")) {
    DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  } else {
    DESeq2::vst(dds, blind = FALSE)
  }
  X <- assay(vsd)
  batch <- coldata_use$batch
  X_resid <- X
  for (i in seq_len(nrow(X))) {
    fit <- lm(X[i, ] ~ batch)
    X_resid[i, ] <- residuals(fit)
  }
  write_csv_safe(as.data.frame(X), file.path(DIR_MAT, "VST_matrix_gene_x_sample.csv"), row.names = TRUE)
  write_csv_safe(as.data.frame(X_resid), file.path(DIR_MAT, "VST_residual_batch_gene_x_sample.csv"), row.names = TRUE)

  pca_batch_shapes  <- c(Run1 = 16, Run2 = 17, Run3 = 15)
  pca_status_colors <- c(control = "#FF725CFF", patient = "#6CC5B0FF")
  fix_batch_to_run <- function(x) { x0 <- trimws(as.character(x)); if (all(grepl("^[0-9]+$", x0))) x0 <- paste0("Run", x0); x0 }
  pca_run <- function(expr_gxs, meta) {
    pc <- prcomp(t(expr_gxs), center = TRUE, scale. = FALSE)
    ve <- (pc$sdev^2) / sum(pc$sdev^2) * 100
    df <- data.frame(sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                     status = meta[rownames(pc$x), "status"], batch = meta[rownames(pc$x), "batch"],
                     sex = meta[rownames(pc$x), "sex"], stringsAsFactors = FALSE)
    df$status <- factor(as.character(df$status), levels = names(pca_status_colors))
    df$batch  <- factor(fix_batch_to_run(df$batch), levels = names(pca_batch_shapes))
    list(df = df, ve = ve, pc = pc)
  }
  pca_plot <- function(pca_obj, title) {
    df <- pca_obj$df; ve <- pca_obj$ve
    ggplot(df, aes(PC1, PC2, shape = batch, colour = status)) +
      geom_point(size = 3, alpha = 0.9) +
      ggrepel::geom_text_repel(aes(label = sample), size = 2.5, show.legend = FALSE) +
      scale_shape_manual(values = pca_batch_shapes, drop = FALSE) +
      scale_color_manual(values = pca_status_colors, drop = FALSE) +
      theme_bw() +
      labs(title = title, x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
      guides(shape = guide_legend(order = 1), colour = guide_legend(order = 2))
  }
  meta0 <- coldata_use[colnames(X), , drop = FALSE]
  pca_pre  <- pca_run(X, meta0)
  pca_post <- pca_run(X_resid, meta0)
  write_csv_safe(pca_pre$df,  file.path(DIR_PCA, "PCA_scores_VST_raw.csv"))
  write_csv_safe(pca_post$df, file.path(DIR_PCA, "PCA_scores_VST_residual_batch.csv"))
  p_both <- pca_plot(pca_pre,  paste0("[", run_tag, "] PCA on VST (before batch correction)")) +
            pca_plot(pca_post, paste0("[", run_tag, "] PCA on VST residuals (batch-regressed)")) +
            patchwork::plot_layout(ncol = 2)
  save_pdf(p_both, "PCA_VST_pre_post_residual_batch", w = 12, h = 6)

  # -----------------------------
  # Volcano (shrunken LFC) with chrX/chrY and extreme-LFC labels (Reviewer 2)
  # -----------------------------
  make_volcano_df <- function(res, cts_all_df, lfc_cut = LFC_CUT, padj_cut = PADJ_CUT) {
    df <- as.data.frame(res)
    df$ENSEMBL <- rownames(df)
    df$SYMBOL  <- cts_all_df[df$ENSEMBL, "SYMBOL"]
    df$Chr     <- cts_all_df[df$ENSEMBL, "Chr"]
    df$neglog10_padj <- -log10(df$padj)
    df$signif <- "NS"
    df$signif[df$log2FoldChange >=  lfc_cut & df$padj < padj_cut] <- "Up"
    df$signif[df$log2FoldChange <= -lfc_cut & df$padj < padj_cut] <- "Down"
    df$sex_chromosome <- df$Chr %in% c("chrX", "chrY")
    df$Label <- ifelse(!is.na(df$SYMBOL) & df$SYMBOL != "", df$SYMBOL, df$ENSEMBL)
    df$extreme <- is.finite(df$log2FoldChange) & abs(df$log2FoldChange) >= EXTREME_LFC
    df
  }
  df_vol <- make_volcano_df(res_pvc_shr, cts_all_use)
  plot_df_vol <- subset(df_vol, is.finite(neglog10_padj) & is.finite(log2FoldChange))
  write_csv_safe(plot_df_vol, file.path(DIR_VOLMA, "Volcano_plotdata.csv"))
  # extreme and sex-chromosome genes among DEGs (Reviewer 2: "label extreme genes")
  extreme_tab <- df_vol[df_vol$signif != "NS" & (df_vol$extreme | df_vol$sex_chromosome), , drop = FALSE]
  extreme_tab <- extreme_tab[order(-abs(extreme_tab$log2FoldChange)), , drop = FALSE]
  write_csv_safe(extreme_tab, file.path(DIR_VOLMA, "DEG_extreme_or_sex_chromosome_genes.csv"))
  mode_log("[Volcano] DEGs with |LFC| >= ", EXTREME_LFC, " or on chrX/chrY: ", nrow(extreme_tab),
           " (chrX/chrY: ", sum(extreme_tab$sex_chromosome), ")")

  lab_df <- plot_df_vol[plot_df_vol$signif != "NS" & (plot_df_vol$extreme | plot_df_vol$sex_chromosome), , drop = FALSE]
  p_vol <- ggplot(plot_df_vol, aes(x = log2FoldChange, y = neglog10_padj)) +
    ggrastr::geom_point_rast(aes(color = signif), alpha = 0.7, size = 1.8, raster.dpi = 300) +
    geom_point(data = subset(plot_df_vol, sex_chromosome & signif != "NS"), shape = 21, size = 3, colour = "black", fill = NA) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed") +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed") +
    ggrepel::geom_text_repel(data = lab_df, aes(label = paste0(Label, " (", Chr, ")")), size = 2.5, max.overlaps = 40) +
    scale_color_manual(values = c(Up = "red", Down = "blue", NS = "black")) +
    theme_bw() +
    labs(title = fig_title(run_tag, "Volcano: patient vs control"),
         subtitle = paste0("LFC shrinkage: ", shr_type, "; circled = chrX/chrY DEGs; labelled = |LFC| >= ", EXTREME_LFC, " or chrX/chrY"),
         x = "log2 fold change", y = "-log10(adjusted p-value)", color = "DEG")
  save_pdf(p_vol, "Volcano", w = 7, h = 6)

  slc_row_v <- plot_df_vol[plot_df_vol$SYMBOL %in% "SLC25A13", , drop = FALSE]
  slc_note <- NULL
  if (nrow(slc_row_v) == 0 && length(slc_id) > 0) {
    # SLC25A13 dropped from the plot data because default Cook's filtering set its padj to NA:
    # place the point using the shrunken LFC and the padj computed without the Cook's filter.
    slc_row_v <- data.frame(log2FoldChange = as.data.frame(res_pvc_shr)[slc_id, "log2FoldChange"],
                            neglog10_padj = -log10(res_pvc_nocook[slc_id, "padj"]))
    slc_note <- "\npositioned with FDR computed without Cook's outlier filtering"
  }
  if (nrow(slc_row_v) > 0) {
    p_vol_slc <- p_vol + theme(legend.position = "none") +
      labs(subtitle = paste(c(sprintf("SLC25A13: log2FC = %.2f, P = %.2g, FDR = %.2g",
                                      res_pvc_nocook[slc_id, "log2FoldChange"], res_pvc_nocook[slc_id, "pvalue"],
                                      res_pvc_nocook[slc_id, "padj"]), slc_note), collapse = "")) +
      theme(plot.subtitle = element_text(size = 8)) +
      geom_point(data = slc_row_v, inherit.aes = FALSE, aes(x = log2FoldChange, y = neglog10_padj),
                 shape = 21, fill = "gold", color = "black", size = 4, stroke = 1.1) +
      geom_text(data = slc_row_v, inherit.aes = FALSE, aes(x = log2FoldChange, y = neglog10_padj, label = "SLC25A13"), vjust = -1.0, size = 4)
    save_pdf(p_vol_slc, "Volcano_SLC25A13_gold", w = 7, h = 6)
  }

  # -----------------------------
  # MA plot (shrunken LFC)
  # -----------------------------
  ma_df <- as.data.frame(res_pvc_shr) %>% tibble::rownames_to_column("ENSEMBL") %>%
    mutate(Amean = log2(baseMean + 1), logFC = log2FoldChange)
  ma_df$SYMBOL <- cts_all_use[ma_df$ENSEMBL, "SYMBOL"]
  ma_df$Chr <- cts_all_use[ma_df$ENSEMBL, "Chr"]
  ma_df <- ma_df %>%
    mutate(Label = if_else(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, ENSEMBL),
           Significance = case_when(
             is.finite(padj) & padj < PADJ_CUT & logFC >=  LFC_CUT ~ "Up (sig)",
             is.finite(padj) & padj < PADJ_CUT & logFC <= -LFC_CUT ~ "Down (sig)",
             TRUE ~ "Not Significant"),
           Significance = factor(Significance, levels = c("Up (sig)", "Down (sig)", "Not Significant")),
           sex_chromosome = Chr %in% c("chrX", "chrY"))
  write_csv_safe(ma_df, file.path(DIR_VOLMA, "MA_plotdata.csv"))
  ma_cols <- c("Up (sig)" = "red", "Down (sig)" = "blue", "Not Significant" = "grey")
  ma_lab <- ma_df %>% dplyr::filter(Significance != "Not Significant" & (abs(logFC) >= EXTREME_LFC | sex_chromosome))
  p_ma_A <- ggplot(ma_df, aes(x = Amean, y = logFC)) +
    ggrastr::geom_point_rast(aes(color = Significance), size = 1.8, alpha = 0.7, raster.dpi = 300) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = ma_cols, drop = FALSE) +
    labs(title = fig_title(run_tag, paste0("MA plot (", shr_type, " shrinkage)")), x = "Average log2 expression (A)", y = "log2 fold change") +
    theme_minimal() + theme(legend.position = "none")
  p_ma_B <- ggplot(ma_df, aes(x = Amean, y = logFC)) +
    ggrastr::geom_point_rast(aes(color = Significance), size = 1.8, alpha = 0.7, raster.dpi = 300) +
    geom_point(data = subset(ma_df, sex_chromosome & Significance != "Not Significant"), shape = 21, size = 3, colour = "black", fill = NA) +
    ggrepel::geom_text_repel(data = ma_lab, aes(label = paste0(Label, " (", Chr, ")")), size = 2.5, max.overlaps = 40) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = ma_cols, drop = FALSE) +
    labs(title = expression("Labeled: |LFC|" >= "5 or chrX/chrY DEGs"), x = "Average log2 expression (A)", y = "log2 fold change") +
    theme_minimal() + theme(legend.position = "none")
  slc_row_m <- ma_df %>% dplyr::filter(SYMBOL %in% "SLC25A13")
  if (nrow(slc_row_m) > 0) {
    hl <- function(p) p + geom_point(data = slc_row_m, inherit.aes = FALSE, aes(x = Amean, y = logFC),
                                     shape = 21, fill = "gold", color = "black", size = 3.5, stroke = 1.1) +
      geom_text(data = slc_row_m, inherit.aes = FALSE, aes(x = Amean, y = logFC, label = "SLC25A13"), vjust = -1.0, size = 4)
    p_ma_A <- hl(p_ma_A); p_ma_B <- hl(p_ma_B)
  }
  save_pdf((p_ma_A + p_ma_B) + patchwork::plot_annotation(title = fig_title(run_tag, "patient vs control: MA plot (DESeq2)")),
           "MA_pair_SLC25A13_highlight", w = 12, h = 6)

  # -----------------------------
  # GSEA (all terms exported with NES and FDR)
  # -----------------------------
  geneList <- make_geneList_from_deseq(res_pvc)
  write_csv_safe(data.frame(ENTREZID = names(geneList), stat = geneList,
                            SYMBOL = safe_mapIds(names(geneList), "SYMBOL", "ENTREZID", "first")),
                 file.path(DIR_GSEA, "ranked_geneList_Wald_stat.csv"))
  mode_log("[GSEA] ranked genes = ", length(geneList),
           if (mode_tag == "ENRICHED") " (captured genes = GSEA background)" else " (whole transcriptome)")

  gsea_objs <- list(
    GO_BP = safe_gseGO(geneList, "BP", run_tag),
    GO_MF = safe_gseGO(geneList, "MF", run_tag),
    GO_CC = safe_gseGO(geneList, "CC", run_tag),
    Reactome = safe_gsePathway(geneList, run_tag)
  )
  gsea_tabs <- list()
  for (nm in names(gsea_objs)) {
    gsea_tabs[[nm]] <- save_gsea_table(gsea_objs[[nm]], paste0("GSEA_", nm), nm)
    save_gsea_dot_ridge(gsea_objs[[nm]], paste0("GSEA_", nm),
                        subtitle = ifelse(mode_tag == "ENRICHED", "background = captured genes", "whole-transcriptome ranking"))
    n_t <- if (is.null(gsea_tabs[[nm]])) 0 else nrow(gsea_tabs[[nm]])
    n_s <- if (is.null(gsea_tabs[[nm]])) 0 else sum(gsea_tabs[[nm]]$p.adjust < PADJ_CUT)
    mode_log("[GSEA] ", nm, ": tested = ", n_t, ", padj < ", PADJ_CUT, " = ", n_s)
  }
  gsea_all <- do.call(rbind, gsea_tabs[!vapply(gsea_tabs, is.null, logical(1))])
  if (!is.null(gsea_all) && nrow(gsea_all)) {
    gsea_all <- gsea_all[order(gsea_all$p.adjust, gsea_all$pvalue), , drop = FALSE]
    write_csv_safe(gsea_all[, setdiff(colnames(gsea_all), "core_enrichment")], file.path(DIR_GSEA, "GSEA_ALL_ONTOLOGIES_results.csv"))
    save_nes_barplot(gsea_all, "GSEA_top_terms_NES_barplot")
    for (nm in names(gsea_tabs)) if (!is.null(gsea_tabs[[nm]])) save_nes_barplot(gsea_tabs[[nm]], paste0("GSEA_", nm, "_top_terms_NES_barplot"))
  }
  save_rds_safe(gsea_objs, file.path(DIR_RDATA, "gsea_objects.rds"))

  # Previously reported terms (Fig. 3C/D, Table S12, text) in this run
  prev <- data.frame(ID = PREVIOUSLY_REPORTED_TERMS, stringsAsFactors = FALSE)
  prev$Description <- NA_character_
  if (!is.null(gsea_all) && nrow(gsea_all)) {
    ix <- match(prev$ID, gsea_all$ID)
    prev$Description <- gsea_all$Description[ix]
    # fill names of terms that were not testable in this background (set size outside 10-500)
    miss <- is.na(prev$Description)
    if (any(miss)) {
      prev$Description[miss] <- tryCatch(
        unname(AnnotationDbi::Term(GO.db::GOTERM[prev$ID[miss]])), error = function(e) NA_character_)
    }
    prev$Ontology <- gsea_all$Ontology[ix]
    prev$setSize <- gsea_all$setSize[ix]
    prev$NES <- gsea_all$NES[ix]
    prev$pvalue <- gsea_all$pvalue[ix]
    prev$p.adjust <- gsea_all$p.adjust[ix]
    prev$status <- ifelse(is.na(ix), "not_testable_in_this_background",
                          ifelse(prev$p.adjust < PADJ_CUT, "FDR<0.05", "not_significant"))
  }
  prev$mode <- mode_tag; prev$design <- design_tag
  write_csv_safe(prev, file.path(DIR_GSEA, "previously_reported_terms_in_this_run.csv"))

  # -----------------------------
  # STRICT GSEA (padj < 0.05): subsets of the same objects
  # -----------------------------
  for (nm in names(gsea_objs)) {
    g2 <- strict_subset(gsea_objs[[nm]])
    if (is.null(g2)) { mode_log("[GSEA_STRICT] ", nm, ": no term with padj < ", PADJ_CUT); next }
    save_gsea_table(g2, paste0("GSEA_", nm, "_padj", PADJ_CUT), nm)
    save_gsea_dot_ridge(g2, paste0("GSEA_", nm, "_padj", PADJ_CUT),
                        subtitle = ifelse(mode_tag == "ENRICHED", "background = captured genes", "whole-transcriptome ranking"))
  }
  mode_log("[GSEA_STRICT] padj<", PADJ_CUT, " exports completed.")

  # -----------------------------
  # PANEL-AWARE ORA (ENRICHED only; universe = captured genes)
  # -----------------------------
  ora_summary <- NULL
  if (mode_tag == "ENRICHED" && !is.null(panel_info_use) && RUN_PANEL_ORA) {
    universe_entrez <- unique(panel_info_use$panel_entrez)
    universe_entrez <- universe_entrez[!is.na(universe_entrez) & universe_entrez != ""]
    sig_entrez <- make_sig_entrez_from_deseq(res_pvc)
    mode_log("[ORA] DEG Entrez = ", length(sig_entrez), "; universe (captured) Entrez = ", length(universe_entrez))
    if (length(sig_entrez) > 0 && length(universe_entrez) > 0) {
      ora_objs <- list(
        GO_BP = safe_enrichGO_panelORA(sig_entrez, universe_entrez, "BP", run_tag),
        GO_MF = safe_enrichGO_panelORA(sig_entrez, universe_entrez, "MF", run_tag),
        GO_CC = safe_enrichGO_panelORA(sig_entrez, universe_entrez, "CC", run_tag),
        Reactome = safe_enrichPathway_panelORA(sig_entrez, universe_entrez, run_tag)
      )
      for (nm in names(ora_objs)) {
        d <- save_ora_table(ora_objs[[nm]], paste0("ORApanel_", nm), nm)
        n_t <- if (is.null(d)) 0 else nrow(d); n_s <- if (is.null(d)) 0 else sum(d$p.adjust < PADJ_CUT)
        ora_summary <- rbind(ora_summary, data.frame(ontology = nm, tested = n_t, padj_lt_0.05 = n_s))
        mode_log("[ORA] ", nm, ": tested = ", n_t, ", padj < ", PADJ_CUT, " = ", n_s)
        if (!is.null(ora_objs[[nm]]) && is_gsea_nonempty(ora_objs[[nm]])) {
          p <- tryCatch(enrichplot::dotplot(ora_objs[[nm]], showCategory = SHOW_N) +
                          ggtitle(paste0("[", run_tag, "] ORA ", nm, " (universe = captured genes)")) + theme_bw(),
                        error = function(e) NULL)
          if (!is.null(p)) save_pdf(p, paste0("ORApanel_", nm, "_dotplot"), w = 9, h = 6)
        }
      }
      write_csv_safe(data.frame(ENTREZID = sig_entrez, SYMBOL = safe_mapIds(sig_entrez, "SYMBOL", "ENTREZID", "first")),
                     file.path(DIR_GSEA, "ORApanel_input_DEG_entrez.csv"))
    }
    mode_log("[ORA] PANEL-AWARE ORA completed for ENRICHED mode.")
  }

  # -----------------------------
  # Reproducibility snapshot + bundle
  # -----------------------------
  write_csv_safe(make_inputs_manifest(c(count_file, info_xlsx, panel_csv)), file.path(DIR_RDATA, "inputs_hashes.csv"))
  capture.output(sessionInfo(), file = file.path(DIR_RDATA, "sessionInfo.txt"))
  file.copy(file.path(DIR_RDATA, "sessionInfo.txt"), DIR_BUNDLE, overwrite = TRUE)
  file.copy(file.path(DIR_RDATA, "inputs_hashes.csv"), DIR_BUNDLE, overwrite = TRUE)
  if (file.exists(script_self_path)) file.copy(script_self_path, DIR_BUNDLE, overwrite = TRUE)

  mode_log("[OK] Finished run: ", run_tag, " -> ", dirs$ROOT)

  # Return a compact summary row for the cross-run table
  n_sig <- function(nm) if (is.null(gsea_tabs[[nm]])) NA_integer_ else sum(gsea_tabs[[nm]]$p.adjust < PADJ_CUT)
  n_tst <- function(nm) if (is.null(gsea_tabs[[nm]])) NA_integer_ else nrow(gsea_tabs[[nm]])
  slc_lfc <- if (length(slc_id) > 0) res_pvc[slc_id, "log2FoldChange"] else NA
  slc_p   <- if (length(slc_id) > 0) res_pvc_nocook[slc_id, "pvalue"] else NA
  slc_padj <- if (length(slc_id) > 0) res_pvc[slc_id, "padj"] else NA
  list(
    summary = data.frame(
      mode = mode_tag, design = design_tag, formula = deparse(design_formula),
      n_genes_in_matrix = nrow(cts_use), n_genes_ranked = length(geneList),
      n_DEG = nrow(de$sig), n_DEG_shrunkenLFC = nrow(de_shr), n_DEG_chrXY = sum(de$sig$Chr %in% c("chrX", "chrY")),
      n_DEG_extreme_LFC = sum(abs(de$sig$log2FoldChange) >= EXTREME_LFC),
      SLC25A13_log2FC = round(slc_lfc, 3), SLC25A13_p = signif(slc_p, 3), SLC25A13_padj = signif(slc_padj, 3),
      GO_BP_tested = n_tst("GO_BP"), GO_BP_FDR05 = n_sig("GO_BP"),
      GO_MF_tested = n_tst("GO_MF"), GO_MF_FDR05 = n_sig("GO_MF"),
      GO_CC_tested = n_tst("GO_CC"), GO_CC_FDR05 = n_sig("GO_CC"),
      Reactome_tested = n_tst("Reactome"), Reactome_FDR05 = n_sig("Reactome"),
      ORA_panel_terms_FDR05 = if (is.null(ora_summary)) NA_integer_ else sum(ora_summary$padj_lt_0.05),
      stringsAsFactors = FALSE),
    previously_reported = prev,
    gsea_all = gsea_all
  )
}

# ============================================================
# 3) RUN ALL MODES x DESIGNS
# ============================================================
run_summaries <- list(); prev_tabs <- list(); gsea_all_runs <- list()
for (mode_tag in RUN_MODES) {
  for (design_tag in names(DESIGNS)) {
    dirs <- make_run_dirs(OUTDIR, mode_tag, design_tag)
    res_run <- tryCatch({
      run_one_pipeline(
        mode_tag = mode_tag, design_tag = design_tag, design_formula = DESIGNS[[design_tag]],
        cts_use = analysis_inputs[[mode_tag]]$cts,
        cts_all_use = analysis_inputs[[mode_tag]]$cts_all,
        coldata_use = analysis_inputs[[mode_tag]]$coldata,
        dirs = dirs,
        panel_info_use = analysis_inputs[[mode_tag]]$panel_info
      )
    }, error = function(e) {
      errfile <- file.path(dirs$LOGS, paste0("ERROR_", mode_tag, "_", design_tag, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
      msg <- paste0("[FATAL ERROR] mode=", mode_tag, " design=", design_tag, "\ntime=", format(Sys.time()),
                    "\nmessage=", conditionMessage(e), "\n")
      message(msg); cat(msg, file = errfile, append = FALSE)
      NULL
    })
    if (is.null(res_run)) { master_log("[WARN] run failed and was skipped: ", mode_tag, "/", design_tag); next }
    key <- paste0(mode_tag, "_", design_tag)
    run_summaries[[key]] <- res_run$summary
    prev_tabs[[key]] <- res_run$previously_reported
    if (!is.null(res_run$gsea_all)) {
      g <- res_run$gsea_all; g$run <- key
      gsea_all_runs[[key]] <- g[, c("run", "mode", "design", "Ontology", "ID", "Description", "setSize", "NES", "pvalue", "p.adjust")]
    }
  }
}

# ============================================================
# 4) SUMMARY (cross-run tables for Methods / Results / Table S12)
# ============================================================
if (length(run_summaries)) {
  S <- do.call(rbind, run_summaries)
  write_csv_safe(S, file.path(SUMMARY_DIR, "summary_all_runs.csv"))
  master_log("[SUMMARY] runs = ", nrow(S))
  print(S[, c("mode", "design", "n_genes_ranked", "n_DEG", "SLC25A13_log2FC", "SLC25A13_p",
              "GO_BP_FDR05", "GO_MF_FDR05", "GO_CC_FDR05", "Reactome_FDR05")], row.names = FALSE)
}
if (length(prev_tabs)) {
  P <- do.call(rbind, prev_tabs)
  write_csv_safe(P, file.path(SUMMARY_DIR, "previously_reported_terms_all_runs.csv"))
  # Wide comparison: one row per term, NES/FDR per run
  term_names <- P %>% filter(!is.na(Description)) %>% group_by(ID) %>%
    summarise(Description = first(Description), Ontology = first(na.omit(Ontology)), .groups = "drop")
  Pw <- P %>% mutate(run = paste0(mode, "_", design)) %>%
    select(ID, run, setSize, NES, p.adjust) %>%
    pivot_wider(names_from = run, values_from = c(setSize, NES, p.adjust), names_glue = "{run}_{.value}") %>%
    left_join(term_names, by = "ID") %>%
    select(ID, Description, Ontology, everything())
  Pw <- Pw[match(PREVIOUSLY_REPORTED_TERMS, Pw$ID), , drop = FALSE]
  write_csv_safe(as.data.frame(Pw), file.path(SUMMARY_DIR, "previously_reported_terms_wide.csv"))
}
if (length(gsea_all_runs)) {
  G <- do.call(rbind, gsea_all_runs)
  write_csv_safe(G, file.path(SUMMARY_DIR, "GSEA_all_runs_long.csv"))
  # Table S12 candidate: ENRICHED/PRIMARY, all ontologies, ordered by FDR, with NES
  tS12 <- G[G$run == "ENRICHED_PRIMARY", , drop = FALSE]
  tS12 <- tS12[order(tS12$p.adjust, tS12$pvalue), , drop = FALSE]
  write_csv_safe(tS12, file.path(SUMMARY_DIR, "TableS12_candidate_ENRICHED_PRIMARY_GSEA.csv"))
  # Top terms per run (for a compact supplementary comparison)
  topk <- G %>% group_by(run) %>% arrange(p.adjust, pvalue, .by_group = TRUE) %>% slice_head(n = SHOW_N) %>% ungroup()
  write_csv_safe(as.data.frame(topk), file.path(SUMMARY_DIR, "GSEA_top15_per_run.csv"))
}
capture.output(sessionInfo(), file = file.path(SUMMARY_DIR, "sessionInfo.txt"))
if (file.exists(script_self_path)) file.copy(script_self_path, SUMMARY_DIR, overwrite = TRUE)

master_log("[OK] Finished. Modes x designs loop completed.")
master_log("[OK] Master output root: ", OUTDIR)
