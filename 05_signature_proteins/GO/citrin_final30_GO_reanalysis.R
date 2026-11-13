## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Original 15-protein GO enrichment reanalysis, citrin DBS proteomics rebuttal
# Run inside 20261029_rebuttal.Rproj:
# source("script/citrin_original15_GO_reanalysis.R", encoding="UTF-8")
# Inputs (project ROOT):
#   data/citrin_original15_GO_universe.csv  (derived from original raw_data.xlsx)
#   data/citrin_original15_QC_summary.csv   (original full-data refit 15 protein IDs)
# Outputs: Figures/GO_original15/ (static directory, overwritten on rerun)
# Requires: BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "GO.db"))
# Optional: install.packages(c("ggplot2", "readr"))  # base read.csv is used below
# GO BP is the PREDEFINED primary analysis, MF and CC descriptive secondary analyses.
# Important: this is a discovery-list enrichment, not classifier validation.

options(stringsAsFactors = FALSE)
required <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "GO.db", "ggplot2")
missing_pkg <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkg)) stop(
  "Missing R packages: ", paste(missing_pkg, collapse=", "),
  "\nInstall using BiocManager::install(c('clusterProfiler','org.Hs.eg.db','AnnotationDbi','GO.db')); install.packages('ggplot2')")

# Project-root detection without changing the working directory.
find_root <- function() {
  candidates <- unique(c(getwd(), dirname(getwd()), file.path(getwd(), "20261029_rebuttal")))
  good <- candidates[file.exists(file.path(candidates, "data", "citrin_original15_GO_universe.csv"))]
  if (!length(good)) stop("Cannot find data/citrin_original15_GO_universe.csv from getwd()=", getwd(),
                          "; open 20261029_rebuttal.Rproj and copy the supplied universe CSV into data/.")
  good[[1]]
}
root <- file.path(REPO, "05_signature_proteins", "GO")   # repository layout
input_universe <- file.path(root, "measured_proteome_GO_universe.csv")
input_15 <- file.path(root, "citrin_final30_panel.csv")
if (!file.exists(input_15)) stop("Missing: ", input_15)
out <- file.path(root, "results")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

raw <- read.csv(input_universe, check.names=FALSE)
qc <- read.csv(input_15, check.names=FALSE)
stopifnot(all(c("Gene", "UniquePeptides", "Eligible_3unique") %in% colnames(raw)))
stopifnot(all(c("Gene", "FeatureID") %in% colnames(qc)))
raw$Gene <- trimws(raw$Gene)
qc$Gene <- trimws(qc$Gene)
background_symbols <- unique(raw$Gene[!is.na(raw$Eligible_3unique) & raw$Eligible_3unique == 1L & nzchar(raw$Gene)])
core_symbols <- unique(qc$Gene[nzchar(qc$Gene)])
if (length(core_symbols) != 30L) stop("Expected 30 panel genes; got ", length(core_symbols))
if (!all(core_symbols %in% background_symbols)) stop("Original signature absent from the defined measurable universe: ",
    paste(setdiff(core_symbols, background_symbols), collapse=", "))
stopifnot(nrow(raw) == 7474L, sum(raw$Eligible_3unique == 1L,na.rm=TRUE) == 5199L)
message("Original 15 genes: ", paste(core_symbols, collapse=", "))
message("Original 5,199 eligible proteins (>=3 unique peptides) correspond to ", length(background_symbols), " unique gene symbols")

# Map SYMBOL to ENTREZ ID, used by GO annotations. Ambiguous symbols are EXCLUDED,
# not silently assigned to the first Entrez ID. Report this mapping explicitly.
annotation <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
  keys=background_symbols, columns="ENTREZID", keytype="SYMBOL")
annotation <- unique(annotation[, c("SYMBOL", "ENTREZID")])
annotation <- annotation[!is.na(annotation$ENTREZID) & nzchar(annotation$ENTREZID),,drop=FALSE]
counts <- table(annotation$SYMBOL)
ambig <- names(counts)[counts > 1]
if (length(ambig) > 0) warning("Gene symbols with multiple Entrez IDs excluded: ", paste(ambig, collapse=", "))
map <- annotation[!annotation$SYMBOL %in% ambig,,drop=FALSE]
if (anyDuplicated(map$SYMBOL)) stop("Nonunique gene-to-Entrez mapping")
write.csv(map, file.path(out,"GO_Entrez_mapping.csv"), row.names=FALSE, na="")
write.csv(data.frame(Gene=sort(background_symbols),
    mapped=sort(background_symbols) %in% map$SYMBOL,
    ambiguous=sort(background_symbols) %in% ambig,
    selected_in_final30=sort(background_symbols) %in% core_symbols),
    file.path(out,"GO_universe_mapping_audit.csv"), row.names=FALSE)
if (!all(core_symbols %in% map$SYMBOL)) stop("Mapping failed for core proteins: ", paste(setdiff(core_symbols,map$SYMBOL),collapse=", "))
background_entrez <- unique(map$ENTREZID)
core_entrez <- unique(map$ENTREZID[match(core_symbols,map$SYMBOL)])

# enrichGO uses hypergeometric over-representation testing.
# minGSSize/maxGSSize match the study's stated gene-set size bounds for RNA enrichment.
# Primary BP test: BH FDR within all tested BP categories (not only those reported).
# MF and CC are SECONDARY analyses; each adjusted within its own ontology.
result_list <- list()
for (ontology in c("BP","MF","CC")) {
  ego <- clusterProfiler::enrichGO(
    gene=core_entrez,
    universe=background_entrez,
    OrgDb=org.Hs.eg.db::org.Hs.eg.db,
    keyType="ENTREZID",
    ont=ontology,
    pAdjustMethod="BH",
    pvalueCutoff=1,
    qvalueCutoff=1,
    minGSSize=10,
    maxGSSize=500,
    readable=TRUE)
  res <- as.data.frame(ego)
  if (nrow(res)) {
    res$Ontology <- ontology
    res$InputN <- length(core_entrez)
    res$MappedUniverseN <- length(background_entrez)
    res <- res[order(res$p.adjust,res$pvalue),,drop=FALSE]
  }
  result_list[[ontology]] <- res
  write.csv(res, file.path(out,paste0("GO_",ontology,"_all_tested.csv")),row.names=FALSE,na="")
  sig <- if (nrow(res)) res[!is.na(res$p.adjust) & res$p.adjust < .05,,drop=FALSE] else res
  write.csv(sig, file.path(out,paste0("GO_",ontology,"_FDR_lt_0.05.csv")),row.names=FALSE,na="")
  message(ontology, ": tested=",nrow(res),", FDR<0.05=",nrow(sig))
}

# Focused audit of claims made in original Table S11 (original background 18,805)
focus <- c("GO:0051592", "GO:0009152", "GO:0009260", "GO:0046390")
ref <- data.frame(ID=focus,
    Original_BgRatio=c("139/18805","157/18805","173/18805","180/18805"),
    Original_pvalue=c(1.684908426548779e-4,2.413067579437091e-4,
                      3.209606707326698e-4,3.605525922819778e-4),
    Original_FDR=rep(.047502804033150575,4), stringsAsFactors=FALSE)
bp <- result_list$BP
if (nrow(bp)) {
  ix <- match(ref$ID,bp$ID)
  ref$Reanalysis_Description <- bp$Description[ix]
  ref$Reanalysis_GeneRatio <- bp$GeneRatio[ix]
  ref$Reanalysis_BgRatio <- bp$BgRatio[ix]
  ref$Reanalysis_Count <- bp$Count[ix]
  ref$Reanalysis_Genes <- bp$geneID[ix]
  ref$Reanalysis_pvalue <- bp$pvalue[ix]
  ref$Reanalysis_FDR <- bp$p.adjust[ix]
  ref$Pass_FDR_0.05 <- !is.na(ref$Reanalysis_FDR) & ref$Reanalysis_FDR < .05
  ref$Status <- ifelse(is.na(ref$Reanalysis_FDR),
      "NOT_TESTED_or_NOT_ENRICHED", ifelse(ref$Pass_FDR_0.05,"SIGNIFICANT","NOT_SIGNIFICANT"))
} else {
  ref$Status <- "NO_TESTABLE_GO_TERMS"
}
write.csv(ref,file.path(out,"GO_original_TableS11_claims_audit.csv"),row.names=FALSE,na="")

# -----------------------------------------------------------------------------
# ADDED: Post-hoc SLC25A13-annotation-guided targeted BP analyses.
# IMPORTANT: These are exploratory, hypothesis-guided analyses performed after
# examining the original enrichment results, NOT preregistered confirmatory tests.
# GOALL = direct GO annotations plus propagated ancestor terms (BP only).
# GO = directly annotated GO terms (BP only).
# Targeted BH correction is performed separately over the full annotated term set,
# assigning a p-value of 1 to unreturned/zero-overlap terms.
# No change is made to the comprehensive BP/MF/CC analysis above.
# -----------------------------------------------------------------------------
slc_map <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db, keys="SLC25A13", keytype="SYMBOL",
  columns="ENTREZID")
slc_entrez <- unique(stats::na.omit(slc_map$ENTREZID))
if (length(slc_entrez) != 1L)
  stop("SLC25A13 must resolve to exactly one Entrez ID; got: ",
       paste(slc_entrez, collapse=", "))

get_slc_go_bp <- function(all = FALSE) {
  col <- if (all) "GOALL" else "GO"
  ont <- if (all) "ONTOLOGYALL" else "ONTOLOGY"
  tab <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys=slc_entrez, keytype="ENTREZID", columns=c(col, ont))
  keep <- !is.na(tab[[col]]) & !is.na(tab[[ont]]) & tab[[ont]] == "BP"
  ids <- unique(as.character(tab[[col]][keep]))
  sort(ids[nzchar(ids)])
}
slc_goall_bp <- get_slc_go_bp(all=TRUE)
slc_direct_bp <- get_slc_go_bp(all=FALSE)
if (!length(slc_goall_bp) || !length(slc_direct_bp))
  stop("No SLC25A13 GO annotations found: inspect org.Hs.eg.db version")
if (length(slc_direct_bp) != 9L)
  warning("Direct BP GO annotations for SLC25A13: ",length(slc_direct_bp),
          " (Claude reported 9); check annotation-version / definition")
if (nrow(bp) != 657L)
  warning("Comprehensive BP returned ",nrow(bp),
          " terms (Claude reported 657); restricted BH values may differ")
if (!all(slc_direct_bp %in% slc_goall_bp))
  warning("Direct SLC25A13 GO:BP terms not fully contained in GOALL: check database version")

# FIXED DENOMINATOR FOR TARGETED MULTIPLICITY:
# Use ALL GO:BP IDs in the specified SLC25A13 annotation set, regardless of
# overlap with the 15 proteins. enrichGO omits zero-overlap GO terms even when
# pvalueCutoff=1, so limiting BH to rows returned by enrichGO is anti-conservative.
# Terms absent from the comprehensive result are assigned p=1. This also
# conservatively includes terms excluded by the pre-defined 10-500 size filter.
# The comprehensive BP p-values for overlapping terms remain unchanged.
# As the hypothesis sets were chosen after inspecting the original results,
# the corrected FDR values remain post-hoc exploratory, not confirmatory.
run_slc_targeted <- function(hypothesis_ids, analysis_name) {
  ids <- sort(unique(hypothesis_ids))
  ix <- match(ids, bp$ID)
  present <- !is.na(ix)
  tbl <- bp[ix, , drop=FALSE]
  tbl$ID <- ids
  tbl$Present_in_full_BP_results <- present
  tbl$pvalue[!present] <- 1
  tbl$Targeted_BH_FDR <- stats::p.adjust(tbl$pvalue, method="BH")
  tbl$Targeted_Significant <- !is.na(tbl$Targeted_BH_FDR) & tbl$Targeted_BH_FDR < .05
  tbl$HypothesisSet <- analysis_name
  tbl$TotalAnnotationTerms <- length(ids)
  tbl$HypothesesWithinComprehensiveTestedBP <- sum(present)
  tbl <- tbl[order(tbl$Targeted_BH_FDR, tbl$pvalue), , drop=FALSE]
  write.csv(tbl, file.path(out, paste0("GO_SLC25A13_",analysis_name,"_targeted.csv")),
            row.names=FALSE, na="")
  message("SLC25A13 ",analysis_name,": BP annotation terms=",length(ids),
          "; terms with a positive overlap returned by comprehensive BP=",sum(present),
          "; within-set FDR<0.05=",sum(tbl$Targeted_Significant))
  tbl
}
target_all <- run_slc_targeted(slc_goall_bp,"GOALL")
target_direct <- run_slc_targeted(slc_direct_bp,"DirectGO")

# Explicitly audit the ORIGINAL FOUR terms against both full and restricted
# analyses. Terms absent from a tested subset remain NA, never p=0.
append_target_to_ref <- function(ref,tab,prefix) {
  idx <- match(ref$ID, tab$ID)
  ref[[paste0(prefix,"_Included")]] <- !is.na(idx)
  ref[[paste0(prefix,"_FDR")]] <- tab$Targeted_BH_FDR[idx]
  ref[[paste0(prefix,"_FDR_lt_0.05")]] <- !is.na(ref[[paste0(prefix,"_FDR")]]) &
    ref[[paste0(prefix,"_FDR")]] < .05
  ref
}
ref <- append_target_to_ref(ref,target_all,"GOALL")
ref <- append_target_to_ref(ref,target_direct,"DirectGO")
write.csv(ref,file.path(out,"GO_original_TableS11_claims_audit.csv"),
          row.names=FALSE,na="")

# Make the annotation and multiplicity denominator explicit and inspectable.
all_ids <- sort(unique(c(slc_goall_bp,slc_direct_bp)))
annotation_audit <- data.frame(
  GO_ID=all_ids,
  Included_in_SLC25A13_GOALL_BP=all_ids %in% slc_goall_bp,
  Included_in_SLC25A13_direct_BP=all_ids %in% slc_direct_bp,
  Present_in_comprehensive_BP_results=all_ids %in% bp$ID,
  stringsAsFactors=FALSE)
write.csv(annotation_audit,
          file.path(out,"GO_SLC25A13_BP_hypothesis_set_audit.csv"),
          row.names=FALSE,na="")

# Keep the existing audit filenames. The all-annotation FDR should now match
# the corrected targeted analysis (it was previously only a sensitivity check).
for (name in c("GOALL","DirectGO")) {
  ids <- if (name == "GOALL") slc_goall_bp else slc_direct_bp
  tab <- if (name == "GOALL") target_all else target_direct
  p_all <- rep(1, length(ids))
  names(p_all) <- ids
  if (nrow(tab)) p_all[tab$ID] <- tab$pvalue
  fdr_all <- stats::p.adjust(p_all,"BH")
  sensitivity <- data.frame(GO_ID=ids,
      Full_BP_pvalue=p_all,
      FDR_including_absent_terms=as.numeric(fdr_all),
      Present_in_full_BP_results=ids %in% bp$ID,
      stringsAsFactors=FALSE)
  write.csv(sensitivity,
      file.path(out,paste0("GO_SLC25A13_",name,"_all_annotation_terms_sensitivity.csv")),
      row.names=FALSE,na="")
}

# Require an explicit sanity check that the original four reported GO terms
# belong to the hypothesis set if that is being claimed in the rebuttal.
if (!all(focus %in% slc_goall_bp)) {
  warning("Not all original Table S11 GO IDs are present in SLC25A13 GOALL BP: ",
          paste(setdiff(focus,slc_goall_bp),collapse=", "))
}

# Supplement-ready figure: show top 15 GO:BP terms; hollow labeling denotes non-significant.
if (nrow(bp)) {
  shown <- head(bp,15)
  shown$Description <- factor(shown$Description,levels=rev(shown$Description))
  shown$FDR_plot <- pmax(shown$p.adjust, .Machine$double.xmin)
  p <- ggplot2::ggplot(shown,ggplot2::aes(x= -log10(FDR_plot), y=Description)) +
    ggplot2::geom_point(ggplot2::aes(size=Count,fill=p.adjust<0.05),shape=21,stroke=.55,color="#273344") +
    ggplot2::scale_fill_manual(values=c("TRUE"="#4269D0","FALSE"="white"),name="BH FDR < 0.05") +
    ggplot2::geom_vline(xintercept=-log10(.05), linetype="dashed", color="#737D87") +
    ggplot2::labs(title="GO:BP enrichment of the final consensus-30 panel",
         subtitle="Background: proteins with >=3 unique peptides; BH-adjusted hypergeometric test",
         x="-log10(BH-adjusted P)",y=NULL,size="Overlap") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),
                  axis.text.y=ggplot2::element_text(color="#17212B"))
  ggplot2::ggsave(file.path(out,"FigureS_final30_GO_BP.pdf"),p,width=10,height=6.5)
  ggplot2::ggsave(file.path(out,"FigureS_final30_GO_BP.png"),p,width=10,height=6.5,dpi=300)
}

# Plain-text, provenance-aware report: no silent claims about significance.
report <- c(
  "Final consensus-30 panel GO ORA (2026-11-08)",
  paste0("Date: ", Sys.Date()),
  paste0("Input: ", input_15),
  paste0("Background source: ",input_universe),
  paste0("Proteins in original raw data: ",nrow(raw)),
  paste0("Unique-peptide-eligible proteins: ",sum(raw$Eligible_3unique==1L)),
  paste0("Distinct eligible gene symbols: ",length(background_symbols)),
  paste0("Unambiguously mapped Entrez gene IDs in eligible universe: ",length(background_entrez)),
  paste0("Selected original signature genes: ",length(core_symbols)),
  paste0("Mapped signature Entrez IDs: ",length(core_entrez)),
  paste0("Ambiguous gene symbols excluded: ",length(ambig)),
  paste0("GO primary ontology: BP; secondary: MF, CC"),
  "Method: clusterProfiler enrichGO; one-sided hypergeometric ORA; background=eligible measured proteins; BH within ontology",
  "Gene set size filter: minGSSize=10; maxGSSize=500",
  paste0("R version: ",R.version.string),
  paste0("clusterProfiler: ",as.character(utils::packageVersion("clusterProfiler"))),
  paste0("org.Hs.eg.db: ",as.character(utils::packageVersion("org.Hs.eg.db"))),
  paste0("GO.db: ",as.character(utils::packageVersion("GO.db"))),
  "Original Table S11 used an 18,805-gene background and its FDR is NOT reused here.",
  "Original 15 proteins are selected from the discovery cohort. Enrichment is exploratory, not independent validation.",
  "Added exploratory post-hoc SLC25A13 GOALL and direct-GO BP hypothesis-set analyses.",
  "Targeted FDR = BH over ALL SLC25A13-annotated GO:BP terms (unreturned terms assigned p=1), NOT only the positively overlapping enriched terms.",
  paste0("SLC25A13 Entrez ID: ",slc_entrez),
  paste0("SLC25A13 GOALL BP annotations: ",length(slc_goall_bp),"; intersected tested GO:BP: ",nrow(target_all)),
  paste0("SLC25A13 directly annotated GO:BP terms: ",length(slc_direct_bp),"; intersected tested GO:BP: ",nrow(target_direct)),
  "See GO_SLC25A13_*_all_annotation_terms_sensitivity.csv for audit of all annotated terms, including non-overlapping terms with p=1.",
  "See GO_original_TableS11_claims_audit.csv to determine whether 'response to calcium ion' remains significant.")
writeLines(report,file.path(out,"GO_analysis_methods_and_audit.txt"))
message("Completed GO reanalysis: ",out)
