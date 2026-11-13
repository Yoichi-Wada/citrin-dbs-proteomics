## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
# Exploratory abundance-gated extension of the two-step procedure (Supplementary Table S29, sections A and B, rows
# "gate + reported classifier"). Post hoc on the same 60 specimens; not an independent validation.
#   Step 1  SLC25A13 not quantified                      -> LoF/LoF
#   Gate    SLC25A13/ACTB above a family-held-out cut    -> negative branch (no classifier call)
#           cut = midpoint between the highest ratio among training carriers + LoF/Mis and the lowest ratio among
#                 training disease + healthy controls, estimated with the specimen's family held out
#   Step 2  remaining specimens: out-of-fold call of the reported ensemble (01_nested_ensemble/results/outer_predictions.csv)
# Also writes the 17-specimen input folder used by run_carrier_vs_LoFMis_retraining.sh (not distributed: contains family_map.csv).
suppressPackageStartupMessages(library(openxlsx))
OUT <- file.path(REPO, "03_sensitivity", "abundance_gate", "results"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
si <- add_family(read.csv(file.path(DATA, "sample_info.csv"), stringsAsFactors = FALSE))
rw <- read.xlsx(file.path(DATA, "raw_data.xlsx"), sheet = "data", check.names = FALSE, sep.names = " ")
M <- as.matrix(rw[, si$SampleID]); storage.mode(M) <- "numeric"; sym <- rw[["Master Gene Symbol"]]
slc <- M[which(sym == "SLC25A13")[1], ]; actb <- M[which(sym == "ACTB")[1], ]
d <- data.frame(SampleID = si$SampleID, group = si$group, batch = si$batch, family_id = si$family_id, stringsAsFactors = FALSE)
d$quantified <- is.finite(slc) & slc > 0
d$ratio_x1e4 <- ifelse(d$quantified, round(slc / actb * 1e4, 3), 0)
ratio <- ifelse(d$quantified, slc / actb * 1e4, NA_real_)
low <- d$group %in% c("WS", "SM"); high <- d$group %in% c("DD", "WW")
cut <- t(sapply(seq_len(nrow(d)), function(i) { tr <- d$family_id != d$family_id[i] & d$quantified
  c(maxlow = max(ratio[tr & low]), mid = (max(ratio[tr & low]) + min(ratio[tr & high])) / 2) }))
d$cut_maxlow_x1e4 <- round(cut[, "maxlow"], 3); d$cut_midpoint_x1e4 <- round(cut[, "mid"], 3)
d$gate_call_midpoint <- ifelse(!d$quantified, "LoF/LoF (SLC25A13 absent)", ifelse(ratio > cut[, "mid"], "Negative branch (high SLC25A13/ACTB)", "to stage 2"))
d$kept_maxlow <- d$quantified & ratio <= cut[, "maxlow"]
cat("gate thresholds (x1e-4) across held-out families:", paste(range(d$cut_midpoint_x1e4), collapse = " - "), "\n"); print(table(d$group, d$gate_call_midpoint))
write.csv(d[, setdiff(names(d), "family_id")], file.path(OUT, "stage1_ratio_gate_familyLOFO.csv"), row.names = FALSE)   # family identifiers are not distributed

## ---- gate + reported classifier (out-of-fold calls), and the reported two-stage procedure for comparison
cp <- function(k, n) { ci <- binom.test(k, n)$conf.int; sprintf("%d/%d (%.3f-%.3f)", k, n, ci[1], ci[2]) }
pr <- read.csv(file.path(REPO, "01_nested_ensemble", "results", "outer_predictions.csv"), stringsAsFactors = FALSE)
truth3 <- ifelse(d$group == "SS", "LoF/LoF", ifelse(d$group == "SM", "LoF/Mis", "other"))
truth4 <- ifelse(d$group == "SS", "LoF/LoF", ifelse(d$group == "SM", "LoF/Mis", ifelse(d$group == "WS", "carrier", "Ref/Ref")))
rows <- list()
for (st in c("consensus30", "diverse30", "adaptive")) for (po in c("full_inner_recall", "inner_spec90")) {
  q <- pr[pr$strategy == st & pr$policy == po, ]; pc <- q$predicted_class[match(d$SampleID, q$SampleID)]
  for (design in c("SLC25A13 gate + existing 44-classifier OOF", "current two-stage (absence + 44-classifier)")) {
    gated <- design == "SLC25A13 gate + existing 44-classifier OOF"
    call <- ifelse(!d$quantified, "LoF/LoF", ifelse(gated & d$gate_call_midpoint != "to stage 2", "Ref/Ref", ifelse(pc == 1, "LoF/Mis", "carrier")))
    fp <- d$group %in% c("WS", "DD", "WW") & call == "LoF/Mis"
    c3 <- ifelse(call %in% c("carrier", "Ref/Ref"), "other", call)
    rows[[length(rows) + 1]] <- data.frame(design = design, strategy = st, policy = po,
      LoFMis_sens = cp(sum(d$group == "SM" & call == "LoF/Mis"), 7), LoFMis_spec_vs37 = cp(37 - sum(fp), 37), FP = paste(d$SampleID[fp], collapse = ";"),
      acc3 = cp(sum(c3 == truth3), 60), acc4 = if (gated) cp(sum(call == truth4), 60) else "n/a", stringsAsFactors = FALSE)
  }
}
hyb <- do.call(rbind, rows); write.csv(hyb, file.path(OUT, "hybrid_gate_plus_existing_classifier_summary.csv"), row.names = FALSE)
print(hyb[hyb$strategy == "consensus30", ], row.names = FALSE)

## ---- input folder for the carrier-vs-LoF/Mis retraining (17 specimens); written outside results/ because it holds family_map.csv
IN <- file.path(REPO, "03_sensitivity", "abundance_gate", "results_rerun", "carrier_vs_LoFMis_retraining_input"); dir.create(IN, recursive = TRUE, showWarnings = FALSE)
k <- d$group %in% c("WS", "SM")
write.csv(read.csv(file.path(DATA, "sample_info.csv"), stringsAsFactors = FALSE)[k, ], file.path(IN, "sample_info.csv"), row.names = FALSE, quote = FALSE)
fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); write.csv(fm[fm$SampleID %in% d$SampleID[k], ], file.path(IN, "family_map.csv"), row.names = FALSE, quote = FALSE)
if (!file.exists(file.path(IN, "raw_data.xlsx"))) file.copy(file.path(DATA, "raw_data.xlsx"), file.path(IN, "raw_data.xlsx"))
cat("retraining input written to", IN, "(", sum(k), "specimens )\n")
