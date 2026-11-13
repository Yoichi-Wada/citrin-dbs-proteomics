## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
# Supplementary Table S26: quantities of the SLC25A13-assigned precursors in each of the 60 specimens, taken without
# modification from the DIA-NN precursor matrix (DBS-Tohoku.pr_matrix.tsv; not distributed here, part of the jPOST deposit).
# Empty cells of the precursor matrix (precursor not identified in the run) are written as 0. The table printed in the
# supplement rounds most values to two significant figures.
pr <- read.delim(file.path(EXT, "DBS-Tohoku.pr_matrix.tsv"), check.names = FALSE, quote = "", na.strings = c("", "NA"))
si <- read.csv(file.path(DATA, "sample_info.csv"), stringsAsFactors = FALSE)
ord <- c(grep("^SS_", si$SampleID, value = TRUE), grep("^WS_", si$SampleID, value = TRUE), grep("^SM_", si$SampleID, value = TRUE),
         grep("^WW_", si$SampleID, value = TRUE), grep("^DD_", si$SampleID, value = TRUE))
num <- function(x) as.integer(sub(".*_", "", x)); ord <- unlist(lapply(c("SS", "WS", "SM", "WW", "DD"), function(g) { s <- ord[startsWith(ord, paste0(g, "_"))]; s[order(num(s))] }))
stopifnot(length(ord) == 60, all(ord %in% names(pr)))
x <- pr[!is.na(pr$Genes) & pr$Genes == "SLC25A13", ]
out <- data.frame(Stripped_sequence = x$Stripped.Sequence, Charge = x$Precursor.Charge, x[, ord], check.names = FALSE)
out[is.na(out)] <- 0; out <- out[order(out$Stripped_sequence), ]
cat("SLC25A13-assigned precursors:", nrow(out), "; specimens with at least one quantified precursor:", sum(colSums(out[, ord] > 0) > 0), "of 60\n")
cat("precursors quantified per specimen, by group (median, range):\n"); n <- colSums(out[, ord] > 0)
print(do.call(rbind, lapply(split(n, substr(names(n), 1, 2)), function(v) c(median = median(v), min = min(v), max = max(v)))))
write.csv(out, file.path(REPO, "04_SLC25A13_measurement", "results", "TableS26_SLC25A13_precursor_quantities_60specimens.csv"), row.names = FALSE)
