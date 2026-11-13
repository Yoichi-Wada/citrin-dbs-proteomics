## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(openxlsx)})
D <- DATA
pr <- fread(file.path(EXT,"DBS-Tohoku.pr_matrix.tsv")); si <- fread(file.path(D,"sample_info.csv"))
smp <- setdiff(names(pr)[sapply(pr, is.numeric)], c("Precursor.Charge","Proteotypic"))
cn <- function(x) sub(".*[\\\\/]","", sub("\\.raw$|\\.d$","",x))
cat("columns (runs) in pr_matrix:", length(smp), "; replicate injections present:", any(duplicated(cn(smp))), "\n")
sel <- pr[Protein.Group=="Q9UJS0"]; M <- as.matrix(sel[, ..smp]); colnames(M) <- cn(colnames(M)); M[is.na(M)] <- 0; rownames(M) <- sel$Stripped.Sequence
grp <- setNames(si$group, si$SampleID)[colnames(M)]
non <- colnames(M)[grp!="SS"]; L <- log2(M[, non]); L[!is.finite(L)] <- NA
## inter-precursor consistency across the 44 specimens where SLC25A13 quantified
ok <- rowSums(!is.na(L)) >= 15
C <- cor(t(L[ok,]), use="pairwise.complete.obs", method="spearman")
cat("\nSpearman correlation between SLC25A13 precursor log-intensities across non-LoF/LoF specimens (precursors quantified in >=15):\n"); print(round(C,2))
cat("median pairwise rho:", round(median(C[upper.tri(C)]),2), " range:", round(range(C[upper.tri(C)]),2), "\n")
## protein-level CV within control groups (ACTB-normalised, linear scale)
raw <- read.xlsx(file.path(D,"raw_data.xlsx"), sheet="data", sep.names=" ")
sc <- grep("^(SS|SM|WS|WW|DD)_[0-9]+$", names(raw), value=TRUE); acc <- raw[["Master accession number"]]
x13 <- as.numeric(raw[acc=="Q9UJS0", sc]); actb <- as.numeric(raw[acc=="P60709", sc]); g <- setNames(si$group, si$SampleID)[sc]
r <- x13/actb
cv <- function(v) sd(v)/mean(v)
cat("\nSLC25A13/ACTB CV within group: ", paste(names(tapply(r[g!="SS"], g[g!="SS"], cv)), round(tapply(r[g!="SS"], g[g!="SS"], cv),2), collapse="; "), "\n")
cat("raw SLC25A13 CV within group: ", paste(names(tapply(x13[g!="SS"], g[g!="SS"], cv)), round(tapply(x13[g!="SS"], g[g!="SS"], cv),2), collapse="; "), "\n")
cat("ACTB raw CV within group: ", paste(names(tapply(actb, g, cv)), round(tapply(actb, g, cv),2), collapse="; "), "\n")
## MaxLFQ (pg_matrix) vs summed precursor vs raw_data protein for SLC25A13: consistency of quantification routes
pg <- fread(file.path(D,"DBS-Tohoku.pg_matrix.tsv")); pgs <- setdiff(names(pg)[sapply(pg, is.numeric)], character(0))
v <- as.numeric(pg[Protein.Group=="Q9UJS0", ..pgs]); names(v) <- cn(pgs); v[is.na(v)] <- 0
common <- intersect(names(v)[v>0], sc[x13>0])
cat("\nSpearman rho MaxLFQ vs raw_data protein intensity (n=", length(common), "):", round(cor(v[common], x13[match(common, sc)], method="spearman"),3), "\n")
cat("MaxLFQ zero in SS:", sum(v[sc[g=="SS"]]==0), "/16\n")
