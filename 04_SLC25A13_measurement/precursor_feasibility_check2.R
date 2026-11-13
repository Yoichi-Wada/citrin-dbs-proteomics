## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(openxlsx)})
D <- EXT   # FASTA (raw/) and pr_matrix are not distributed
out <- file.path(REPO, "04_SLC25A13_measurement", "results")
fa <- readLines(file.path(D,"raw/250404_GC_Human_UP000005640_E20644.fasta"))
getseq <- function(acc){ i <- grep(paste0("^>sp\\|",acc,"\\|"), fa); j <- grep("^>", fa); nxt <- j[j>i][1]; if(is.na(nxt)) nxt <- length(fa)+1; paste(fa[(i+1):(nxt-1)], collapse="") }
s13 <- getseq("Q9UJS0"); s12 <- getseq("O75746"); cat("len13",nchar(s13),"len12",nchar(s12),"\n")
pr <- fread(file.path(D,"DBS-Tohoku.pr_matrix.tsv"))
smp <- names(pr)[sapply(pr, is.numeric)]; smp <- setdiff(smp, c("Precursor.Charge","Proteotypic"))
sel <- pr[grepl("Q9UJS0|O75746", Protein.Group)]
sel[, n_det := rowSums(!is.na(.SD) & .SD>0), .SDcols=smp]
sel[, pos13 := mapply(function(p) regexpr(p, s13, fixed=TRUE), Stripped.Sequence)]
sel[, pos12 := mapply(function(p) regexpr(p, s12, fixed=TRUE), Stripped.Sequence)]
sel[, class := ifelse(pos13>0 & pos12>0, "shared", ifelse(pos13>0, "SLC25A13", ifelse(pos12>0,"SLC25A12","none")))]
res <- sel[, .(Protein.Group, Genes, Stripped.Sequence, z=Precursor.Charge, class, pos13, end13=ifelse(pos13>0,pos13+nchar(Stripped.Sequence)-1,NA), pos12, n_det)][order(class,pos13)]
print(res, nrows=40)
fwrite(res, file.path(out,"slc25a13_precursors_positions.csv"))
cov <- rep(FALSE,nchar(s13)); for(i in which(sel$pos13>0)) cov[sel$pos13[i]:(sel$pos13[i]+nchar(sel$Stripped.Sequence[i])-1)] <- TRUE
cat("SLC25A13 residues covered by any detected precursor:", sum(cov), "/", nchar(s13), sprintf("(%.1f%%)\n", 100*mean(cov)))
## per-specimen number of SLC25A13-specific precursors detected, by group
si <- fread(file.path(DATA,"sample_info.csv"))
spec <- sel[class=="SLC25A13"]
cn <- function(x) sub(".*[\\\\/]","", sub("\\.raw$|\\.d$","",x))
m <- as.matrix(spec[, ..smp]); colnames(m) <- cn(colnames(m)); m[is.na(m)] <- 0
ndet <- colSums(m>0); sumI <- colSums(m)
ps <- data.table(SampleID=names(ndet), n_specific_precursors_detected=ndet, sum_specific_intensity=sumI)
ps <- merge(ps, si, by="SampleID")
print(ps[, .(n=.N, median_nprec=median(n_specific_precursors_detected), min_nprec=min(n_specific_precursors_detected), max_nprec=max(n_specific_precursors_detected), n_zero=sum(n_specific_precursors_detected==0)), by=group][order(group)])
fwrite(ps, file.path(out,"slc25a13_specific_precursors_per_specimen.csv"))
## 4. empirical detection floor: SLC25A13 protein raw intensity vs dataset-wide quantified intensities
raw <- read.xlsx(file.path(DATA,"raw_data.xlsx"), sheet="data", sep.names=" ")
cat("raw_data.xlsx dims", dim(raw), " first cols:", paste(head(names(raw),4), collapse=","), "\n")
sc <- grep("^(SS|SM|WS|WW|DD)_[0-9]+$", names(raw), value=TRUE); stopifnot(length(sc)==60)
acc <- raw[["Master accession number"]]
M <- as.matrix(raw[, sc]); storage.mode(M) <- "numeric"; M[is.na(M)] <- 0; rownames(M) <- acc
npep <- raw[["Number"]]  # unique peptides? check name
cat("columns of raw_data.xlsx:", paste(setdiff(names(raw), sc), collapse=" | "), "\n")
grp <- setNames(si$group, si$SampleID)[sc]
x13 <- M["Q9UJS0", ]; actb <- M["P60709", ]
cat("\nSLC25A13 raw intensity detected (>0) by group:\n"); print(table(grp, x13>0))
## detection floor: per specimen, the smallest non-zero protein intensity and its 1st/5th percentile; where does SLC25A13 sit?
floor_tab <- data.table(SampleID=sc, group=grp,
  min_quant = apply(M, 2, function(v) min(v[v>0])),
  p01_quant = apply(M, 2, function(v) quantile(v[v>0], 0.01)),
  p05_quant = apply(M, 2, function(v) quantile(v[v>0], 0.05)),
  n_quant = colSums(M>0),
  SLC25A13 = x13,
  SLC25A13_pct_rank = sapply(seq_along(sc), function(j) if(x13[j]>0) mean(M[M[,j]>0, j] <= x13[j]) else NA_real_),
  SLC25A13_over_ACTB = ifelse(x13>0, x13/actb, NA))
print(floor_tab[, .(n=.N, median_min_quant=median(min_quant), median_p01=median(p01_quant),
  SLC25A13_median=median(SLC25A13[SLC25A13>0]), SLC25A13_min=suppressWarnings(min(SLC25A13[SLC25A13>0])),
  SLC25A13_pctrank_median=median(SLC25A13_pct_rank, na.rm=TRUE), SLC25A13_pctrank_min=suppressWarnings(min(SLC25A13_pct_rank, na.rm=TRUE))), by=group][order(group)])
fwrite(floor_tab, file.path(out,"detection_floor_per_specimen.csv"))
## ratio of the lowest detected SLC25A13 (any non-LoF/LoF specimen) to control median, on ACTB-normalised and raw scales
ctrl <- grp %in% c("DD","WW")
r_actb <- (x13/actb); r_actb[x13==0] <- NA
cat("\nLowest detected SLC25A13/ACTB relative to control median:", round(min(r_actb, na.rm=TRUE)/median(r_actb[ctrl], na.rm=TRUE),3),
    " which specimen:", names(which.min(r_actb)), "\n")
cat("Lowest detected raw SLC25A13 relative to control median raw:", round(min(x13[x13>0])/median(x13[ctrl]),3), "\n")
## SS specimens: is the raw ACTB and the number of quantified proteins comparable (sample quality)?
print(floor_tab[, .(median_nquant=median(n_quant), min_nquant=min(n_quant)), by=group][order(group)])
