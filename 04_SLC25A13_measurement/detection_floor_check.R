## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table)})
D <- EXT   # FASTA (raw/) and pr_matrix are not distributed
out <- file.path(REPO, "04_SLC25A13_measurement", "results")
## 1. SLC25A13 (Q9UJS0) and SLC25A12 (O75746) sequences from the search FASTA
fa <- readLines(file.path(D,"raw/250404_GC_Human_UP000005640_E20644.fasta"))
getseq <- function(acc){ i <- grep(paste0("^>sp\\|",acc,"\\|"), fa); j <- grep("^>", fa); nxt <- j[j>i][1]; if(is.na(nxt)) nxt <- length(fa)+1
  paste(fa[(i+1):(nxt-1)], collapse="") }
s13 <- getseq("Q9UJS0"); s12 <- getseq("O75746")
cat("SLC25A13 length", nchar(s13), " SLC25A12 length", nchar(s12), "\n")
## 2. precursors in pr_matrix assigned to Q9UJS0 / O75746 (incl. shared groups)
pr <- fread(file.path(D,"DBS-Tohoku.pr_matrix.tsv"))
sel <- pr[grepl("Q9UJS0|O75746", Protein.Group) | grepl("Q9UJS0|O75746", Protein.Ids)]
smp <- grep("\\.raw$|\\.d$|^DD_|^SS_|^SM_|^WS_|^WW_", names(pr), value=TRUE)
smp <- names(pr)[sapply(pr, is.numeric)]; smp <- smp[!smp %in% c("Precursor.Charge","Proteotypic")]
sel[, n_det := rowSums(!is.na(.SD) & .SD>0), .SDcols=smp]
sel[, pos13 := regexpr(Stripped.Sequence, s13, fixed=TRUE)]
sel[, pos12 := regexpr(Stripped.Sequence, s12, fixed=TRUE)]
print(sel[, .(Protein.Group, Genes, Stripped.Sequence, Precursor.Charge, pos13, end13=ifelse(pos13>0,pos13+nchar(Stripped.Sequence)-1,NA), pos12, n_det)][order(pos13)])
# (the richer precursor table with class/z/end13 used for Table S25 is written by precursor_feasibility_check2.R)
fwrite(sel[, .(Protein.Group, Genes, Stripped.Sequence, Precursor.Charge, pos13, pos12, n_det)], file.path(out,"slc25a13_precursors_positions_detection_floor_check.csv"))
## 3. in-silico tryptic peptides of SLC25A13 (missed cleavage 0/1), with mass and m/z for z=2..4; library window 7-45 aa, m/z 495-705
mono <- c(G=57.02146,A=71.03711,S=87.03203,P=97.05276,V=99.06841,T=101.04768,C=103.00919,L=113.08406,I=113.08406,N=114.04293,D=115.02694,Q=128.05858,K=128.09496,E=129.04259,M=131.04049,H=137.05891,F=147.06841,R=156.10111,Y=163.06333,W=186.07931)
pepmass <- function(p) sum(mono[strsplit(p,"")[[1]]]) + 18.010565
digest <- function(s, mc=1){ aa <- strsplit(s,"")[[1]]; cuts <- which(aa %in% c("K","R") & c(aa[-1],"") != "P"); st <- c(1, cuts+1); en <- c(cuts, length(aa))
  res <- list(); for(i in seq_along(st)) for(m in 0:mc){ if(i+m>length(st)) next; res[[length(res)+1]] <- data.table(start=st[i], end=en[i+m], mc=m, seq=substr(s, st[i], en[i+m])) }
  rbindlist(res) }
tp <- digest(s13); tp[, len := nchar(seq)]; tp[, mass := sapply(seq, pepmass)]
tp[, mz2 := (mass+2*1.007276)/2]; tp[, mz3 := (mass+3*1.007276)/3]; tp[, mz4 := (mass+4*1.007276)/4]
tp[, in_lib_window := len>=7 & len<=45 & ((mz2>=495 & mz2<=705)|(mz3>=495 & mz3<=705)|(mz4>=495 & mz4<=705))]
tp[, in_SLC25A12 := sapply(seq, function(p) grepl(p, s12, fixed=TRUE))]
tp[, detected := seq %in% sel$Stripped.Sequence]
fwrite(tp, file.path(out,"slc25a13_insilico_tryptic_peptides.csv"))
cat("\nTryptic peptides (mc<=1):", nrow(tp), " in library window:", sum(tp$in_lib_window), " of which detected:", sum(tp$in_lib_window & tp$detected), "\n")
cat("Sequence coverage by detected precursors (residues):", {cov <- rep(FALSE,nchar(s13)); for(i in which(sel$pos13>0)) cov[sel$pos13[i]:(sel$pos13[i]+nchar(sel$Stripped.Sequence[i])-1)] <- TRUE; sum(cov)}, "/", nchar(s13), "\n")
## Known Japanese SLC25A13 missense/in-frame variants (for feasibility only; actual genotypes to be taken from the clinical table)
known <- data.table(variant=c("p.Arg184His? (c.550C>T)","p.Arg360Ter(nonsense, ref only)","p.Leu598Arg (c.1792T>G)","p.Gly437Glu (c.1310G>A)","p.Arg588Gln (c.1763G>A)","p.Ala312Thr? ","p.Glu601Ter"), pos=c(184,360,598,437,588,312,601))
## report which tryptic peptide covers each candidate position and whether it is in the library window
cov_tab <- known[, {p <- pos; r <- tp[start<=p & end>=p & mc==0]; .(peptide=r$seq, len=r$len, mz2=round(r$mz2,2), mz3=round(r$mz3,2), in_window=r$in_lib_window, wt_detected=r$detected)}, by=.(variant,pos)]
print(cov_tab)
fwrite(cov_tab, file.path(out,"candidate_missense_positions_wt_peptides.csv"))
