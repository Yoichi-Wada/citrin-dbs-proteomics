## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table)})
D <- EXT   # FASTA and pr_matrix are not distributed
fa <- readLines(file.path(D,"raw/250404_GC_Human_UP000005640_E20644.fasta"))
getseq <- function(acc){ i <- grep(paste0("^>sp\\|",acc,"\\|"), fa); j <- grep("^>", fa); nxt <- j[j>i][1]; if(is.na(nxt)) nxt <- length(fa)+1; paste(fa[(i+1):(nxt-1)], collapse="") }
s13 <- getseq("Q9UJS0"); s12 <- getseq("O75746")
mono <- c(G=57.02146,A=71.03711,S=87.03203,P=97.05276,V=99.06841,T=101.04768,C=103.00919,L=113.08406,I=113.08406,N=114.04293,D=115.02694,Q=128.05858,K=128.09496,E=129.04259,M=131.04049,H=137.05891,F=147.06841,R=156.10111,Y=163.06333,W=186.07931)
pepmass <- function(p) sum(mono[strsplit(p,"")[[1]]]) + 18.010565
digest <- function(s, mc=1){ aa <- strsplit(s,"")[[1]]; cuts <- which(aa %in% c("K","R") & c(aa[-1],"") != "P"); st <- c(1, cuts+1); en <- c(cuts, length(aa))
  res <- list(); for(i in seq_along(st)) for(m in 0:mc){ if(i+m>length(st)) next; res[[length(res)+1]] <- data.table(start=st[i], end=en[i+m], mc=m, seq=substr(s, st[i], en[i+m])) }
  rbindlist(res) }
annot <- function(tp){ tp[, len := nchar(seq)]; tp[, mass := sapply(seq, pepmass)]
  tp[, mz2 := (mass+2*1.007276)/2]; tp[, mz3 := (mass+3*1.007276)/3]; tp[, mz4 := (mass+4*1.007276)/4]
  tp[, in_lib_window := len>=7 & len<=45 & ((mz2>=495 & mz2<=705)|(mz3>=495 & mz3<=705)|(mz4>=495 & mz4<=705))]; tp }
## missense variants from Table S[XXX]
mis <- data.table(variant=c("c.701G>C","c.1178G>T","c.1478A>G","c.1592G>A","c.1801G>A"),
  prot=c("p.Arg234Thr","p.Gly393Val","p.Asp493Gly","p.Gly531Asp","p.Glu601Lys"),
  pos=c(234,393,493,531,601), ref=c("R","G","D","G","E"), alt=c("T","V","G","D","K"),
  individuals=c("SM_2, SM_3","SM_1","SM_7","SM_4, SM_5","SM_6"))
stopifnot(all(sapply(mis$pos, function(p) substr(s13,p,p))==mis$ref))
## detected precursors
pr <- fread(file.path(D,"DBS-Tohoku.pr_matrix.tsv"))
smp <- setdiff(names(pr)[sapply(pr, is.numeric)], c("Precursor.Charge","Proteotypic"))
cn <- function(x) sub(".*[\\\\/]","", sub("\\.raw$|\\.d$","",x))
sel <- pr[grepl("Q9UJS0", Protein.Group)]
sel[, pos13 := mapply(function(p) regexpr(p, s13, fixed=TRUE), Stripped.Sequence)]
sel <- sel[pos13>0]
M <- as.matrix(sel[, ..smp]); colnames(M) <- cn(colnames(M)); M[is.na(M)] <- 0
rownames(M) <- paste0(sel$Stripped.Sequence, "/z", sel$Precursor.Charge)
## per-variant: WT tryptic peptide(s) covering the position (mc 0) and the variant peptide(s)
tpWT <- annot(digest(s13, 1))
rows <- list()
for(i in seq_len(nrow(mis))){
  p <- mis$pos[i]; s_alt <- s13; substr(s_alt, p, p) <- mis$alt[i]
  tpV <- annot(digest(s_alt, 1))
  wt <- tpWT[start<=p & end>=p]; vt <- tpV[start<=p & end>=p]
  vt <- vt[!seq %in% tpWT$seq]  # peptides that differ from WT
  wt[, type := "WT"]; vt[, type := "variant"]
  r <- rbind(wt, vt); r[, variant := mis$prot[i]]; r[, individuals := mis$individuals[i]]
  r[, in_SLC25A12 := sapply(seq, function(q) grepl(q, s12, fixed=TRUE))]
  r[, identified_in_dataset := seq %in% sel$Stripped.Sequence]
  r[, n_specimens_detected := sapply(seq, function(q){ k <- which(sel$Stripped.Sequence==q); if(length(k)) sum(colSums(M[k,,drop=FALSE])>0) else 0L })]
  rows[[i]] <- r }
tab <- rbindlist(rows)[, .(variant, individuals, type, start, end, mc, seq, len, mz2=round(mz2,2), mz3=round(mz3,2), in_lib_window, in_SLC25A12, identified_in_dataset, n_specimens_detected)]
options(width=250); print(tab[mc==0]); cat("\n(with one missed cleavage)\n"); print(tab[mc==1 & in_lib_window==TRUE])
fwrite(tab, file.path(REPO, "04_SLC25A13_measurement", "results", "missense_variant_peptides_insilico.csv"))
## detection of each SLC25A13-specific precursor in SM specimens (and group summary)
si <- fread(file.path(DATA,"sample_info.csv")); grp <- setNames(si$group, si$SampleID)
sm <- paste0("SM_",1:7)
det <- (M[, sm] > 0) * 1
ord <- order(sel$pos13); det <- det[ord, ]
cat("\nSLC25A13-specific precursor detection (1 = quantified) in LoF/Mis specimens:\n"); print(det)
cat("\nDetection frequency by group (n specimens quantified / n):\n")
freq <- sapply(split(colnames(M), grp[colnames(M)]), function(cols) rowSums(M[ord, cols, drop=FALSE] > 0))
print(freq)
out <- data.table(precursor=rownames(det), start=sel$pos13[ord], det, freq)
fwrite(out, file.path(REPO, "04_SLC25A13_measurement", "results", "slc25a13_specific_precursor_detection_by_group.csv"))
## intensities of the two covered WT precursors in SM vs others
for(q in c("GLLPQLLGVAPEK","SSPQFGVTLLTYELLQR")){ k <- which(sel$Stripped.Sequence==q); v <- colSums(M[k,,drop=FALSE])
  cat("\n", q, ": quantified in", sum(v>0), "specimens; by group:", paste(names(table(grp[names(v)[v>0]])), table(grp[names(v)[v>0]]), collapse=" "), "\n")
  cat("  SM values:", paste(sm, signif(v[sm],3), collapse="; "), "\n") }
