## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
D <- EXT   # search FASTA (raw/) is not distributed
fa <- readLines(file.path(D,"raw/250404_GC_Human_UP000005640_E20644.fasta"))
hdr <- grep("^>", fa); ends <- c(hdr[-1]-1, length(fa))
seqs <- mapply(function(a,b) paste(fa[a:b], collapse=""), hdr+1, ends)
names(seqs) <- sub("^>(\\S+).*","\\1", fa[hdr])
cat("FASTA entries:", length(seqs), "\n")
seqsIL <- gsub("I","L", seqs)
peps <- c("TVELLSGVVDQTK","EEFVLAAQK","STGSFVGELMYK","GLLPQLLGVAPEK","SSPQFGVTLLTYELLQR","INLPAPNPDHVGGYK","LAVATFAGIENK","FGLYLPLFKPSVSTSK",
          "KIYSTLAGTR","YEGFFGLYR","LQVAGEITTGPR")
res <- rbindlist(lapply(peps, function(p){
  hit <- names(seqs)[grepl(p, seqs, fixed=TRUE)]
  hitIL <- names(seqs)[grepl(gsub("I","L",p), seqsIL, fixed=TRUE)]
  data.table(peptide=p, exact_hits=paste(sub("^sp\\|[^|]+\\|","",hit), collapse=";"), n_exact=length(hit),
             IL_collapsed_hits=paste(sub("^sp\\|[^|]+\\|","",hitIL), collapse=";"), n_IL=length(hitIL))}))
options(width=200); print(res)
fwrite(res, file.path(REPO, "04_SLC25A13_measurement", "results", "slc25a13_peptide_uniqueness_fasta.csv"))
s12 <- seqs[grep("O75746", names(seqs))]; m <- regmatches(s12, regexpr("GL.PQL.GVAPEK", s12)); cat("SLC25A12 homologous segment:", m, "\n")
cat("Fisher exact, detection 0/16 LoF/LoF vs 44/44 others: P =", format(fisher.test(matrix(c(0,16,44,0),2))$p.value, digits=3), "\n")
