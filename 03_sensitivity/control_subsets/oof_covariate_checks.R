## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
## -------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(openxlsx)})
D <- DATA
P <- file.path(REPO, "01_nested_ensemble", "results")
f <- Sys.getenv("CITRIN_TABLE_S_WORKBOOK", "")   # supplementary-table workbook (ages/sex); if absent, the distributed covariates CSV is used
si <- read.csv(file.path(D,"sample_info.csv")); cat("sample_info group x batch:\n"); print(table(si$group, si$batch))
p <- add_family(read.csv(file.path(P,"outer_primary_predictions.csv")))
## ages/sex: rebuilt from the Supplementary workbook when present; otherwise the previously extracted table is reused (ages/sex unchanged)
if (file.exists(f)) {
s1 <- read.xlsx(f, sheet=4, sep.names=" ", colNames=FALSE); s3 <- read.xlsx(f, sheet=6, sep.names=" ", colNames=FALSE); sx <- read.xlsx(f, sheet=7, sep.names=" ", colNames=FALSE)
age_num <- function(a){ a <- as.character(a); if(is.na(a)||a=="") return(NA); y <- regmatches(a, regexpr("[0-9.]+ *y", a)); m <- regmatches(a, regexpr("[0-9.]+ *m", a)); d <- regmatches(a, regexpr("[0-9.]+ *d", a))
  v <- 0; if(length(y)) v <- v + as.numeric(sub(" *y","",y)); if(length(m)) v <- v + as.numeric(sub(" *m","",m))/12; if(length(d)) v <- v + as.numeric(sub(" *d","",d))/365; v }
sexnorm <- function(s){ s <- tolower(as.character(s)); ifelse(s %in% c("woman","girl"), "F", ifelse(s %in% c("man","boy"), "M", NA)) }
cov <- rbind(
  data.frame(SampleID=as.character(unlist(s1[2, 3:ncol(s1)])), age=sapply(unlist(s1[3, 3:ncol(s1)]), age_num), sex=sexnorm(unlist(s1[4, 3:ncol(s1)]))),
  data.frame(SampleID=as.character(unlist(s3[2, 2:ncol(s3)])), age=sapply(unlist(s3[3, 2:ncol(s3)]), age_num), sex=sexnorm(unlist(s3[4, 2:ncol(s3)]))),
  data.frame(SampleID=as.character(sx[3:18,2]), age=sapply(sx[3:18,4], age_num), sex=sexnorm(sx[3:18,5])))
cov <- cov[!is.na(cov$SampleID) & cov$SampleID!="",]
write.csv(cov, file.path(REPO, "03_sensitivity", "control_subsets", "results", "covariates_from_tables.csv"), row.names=FALSE)
} else cov <- read.csv(file.path(REPO, "03_sensitivity", "control_subsets", "results", "covariates_from_tables.csv"), stringsAsFactors=FALSE)
p <- merge(p, cov, by="SampleID", all.x=TRUE)
cat("\nage available:", sum(!is.na(p$age)), "of 44; missing:", paste(p$SampleID[is.na(p$age)], collapse=","), "\n")
cat("\nAge by group (44):\n"); print(do.call(rbind, lapply(split(p$age, p$group), function(a) c(n=sum(!is.na(a)), median=median(a,na.rm=TRUE), min=min(a,na.rm=TRUE), max=max(a,na.rm=TRUE)))))
cat("\nSex by group (44):\n"); print(table(p$group, p$sex, useNA="ifany"))
auc <- function(s, y){ r <- rank(s); n1 <- sum(y==1); n0 <- sum(y==0); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
## subgroup re-aggregations of OOF predictions
SUB <- list(); sub <- function(keep, lab){ d <- p[p$truth==1 | keep,]; SUB[[length(SUB) + 1]] <<- data.frame(evaluation_set = lab, specimens = nrow(d), controls = sum(d$truth==0), AUC = round(auc(d$probability, d$truth), 3), sensitivity = sprintf("%d/7", sum(d$predicted_class[d$truth==1]==1)), specificity = sprintf("%d/%d", sum(d$predicted_class[d$truth==0]==0), sum(d$truth==0)), stringsAsFactors = FALSE); cat(sprintf("%-40s n=%d  AUC=%.3f  sens=%d/7  spec=%d/%d\n", lab, nrow(d), auc(d$probability, d$truth), sum(d$predicted_class[d$truth==1]==1), sum(d$predicted_class[d$truth==0]==0), sum(d$truth==0))) }
cat("\nOOF re-aggregations (consensus30, full_inner_recall):\n")
sub(p$truth==0, "all 37 controls")
sub(p$group %in% c("DD","WS"), "excluding healthy adults (DD+WS)")
sub(p$group=="DD", "disease controls only")
sub(p$group=="WS", "carriers only")
sub(p$group=="WW", "healthy adults only")
sub(p$truth==0 & !is.na(p$age) & p$age < 18, "controls aged <18 y (DD/WW only have age)")
sub(p$truth==0 & p$batch=="Run3", "Run3 controls only")
sub(p$truth==0 & p$batch=="Run2", "Run2 controls only")
## covariate association among controls
ctrl <- p[p$truth==0,]
cat("\nControls: Spearman score vs age (n=", sum(!is.na(ctrl$age)), "): rho=", round(cor(ctrl$probability, ctrl$age, method="spearman", use="complete"),3), " p=", signif(cor.test(ctrl$probability, ctrl$age, method="spearman", exact=FALSE)$p.value,2), "\n")
cat("Controls: score by sex (Wilcoxon): "); print(tapply(ctrl$probability, ctrl$sex, median)); cat(" p=", signif(wilcox.test(probability ~ sex, ctrl)$p.value,2), "\n")
cat("Controls: score by run: "); print(tapply(ctrl$probability, ctrl$batch, median)); cat(" p=", signif(wilcox.test(probability ~ batch, ctrl)$p.value,2), "\n")
cat("Controls: FP rate by run: "); print(tapply(ctrl$predicted_class, ctrl$batch, function(x) paste0(sum(x), "/", length(x))))
cat("\nFalse positives with age/sex:\n"); print(ctrl[ctrl$predicted_class==1, c("SampleID","group","batch","age","sex","probability","threshold")], row.names=FALSE)
cat("\nLoF/Mis with age/sex/score:\n"); print(p[p$truth==1, c("SampleID","family_id","age","sex","probability","threshold","predicted_class")], row.names=FALSE)
## score vs age within SM (7) and overall logistic-type check: AUC of age alone for SM vs controls with age
d <- p[!is.na(p$age),]; cat("\nAge alone as a classifier (SM vs controls with known age, n=", nrow(d), "): AUC=", round(auc(-d$age, d$truth),3), " (younger=positive)\n")
write.csv(p[, setdiff(names(p), "family_id")], file.path(REPO, "03_sensitivity", "control_subsets", "results", "oof_predictions_with_covariates.csv"), row.names=FALSE)
## tables of the re-aggregations and of the covariate checks (Supplementary Table S13)
write.csv(do.call(rbind, SUB), file.path(REPO, "03_sensitivity", "control_subsets", "results", "control_subgroup_performance.csv"), row.names = FALSE)
ct <- cor.test(ctrl$probability, ctrl$age, method = "spearman", exact = FALSE)
write.csv(data.frame(statistic = c("controls_with_age", "spearman_rho_score_vs_age", "spearman_P", "AUC_age_alone_LoFMis_vs_controls_with_age"),
  value = c(sum(!is.na(ctrl$age)), round(unname(ct$estimate), 3), signif(ct$p.value, 2), round(auc(-d$age, d$truth), 3))),
  file.path(REPO, "03_sensitivity", "control_subsets", "results", "control_covariate_association.csv"), row.names = FALSE)
