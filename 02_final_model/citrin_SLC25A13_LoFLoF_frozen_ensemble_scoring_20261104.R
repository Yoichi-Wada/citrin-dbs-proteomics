#!/usr/bin/env Rscript
# ===========================================================================
# Citrin DBS proteomics rebuttal: score the 16 LoF/LoF (SS) specimens using
# the FROZEN full-training ensemble fitted on 44 non-SS specimens.
#
# REQUIRED: final_training_model_NOT_VALIDATED.rds from
#   citrin_nested_ensemble_familyLOFO_v1_1.R (model version 1.1).
#
# NO MODEL RE-FITTING, feature selection, training-set normalization fitting,
# batch-parameter estimation, or threshold re-optimization is performed here.
# Only the held-out/test-side MinProb transform from the v2.2 engine is used.
# The RDS fitted model is never modified.
#
# IMPORTANT:
#   - The model was trained on Run2/Run3, not Run1.
#   - no_combat: Run1 can be numerically scored but is an UNSEEN BATCH;
#                report Run1 separately, and do not claim validation.
#   - combat_mean: no frozen correction exists for Run1; no Run1 scores.
#   - These predictions are exploratory out-of-training-PHENOTYPE scores,
#     not an externally validated diagnostic sensitivity estimate.
#   - The final-model threshold is from training-set inner validation, NOT
#     the fold-specific thresholds used in the 44-sample OOF evaluation.
#
# Output directory (automatically created):
#   <project root>/result/LoFLoF_frozen_ensemble/
# In RStudio open 20261029_rebuttal.Rproj and run:
#   source("script/citrin_SLC25A13_LoFLoF_frozen_ensemble_scoring_20261103.R")
# ===========================================================================

options(stringsAsFactors = FALSE)

# USER OVERRIDES: normally leave blank.
ROOT_OVERRIDE <- normalizePath(Sys.getenv("CITRIN_REPO", "."))   # repository root (public release)
MODEL_RDS_OVERRIDE <- file.path(ROOT_OVERRIDE, "02_final_model", "FinalModel_reproducibility", "final_training_model_FROZEN_NOT_VALIDATED.rds")
DATA_DIR_OVERRIDE <- Sys.getenv("CITRIN_DATA_DIR", file.path(ROOT_OVERRIDE, "data"))
OUTPUT_DIR_NAME <- file.path("02_final_model", "results", "LoFLoF_frozen_ensemble")

# Preserve the full-training context and random-seed implementation of v1.1.
FULL_TRAINING_CONTEXT <- "FINAL_ALL44/fit"
ENGINE_BASE_SEED <- 2025L
ENGINE_SEED_SALT <- 11111L

needed <- c("openxlsx", "glmnet")
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1L), quietly=TRUE)]
if (length(missing_pkgs)) {
  stop("Required R packages missing: ", paste(missing_pkgs, collapse=", "),
       ". Install with install.packages(c('openxlsx','glmnet')).")
}

norm <- function(x) toupper(trimws(as.character(x)))
canon <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
find_column <- function(z, options) {
  i <- match(canon(options), canon(names(z)))
  i <- i[!is.na(i)]
  if (!length(i)) {
    stop("Required column not found (", paste(options, collapse=" / "),
         "). Columns: ", paste(names(z), collapse=", "))
  }
  names(z)[i[1L]]
}

find_project <- function() {
  if (nzchar(ROOT_OVERRIDE)) {
    if (!dir.exists(ROOT_OVERRIDE)) stop("ROOT_OVERRIDE directory missing")
    return(normalizePath(ROOT_OVERRIDE, winslash="/", mustWork=TRUE))
  }
  candidates <- c(getwd())
  if (requireNamespace("rstudioapi", quietly=TRUE)) {
    p <- tryCatch(rstudioapi::getActiveProject(), error=function(e) NULL)
    if (is.character(p) && length(p)==1L && nzchar(p)) candidates <- c(p, candidates)
  }
  cur <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
  for (i in seq_len(6L)) {
    candidates <- c(candidates, cur)
    if (identical(cur, dirname(cur))) break
    cur <- dirname(cur)
  }
  candidates <- unique(candidates[dir.exists(candidates)])
  matches <- candidates[file.exists(file.path(candidates, "20261029_rebuttal.Rproj"))]
  if (!length(matches)) {
    stop("Project root not found. Open 20261029_rebuttal.Rproj in RStudio, ",
         "or set ROOT_OVERRIDE at the beginning of this script.")
  }
  normalizePath(matches[1L], winslash="/", mustWork=TRUE)
}

find_data_dir <- function(root) {
  candidates <- c(DATA_DIR_OVERRIDE, Sys.getenv("CITRIN_DATA_DIR", ""),
                  file.path(dirname(root), "data"), file.path(root,"data"),
                  file.path(root,"R","data"), file.path(dirname(dirname(root)),"data"))
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[vapply(candidates, function(d)
    file.exists(file.path(d, "raw_data.xlsx")) &&
    file.exists(file.path(d, "sample_info.csv")), logical(1L))]
  if (!length(hit)) {
    stop("Cannot find raw_data.xlsx + sample_info.csv. Set DATA_DIR_OVERRIDE ",
         "(the same input directory used for the 44-sample ensemble).")
  }
  normalizePath(hit[1L], winslash="/", mustWork=TRUE)
}

find_model <- function(root) {
  candidates <- c(MODEL_RDS_OVERRIDE, Sys.getenv("CITRIN_SS_MODEL_RDS", ""),
    file.path(root,"citrin_nested_ensemble_v1_1","final_training_model_NOT_VALIDATED.rds"),
    file.path(root,"Claude","citrin_nested_ensemble_v1_1_lofo",
              "final_training_model_NOT_VALIDATED.rds"),
    file.path(root,"result","citrin_nested_ensemble_v1_1",
              "final_training_model_NOT_VALIDATED.rds"),
    file.path(root,"final_training_model_NOT_VALIDATED.rds"))
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Frozen full-training RDS is missing. Look for ",
         "final_training_model_NOT_VALIDATED.rds in the ensemble results. ",
         "This script will NOT train a replacement model. Set ",
         "MODEL_RDS_OVERRIDE to the original final-fit RDS path if needed.")
  }
  normalizePath(hit[1L], winslash="/", mustWork=TRUE)
}

# Exact v2.2 / v1.1 seed algorithm for a final-training member.
# v1.1: key_seed("fit", "FINAL_ALL44/fit", member_seed)
# v2.2: family_seed(paste(...,sep="|"), salt=11111L), SEED=2025L.
model_scoring_seed <- function(member_seed) {
  key <- paste("fit", FULL_TRAINING_CONTEXT, as.integer(member_seed), sep="|")
  z <- tempfile(); on.exit(unlink(z), add=TRUE)
  writeLines(paste0(ENGINE_BASE_SEED, "|", ENGINE_SEED_SALT, "|", key),
             z, useBytes=TRUE)
  h <- unname(tools::md5sum(z))
  seed_fit <- as.integer((strtoi(substr(h,1L,7L),base=16L) +
                         as.double(ENGINE_BASE_SEED) +
                         as.double(ENGINE_SEED_SALT)) %% 2000000000) + 1L
  seed_fit + 4L  # v1.1 passes seed+4 to score_preprocessor()
}

validate_model <- function(m) {
  if (!is.list(m) || !identical(m$version, "citrin_nested_ensemble_familyLOFO_1.1") ||
      !identical(m$engine_version, "two_run_fixed_panel_2.2")) {
    stop("Unexpected RDS schema/version. Expected frozen v1.1 ensemble + v2.2 engine.")
  }
  if (!is.character(m$fit_mode) || length(m$fit_mode)!=1L ||
      !m$fit_mode %in% c("no_combat","combat_mean"))
    stop("Invalid model fit_mode; refusing to infer batch handling")
  if (!is.character(m$optimization_status) ||
      !identical(m$optimization_status, "TRAINING-DERIVED NOT EXTERNALLY VALIDATED"))
    stop("Model optimization status is not the expected full-training specification")
  if (length(m$feature_set)!=30L || anyDuplicated(m$feature_set) ||
      length(m$member_seeds)!=10L || length(m$models)!=10L ||
      anyNA(m$member_seeds) || anyDuplicated(m$member_seeds))
    stop("Model does not have a unique consensus30 + 10-member architecture")
  if (!identical(as.integer(m$member_seeds), 2025:2034))
    stop("Unexpected member seeds: inspect before extrapolating to LoF/LoF")
  req <- c("strategy","full_threshold","spec90_threshold")
  if (!is.data.frame(m$inner_thresholds) ||
      !all(req %in% names(m$inner_thresholds)) ||
      !setequal(as.character(m$inner_thresholds$strategy), c("consensus30","diverse30")) ||
      nrow(m$inner_thresholds)!=2L)
    stop("Missing/unexpected final-training inner-derived threshold table")
  for (nm in c("full_threshold", "spec90_threshold")) {
    v <- as.numeric(m$inner_thresholds[[nm]])
    if (any(!is.finite(v) | v < 0 | v > (1 + 1e-10)))
      stop("Invalid final training threshold column: ", nm)
  }
  for (k in seq_along(m$models)) {
    mm <- m$models[[k]]
    pre <- mm$preprocessor
    if (!is.list(pre) || !identical(pre$mode,m$fit_mode) ||
        !is.character(pre$kept) || anyDuplicated(pre$kept) ||
        !all(c("mu","sd","minprob_sigma") %in% names(pre)))
      stop("Invalid saved preprocessor for member ",k)
    if (!all(m$feature_set %in% pre$kept) ||
        !all(mm$diverse_ids %in% pre$kept) || length(mm$diverse_ids)!=30L)
      stop("Member ",k, " lacks required frozen protein features")
    if (!all(pre$kept %in% names(pre$mu)) || !all(pre$kept %in% names(pre$sd)) ||
        any(!is.finite(pre$mu[pre$kept])) ||
        any(!is.finite(pre$sd[pre$kept]) | pre$sd[pre$kept] <= 0) ||
        !is.finite(pre$minprob_sigma) || pre$minprob_sigma <= 0)
      stop("Invalid frozen means/SD/imputation scale in member ",k)
    if (m$fit_mode == "combat_mean") {
      if (!is.matrix(pre$shifts) ||
          !all(c("Run2","Run3") %in% colnames(pre$shifts)))
        stop("Missing frozen Run2/Run3 ComBat shifts in member ",k)
    }
    for (a in c("consensus_fit","diverse_fit")) {
      fit <- mm[[a]]
      if (!is.list(fit) || is.null(fit$fit) ||
          length(fit$lambda)!=1L || !is.finite(fit$lambda))
        stop("Missing fitted glmnet object/lambda for member ",k," (",a,")")
    }
  }
  invisible(TRUE)
}


verify_analysis_plan_if_present <- function(model, model_file, data_dir) {
  plan_file <- file.path(dirname(model_file), "immutable_analysis_plan.csv")
  if (!file.exists(plan_file)) {
    warning("immutable_analysis_plan.csv not adjacent to fitted model RDS; ",
            "input-file hash concordance cannot be independently verified.")
    return(invisible(FALSE))
  }
  plan <- utils::read.csv(plan_file, stringsAsFactors=FALSE,
                          check.names=FALSE)
  if (!all(c("name","value") %in% names(plan)))
    stop("Invalid immutable_analysis_plan.csv: missing name/value columns")
  kv <- setNames(as.character(plan$value),as.character(plan$name))
  for (item in c("raw_md5","sample_info_md5")) {
    if (!item %in% names(kv)) stop("Missing frozen input checksum: ",item)
    f <- if (item=="raw_md5") "raw_data.xlsx" else "sample_info.csv"
    actual <- unname(tools::md5sum(file.path(data_dir,f)))
    if (!identical(tolower(kv[[item]]),tolower(actual)))
      stop("Input MD5 mismatch for ",f,
           ": model was trained from a different version of the file")
  }
  if (!"mode" %in% names(kv) || !identical(kv[["mode"]],model$fit_mode))
    stop("Saved ensemble fit_mode differs from immutable training plan")
  if (!"seeds" %in% names(kv) ||
      !identical(kv[["seeds"]],paste(model$member_seeds,collapse=",")))
    stop("Saved member seeds differ from immutable training plan")
  message("[CHECK] Frozen RDS training-plan hashes, mode and seeds are consistent")
  invisible(TRUE)
}

load_ss_features <- function(data_dir) {
  meta <- utils::read.csv(file.path(data_dir,"sample_info.csv"),
                          stringsAsFactors=FALSE, check.names=FALSE,
                          fileEncoding="UTF-8-BOM")
  sid_col <- find_column(meta, c("SampleID","sample_id","sample id"))
  grp_col <- find_column(meta, "group")
  bat_col <- find_column(meta, "batch")
  meta$SampleID <- norm(meta[[sid_col]])
  meta$group <- norm(meta[[grp_col]])
  meta$batch <- trimws(as.character(meta[[bat_col]]))
  if (anyDuplicated(meta$SampleID) || anyNA(meta$SampleID) ||
      anyNA(meta$group) || anyNA(meta$batch))
    stop("Invalid SampleID/group/batch in sample_info.csv")
  ss <- meta[meta$group == "SS",,drop=FALSE]
  if (nrow(ss) != 16L || !setequal(ss$SampleID, paste0("SS_",1:16)))
    stop("Expected SS_1..SS_16 (16 LoF/LoF specimens), got: ",
         paste(ss$SampleID,collapse=", "))
  ss <- ss[order(as.integer(sub("^SS_", "", ss$SampleID))),,drop=FALSE]
  if (!all(ss$batch %in% c("Run1","Run2","Run3")))
    stop("Unexpected SS acquisition batch label")
  if (!identical(as.integer(table(factor(ss$batch, levels=c("Run1","Run2","Run3")))),
                 c(8L,4L,4L)))
    stop("Unexpected SS batch composition: expected Run1=8, Run2=4, Run3=4")

  raw <- openxlsx::read.xlsx(file.path(data_dir,"raw_data.xlsx"),
                             check.names=FALSE)
  col_acc <- find_column(raw, c("Master Accession Number", "master_accession_number"))
  col_pep <- find_column(raw, c("Unique peptides", "unique_peptides"))
  acc <- norm(raw[[col_acc]])
  if (anyNA(acc) || any(!nzchar(acc))) stop("Invalid protein accession identifiers")
  if (anyDuplicated(acc)) {
    # Same deduplication order as v2.2 read_44().
    keep <- !duplicated(acc)
    raw <- raw[keep,,drop=FALSE]
    acc <- acc[keep]
  }
  if (!"P60709" %in% acc) stop("ACTB P60709 not present in raw_data.xlsx")
  col_idx <- match(ss$SampleID, norm(names(raw)))
  if (anyNA(col_idx)) stop("Missing SS sample columns in raw_data.xlsx")
  vv <- as.matrix(raw[,col_idx,drop=FALSE])
  suppressWarnings(storage.mode(vv) <- "numeric")
  colnames(vv) <- ss$SampleID
  rownames(vv) <- acc
  actb <- as.numeric(vv["P60709",])
  if (any(!is.finite(actb) | actb <= 0))
    stop("Missing/nonpositive ACTB in one or more SS specimens")
  vv[!is.finite(vv) | vv <= 0] <- NA_real_
  slc <- if ("Q9UJS0" %in% rownames(vv))
    as.numeric(vv["Q9UJS0",]) / actb else rep(NA_real_,ncol(vv))
  slc_detected <- if ("Q9UJS0" %in% rownames(vv))
    is.finite(vv["Q9UJS0",]) else rep(NA,ncol(vv))

  x <- sweep(log2(vv),2L,log2(actb),"-")
  pep <- suppressWarnings(as.numeric(raw[[col_pep]]))
  eligible <- is.finite(pep) & pep >= 3L
  if (sum(eligible) != 5199L)
    stop("Expected 5,199 initial protein features (>=3 unique peptides), found ",
         sum(eligible),". Check raw-data revision/model provenance.")
  x <- t(x[eligible,,drop=FALSE])  # SS specimens x protein accessions
  stopifnot(identical(rownames(x),ss$SampleID), !anyDuplicated(colnames(x)))
  list(X=x, ss=ss, slc_ratio=slc, slc_detected=slc_detected,
       actb=actb, n_features=ncol(x))
}

# Test-side MinProb transformation, identical to v2.2 score_preprocessor():
# q=0.01 from each SS specimen, sigma frozen from model training; rnorm() for
# left-censored features. The RNG is fixed per member, with ascending SS_ID
# sample ordering. No SS class labels or model refitting are used.
score_preprocessor_frozen <- function(pre, X_test, test_batch, scoring_seed) {
  if (nrow(X_test)!=length(test_batch)) stop("Batch labels are misaligned")
  z <- X_test
  if (!is.matrix(z)) z <- as.matrix(z)
  set.seed(scoring_seed)
  for (i in seq_len(nrow(z))) {
    observed <- as.numeric(z[i,is.finite(z[i,])])
    if (length(observed) < 10L)
      stop("<10 measured proteins in SS specimen ",rownames(z)[i])
    q <- unname(stats::quantile(observed,probs=0.01,
                               na.rm=TRUE,names=FALSE))
    idx <- which(!is.finite(z[i,]))
    if (length(idx)) z[i,idx] <- stats::rnorm(length(idx),q,pre$minprob_sigma)
  }
  if (any(!is.finite(z))) stop("Nonfinite SS imputed protein abundance")
  if (pre$mode == "combat_mean") {
    if (!all(test_batch %in% colnames(pre$shifts)))
      stop("Attempted unseen batch under frozen ComBat. No extrapolation allowed.")
    for (b in unique(test_batch)) {
      ix <- which(test_batch == b)
      if (!all(colnames(z) %in% rownames(pre$shifts)))
        stop("Frozen ComBat shift feature mismatch")
      z[ix,] <- sweep(z[ix,,drop=FALSE],2L,
                     pre$shifts[colnames(z),b],"-")
    }
  }
  if (!all(pre$kept %in% colnames(z)))
    stop("Model expects protein features missing from raw data")
  z <- z[,pre$kept,drop=FALSE]
  z <- sweep(sweep(z,2L,pre$mu[pre$kept],"-"),
             2L,pre$sd[pre$kept],"/")
  if (any(!is.finite(z))) stop("Nonfinite SS standardized protein abundance")
  z
}

score_frozen_model <- function(model, input) {
  ss <- input$ss
  mode <- model$fit_mode
  known_batches <- c("Run2", "Run3")
  in_training_batch <- ss$batch %in% known_batches
  # For no_combat, RUN1 numerical scoring is allowed but clearly flagged.
  # For combat_mean, no Run1 correction was learned, so retain Run1 as NA.
  scoreable <- if (mode == "no_combat") rep(TRUE,nrow(ss)) else in_training_batch
  if (!any(scoreable)) stop("No SS specimens can be scored by this model")

  X <- input$X[scoreable,,drop=FALSE]
  batch <- ss$batch[scoreable]
  n <- nrow(X)
  K <- length(model$models)
  con_member <- matrix(NA_real_,n,K,
                       dimnames=list(rownames(X),as.character(model$member_seeds)))
  div_member <- con_member
  member_logs <- vector("list",2L*K)
  j <- 0L
  for (k in seq_len(K)) {
    m <- model$models[[k]]
    pre <- m$preprocessor
    seed <- model_scoring_seed(model$member_seeds[k])
    z <- score_preprocessor_frozen(pre,X,batch,seed)
    cfit <- m$consensus_fit
    dfit <- m$diverse_fit
    cprob <- as.numeric(stats::predict(cfit$fit,
      newx=z[,model$feature_set,drop=FALSE], s=cfit$lambda, type="response"))
    dprob <- as.numeric(stats::predict(dfit$fit,
      newx=z[,m$diverse_ids,drop=FALSE], s=dfit$lambda, type="response"))
    if (length(cprob)!=n || length(dprob)!=n ||
        any(!is.finite(cprob) | cprob < 0 | cprob > 1) ||
        any(!is.finite(dprob) | dprob < 0 | dprob > 1))
      stop("Invalid predicted probabilities from frozen member ",k)
    con_member[,k] <- cprob
    div_member[,k] <- dprob
    j <- j+1L
    member_logs[[j]] <- data.frame(SampleID=rownames(X),
      strategy="consensus30",member_seed=model$member_seeds[k],
      test_MinProb_seed=seed,member_probability=cprob)
    j <- j+1L
    member_logs[[j]] <- data.frame(SampleID=rownames(X),
      strategy="diverse30",member_seed=model$member_seeds[k],
      test_MinProb_seed=seed,member_probability=dprob)
  }
  probs <- list(consensus30=rowMeans(con_member),diverse30=rowMeans(div_member))
  tbl <- model$inner_thresholds
  rows <- lapply(c("consensus30","diverse30"),function(strategy) {
    ix <- match(strategy,as.character(tbl$strategy))
    tau1 <- as.numeric(tbl$full_threshold[ix])
    tau2 <- as.numeric(tbl$spec90_threshold[ix])
    p <- rep(NA_real_,nrow(ss))
    p[scoreable] <- probs[[strategy]]
    status <- ifelse(scoreable,
      ifelse(in_training_batch,"SCORED_KNOWN_BATCH",
             "SCORED_UNSEEN_RUN1_NO_COMBAT"),
      "NOT_SCORED_UNSEEN_RUN1_COMBAT")
    data.frame(
      SampleID=ss$SampleID, genotype_group="LoF/LoF", batch=ss$batch,
      training_batch_seen=in_training_batch, scoring_status=status,
      strategy=strategy, number_of_models=K,
      SLC25A13_ACTB_raw_ratio=input$slc_ratio,
      SLC25A13_raw_quantified=input$slc_detected,
      number_of_raw_eligible_features=input$n_features,
      number_of_quantified_eligible_features=rowSums(is.finite(input$X)),
      probability=p,
      threshold_full_inner_recall=tau1,
      above_full_inner_recall=ifelse(is.finite(p),as.integer(p>=tau1),NA_integer_),
      threshold_inner_spec90=tau2,
      above_inner_spec90=ifelse(is.finite(p),as.integer(p>=tau2),NA_integer_),
      stringsAsFactors=FALSE)
  })
  list(scores=do.call(rbind,rows), member_scores=do.call(rbind,member_logs))
}

summarize_ss <- function(scores) {
  out <- list()
  for (s in c("consensus30","diverse30")) {
    for (b in c("ALL","Run1","Run2","Run3")) {
      x <- scores[scores$strategy==s & (b=="ALL" | scores$batch==b),,drop=FALSE]
      z <- x$probability[is.finite(x$probability)]
      out[[length(out)+1L]] <- data.frame(
        strategy=s, batch=b, n_SS_specimens=nrow(x),
        n_scored=length(z), n_unscored=sum(!is.finite(x$probability)),
        median_probability=if(length(z)) stats::median(z) else NA_real_,
        range_min=if(length(z)) min(z) else NA_real_,
        range_max=if(length(z)) max(z) else NA_real_,
        n_above_full_inner_recall=sum(x$above_full_inner_recall==1L,na.rm=TRUE),
        n_above_inner_spec90=sum(x$above_inner_spec90==1L,na.rm=TRUE),
        stringsAsFactors=FALSE)
    }
  }
  do.call(rbind,out)
}

main <- function() {
  root <- find_project()
  data_dir <- find_data_dir(root)
  model_file <- find_model(root)
  outdir <- file.path(root,OUTPUT_DIR_NAME)
  # root is the directory containing 20261029_rebuttal.Rproj.
  # Never create an output folder at the project root.
  dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
  cat("[PROJECT] ",root,"\n[DATA] ",data_dir,"\n[MODEL] ",model_file,"\n",sep="")
  model <- readRDS(model_file)
  validate_model(model)
  verify_analysis_plan_if_present(model,model_file,data_dir)
  input <- load_ss_features(data_dir)
  if (!all(model$feature_set %in% colnames(input$X)))
    stop("Original 30 consensus features absent from the eligible matrix")
  if (model$fit_mode=="no_combat") {
    warning("8/16 LoF/LoF specimens are from Run1, an UNSEEN batch. ",
            "They will be scored without any ComBat adjustment and reported ",
            "separately. Do not interpret Run1 scores as batch-validated.")
  } else {
    warning("Frozen ComBat was fitted only in Run2/Run3; Run1 (8 SS) ",
            "cannot be scored. Run1 rows will be NA with explicit status.")
  }
  results <- score_frozen_model(model,input)
  summary <- summarize_ss(results$scores)
  utils::write.csv(results$scores,
    file.path(outdir,"SS16_frozen_ensemble_scores.csv"),row.names=FALSE,na="")
  utils::write.csv(results$member_scores,
    file.path(outdir,"SS16_frozen_member_scores.csv"),row.names=FALSE,na="")
  utils::write.csv(summary,
    file.path(outdir,"SS16_frozen_scores_summary.csv"),row.names=FALSE,na="")
  provenance <- c(
    "LoF/LoF (SS) scoring using frozen final 44-non-SS full-training ensemble",
    paste0("RDS: ",model_file),
    paste0("RDS MD5: ",unname(tools::md5sum(model_file))),
    paste0("raw_data.xlsx MD5: ",unname(tools::md5sum(file.path(data_dir,"raw_data.xlsx")))),
    paste0("sample_info.csv MD5: ",unname(tools::md5sum(file.path(data_dir,"sample_info.csv")))),
    paste0("Model version: ",model$version),
    paste0("Engine version: ",model$engine_version),
    paste0("Model mode: ",model$fit_mode),
    paste0("Member seeds: ",paste(model$member_seeds,collapse=",")),
    "Test MinProb: specimen-specific q=0.01, model training sigma, deterministic member-specific seed",
    "Test MinProb RNG: v1.1 key_seed(fit,FINAL_ALL44/fit,member_seed)+4; SS_1..SS_16 order",
    "ACTB reference: P60709; log2(sample protein intensity / sample ACTB intensity)",
    "Model class labels are not used for scoring; SS group is only used to select test specimens",
    "Model coefficients, member feature lists, preprocessing stats, ComBat shifts and threshold are unchanged",
    "Run1 is unseen in training. If no_combat, Run1 numerical scores are extrapolation; if combat_mean, NA.",
    "Full inner recall threshold was calibrated on the 44 non-SS training samples, not on SS samples.",
    "No SS score is an independent validation of a citrin-severity or clinical diagnostic scale.",
    "Do not compare the 16 full-data-training-fit SS scores directly to 44 OOF scores as if one locked model scored both.",
    paste0("R version: ",R.version.string),
    paste0("glmnet version: ",as.character(utils::packageVersion("glmnet"))))
  writeLines(provenance,file.path(outdir,"SS16_scoring_methods_and_caveats.txt"))
  utils::capture.output(sessionInfo(),file=file.path(outdir,"R_sessionInfo.txt"))
  cat("\n[RESULTS]\n")
  print(summary,row.names=FALSE)
  cat("\nFiles saved under: ",outdir,"\n",sep="")
  invisible(results)
}

main()
