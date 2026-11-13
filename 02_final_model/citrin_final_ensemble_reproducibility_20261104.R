#!/usr/bin/env Rscript
# ============================================================================
# CITRIN rebuttal - frozen 10-member ensemble reproducibility
# Version 1.0 / 2026-11-04
#
# PURPOSE
#   Export the complete final, FULL-TRAINING fitted research ensemble and a
#   worked numerical scoring example. No training/refitting/feature selection,
#   threshold tuning, or new validation is performed.
#
# REQUIRED EXISTING OBJECT
#   citrin_nested_ensemble_v1_1/final_training_model_NOT_VALIDATED.rds
#   created by citrin_nested_ensemble_familyLOFO_v1_1.R, v1.1 + engine v2.2.
#   A separate frozen RDS for each outer fold is NOT the final full-data model.
#
# PROJECT AND OUTPUT (all under project; no new root-level result folder)
#   Open 20261029_rebuttal.Rproj in RStudio and run:
#     source("script/citrin_final_ensemble_reproducibility_20261104.R",
#            encoding="UTF-8")
#   Outputs: <project>/result/FinalModel_reproducibility/
#
# IMPORTANT SCIENTIFIC INTERPRETATION
#   - consensus30: ONE feature set, 10 independently fitted glmnet member
#     models; ensemble score = mean of TEN member logistic probabilities.
#   - diverse30: 10 possibly different feature sets, secondary/post-hoc.
#   - Inner-OOF training thresholds are fixed; different from outer-fold
#     thresholds in the 44-subject nested validation.
#   - Demonstration specimens SM_1 and WS_1 were IN THE TRAINING COHORT.
#     Worked examples are implementation checks, NOT validation or OOF scores.
#   - Test-side MinProb uses sample-specific 1% quantile, a FROZEN training SD,
#     and seeded stochastic draws; sample is scored ALONE (not jointly with
#     another specimen) to eliminate test-batch/order-dependent RNG draws.
#   - Run1 was not represented in model fitting (Run2/Run3 only). Under frozen
#     ComBat, an unseen batch is not scoreable. The examples use known batches.
# ============================================================================

options(stringsAsFactors = FALSE, digits = 17)

# --- User overrides: leave blank when using the standard project structure ---
ROOT_OVERRIDE <- normalizePath(Sys.getenv("CITRIN_REPO", "."))   # repository root (public release)
MODEL_RDS_OVERRIDE <- file.path(ROOT_OVERRIDE, "02_final_model", "FinalModel_reproducibility", "final_training_model_FROZEN_NOT_VALIDATED.rds")
DATA_DIR_OVERRIDE <- Sys.getenv("CITRIN_DATA_DIR", file.path(ROOT_OVERRIDE, "data"))
EXAMPLE_IDS <- c("SM_1", "WS_1")    # in-sample worked calculations only
OUTPUT_SUBDIR <- file.path("02_final_model", "results", "FinalModel_reproducibility_rerun")
COPY_FROZEN_MODEL_TO_OUTPUT <- TRUE  # exact-byte archive; no editing or overwriting
EXPECTED_PROJECT <- "20261029_rebuttal.Rproj"
CONSENSUS_SEEDS <- 2025:2034
SCORING_CONTEXT <- "FINAL_ALL44/fit"
RNG_MASTER_SEED <- 2025L
RNG_SALT <- 11111L
MIN_UNIQUE_PEPTIDES <- 3L
MINPROB_Q <- 0.01
CALC_TOL <- 1e-7

for (pkg in c("openxlsx", "glmnet")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Missing R package '", pkg, "'. Install with install.packages('", pkg, "').")
}

fail <- function(...) stop(paste0(...), call. = FALSE)
check <- function(cond, ...) if (!isTRUE(cond)) fail(...)
canon <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
norm_id <- function(x) toupper(trimws(as.character(x)))
md5 <- function(f) unname(tools::md5sum(f))
write_table <- function(d, out, name) {
  utils::write.csv(d, file.path(out, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
col_lookup <- function(d, choices) {
  hit <- match(canon(choices), canon(names(d)))
  hit <- hit[!is.na(hit)]
  if (!length(hit)) fail("Could not locate column among: ", paste(choices, collapse = ", "),
                        "; available: ", paste(names(d), collapse = ", "))
  names(d)[hit[1L]]
}

find_root <- function() {
  if (nzchar(ROOT_OVERRIDE)) {
    root <- normalizePath(ROOT_OVERRIDE, mustWork = TRUE)
    if (!file.exists(file.path(root, EXPECTED_PROJECT))) message("[NOTE] ", EXPECTED_PROJECT, " not found; running in the public-repository layout.")
    return(root)
  }
  candidates <- character(0)
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    rp <- tryCatch(rstudioapi::getActiveProject(), error = function(e) NULL)
    if (is.character(rp) && length(rp) == 1L && !is.na(rp) && nzchar(rp))
      candidates <- c(candidates, rp)
  }
  cur <- normalizePath(getwd(), mustWork = TRUE)
  for (i in 0:7) {
    candidates <- c(candidates, cur)
    nxt <- dirname(cur)
    if (identical(nxt, cur)) break
    cur <- nxt
  }
  candidates <- unique(candidates[dir.exists(candidates)])
  hits <- candidates[file.exists(file.path(candidates, EXPECTED_PROJECT))]
  if (!length(hits)) fail("Project not found: open ", EXPECTED_PROJECT,
                          " in RStudio or set ROOT_OVERRIDE.")
  normalizePath(hits[1L], winslash = "/", mustWork = TRUE)
}

find_data <- function(root) {
  candidates <- unique(c(DATA_DIR_OVERRIDE, Sys.getenv("CITRIN_DATA_DIR", ""),
    file.path(dirname(root), "data"), file.path(root, "data"),
    file.path(root, "R", "data"), file.path(dirname(dirname(root)), "data")))
  candidates <- candidates[nzchar(candidates)]
  hits <- candidates[vapply(candidates, function(p)
    file.exists(file.path(p, "raw_data.xlsx")) &&
      file.exists(file.path(p, "sample_info.csv")), logical(1))]
  if (!length(hits)) fail("Cannot find raw_data.xlsx and sample_info.csv. Set DATA_DIR_OVERRIDE.")
  normalizePath(hits[1L], winslash = "/", mustWork = TRUE)
}

find_model <- function(root) {
  candidates <- unique(c(MODEL_RDS_OVERRIDE, Sys.getenv("CITRIN_SS_MODEL_RDS", ""),
    file.path(root, "citrin_nested_ensemble_v1_1", "final_training_model_NOT_VALIDATED.rds"),
    file.path(root, "Claude", "citrin_nested_ensemble_v1_1_lofo", "final_training_model_NOT_VALIDATED.rds"),
    file.path(root, "result", "citrin_nested_ensemble_v1_1", "final_training_model_NOT_VALIDATED.rds"),
    file.path(root, "final_training_model_NOT_VALIDATED.rds")))
  candidates <- candidates[nzchar(candidates)]
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) fail(
    "Frozen FULL-TRAINING RDS not found. Expected final_training_model_NOT_VALIDATED.rds\n",
    "This script does not fit a replacement model. Set MODEL_RDS_OVERRIDE to its path,\n",
    "or generate the final fit with the original nested ensemble v1.1 pipeline.")
  normalizePath(hits[1L], winslash = "/", mustWork = TRUE)
}

# Same deterministic algorithm as v1.1 key_seed("fit","FINAL_ALL44/fit",member_seed)
# and engine v2.2 family_seed(..., salt=11111): test-side seed is +4.
scoring_seed <- function(member_seed) {
  key <- paste("fit", SCORING_CONTEXT, member_seed, sep = "|")
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(paste0(RNG_MASTER_SEED, "|", RNG_SALT, "|", key), tmp, useBytes = TRUE)
  h <- md5(tmp)
  parent <- as.integer((strtoi(substr(h, 1L, 7L), base = 16L) +
     as.double(RNG_MASTER_SEED) + as.double(RNG_SALT)) %% 2000000000) + 1L
  parent + 4L
}

extract_coeff <- function(fitted, features, label) {
  check(is.list(fitted) && inherits(fitted$fit, "glmnet") &&
        length(fitted$lambda) == 1L && is.finite(fitted$lambda) && fitted$lambda > 0,
        "Invalid glmnet fit/lambda in ", label)
  cm <- as.matrix(stats::coef(fitted$fit, s = fitted$lambda))
  check(ncol(cm) == 1L && "(Intercept)" %in% rownames(cm),
        "Unexpected coefficient format: ", label)
  check(setequal(setdiff(rownames(cm), "(Intercept)"), features),
        "Coefficient feature IDs differ from frozen panel in ", label)
  check(all(is.finite(cm)), "Nonfinite coefficients in ", label)
  co <- setNames(as.numeric(cm[, 1L]), rownames(cm))
  check(!anyDuplicated(names(co)), "Duplicate feature names in coefficients: ", label)
  co
}

validate_frozen <- function(m) {
  check(is.list(m) && identical(m$version, "citrin_nested_ensemble_familyLOFO_1.1") &&
        identical(m$engine_version, "two_run_fixed_panel_2.2"),
        "Unexpected RDS version: expected ensemble v1.1 and engine v2.2")
  check(identical(m$optimization_status, "TRAINING-DERIVED NOT EXTERNALLY VALIDATED"),
        "This is not the expected fully trained exploratory RDS")
  check(m$fit_mode %in% c("no_combat", "combat_mean") && length(m$fit_mode) == 1L,
        "Unexpected fit_mode")
  check(identical(as.integer(m$member_seeds), as.integer(CONSENSUS_SEEDS)) &&
        length(m$models) == 10L, "Expected exactly ten members with seeds 2025--2034")
  check(is.character(m$feature_set) && length(m$feature_set) == 30L &&
        !anyNA(m$feature_set) && !anyDuplicated(m$feature_set),
        "Invalid consensus-30 feature set")
  required_tau <- c("strategy", "full_threshold", "spec90_threshold")
  check(is.data.frame(m$inner_thresholds) &&
        all(required_tau %in% names(m$inner_thresholds)) &&
        setequal(m$inner_thresholds$strategy, c("consensus30", "diverse30")) &&
        nrow(m$inner_thresholds) == 2L,
        "Missing final inner-derived threshold table")
  for (v in c("full_threshold", "spec90_threshold"))
    check(all(is.finite(m$inner_thresholds[[v]])) &&
          all(m$inner_thresholds[[v]] >= 0) &&
          all(m$inner_thresholds[[v]] <= 1 + 1e-10), "Bad threshold: ", v)
  for (i in seq_along(m$models)) {
    mm <- m$models[[i]]
    pre <- mm$preprocessor
    check(is.list(pre) && identical(pre$mode, m$fit_mode),
          "Preprocessor mode mismatch for member ", i)
    check(is.character(pre$kept) && !anyDuplicated(pre$kept) &&
          length(pre$kept) >= 30L &&
          all(m$feature_set %in% pre$kept) &&
          length(mm$diverse_ids) == 30L && !anyDuplicated(mm$diverse_ids) &&
          all(mm$diverse_ids %in% pre$kept), "Member ", i, " has invalid panel")
    check(all(pre$kept %in% names(pre$mu)) &&
          all(pre$kept %in% names(pre$sd)) &&
          all(is.finite(pre$mu[pre$kept])) &&
          all(is.finite(pre$sd[pre$kept]) & pre$sd[pre$kept] > 0) &&
          length(pre$minprob_sigma) == 1L &&
          is.finite(pre$minprob_sigma) && pre$minprob_sigma > 0,
          "Invalid frozen preprocessing parameters in member ", i)
    if (m$fit_mode == "combat_mean")
      check(is.matrix(pre$shifts) &&
            all(c("Run2", "Run3") %in% colnames(pre$shifts)) &&
            all(pre$kept %in% rownames(pre$shifts)),
            "Frozen ComBat shifts missing in member ", i)
    extract_coeff(mm$consensus_fit, m$feature_set, paste0("consensus member ", i))
    extract_coeff(mm$diverse_fit, mm$diverse_ids, paste0("diverse member ", i))
  }
  invisible(TRUE)
}

verify_training_plan <- function(model, model_path, data_dir) {
  f <- file.path(dirname(model_path), "immutable_analysis_plan.csv")
  if (!file.exists(f)) {
    warning("immutable_analysis_plan.csv not adjacent to the RDS: training-source hashes cannot be verified")
    return(NULL)
  }
  p <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  check(all(c("name", "value") %in% names(p)) && !anyDuplicated(p$name),
        "Malformed immutable_analysis_plan.csv")
  kv <- stats::setNames(as.character(p$value), p$name)
  for (key in c("raw_md5", "sample_info_md5")) {
    check(key %in% names(kv), "Training plan missing checksum: ", key)
    fn <- if (key == "raw_md5") "raw_data.xlsx" else "sample_info.csv"
    check(identical(tolower(kv[[key]]), tolower(md5(file.path(data_dir, fn)))),
          "Training input does not match saved RDS training plan: ", fn)
  }
  check("mode" %in% names(kv) && identical(kv[["mode"]], model$fit_mode),
        "Training plan batch-correction mode mismatch")
  check("seeds" %in% names(kv) &&
        identical(kv[["seeds"]], paste(model$member_seeds, collapse = ",")),
        "Training plan seed mismatch")
  p
}

load_worked_examples <- function(datadir, ids) {
  check(length(ids) >= 1L && !anyDuplicated(ids), "Invalid EXAMPLE_IDS")
  si <- utils::read.csv(file.path(datadir, "sample_info.csv"),
                        check.names = FALSE, fileEncoding = "UTF-8-BOM")
  sid_col <- col_lookup(si, c("SampleID", "sample_id", "sample id"))
  group_col <- col_lookup(si, "group")
  batch_col <- col_lookup(si, "batch")
  si$SampleID <- norm_id(si[[sid_col]])
  si$group <- norm_id(si[[group_col]])
  si$batch <- trimws(as.character(si[[batch_col]]))
  check(!anyDuplicated(si$SampleID) && !anyNA(si$SampleID), "Invalid sample_info.csv IDs")
  ix <- match(norm_id(ids), si$SampleID)
  check(!anyNA(ix), "Example specimen ID absent from sample_info.csv")
  meta <- si[ix, c("SampleID", "group", "batch"), drop = FALSE]
  check(all(meta$group != "SS") && all(meta$batch %in% c("Run2", "Run3")),
        "Worked examples must be from training groups / known batches")
  raw <- openxlsx::read.xlsx(file.path(datadir, "raw_data.xlsx"), check.names = FALSE)
  acc <- norm_id(raw[[col_lookup(raw, c("Master Accession Number", "master_accession_number"))]])
  check(!anyNA(acc) && all(nzchar(acc)), "Invalid accession IDs in raw_data.xlsx")
  if (anyDuplicated(acc)) {
    keep <- !duplicated(acc) # same choice as v2.2 original engine
    raw <- raw[keep, , drop = FALSE]
    acc <- acc[keep]
  }
  check("P60709" %in% acc, "ACTB P60709 not in raw_data.xlsx")
  colidx <- match(meta$SampleID, norm_id(names(raw)))
  check(!anyNA(colidx), "Example specimen not found in raw_data.xlsx")
  vv <- as.matrix(raw[, colidx, drop = FALSE])
  suppressWarnings(storage.mode(vv) <- "numeric")
  colnames(vv) <- meta$SampleID
  rownames(vv) <- acc
  actb <- as.numeric(vv["P60709", ])
  check(all(is.finite(actb) & actb > 0), "Invalid ACTB intensity in examples")
  vv[!is.finite(vv) | vv <= 0] <- NA_real_
  lg <- sweep(log2(vv), 2L, log2(actb), "-")
  peps <- suppressWarnings(as.numeric(raw[[col_lookup(raw, c("Unique peptides", "unique_peptides"))]]))
  elig <- is.finite(peps) & peps >= MIN_UNIQUE_PEPTIDES
  check(sum(elig) == 5199L, "Protein universe differs from original 5,199 features")
  xx <- t(lg[elig, , drop = FALSE])
  check(identical(rownames(xx), meta$SampleID) &&
        !anyDuplicated(colnames(xx)), "Example intensity matrix not aligned")
  list(X = xx, meta = meta, ACTB = stats::setNames(actb, meta$SampleID), n_features = ncol(xx))
}

# One test sample at a time: exact v2.2 apply_minprob() + score_preprocessor().
# Randomness is reproducible for a fixed member and that single sample.
score_one_preprocessor <- function(pre, raw_logratio, batch, seed) {
  check(is.matrix(raw_logratio) && nrow(raw_logratio) == 1L &&
        length(batch) == 1L, "One-sample scoring expects 1 x features matrix")
  obs <- as.numeric(raw_logratio[1L, is.finite(raw_logratio[1L, ])])
  check(length(obs) >= 10L, "Example has <10 observed proteins")
  q <- unname(stats::quantile(obs, probs = MINPROB_Q,
                              na.rm = TRUE, names = FALSE))
  ximp <- raw_logratio
  misses <- which(!is.finite(ximp[1L, ]))
  set.seed(seed)
  if (length(misses)) ximp[1L, misses] <- stats::rnorm(length(misses), q, pre$minprob_sigma)
  check(all(is.finite(ximp)), "Imputation left nonfinite values")
  after_batch <- ximp
  shifts <- stats::setNames(rep(0, ncol(ximp)), colnames(ximp))
  if (pre$mode == "combat_mean") {
    check(batch %in% colnames(pre$shifts),
          "Cannot use frozen ComBat for unseen batch ", batch)
    check(all(colnames(ximp) %in% rownames(pre$shifts)),
          "Frozen shift lacks raw protein features")
    shifts <- pre$shifts[colnames(ximp), batch]
    after_batch <- sweep(after_batch, 2L, shifts, "-")
  }
  check(all(pre$kept %in% colnames(after_batch)), "Frozen preprocessor feature mismatch")
  scaled <- after_batch[, pre$kept, drop = FALSE]
  scaled <- sweep(sweep(scaled, 2L, pre$mu[pre$kept], "-"),
                  2L, pre$sd[pre$kept], "/")
  check(all(is.finite(scaled)), "Nonfinite standardized features")
  list(z = scaled, raw = raw_logratio, imputed = ximp,
       after_batch = after_batch, shifts = shifts,
       q01 = q, missing_count = length(misses), observed_count = length(obs))
}

make_model_tables <- function(m) {
  allco <- list(); allpar <- list(); allpre <- list(); allbat <- list();
  featdiv <- list(); kco <- 0L; kpar <- 0L
  for (k in seq_along(m$models)) {
    seed <- as.integer(m$member_seeds[k]); mm <- m$models[[k]]; pre <- mm$preprocessor
    allpre[[k]] <- data.frame(member_seed = seed, FeatureID = pre$kept,
      center_mu = unname(pre$mu[pre$kept]), scale_sd = unname(pre$sd[pre$kept]),
      minprob_training_sigma = as.numeric(pre$minprob_sigma),
      batch_correction = pre$mode)
    if (m$fit_mode == "combat_mean") {
      allbat[[k]] <- data.frame(member_seed = seed,
        FeatureID = rep(rownames(pre$shifts), times = ncol(pre$shifts)),
        batch = rep(colnames(pre$shifts), each = nrow(pre$shifts)),
        frozen_shift = as.vector(pre$shifts))
    }
    featdiv[[k]] <- data.frame(member_seed = seed, panel_rank = seq_along(mm$diverse_ids),
                              FeatureID = mm$diverse_ids)
    for (strategy in c("consensus30", "diverse30")) {
      ob <- if (strategy == "consensus30") mm$consensus_fit else mm$diverse_fit
      ids <- if (strategy == "consensus30") m$feature_set else mm$diverse_ids
      co <- extract_coeff(ob, ids, paste(strategy, seed))
      kco <- kco + 1L
      allco[[kco]] <- data.frame(member_seed = seed, strategy = strategy,
         FeatureID = names(co), is_intercept = names(co) == "(Intercept)",
         coefficient = unname(co), nonzero = unname(co) != 0,
         lambda = as.numeric(ob$lambda),
         lambda_method = if (!is.null(ob$lambda_method)) as.character(ob$lambda_method) else NA_character_)
      kpar <- kpar + 1L
      allpar[[kpar]] <- data.frame(member_seed = seed, strategy = strategy,
        lambda = as.numeric(ob$lambda),
        lambda_method = if (!is.null(ob$lambda_method)) as.character(ob$lambda_method) else NA_character_,
        intercept = unname(co["(Intercept)"]),
        n_features_in_panel = length(ids),
        n_nonzero_coefficients = sum(co[names(co) != "(Intercept)"] != 0),
        minprob_training_sigma = as.numeric(pre$minprob_sigma),
        training_preprocessor_kept_features = length(pre$kept),
        seed_for_single_specimen_score = scoring_seed(seed),
        fit_mode = m$fit_mode)
    }
  }
  list(coeff = do.call(rbind, allco),
       params = do.call(rbind, allpar),
       preprocess = do.call(rbind, allpre),
       shifts = if (length(allbat)) do.call(rbind, allbat) else NULL,
       diverse_features = do.call(rbind, featdiv))
}

score_worked_examples <- function(m, ex) {
  metrics <- list(); member_rows <- list(); contrib_rows <- list(); summaries <- list()
  cidx <- 0L; midx <- 0L; sidx <- 0L
  for (i in seq_len(nrow(ex$meta))) {
    id <- ex$meta$SampleID[i]
    group <- ex$meta$group[i]
    batch <- ex$meta$batch[i]
    raw_single <- ex$X[id, , drop = FALSE]
    for (k in seq_along(m$models)) {
      member_seed <- as.integer(m$member_seeds[k])
      mm <- m$models[[k]]
      pre <- mm$preprocessor
      seed <- scoring_seed(member_seed)
      prep <- score_one_preprocessor(pre, raw_single, batch, seed)
      midx <- midx + 1L
      metrics[[midx]] <- data.frame(SampleID = id, group = group, batch = batch,
        member_seed = member_seed, seed_for_imputation = seed,
        observed_eligible_features = prep$observed_count,
        imputed_eligible_features = prep$missing_count,
        specimen_log2_1pct_quantile = prep$q01,
        training_minprob_sigma = pre$minprob_sigma)
      for (strategy in c("consensus30", "diverse30")) {
        obj <- if (strategy == "consensus30") mm$consensus_fit else mm$diverse_fit
        ids <- if (strategy == "consensus30") m$feature_set else mm$diverse_ids
        co <- extract_coeff(obj, ids, paste(strategy, member_seed))
        zz <- prep$z[, ids, drop = FALSE]
        beta <- co[ids]
        terms <- as.numeric(zz[1L, ]) * unname(beta)
        eta <- unname(co["(Intercept)"]) + sum(terms)
        p_manual <- stats::plogis(eta)
        p_glmnet <- as.numeric(stats::predict(obj$fit, newx = zz,
                                              s = obj$lambda, type = "response"))
        delta <- abs(p_manual - p_glmnet)
        check(is.finite(delta) && delta < CALC_TOL,
              "Manual score differs from glmnet for sample ", id,
              " strategy ", strategy, " seed ", member_seed,
              ": delta=", format(delta, digits=12))
        sidx <- sidx + 1L
        member_rows[[sidx]] <- data.frame(SampleID=id, group=group,
          batch=batch, strategy=strategy, member_seed=member_seed,
          selected_lambda=obj$lambda,
          intercept=unname(co["(Intercept)"]),
          sum_beta_z=sum(terms), logit_eta=eta,
          probability_manual=p_manual, probability_glmnet=p_glmnet,
          absolute_difference=delta,
          test_minprob_seed=seed)
        cidx <- cidx + 1L
        contrib_rows[[cidx]] <- data.frame(SampleID=id, group=group,
          batch=batch, strategy=strategy, member_seed=member_seed,
          FeatureID=ids,
          raw_log2_protein_over_ACTB=as.numeric(prep$raw[1L, ids]),
          post_MinProb_log2_ratio=as.numeric(prep$imputed[1L, ids]),
          frozen_batch_mean_shift=unname(prep$shifts[ids]),
          post_batch_log2_ratio=as.numeric(prep$after_batch[1L, ids]),
          training_mu=unname(pre$mu[ids]),
          training_sd=unname(pre$sd[ids]),
          standardized_value_z=as.numeric(zz[1L, ]),
          beta=unname(beta), beta_times_z=terms)
      }
    }
  }
  member <- do.call(rbind, member_rows)
  for (id in ex$meta$SampleID) for (st in c("consensus30", "diverse30")) {
    sub <- member[member$SampleID == id & member$strategy == st, , drop = FALSE]
    check(nrow(sub) == 10L, "Expected ten scoring members for ", id, "/", st)
    ix <- match(st, m$inner_thresholds$strategy)
    tau <- as.numeric(m$inner_thresholds$full_threshold[ix])
    alt <- as.numeric(m$inner_thresholds$spec90_threshold[ix])
    p <- mean(sub$probability_glmnet)
    check(abs(p - mean(sub$probability_manual)) < CALC_TOL,
          "Ensemble manual vs glmnet discrepancy for ", id)
    summaries[[length(summaries)+1L]] <- data.frame(SampleID=id,
      genotype_group=ex$meta$group[match(id, ex$meta$SampleID)],
      batch=ex$meta$batch[match(id, ex$meta$SampleID)],
      strategy=st, n_members=nrow(sub), ensemble_probability=p,
      full_inner_recall_threshold=tau,
      predicted_full_inner_recall=as.integer(p >= tau),
      spec90_threshold=alt,
      predicted_spec90=as.integer(p >= alt),
      max_manual_vs_glmnet_member_difference=max(sub$absolute_difference),
      interpretation="IN-SAMPLE WORKED EXAMPLE ONLY - NOT OOF OR VALIDATION")
  }
  list(summary=do.call(rbind, summaries),
       members=member, contributions=do.call(rbind, contrib_rows),
       imputation=do.call(rbind, metrics))
}

verify_written_examples <- function(out, tol = CALC_TOL) {
  m <- utils::read.csv(file.path(out, "scoring_example_member_calculations.csv"),
                       check.names=FALSE)
  contributions <- utils::read.csv(file.path(out, "scoring_example_feature_contributions.csv"),
                       check.names=FALSE)
  keys <- unique(m[, c("SampleID", "strategy", "member_seed")])
  for (j in seq_len(nrow(keys))) {
    k <- keys[j, ]
    r <- m[m$SampleID == k$SampleID & m$strategy == k$strategy &
             m$member_seed == k$member_seed, , drop = FALSE]
    f <- contributions[contributions$SampleID == k$SampleID &
             contributions$strategy == k$strategy &
             contributions$member_seed == k$member_seed, , drop = FALSE]
    check(nrow(r)==1L && nrow(f)==30L, "Exported scoring example incomplete")
    # Independent numeric reconstruction using only the OUTPUT CSV fields.
    z <- (f$post_batch_log2_ratio - f$training_mu) / f$training_sd
    check(max(abs(z - f$standardized_value_z)) < tol,
          "Standardization does not reconstruct from exported parameters")
    reconstructed_p <- stats::plogis(r$intercept + sum(z * f$beta))
    check(abs(reconstructed_p - r$probability_glmnet) < tol,
          "Frozen exported coefficient/scaling calculations disagree")
  }
  invisible(TRUE)
}

main <- function() {
  root <- find_root()
  datadir <- find_data(root)
  modelpath <- find_model(root)
  out <- file.path(root, OUTPUT_SUBDIR)
  # The only new folder is directly underneath the existing project/result/.
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)   # public-repository layout
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  model <- readRDS(modelpath)
  validate_frozen(model)
  plan <- verify_training_plan(model, modelpath, datadir)
  if (COPY_FROZEN_MODEL_TO_OUTPUT) {
    archive <- file.path(out, "final_training_model_FROZEN_NOT_VALIDATED.rds")
    if (file.exists(archive)) {
      check(identical(tolower(md5(archive)), tolower(md5(modelpath))),
            "Existing archived model has a DIFFERENT MD5. Refusing to overwrite: ", archive)
    } else {
      check(isTRUE(file.copy(modelpath, archive, overwrite = FALSE)),
            "Could not copy frozen model into the reproducibility folder")
      check(identical(tolower(md5(archive)), tolower(md5(modelpath))),
            "Copied frozen RDS checksum mismatch")
    }
  }
  ex <- load_worked_examples(datadir, EXAMPLE_IDS)
  check(all(model$feature_set %in% colnames(ex$X)),
        "Consensus feature set differs from original raw protein matrix")
  model_tables <- make_model_tables(model)
  worked <- score_worked_examples(model, ex)
  threshold <- model$inner_thresholds

  # Check frozen RDS versus published final training CSV snapshots, when present.
  f_features <- file.path(dirname(modelpath),
                          "final_consensus30_proteins_NOT_INDEPENDENTLY_VALIDATED.csv")
  if (file.exists(f_features)) {
    frozen_features <- utils::read.csv(f_features, check.names = FALSE)
    check("FeatureID" %in% names(frozen_features) &&
          identical(as.character(frozen_features$FeatureID), model$feature_set),
          "Selected feature CSV does NOT match frozen final model RDS")
  }
  f_threshold <- file.path(dirname(modelpath),
                           "final_training_inner_thresholds_NOT_VALIDATION.csv")
  if (file.exists(f_threshold)) {
    frozen_t <- utils::read.csv(f_threshold, check.names = FALSE)
    check(all(c("strategy", "full_threshold", "spec90_threshold") %in% names(frozen_t)),
          "Malformed existing threshold snapshot")
    for (st in c("consensus30", "diverse30")) {
      a <- threshold[threshold$strategy == st, , drop = FALSE]
      b <- frozen_t[frozen_t$strategy == st, , drop = FALSE]
      check(nrow(b) == 1L && all(abs(as.numeric(unlist(a[c("full_threshold", "spec90_threshold")], use.names=FALSE)) -
                  as.numeric(unlist(b[c("full_threshold", "spec90_threshold")], use.names=FALSE))) < 1e-12),
            "Threshold snapshot mismatch for ", st)
    }
  }

  write_table(data.frame(rank=seq_along(model$feature_set),
                        FeatureID=model$feature_set),out,"consensus30_frozen_protein_list.csv")
  write_table(model_tables$diverse_features,out,"diverse30_per_member_protein_lists.csv")
  write_table(model_tables$coeff,out,"member_glmnet_coefficients_ALL_10.csv")
  write_table(model_tables$params,out,"member_model_parameters.csv")
  write_table(model_tables$preprocess,out,"member_preprocessing_center_scale_ALL.csv")
  if (!is.null(model_tables$shifts))
    write_table(model_tables$shifts,out,"frozen_ComBat_batch_shifts_ALL.csv")
  write_table(threshold,out,"final_training_inner_thresholds.csv")
  if (!is.null(plan)) write_table(plan,out,"immutable_training_plan_snapshot.csv")
  write_table(worked$imputation,out,"scoring_example_imputation_quantiles.csv")
  write_table(worked$members,out,"scoring_example_member_calculations.csv")
  write_table(worked$contributions,out,"scoring_example_feature_contributions.csv")
  write_table(worked$summary,out,"scoring_example_ensemble_summary.csv")
  verify_written_examples(out)

  manifest <- data.frame(key=c("model_rds_path", "model_rds_md5", "raw_data_md5",
      "sample_info_md5", "ensemble_version", "engine_version", "batch_correction",
      "consensus_n_proteins", "members", "member_seeds", "training_n_samples",
      "training_n_positive", "fit_status", "test_minprob_q", "full_training_context",
      "worked_example_ids", "worked_examples_are_oof", "frozen_model_copy_in_output",
      "manual_glmnet_max_difference",
      "training_input_checksums_verified"),
    value=as.character(c(modelpath, md5(modelpath), md5(file.path(datadir,"raw_data.xlsx")),
      md5(file.path(datadir,"sample_info.csv")), model$version, model$engine_version,
      model$fit_mode, length(model$feature_set), length(model$member_seeds),
      paste(model$member_seeds,collapse=","), 44L, 7L,
      model$optimization_status, MINPROB_Q, SCORING_CONTEXT,
      paste(EXAMPLE_IDS,collapse=","), "FALSE", COPY_FROZEN_MODEL_TO_OUTPUT,
      max(worked$members$absolute_difference),
      !is.null(plan))))
  write_table(manifest,out,"model_manifest.csv")
  methods <- c(
    "CITRIN 10-member full-training frozen ensemble: model reproducibility audit",
    paste0("Project root: ",root), paste0("Output folder: ",out),
    paste0("Frozen model path: ",modelpath),
    paste0("Frozen RDS MD5: ",md5(modelpath)),
    paste0("Model version: ",model$version), paste0("Engine version: ",model$engine_version),
    paste0("Fit mode: ",model$fit_mode),
    "The frozen RDS is trained on 44 non-SS specimens (7 LoF/Mis, 37 negatives) from Run2 and Run3.",
    "The output is a FINAL FULL-DATA TRAINING FIT, NOT an independently validated external model.",
    "Original submitted refit: 30 candidates/15 nonzero coefficients; DIFFERENT from revised 10-member ensemble.",
    "Consensus30: shared list of 30 UniProt features; each of ten members has its own fitted glmnet model.",
    "Diverse30: each member uses its own 30 proteins; optional post-hoc comparator.",
    "Scores are the MEAN of 10 member probabilities, NOT the sigmoid of mean coefficients.",
    "One member: raw intensity -> exclude nonpositive -> log2(protein/ACTB P60709).",
    "Protein universe at input: >=3 unique peptides, 5,199 candidates; sample 1% quantile based on ALL eligible observed proteins.",
    "Test MinProb: set.seed(member_seed_hash+4), for each missing eligible feature draw N(sample q01, frozen training minprob_sigma^2).",
    "The hash is v2.2 family_seed('fit|FINAL_ALL44/fit|member_seed',salt=11111); base seed 2025.",
    "A specimen is scored alone so each example's random draws do not depend on other example rows or order.",
    "When combat_mean, subtract saved protein- and batch-specific FROZEN shifts before member center/scale.",
    "ComBat frozen mean shifts are available only for observed training batches Run2 and Run3; Run1 cannot be scored under frozen ComBat.",
    "For each member, standardized z_j = (post_batch_log2_ratio_j - mu_j) / sd_j.",
    "Member logit eta_k = intercept_k + sum_{j in member panel}(beta_jk * z_jk).",
    "Member probability p_k = 1/(1+exp(-eta_k)). Ensemble probability = mean_k(p_k).",
    "The primary fixed cutoff is full_threshold from the 44-subject training's INNER family-LOFO predictions.",
    "The secondary cutoff is spec90_threshold; neither is re-estimated using the worked example.",
    "The final-training thresholds are DIFFERENT from outer-fold thresholds used to evaluate 44 OOF subjects.",
    "SM_1 and WS_1 are in the 44-subject TRAINING COHORT. Example scores are IN-SAMPLE calculation demonstrations.",
    "Do not report scoring_example_ensemble_summary.csv as model accuracy or cross-validated predictions.",
    "Missing-value draws are stochastic under the original pipeline; identical reruns require identical R RNG, raw feature order and seed.",
    paste0("Manual-vs-glmnet max absolute discrepancy: ",
           format(max(worked$members$absolute_difference),digits=12)),
    "All exported feature contributions were re-read and recalculated from their CSVs; checks passed.",
    paste0("R: ",R.version.string),
    paste0("glmnet: ",as.character(utils::packageVersion("glmnet"))),
    paste0("openxlsx: ",as.character(utils::packageVersion("openxlsx")))
  )
  writeLines(methods,file.path(out,"scoring_rule_methods_and_caveats.txt"),useBytes=TRUE)
  utils::capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
  cat("\n[FROZEN MODEL] ",model$version," / ",model$fit_mode,"\n",sep="")
  cat("[SCORING EXAMPLES] in-sample demonstration, NOT validation:\n")
  print(worked$summary,row.names=FALSE)
  cat("[PASS] All ten-member manual calculations agree with glmnet predictions.\n")
  cat("[PASS] Exported CSVs independently reproduce the worked example calculations.\n")
  cat("[OUTPUT] ",out,"\n",sep="")
  invisible(list(model=model,examples=worked,output=out))
}

main()
