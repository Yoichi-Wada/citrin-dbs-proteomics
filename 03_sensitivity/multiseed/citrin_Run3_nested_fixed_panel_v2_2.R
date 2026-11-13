#!/usr/bin/env Rscript
# ============================================================================
# CITRIN / RUN 3 ONLY: NESTED FAMILY-GROUPED VALIDATION (R4, 2026-10-29)
# ============================================================================
# Purpose: validate 7 SM (LoF/Mis) vs 28 non-SS Run3 controls, using the
# original manuscript's model DEVELOPMENT specification as closely as feasible.
#
# PRESERVED: raw_data.xlsx; ACTB log2 ratio; >=3 unique peptides;
# MinProb q=.01; elastic net alpha=.5; inverse-prevalence positive weighting;
# B=30 class-stratified 80% stability subsamples; top-1000 variance shortlist;
# three-fold AUC criterion for lambda; candidate m=5,10,15,20,30,40,50;
# threshold maximizes sensitivity subject to specificity>=.90.
#
# NECESSARY CHANGES: (i) Run3-only -> no cross-run ComBat; (ii) ALL adaptive
# steps must be re-fitted within training folds; (iii) all related participants
# withheld together; (iv) model size/threshold tuned using inner grouped LOOCV.
#
# SMALL-SAMPLE SAFEGUARDS (reported as diagnostics, not concealed):
# - At least three positive observations are retained in a stability subsample
#   when available; with 3 positives this uses 100%, not 80% of positives.
# - Preserve original cv.glmnet(type.measure="auc") for >=30-sample fits
#   where grouped 3-fold AUC is feasible; otherwise use explicit pooled 3-fold
#   out-of-fold AUC rather than silently switching to deviance.
# - If NO valid three-fold family-grouped AUC tuning is possible, lambda is
#   selected by an explicitly disclosed fallback: 0.10 * training lambda_max.
#   Check diagnostic counts; this is not a completely unchanged algorithm.
#
# NO CLASS LABEL IS EVER USED TO PREPROCESS A HELD-OUT FAMILY.
# A held-out specimen's own protein intensities provide its MinProb q=.01
# sample location; imputation SD is estimated ONLY from training specimens.
#
# ANALYTIC LIMITATIONS: upstream DIA-NN cross-run normalization, spectral
# quantification, and the >=3 peptide count from original software are not
# reconstructed within CV; new Run3-only model is not externally validated.
#
# INPUTS: R/data/raw_data.xlsx and R/data/sample_info.csv; family_id column in
# sample_info.csv OR R/data/family_map.csv OR edited v3 template.
# OUTPUTS: R/[current project]/result/run3_original_nested_cv/<run_key>/*
# CACHES: saved per outer family, making interrupted runs restartable.
# ============================================================================

# IMPORTANT INFERENCE NOTE:
# fixed15 evaluates a PROCEDURE that re-selects 15 protein features within
# each training fold; it does NOT externally validate the exact 15 proteins
# originally chosen using all 44 samples. In the submitted manuscript,
# the 15 nonzero coefficients came from refitting a 30-CANDIDATE panel.
# Thus fixed30 is closer to the original model development procedure.
# This analysis was devised after inspecting the outcome data: treat the
# fixed-panel comparison as exploratory, not a prospective confirmation.

SCRIPT_VERSION <- "run3_fixed_panel_original_1.2"
SEED <- as.integer(Sys.getenv("CITRIN_SEED", "2025"))
SUBSAMPLE_UNIT <- tolower(Sys.getenv("CITRIN_SUBSAMPLE_UNIT", "sample"))
if (!SUBSAMPLE_UNIT %in% c("sample", "family")) stop("CITRIN_SUBSAMPLE_UNIT: use sample or family")
POS <- "SM"; EXCLUDE <- "SS"; RUN3 <- Sys.getenv("CITRIN_RUN3", "Run3")
ALPHA <- 0.5; MINPROB_Q <- 0.01; MIN_UNIQUE_PEPTIDES <- 3L
B_STABILITY <- as.integer(Sys.getenv("CITRIN_B", "30"))
STABILITY_FRAC <- 0.80; N_VARIANCE <- 1000L
PANEL_GRID <- c(5L, 10L, 15L, 20L, 30L, 40L, 50L)
SPEC_FLOOR <- 0.90
PANEL_MODE <- tolower(Sys.getenv("CITRIN_PANEL_MODE", "select"))
INNER_MODE <- tolower(Sys.getenv("CITRIN_INNER_MODE", "lofo"))
LAMBDA_FALLBACK_RATIO <- 0.10
N_BOOT <- as.integer(Sys.getenv("CITRIN_N_BOOT", "2000"))
OUTER_LIMIT <- as.integer(Sys.getenv("CITRIN_OUTER_LIMIT", "0"))
REQUIRE_35 <- tolower(Sys.getenv("CITRIN_EXPECT_35", "true")) == "true"

if (!is.finite(SEED) || !is.finite(B_STABILITY) || B_STABILITY < 1L ||
    !is.finite(N_BOOT) || N_BOOT < 0L || !is.finite(OUTER_LIMIT))
  stop("Invalid CITRIN_SEED, CITRIN_B, CITRIN_N_BOOT, or CITRIN_OUTER_LIMIT")
if (!PANEL_MODE %in% c("select", "fixed30", "fixed15"))
  stop("CITRIN_PANEL_MODE must be 'select', 'fixed30' or 'fixed15'")
if (!INNER_MODE %in% c("lofo", "group3"))
  stop("CITRIN_INNER_MODE must be 'lofo' or 'group3'")

options(stringsAsFactors = FALSE)
needed <- c("openxlsx", "imputeLCMD", "glmnet")
absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent)) stop("Install required R packages: install.packages(c(",
                         paste(sprintf("'%s'", absent), collapse = ", "), "))")

# Stable keyed seed: independent of fold positions / numeric family labels.
# Uses MD5 of the complete key (not a weighted character sum).
family_seed <- function(key, salt = 0L) {
  temp <- tempfile(); on.exit(unlink(temp), add = TRUE)
  writeLines(paste0(SEED, "|", salt, "|", as.character(key)), temp, useBytes = TRUE)
  h <- unname(tools::md5sum(temp))
  as.integer((strtoi(substr(h, 1L, 7L), base = 16L) +
                as.double(SEED) + as.double(salt)) %% 2000000000) + 1L
}

# Refit coefficients including intercept; return stable accession names.
extract_outer_coefficients <- function(fit, ids, fam, freq) {
  cm <- as.matrix(stats::coef(fit$fit, s = fit$lambda))
  vals <- as.numeric(cm[, 1L])
  nm <- rownames(cm)
  data.frame(held_family = fam, FeatureID = nm,
             coefficient = vals, nonzero = is.finite(vals) & abs(vals) > 0,
             is_intercept = nm == "(Intercept)",
             stability_frequency = c(NA_real_, as.numeric(freq[nm[-1L]])),
             stringsAsFactors = FALSE)
}

normalize_id <- function(x) toupper(trimws(as.character(x)))
canon <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
find_col <- function(df, candidates, required = TRUE) {
  idx <- which(canon(names(df)) %in% canon(candidates))
  if (!length(idx)) {
    if (!required) return(NULL)
    stop("Missing column among [", paste(candidates, collapse = ", "),
         "]. Available: ", paste(names(df), collapse = ", "))
  }
  names(df)[idx[1L]]
}
null_family <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | toupper(x) %in% c("NA", "NULL", "N/A")] <- NA_character_
  x
}
canon_batch <- function(x) {
  x <- canon(x)
  ifelse(grepl("^[0-9]+$", x), paste0("run", x), x)
}
find_input_dir <- function() {
  env <- Sys.getenv("CITRIN_DATA_DIR", "")
  if (nzchar(env)) return(env)
  base <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  paths <- character(0)
  for (k in 0:4) {
    paths <- c(paths, file.path(base, "data"), file.path(base, "R", "data"))
    nextbase <- dirname(base)
    if (identical(nextbase, base)) break
    base <- nextbase
  }
  paths <- unique(paths)
  hits <- paths[file.exists(file.path(paths, "raw_data.xlsx")) &
                  file.exists(file.path(paths, "sample_info.csv"))]
  if (!length(hits)) stop("Cannot locate raw_data.xlsx + sample_info.csv. Set ",
                           "Sys.setenv(CITRIN_DATA_DIR='/full/path/to/R/data').")
  hits[1L]
}
DATA_DIR <- find_input_dir()
RAW_FILE <- file.path(DATA_DIR, "raw_data.xlsx")
SAMPLE_FILE <- file.path(DATA_DIR, "sample_info.csv")
OUT_BASE <- file.path(getwd(), "result", "run3_original_nested_cv")

# Resolve family-ID sources WITHOUT creating or overwriting ANY input file.
read_metadata <- function() {
  si <- read.csv(SAMPLE_FILE, check.names = FALSE, stringsAsFactors = FALSE,
                 fileEncoding = "UTF-8-BOM")
  col_sid <- find_col(si, c("SampleID", "sample_id", "sample id"))
  col_group <- find_col(si, "group")
  col_batch <- find_col(si, "batch")
  si$SampleID <- normalize_id(si[[col_sid]])
  si$group <- toupper(trimws(as.character(si[[col_group]])))
  si$batch <- trimws(as.character(si[[col_batch]]))
  if (anyNA(si$SampleID) || any(si$SampleID == "") || anyDuplicated(si$SampleID))
    stop("Invalid/duplicated SampleID in sample_info.csv")
  if (anyNA(si$group) || any(si$group == "") ||
      anyNA(si$batch) || any(si$batch == ""))
    stop("Missing group or batch in sample_info.csv")
  fcol <- find_col(si, c("family_id", "familyid", "family"), FALSE)
  fid <- if (is.null(fcol)) rep(NA_character_, nrow(si)) else null_family(si[[fcol]])
  source <- ifelse(is.na(fid), "missing", "sample_info.csv")
  family_candidates <- unique(c(
    Sys.getenv("CITRIN_FAMILY_FILE", ""),
    file.path(DATA_DIR, "family_map.csv"),
    file.path(getwd(), "result", "nested_validation", "no_combat_primary", "family_map_TEMPLATE.csv"),
    file.path(getwd(), "..", "result", "nested_validation", "no_combat_primary", "family_map_TEMPLATE.csv"),
    file.path(getwd(), "result", "nested_validation", "batch_covariate_sensitivity", "family_map_TEMPLATE.csv"),
    file.path(DATA_DIR, "family_map_TEMPLATE.csv")
  ))
  family_candidates <- family_candidates[nzchar(family_candidates)]
  for (fp in family_candidates) {
    if (!file.exists(fp)) next
    m <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE,
                  colClasses = "character", fileEncoding = "UTF-8-BOM")
    isid <- find_col(m, c("SampleID", "sample_id", "sample id"))
    ifam <- find_col(m, c("family_id", "familyid", "family"))
    key <- normalize_id(m[[isid]])
    if (anyDuplicated(key) || anyNA(key) || any(key == ""))
      stop("Invalid/duplicated SampleID in ", fp)
    val <- null_family(m[[ifam]][match(si$SampleID, key)])
    conflict <- !is.na(fid) & !is.na(val) & fid != val
    if (any(conflict)) stop("Conflicting family IDs for: ",
                            paste(si$SampleID[conflict], collapse = ", "),
                            "; sources include ", fp)
    fill <- is.na(fid) & !is.na(val)
    fid[fill] <- val[fill]
    source[fill] <- basename(fp)
  }
  si$family_id <- fid
  rownames(si) <- si$SampleID
  message("[INFO] Read metadata for ", nrow(si), " specimens. Family sources: ",
          paste(names(table(source)), as.integer(table(source)), collapse = ", "))
  si
}

read_run3 <- function() {
  si <- read_metadata()
  b_target <- canon_batch(RUN3)
  is3 <- canon_batch(si$batch) == b_target
  if (!any(is3)) stop("No samples match Run3='", RUN3, "'. Batch values: ",
                      paste(unique(si$batch), collapse = ", "),
                      ". Override with Sys.setenv(CITRIN_RUN3=...).")
  in3 <- si[is3, , drop = FALSE]
  cat("\n[RUN3] Full group counts (before excluding SS):\n")
  print(table(in3$group))
  keep <- in3$group != EXCLUDE
  sm <- in3[keep, , drop = FALSE]
  if (anyNA(sm$family_id))
    stop("Run3 missing family IDs: ", paste(sm$SampleID[is.na(sm$family_id)],
                                            collapse = ", "),
         ". Save family IDs in sample_info.csv or family_map.csv.")
  if (REQUIRE_35 && !(nrow(sm) == 35L && sum(sm$group == POS) == 7L &&
                      sum(sm$group != POS) == 28L))
    stop("Expected Run3 non-SS 35 = 7 SM + 28 others, observed ", nrow(sm),
         " with ", sum(sm$group == POS), " SM and ",
         sum(sm$group != POS), " others. Check batch labels.")
  sm_fam <- table(sm$family_id[sm$group == POS])
  cat("\n[CHECK] SM-positive specimens per family:\n")
  print(sm_fam)
  if (length(sm_fam) != 5L)
    stop("Expected exactly 5 SM-positive families, found ", length(sm_fam))
  if (length(unique(sm$family_id[sm$group == POS])) < 3L)
    stop("<3 independent SM-positive families. A three-fold family-grouped ",
         "inner model is not identifiable; review family relationships.")
  raw <- openxlsx::read.xlsx(RAW_FILE, check.names = FALSE)
  acc_col <- find_col(raw, c("Master Accession Number", "master_accession_number"))
  pep_col <- find_col(raw, c("Unique peptides", "unique_peptides"))
  features <- normalize_id(raw[[acc_col]])
  if (anyNA(features) || any(features == ""))
    stop("Missing Master Accession Number in raw_data.xlsx")
  if (anyDuplicated(features)) {
    warning("Duplicated accession rows: retaining first row, as original code.")
    keepuniq <- !duplicated(features)
    raw <- raw[keepuniq, , drop = FALSE]
    features <- features[keepuniq]
  }
  if (!"P60709" %in% features) stop("ACTB accession P60709 not found")
  ix <- match(sm$SampleID, normalize_id(names(raw)))
  if (anyNA(ix)) stop("Run3 SampleIDs absent from raw_data.xlsx: ",
                      paste(sm$SampleID[is.na(ix)], collapse = ", "))
  X_raw <- as.matrix(raw[, ix, drop = FALSE])
  suppressWarnings(storage.mode(X_raw) <- "numeric")
  colnames(X_raw) <- sm$SampleID
  rownames(X_raw) <- features
  actb <- as.numeric(X_raw["P60709", ])
  if (anyNA(actb) || any(!is.finite(actb) | actb <= 0))
    stop("Missing/zero ACTB intensity in Run3 specimens")
  X_raw[!is.finite(X_raw) | X_raw <= 0] <- NA_real_
  X_log_norm <- sweep(log2(X_raw), 2L, log2(actb), "-")
  peptides <- suppressWarnings(as.numeric(raw[[pep_col]]))
  retain <- is.finite(peptides) & peptides >= MIN_UNIQUE_PEPTIDES
  X_log_norm <- X_log_norm[retain, , drop = FALSE]
  X <- t(X_log_norm)             # 35 specimens x proteins, NA still present
  if (!identical(rownames(X), sm$SampleID)) stop("Metadata/abundance alignment failed")
  y <- as.integer(sm$group == POS)
  cat("[RUN3] Classification: ", nrow(X), " subjects (", sum(y), " SM / ",
      sum(!y), " others), ", ncol(X), " proteins, ",
      length(unique(sm$family_id)), " families.\n", sep = "")
  list(X = X, y = y, metadata = sm)
}

# One-row-at-a-time MinProb transform for prospectively unknown specimens.
# MinProb's location is the sample's own q=.01 abundance quantile; the scatter
# comes from the TRAINING proteomic SD distribution, never the held-out rows.
# Training data are imputed with imputeLCMD::impute.MinProb exactly.
train_minprob <- function(X, seed) {
  stopifnot(is.matrix(X), nrow(X) >= 3L)
  # imputeLCMD works on proteins x samples and uses the sample-wise q-quantile.
  obs_fraction <- colMeans(is.finite(X))
  sds <- apply(X[, obs_fraction > 0.5, drop = FALSE], 2L,
               stats::sd, na.rm = TRUE)
  sigma <- stats::median(sds[is.finite(sds)], na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0)
    stop("Cannot estimate MinProb sigma from >50%-detected proteins")
  if (any(rowSums(is.finite(X)) < 10L))
    stop("A training sample has fewer than 10 observed proteins")
  set.seed(seed)
  imputed <- imputeLCMD::impute.MinProb(t(X), q = MINPROB_Q)
  imputed <- t(imputed)
  if (!identical(dim(imputed), dim(X))) stop("MinProb returned unexpected dimensions")
  dimnames(imputed) <- dimnames(X)
  if (any(!is.finite(imputed))) stop("MinProb produced nonfinite training intensities")
  list(X = imputed, sigma = sigma, n_imputed = sum(!is.finite(X)))
}
apply_minprob <- function(X_test, train_sigma, seed) {
  if (!is.matrix(X_test)) X_test <- as.matrix(X_test)
  ans <- X_test
  set.seed(seed)
  for (i in seq_len(nrow(ans))) {
    observed <- as.numeric(ans[i, is.finite(ans[i, ])])
    if (length(observed) < 10L) stop("Held-out SampleID ", rownames(ans)[i],
                                     " has <10 observed proteins")
    q <- unname(stats::quantile(observed, probs = MINPROB_Q,
                                na.rm = TRUE, names = FALSE))
    missed <- which(!is.finite(ans[i, ]))
    if (length(missed)) ans[i, missed] <- rnorm(length(missed), q, train_sigma)
  }
  if (any(!is.finite(ans))) stop("Nonfinite held-out MinProb values")
  ans
}

make_preprocessor <- function(X_train, seed) {
  imp <- train_minprob(X_train, seed)
  xx <- imp$X
  good <- which(colSums(is.finite(xx)) == nrow(xx) &
                apply(xx, 2L, stats::var) > 0)
  if (length(good) < max(PANEL_GRID))
    stop("Too few training proteins after zero-variance/nonfinite filtering")
  xx <- xx[, good, drop = FALSE]
  mu <- colMeans(xx)
  ss <- apply(xx, 2L, stats::sd)
  scaled <- sweep(sweep(xx, 2L, mu, "-"), 2L, ss, "/")
  list(X = scaled, kept = colnames(xx), mu = mu, sd = ss,
       minprob_sigma = imp$sigma, n_imputed = imp$n_imputed)
}
score_preprocessor <- function(pre, X_test, seed) {
  xx <- apply_minprob(X_test, pre$minprob_sigma, seed)
  xx <- xx[, pre$kept, drop = FALSE]
  z <- sweep(sweep(xx, 2L, pre$mu, "-"), 2L, pre$sd, "/")
  if (any(!is.finite(z))) stop("Nonfinite values in held-out standardization")
  z
}

# AUC via Wilcoxon rank statistic. No ROC package needed.
auc <- function(y, prob) {
  np <- sum(y == 1L); nn <- sum(y == 0L)
  if (!np || !nn || any(!is.finite(prob))) return(NA_real_)
  ranks <- rank(prob, ties.method = "average")
  (sum(ranks[y == 1L]) - np * (np + 1) / 2) / (np * nn)
}
weights_original <- function(y) {
  npos <- sum(y == 1L); nneg <- sum(y == 0L)
  if (npos < 2L || nneg < 2L) stop("glmnet requires >=2 of each class")
  ifelse(y == 1L, nneg / npos, 1)
}

# Family-stratified 3-fold allocation. Every test and training fold must have
# both classes; a grouped allocation is not always mathematically possible.
make_group_3fold <- function(y, family) {
  fam <- as.character(family)
  u <- unique(fam)
  npos <- vapply(u, function(f) sum(y[fam == f] == 1L), integer(1))
  nneg <- vapply(u, function(f) sum(y[fam == f] == 0L), integer(1))
  positive_fams <- u[npos > 0L]
  negative_only <- u[npos == 0L]
  if (length(positive_fams) < 3L || length(negative_only) < 3L)
    return(NULL)
  for (attempt in seq_len(30L)) {
    assignment <- setNames(rep(NA_integer_, length(u)), u)
    # Positive families: distribute round-robin so every fold has >=1.
    p_order <- sample(positive_fams)
    for (i in seq_along(p_order)) assignment[p_order[i]] <- 1L + (i - 1L) %% 3L
    # Negative-only families: distribute round-robin to guarantee negatives.
    n_order <- sample(negative_only)
    for (i in seq_along(n_order)) assignment[n_order[i]] <- 1L + (i - 1L) %% 3L
    folds <- unname(assignment[fam])
    ok <- all(vapply(1:3, function(f) {
      train_y <- y[folds != f]; test_y <- y[folds == f]
      sum(train_y == 1L) >= 2L && sum(train_y == 0L) >= 2L &&
        sum(test_y == 1L) >= 1L && sum(test_y == 0L) >= 1L
    }, logical(1)))
    if (ok) return(folds)
  }
  NULL
}

# Explicit 3-fold cross-fitted AUC lambda selection. cv.glmnet(AUC) may
# automatically switch to deviance if very few observations per fold; by
# computing pooled out-of-fold AUC we preserve the AUC criterion here.
fit_en_original <- function(X_train, y, family, w = NULL, seed = 1L) {
  stopifnot(nrow(X_train) == length(y), length(y) == length(family))
  if (is.null(w)) w <- weights_original(y)
  if (sum(y == 1L) < 2L || sum(y == 0L) < 2L)
    stop("<2 observations from a class in model training")
  set.seed(seed)
  folds <- make_group_3fold(y, family)
  # When the mean fold size is >=10, run the ORIGINAL cv.glmnet(AUC)
  # procedure, but with family-grouped foldid. For tiny folds, glmnet may
  # internally replace AUC by deviance; in that setting explicitly compute
  # pooled out-of-fold AUC below to retain the original optimization metric.
  if (!is.null(folds) && nrow(X_train) >= 30L) {
    warnings_cv <- character(0)
    trial_cv <- tryCatch(withCallingHandlers(
      glmnet::cv.glmnet(X_train, y, family = "binomial", alpha = ALPHA,
                         weights = w, standardize = FALSE,
                         type.measure = "auc", foldid = folds,
                         maxit = 100000L),
      warning = function(ww) {
        warnings_cv <<- c(warnings_cv, conditionMessage(ww))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    switches_metric <- any(grepl("deviance|too few observations|auc.*cannot|auc.*not",
                                  warnings_cv, ignore.case = TRUE))
    if (!inherits(trial_cv, "error") && !switches_metric &&
        is.finite(trial_cv$lambda.min) && trial_cv$lambda.min > 0) {
      return(list(fit = trial_cv$glmnet.fit,
                  lambda = as.numeric(trial_cv$lambda.min),
                  lambda_method = "original_cv.glmnet_grouped_auc",
                  lambda_auc = NA_real_, fallback_reason = ""))
    }
  }
  fit <- glmnet::glmnet(X_train, y, family = "binomial", alpha = ALPHA,
                        weights = w, standardize = FALSE, maxit = 100000L)
  lambda_path <- as.numeric(fit$lambda)
  reason <- ""
  method <- "grouped_3fold_pooled_auc"
  if (!is.null(folds) && length(lambda_path) >= 2L) {
    oof <- matrix(NA_real_, nrow = length(y), ncol = length(lambda_path))
    err <- NULL
    for (f in 1:3) {
      tr <- which(folds != f); te <- which(folds == f)
      trial <- tryCatch({
        # Keep the training-only master lambda sequence across folds.
        ff <- glmnet::glmnet(X_train[tr, , drop = FALSE], y[tr],
                             family = "binomial", alpha = ALPHA,
                             weights = w[tr], standardize = FALSE,
                             lambda = lambda_path, maxit = 100000L)
        as.matrix(predict(ff, newx = X_train[te, , drop = FALSE],
                          s = lambda_path, type = "response"))
      }, error = function(e) e)
      if (inherits(trial, "error")) { err <- conditionMessage(trial); break }
      if (nrow(trial) != length(te) || ncol(trial) != length(lambda_path)) {
        err <- "Unexpected glmnet prediction dimensions"; break
      }
      oof[te, ] <- trial
    }
    if (is.null(err) && all(is.finite(oof))) {
      scores <- vapply(seq_along(lambda_path),
                       function(j) auc(y, oof[, j]), numeric(1))
      if (any(is.finite(scores))) {
        # glmnet path goes from larger to smaller lambda. Favor the larger
        # lambda in an AUC tie, without referring to any held-out family.
        idx <- which(scores == max(scores, na.rm = TRUE))[1L]
        return(list(fit = fit, lambda = lambda_path[idx],
                    lambda_method = method, lambda_auc = scores[idx],
                    fallback_reason = ""))
      }
      err <- "All pooled AUCs nonfinite"
    }
    reason <- if (is.null(err)) "Nonfinite CV predictions" else err
  } else {
    reason <- "Insufficient independent positive/negative families for grouped 3-fold AUC CV"
  }
  lambda_fixed <- max(lambda_path) * LAMBDA_FALLBACK_RATIO
  list(fit = fit, lambda = lambda_fixed, lambda_method = "fallback_lambda0.1max",
       lambda_auc = NA_real_, fallback_reason = reason)
}

# 30 iterations of ORIGINAL class-stratified, sample-level 80% subsampling.
# Parent training data are already standardized as in the original script.
# If 3 positives exist, we retain all three (instead of original 2) because
# otherwise every 3-fold logistic tuning split would have <=1 positive.
stability_rank <- function(X_train, y, family, seed) {
  p <- ncol(X_train)
  freq <- integer(p); names(freq) <- colnames(X_train)
  n_success <- 0L; n_fallback <- 0L
  n_raised_pos <- 0L; failures <- character(0)
  idxp <- which(y == 1L); idxn <- which(y == 0L)
  if (length(idxp) < 2L || length(idxn) < 2L)
    stop("Stability ranking requires >=2 observations of both classes")
  # Class weights are those from the full CURRENT training set, as original.
  parent_w <- weights_original(y)
  set.seed(seed)
  for (b in seq_len(B_STABILITY)) {
    set.seed(family_seed(paste0("stability|", seed, "|", b), 70000L))
    p_count_original <- max(2L, floor(length(idxp) * STABILITY_FRAC))
    p_count <- min(length(idxp), max(3L, p_count_original))
    n_count <- min(length(idxn), max(2L, floor(length(idxn) * STABILITY_FRAC)))
    if (SUBSAMPLE_UNIT == "family") {
      # Include ALL subjects (positive and negative) from selected families.
      # Sample the positive-containing families and negative-only families separately.
      pf <- sort(unique(as.character(family[idxp])))
      nf <- sort(setdiff(unique(as.character(family)), pf))
      kp <- min(length(pf), max(3L, floor(length(pf) * STABILITY_FRAC)))
      kn <- min(length(nf), max(2L, floor(length(nf) * STABILITY_FRAC)))
      chosen_fam <- c(sample(pf, kp), sample(nf, kn))
      tr <- which(as.character(family) %in% chosen_fam)
      if (length(pf) > 0L && length(unique(as.character(family[tr][y[tr] == 1L]))) < 3L)
        stop("Family subsampling left fewer than 3 positive families")
    } else {
      if (p_count > p_count_original) n_raised_pos <- n_raised_pos + 1L
      tr <- c(sample(idxp, p_count), sample(idxn, n_count))
    }
    Xi <- X_train[tr, , drop = FALSE]
    yi <- y[tr]; fi <- family[tr]; wi <- parent_w[tr]
    # Original top-1000 rule: variance computed AFTER parent standardization
    # but WITHIN the 80%-sample subset (not variance of raw intensities).
    if (ncol(Xi) > N_VARIANCE) {
      variances <- apply(Xi, 2L, stats::var)
      variances[!is.finite(variances)] <- -Inf
      selected <- order(variances, decreasing = TRUE)[seq_len(N_VARIANCE)]
      Xi <- Xi[, selected, drop = FALSE]
    }
    # Do NOT rescale Xi here: original pipeline standardize=FALSE and scales
    # X_clean BEFORE stability subsampling.
    model <- tryCatch(fit_en_original(Xi, yi, fi, w = wi,
                                      seed = seed + 1000L + b),
                      error = function(e) e)
    if (inherits(model, "error")) {
      failures <- c(failures, conditionMessage(model)); next
    }
    co <- as.matrix(stats::coef(model$fit, s = model$lambda))[-1L, 1L]
    chosen <- names(co)[is.finite(co) & abs(co) > 0]
    if (length(chosen)) freq[chosen] <- freq[chosen] + 1L
    n_success <- n_success + 1L
    if (model$lambda_method == "fallback_lambda0.1max")
      n_fallback <- n_fallback + 1L
  }
  min_success <- if (B_STABILITY < 5L) 1L else max(3L, ceiling(B_STABILITY / 2))
  if (n_success < min_success)
    stop("Stability selection failed: ", n_success, "/", B_STABILITY,
         " successful fits. First failure: ",
         if (length(failures)) failures[1L] else "no eligible fits")
  # Original stability frequency sorted descending. Original feature-order
  # ties are retained, rather than inventing another ranking criterion.
  ord <- order(-freq, seq_along(freq))
  list(ids = names(freq)[ord], frequency = freq[ord] / n_success,
       n_success = n_success, n_fallback = n_fallback,
       n_raised_pos = n_raised_pos,
       n_failed = B_STABILITY - n_success,
       first_error = if (length(failures)) failures[1L] else "")
}

# Same threshold rule as original source, including original tie-breaking.
select_threshold <- function(prob, y, floor_spec = SPEC_FLOOR, grid_len = 200L) {
  grid <- as.numeric(stats::quantile(prob,
                     probs = seq(0, 1, length.out = grid_len),
                     na.rm = TRUE))
  grid <- sort(unique(c(0, grid, 1)))
  best <- NULL
  for (t in grid) {
    pr <- as.integer(prob >= t)
    se <- sum(pr == 1L & y == 1L) / sum(y == 1L)
    sp <- sum(pr == 0L & y == 0L) / sum(y == 0L)
    if (is.finite(sp) && sp >= floor_spec) {
      if (is.null(best) || se > best$sens ||
          (se == best$sens && sp > best$spec) ||
          (se == best$sens && sp == best$spec && t > best$tau)) {
        best <- list(tau = t, sens = se, spec = sp, policy = "spec_floor")
      }
    }
  }
  if (is.null(best)) {
    youden <- vapply(grid, function(t) {
      pr <- as.integer(prob >= t)
      sum(pr == 1L & y == 1L) / sum(y == 1L) +
        sum(pr == 0L & y == 0L) / sum(y == 0L) - 1
    }, numeric(1))
    t <- grid[which.max(youden)]
    pr <- as.integer(prob >= t)
    best <- list(tau = t,
                 sens = sum(pr == 1L & y == 1L) / sum(y == 1L),
                 spec = sum(pr == 0L & y == 0L) / sum(y == 0L),
                 policy = "youden_fallback")
  }
  best
}

predict_panel <- function(fit, pre, X_test, selected, seed) {
  zz <- score_preprocessor(pre, X_test, seed)
  as.numeric(predict(fit$fit, newx = zz[, selected, drop = FALSE],
                     s = fit$lambda, type = "response"))
}

# Nested inner procedure: ranked features re-selected for EACH inner fold,
# MinProb/centering/scale re-fitted on EACH inner training split, then each
# panel size gets its own lambda and out-of-fold predictions.
inner_select <- function(X, y, family, outer_fam) {
  families <- sort(unique(as.character(family)))
  if (INNER_MODE == "lofo") {
    inner_fold <- match(as.character(family), families)
  } else {
    set.seed(family_seed(outer_fam, 10000L))
    inner_fold <- make_group_3fold(y, family)
    if (is.null(inner_fold)) stop("Unable to construct grouped 3-fold inner split")
  }
  m_grid <- if (PANEL_MODE == "fixed15") 15L else if (PANEL_MODE == "fixed30") 30L else PANEL_GRID
  predictions <- matrix(NA_real_, nrow = nrow(X), ncol = length(m_grid),
                        dimnames = list(rownames(X), as.character(m_grid)))
  audit <- vector("list", length(unique(inner_fold)))
  for (j in seq_along(sort(unique(inner_fold)))) {
    k <- sort(unique(inner_fold))[j]
    tr <- which(inner_fold != k); te <- which(inner_fold == k)
    if (sum(y[tr] == 1L) < 2L || sum(y[tr] == 0L) < 2L)
      stop("Inner training split has <2 cases of a class, fold=", k,
           ". This family structure cannot support the original model.")
    seed <- family_seed(paste(outer_fam, paste(sort(unique(family[te])), collapse = ","), sep = "|"), 5000000L)
    pre <- make_preprocessor(X[tr, , drop = FALSE], seed + 1L)
    ranking <- stability_rank(pre$X, y[tr], family[tr], seed + 2L)
    zd <- score_preprocessor(pre, X[te, , drop = FALSE], seed + 50L)
    measures <- vector("list", length(m_grid))
    for (mi in seq_along(m_grid)) {
      m <- m_grid[mi]
      ids <- head(ranking$ids, m)
      fit <- fit_en_original(pre$X[, ids, drop = FALSE], y[tr], family[tr],
                             seed = seed + 1000L + mi)
      predictions[te, as.character(m)] <- as.numeric(predict(
        fit$fit, newx = zd[, ids, drop = FALSE],
        s = fit$lambda, type = "response"))
      measures[[mi]] <- data.frame(inner_fold = k, n_train = length(tr),
                                  pos_train = sum(y[tr] == 1L),
                                  n_held = length(te), m = m,
                                  lambda_method = fit$lambda_method,
                                  lambda = fit$lambda,
                                  reason = fit$fallback_reason,
                                  ranking_fallbacks = ranking$n_fallback,
                                  ranking_success = ranking$n_success,
                                  raised_pos_draws = ranking$n_raised_pos)
    }
    audit[[j]] <- do.call(rbind, measures)
  }
  if (any(!is.finite(predictions))) stop("Nonfinite/missing inner cross-fitted predictions")
  grid <- do.call(rbind, lapply(m_grid, function(m) {
    p <- predictions[, as.character(m)]
    tau <- select_threshold(p, y)
    data.frame(m = m, AUC = auc(y, p),
               sensitivity = tau$sens, specificity = tau$spec,
               threshold = tau$tau, threshold_policy = tau$policy)
  }))
  # Original candidate selection policy: max AUC, then sensitivity, then
  # specificity, then smaller panel. ALL chosen from outer-training only.
  best <- order(-grid$AUC, -grid$sensitivity,
                -grid$specificity, grid$m)[1L]
  list(m = as.integer(grid$m[best]), threshold = grid$threshold[best],
       threshold_policy = grid$threshold_policy[best],
       metrics = grid, audit = do.call(rbind, audit))
}

metrics <- function(truth, prob, pred) {
  tp <- sum(truth == 1L & pred == 1L)
  tn <- sum(truth == 0L & pred == 0L)
  fp <- sum(truth == 0L & pred == 1L)
  fn <- sum(truth == 1L & pred == 0L)
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  c(AUC = auc(truth, prob), sensitivity = sens, specificity = spec,
    balanced_accuracy = (sens + spec) / 2,
    Brier = mean((prob - truth)^2), accuracy = mean(pred == truth),
    TP = tp, TN = tn, FP = fp, FN = fn)
}

family_bootstrap_ci <- function(df, reps, seed = 2025L) {
  point <- metrics(df$truth, df$probability, df$predicted_class)
  metric_names <- c("AUC", "sensitivity", "specificity",
                    "balanced_accuracy", "Brier", "accuracy")
  lower <- upper <- rep(NA_real_, length(metric_names))
  names(lower) <- names(upper) <- metric_names
  if (reps > 0L) {
    fams <- unique(df$family_id)
    by_family <- split(seq_len(nrow(df)), df$family_id)
    set.seed(seed)
    vals <- matrix(NA_real_, nrow = reps, ncol = length(metric_names))
    colnames(vals) <- metric_names
    for (b in seq_len(reps)) {
      chosen <- sample(fams, size = length(fams), replace = TRUE)
      row_idx <- unlist(by_family[chosen], use.names = FALSE)
      dd <- df[row_idx, , drop = FALSE]
      if (any(table(factor(dd$truth, levels = c(0, 1))) == 0L)) next
      vals[b, ] <- metrics(dd$truth, dd$probability, dd$predicted_class)[metric_names]
    }
    for (m in metric_names) {
      vv <- vals[, m]
      vv <- vv[is.finite(vv)]
      if (length(vv) > 20L) {
        lower[m] <- unname(quantile(vv, .025))
        upper[m] <- unname(quantile(vv, .975))
      }
    }
  }
  exact_sens <- binom.test(point["TP"], point["TP"] + point["FN"])$conf.int
  exact_spec <- binom.test(point["TN"], point["TN"] + point["FP"])$conf.int
  out <- data.frame(metric = metric_names,
                    estimate = as.numeric(point[metric_names]),
                    cluster_bootstrap_low = lower[metric_names],
                    cluster_bootstrap_high = upper[metric_names],
                    exact_binomial_low = NA_real_,
                    exact_binomial_high = NA_real_)
  out$exact_binomial_low[out$metric == "sensitivity"] <- exact_sens[1L]
  out$exact_binomial_high[out$metric == "sensitivity"] <- exact_sens[2L]
  out$exact_binomial_low[out$metric == "specificity"] <- exact_spec[1L]
  out$exact_binomial_high[out$metric == "specificity"] <- exact_spec[2L]
  list(point = point, ci = out)
}

# Resumable full nested outer leave-one-family-out run.
run_outer <- function(dat, outdir) {
  X <- dat$X; y <- dat$y; sm <- dat$metadata
  family <- as.character(sm$family_id)
  families <- sort(unique(family))
  if (OUTER_LIMIT > 0L) families <- head(families, OUTER_LIMIT)
  dir.create(file.path(outdir, "outer_cache"), recursive = TRUE, showWarnings = FALSE)
  results <- vector("list", length(families))
  for (i in seq_along(families)) {
    fam <- families[i]
    fn <- file.path(outdir, "outer_cache", paste0("family_", substr(sprintf("%08x", family_seed(fam, 0L)), 1, 8), ".rds"))
    if (file.exists(fn)) {
      one <- readRDS(fn)
      if (!identical(as.character(one$held_family), fam))
        stop("Cache family mismatch at ", fn, ". Do not mix run configurations.")
      results[[i]] <- one
      cat(sprintf("Outer %d/%d: family %s [cached]\n", i, length(families), fam))
      next
    }
    tr <- which(family != fam); te <- which(family == fam)
    cat(sprintf("Outer %d/%d: held family %s, n=%d (SM=%d)\n",
                i, length(families), fam, length(te), sum(y[te] == 1L)))
    inner <- inner_select(X[tr, , drop = FALSE], y[tr], family[tr], outer_fam=fam)
    seed <- family_seed(fam, 10000000L)
    pre <- make_preprocessor(X[tr, , drop = FALSE], seed + 1L)
    rank <- stability_rank(pre$X, y[tr], family[tr], seed + 2L)
    ids <- head(rank$ids, inner$m)
    fit <- fit_en_original(pre$X[, ids, drop = FALSE], y[tr], family[tr],
                           seed = seed + 3L)
    prob <- predict_panel(fit, pre, X[te, , drop = FALSE], ids, seed + 4L)
    obs <- data.frame(SampleID = rownames(X)[te], family_id = family[te],
                      batch = sm$batch[te], group = sm$group[te],
                      truth = y[te], probability = prob,
                      threshold = inner$threshold,
                      predicted_class = as.integer(prob >= inner$threshold),
                      selected_m = inner$m, selected_lambda = fit$lambda,
                      lambda_method = fit$lambda_method,
                      lambda_fallback_reason = fit$fallback_reason,
                      inner_threshold_policy = inner$threshold_policy)
    coefs <- extract_outer_coefficients(fit, ids, fam, rank$frequency)
    one <- list(held_family = fam, predictions = obs, coefficients = coefs,
                inner_grid = cbind(held_family = fam, inner$metrics),
                lambda_audit = cbind(held_family = fam, inner$audit),
                features = data.frame(held_family = fam, rank = seq_along(ids),
                                      FeatureID = ids,
                                      stability_frequency = unname(rank$frequency[ids])),
                stability = data.frame(held_family = fam,
                                       successful_refits = rank$n_success,
                                       fallback_lambda_refits = rank$n_fallback,
                                       positive_draw_adjustments = rank$n_raised_pos,
                                       failed_refits = rank$n_failed))
    saveRDS(one, fn)
    results[[i]] <- one
    # Small checkpoint after each finished fold; safe if interrupted later.
    write.csv(do.call(rbind, lapply(results[!vapply(results, is.null, logical(1))],
                                       `[[`, "predictions")),
              file.path(outdir, "outer_predictions_PARTIAL.csv"), row.names = FALSE)
  }
  if (OUTER_LIMIT > 0L && OUTER_LIMIT < length(unique(family))) {
    message("[SMOKE TEST ONLY] Only ", OUTER_LIMIT, " outer families processed. ",
            "No valid CV performance estimates will be written.")
    return(invisible(NULL))
  }
  pred <- do.call(rbind, lapply(results, `[[`, "predictions"))
  pred <- pred[match(rownames(X), pred$SampleID), , drop = FALSE]
  if (anyNA(pred$probability) || anyDuplicated(pred$SampleID))
    stop("Missing or duplicate outer predictions")
  write.csv(pred, file.path(outdir, "outer_predictions.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "inner_grid")),
            file.path(outdir, "outer_inner_grid.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "lambda_audit")),
            file.path(outdir, "outer_inner_lambda_audit.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "features")),
            file.path(outdir, "outer_selected_features.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "coefficients")),
            file.path(outdir, "outer_coefficients.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "stability")),
            file.path(outdir, "outer_stability_diagnostics.csv"), row.names = FALSE)
  summarised <- family_bootstrap_ci(pred, N_BOOT, SEED + 77777L)
  write.csv(data.frame(metric = names(summarised$point),
                       estimate = as.numeric(summarised$point)),
            file.path(outdir, "outer_performance.csv"), row.names = FALSE)
  write.csv(summarised$ci, file.path(outdir, "outer_95CI.csv"), row.names = FALSE)
  cat("\n[RUN3 NESTED OUTER RESULTS]\n")
  print(data.frame(metric = names(summarised$point),
                   estimate = as.numeric(summarised$point)), row.names = FALSE)
  cat("\nSelected panel sizes (outer-fold family counts):\n")
  print(table(vapply(results, function(z) z$predictions$selected_m[1L], numeric(1))))
  cat("\nFallbacks in inner stability and lambda fits are in *_audit.csv.\n")
  invisible(summarised)
}

# Optional final discovery fit on all 35 Run3 samples. This is NOT the tested
# fixed model: out-of-fold metrics estimate the entire adaptive algorithm.
final_discovery_fit <- function(dat, outdir) {
  cat("\nFinal Run3 full-sample discovery fit (not independently validated) ...\n")
  X <- dat$X; y <- dat$y; fam <- as.character(dat$metadata$family_id)
  inner <- inner_select(X, y, fam, outer_fam = "FINAL_DISCOVERY")
  pre <- make_preprocessor(X, SEED + 500001L)
  rank <- stability_rank(pre$X, y, fam, SEED + 500002L)
  ids <- head(rank$ids, inner$m)
  fit <- fit_en_original(pre$X[, ids, drop = FALSE], y, fam,
                         seed = SEED + 500003L)
  co <- as.matrix(stats::coef(fit$fit, s = fit$lambda))[, 1L]
  coefs <- data.frame(FeatureID = rownames(as.matrix(stats::coef(fit$fit, s = fit$lambda))), coefficient = as.numeric(co),
                      selected_panel_size = inner$m,
                      lambda = fit$lambda, threshold = inner$threshold)
  write.csv(coefs, file.path(outdir, "final_discovery_coefficients.csv"),
            row.names = FALSE)
  write.csv(data.frame(FeatureID = ids,
                       center = unname(pre$mu[ids]),
                       scale = unname(pre$sd[ids]),
                       stability_frequency = unname(rank$frequency[ids])),
            file.path(outdir, "final_discovery_feature_scaling.csv"),
            row.names = FALSE)
  saveRDS(list(preprocessing = pre, features = ids, glmnet_fit = fit$fit,
               lambda = fit$lambda, threshold = inner$threshold,
               minprob_q = MINPROB_Q, trained_run = RUN3),
          file.path(outdir, "final_discovery_model.rds"))
  write.csv(inner$metrics, file.path(outdir, "final_discovery_inner_grid.csv"),
            row.names = FALSE)
  message("[INFO] Final exploratory model saved. Do not report its training ",
          "accuracy as validated performance.")
}

hash_config <- function(dat) {
  files <- c(RAW_FILE, SAMPLE_FILE)
  input_hash <- unname(tools::md5sum(files))
  config <- paste(SCRIPT_VERSION, paste(input_hash, collapse = "|"),
                  paste(dat$metadata$SampleID, dat$metadata$family_id,
                        dat$metadata$group, sep = ":", collapse = "|"),
                  SEED, B_STABILITY, PANEL_MODE, INNER_MODE, SUBSAMPLE_UNIT,
                  SPEC_FLOOR, MINPROB_Q, LAMBDA_FALLBACK_RATIO, RUN3, sep = "|")
  tmp <- tempfile(); writeLines(config, tmp, useBytes = TRUE)
  on.exit(unlink(tmp), add = TRUE)
  substr(unname(tools::md5sum(tmp)), 1L, 12L)
}

main <- function() {
  dat <- read_run3()
  dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)
  key <- hash_config(dat)
  label <- paste0("B", B_STABILITY, "_", PANEL_MODE, "_", INNER_MODE, "_", key)
  outdir <- file.path(OUT_BASE, label)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  params <- data.frame(name = c("script_version", "data_dir", "run3", "n_samples",
                                "n_SM", "n_controls", "n_proteins", "n_families",
                                "B_stability", "panel_mode", "inner_mode", "subsample_unit",
                                "seed", "lambda_fallback", "unique_peptides_min",
                                "minprob_q", "grouped_outer", "ComBat"),
                       value = c(SCRIPT_VERSION, DATA_DIR, RUN3, nrow(dat$X),
                                 sum(dat$y), sum(dat$y == 0L), ncol(dat$X),
                                 length(unique(dat$metadata$family_id)),
                                 B_STABILITY, PANEL_MODE, INNER_MODE, SUBSAMPLE_UNIT, SEED,
                                 LAMBDA_FALLBACK_RATIO, MIN_UNIQUE_PEPTIDES,
                                 MINPROB_Q, "leave_one_family_out", "NOT_USED_single_batch"))
  write.csv(params, file.path(outdir, "analysis_parameters.csv"), row.names = FALSE)
  write.csv(dat$metadata[, c("SampleID", "group", "batch", "family_id")],
            file.path(outdir, "run3_cohort_used.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
  cat("\nOutput directory: ", normalizePath(outdir, mustWork = FALSE), "\n", sep = "")
  ans <- run_outer(dat, outdir)
  if (is.null(ans)) return(invisible(NULL))
  if (tolower(Sys.getenv("CITRIN_FINAL_MODEL", "false")) == "true")
    final_discovery_fit(dat, outdir)
  cat("\n[COMPLETE] Run3 nested validation; files in ", outdir, "\n", sep = "")
  invisible(ans)
}

# Sourcing the file executes main(), as prior v3. Re-running main() uses
# cached successful outer folds. No existing family-ID file is ever edited.
main()
