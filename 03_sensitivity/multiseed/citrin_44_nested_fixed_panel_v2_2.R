#!/usr/bin/env Rscript
# ============================================================================
# CITRIN / 44-SAMPLE TWO-RUN FAMILY-GROUPED NESTED CROSS-VALIDATION
# Stage 1: 44 samples (7 SM + 37 controls) WITHOUT ComBat.
# Stage 2: SAME splits, seeds, panel grid, thresholds, training-only
#          sva::ComBat(mod=~group, mean.only=TRUE); frozen batch shifts on tests.
# Stage 3: matched comparison and optional Run3-only 35-sample reference.
# ============================================================================
# Manuscript settings retained as closely as possible:
#  ACTB P60709 ratio; >=3 unique peptides; MinProb(q=.01); elastic net alpha=.5;
#  class-inverse weights; B=30 80% class-stratified stability subsamples;
#  top 1000 within-subsample variance; 3-fold grouped AUC lambda selection;
#  panels m=5,10,15,20,30,40,50; threshold maximizes sensitivity with
#  training specificity >=.90; family-grouped nested outer and inner CV.
# Changed for validation:
#  - ALL adaptive preprocessing and model decisions occur within training folds.
#  - Frozen mean-only ComBat replaces original full-data mean+variance ComBat.
#  - Test MinProb uses its own observed q=.01 location and training-derived SD.
#  - Grouped 3-fold AUC and small-sample fallback logged in audit CSVs.
# NOTE: Upstream DIA-NN processing is NOT re-estimated within CV.
# NO original input or family-ID file is created, altered, or overwritten.
# Outer/inner seeds keyed by held family rather than its sorted position.
# Bootstrap intervals are descriptive and conditional on realized OOF scores.
# ============================================================================

# IMPORTANT INFERENCE NOTE:
# fixed15 evaluates a PROCEDURE that re-selects 15 protein features within
# each training fold; it does NOT externally validate the exact 15 proteins
# originally chosen using all 44 samples. In the submitted manuscript,
# the 15 nonzero coefficients came from refitting a 30-CANDIDATE panel.
# Thus fixed30 is closer to the original model development procedure.
# This analysis was devised after inspecting the outcome data: treat the
# fixed-panel comparison as exploratory, not a prospective confirmation.

SCRIPT_VERSION <- "two_run_fixed_panel_2.2"
SEED <- as.integer(Sys.getenv("CITRIN_SEED", "2025"))
SUBSAMPLE_UNIT <- tolower(Sys.getenv("CITRIN_SUBSAMPLE_UNIT", "sample"))
if (!SUBSAMPLE_UNIT %in% c("sample", "family")) stop("CITRIN_SUBSAMPLE_UNIT: use sample or family")
POS <- "SM"; EXCLUDE <- "SS"
RUN2 <- Sys.getenv("CITRIN_RUN2", "Run2")
RUN3 <- Sys.getenv("CITRIN_RUN3", "Run3")
ALPHA <- 0.5; MINPROB_Q <- 0.01; MIN_UNIQUE_PEPTIDES <- 3L
B_STABILITY <- as.integer(Sys.getenv("CITRIN_B", "30"))
STABILITY_FRAC <- 0.80; N_VARIANCE <- 1000L
PANEL_GRID <- c(5L, 10L, 15L, 20L, 30L, 40L, 50L)
SPEC_FLOOR <- 0.90
PANEL_MODE <- tolower(Sys.getenv("CITRIN_PANEL_MODE", "select"))
INNER_MODE <- tolower(Sys.getenv("CITRIN_INNER_MODE", "lofo"))
STAGE <- tolower(Sys.getenv("CITRIN_STAGE", "both"))
LAMBDA_FALLBACK_RATIO <- 0.10
N_BOOT <- as.integer(Sys.getenv("CITRIN_N_BOOT", "2000"))
OUTER_LIMIT <- as.integer(Sys.getenv("CITRIN_OUTER_LIMIT", "0"))
REQUIRE_44 <- tolower(Sys.getenv("CITRIN_EXPECT_44", "true")) == "true"

if (!is.finite(SEED) || !is.finite(B_STABILITY) || B_STABILITY < 1L ||
    !is.finite(N_BOOT) || N_BOOT < 0L || !is.finite(OUTER_LIMIT) || OUTER_LIMIT < 0L)
  stop("Invalid CITRIN_SEED, CITRIN_B, CITRIN_N_BOOT, or CITRIN_OUTER_LIMIT")
if (!PANEL_MODE %in% c("select", "fixed30", "fixed15"))
  stop("CITRIN_PANEL_MODE must be 'select', 'fixed30' or 'fixed15'")
if (!INNER_MODE %in% c("lofo", "group3"))
  stop("CITRIN_INNER_MODE must be 'lofo' or 'group3'")
if (!STAGE %in% c("both", "stage1", "stage2", "compare"))
  stop("CITRIN_STAGE must be 'both', 'stage1', 'stage2', or 'compare'")
if (tolower(RUN2) == tolower(RUN3))
  stop("CITRIN_RUN2 and CITRIN_RUN3 must identify different runs")

options(stringsAsFactors = FALSE)
needed <- c("openxlsx", "imputeLCMD", "glmnet")
if (STAGE %in% c("both", "stage2")) needed <- c(needed, "sva")
absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent)) {
  stop("Missing R package(s): ", paste(absent, collapse = ", "),
       ". CRAN: install.packages(c('openxlsx','glmnet')); ",
       "Bioconductor: BiocManager::install(c('imputeLCMD','sva'))")
}

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
OUT_BASE <- Sys.getenv("CITRIN_OUT_DIR", file.path(getwd(), "result", "two_run_nested_cv"))

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

# Explicitly retain ONLY the 44 non-SS subjects from Run2 and Run3.
read_44 <- function() {
  si <- read_metadata()
  r2 <- canon_batch(RUN2)
  r3 <- canon_batch(RUN3)
  is2 <- canon_batch(si$batch) == r2
  is3 <- canon_batch(si$batch) == r3
  if (!any(is2) || !any(is3))
    stop("Both Run2 and Run3 must be present. Available: ",
         paste(unique(si$batch), collapse=", "))
  in_both <- si[is2 | is3, , drop=FALSE]
  cat("\n[TWO RUNS] Group x batch counts before excluding SS:\n")
  print(table(in_both$batch, in_both$group))
  sm <- in_both[in_both$group != EXCLUDE, , drop=FALSE]
  sm$batch <- ifelse(canon_batch(sm$batch) == r2, "Run2", "Run3")
  if (anyNA(sm$family_id))
    stop("Missing family IDs for these subjects: ",
         paste(sm$SampleID[is.na(sm$family_id)], collapse=", "),
         ". Set Sys.setenv(CITRIN_FAMILY_FILE='/full/path/to/family_map_TEMPLATE.csv').")
  if (REQUIRE_44 && !(nrow(sm) == 44L && sum(sm$group == POS) == 7L &&
                      sum(sm$group != POS) == 37L &&
                      sum(sm$batch == "Run2") == 9L &&
                      sum(sm$batch == "Run3") == 35L))
    stop("Expected 44 subjects (SM=7, Others=37, Run2=9, Run3=35). Found ",
         nrow(sm), "; SM=", sum(sm$group == POS),
         "; Run2=", sum(sm$batch == "Run2"),
         "; Run3=", sum(sm$batch == "Run3"))
  if (any(sm$group[sm$batch == "Run2"] == POS))
    warning("Run2 has SM positives: verify the intended cohort!")
  sm_fam <- table(sm$family_id[sm$group == POS])
  cat("\n[CHECK] SM-positive specimens per family:\n")
  print(sm_fam)
  if (length(sm_fam) != 5L)
    stop("Expected exactly 5 SM-positive families, found ", length(sm_fam))
  if (length(unique(sm$family_id[sm$group == POS])) < 3L)
    stop("<3 independent SM-positive families. Grouped 3-fold AUC not identifiable.")
  raw <- openxlsx::read.xlsx(RAW_FILE, check.names=FALSE)
  acc_col <- find_col(raw, c("Master Accession Number", "master_accession_number"))
  pep_col <- find_col(raw, c("Unique peptides", "unique_peptides"))
  features <- normalize_id(raw[[acc_col]])
  if (anyNA(features) || any(features == ""))
    stop("Missing accession in raw_data.xlsx")
  if (anyDuplicated(features)) {
    warning("Duplicated accessions: retaining first row, as in original code.")
    retain_unique <- !duplicated(features)
    raw <- raw[retain_unique, , drop=FALSE]
    features <- features[retain_unique]
  }
  if (!"P60709" %in% features) stop("ACTB P60709 accession not found")
  ix <- match(sm$SampleID, normalize_id(names(raw)))
  if (anyNA(ix)) stop("Sample IDs absent from raw_data.xlsx: ",
                      paste(sm$SampleID[is.na(ix)], collapse=", "))
  X_raw <- as.matrix(raw[, ix, drop=FALSE])
  suppressWarnings(storage.mode(X_raw) <- "numeric")
  colnames(X_raw) <- sm$SampleID
  rownames(X_raw) <- features
  actb <- as.numeric(X_raw["P60709", ])
  if (anyNA(actb) || any(!is.finite(actb) | actb <= 0))
    stop("Missing, zero, or nonfinite ACTB abundance")
  X_raw[!is.finite(X_raw) | X_raw <= 0] <- NA_real_
  X_log_norm <- sweep(log2(X_raw), 2L, log2(actb), "-")
  peptides <- suppressWarnings(as.numeric(raw[[pep_col]]))
  retain <- is.finite(peptides) & peptides >= MIN_UNIQUE_PEPTIDES
  X_log_norm <- X_log_norm[retain, , drop=FALSE]
  X <- t(X_log_norm)  # subjects x proteins; missing entries retained
  if (!identical(rownames(X), sm$SampleID))
    stop("Metadata/intensity misalignment")
  y <- as.integer(sm$group == POS)
  cat(sprintf("[TWO RUNS] n=%d (SM=%d, others=%d), features=%d, families=%d.\n",
              nrow(X), sum(y), sum(y == 0L), ncol(X),
              length(unique(sm$family_id))))
  print(table(batch=sm$batch, group=sm$group))
  list(X=X, y=y, metadata=sm)
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

# Training-fitted, label-free-at-test batch correction.
# For mean.only=TRUE, x_train - ComBat(x_train) is batch-specific constant.
# We verify this identity numerically and freeze that shift for unseen tests.
fit_meanonly_combat <- function(X_train, batch_train, group_train) {
  stopifnot(is.matrix(X_train), nrow(X_train) == length(batch_train),
            length(batch_train) == length(group_train),
            all(is.finite(X_train)))
  bb <- factor(as.character(batch_train), levels=c("Run2", "Run3"))
  if (anyNA(bb) || length(unique(bb)) != 2L)
    stop("Training-only ComBat requires both Run2 and Run3 in every fit")
  if (any(table(bb) < 2L))
    stop("At least two training specimens per Run are required for ComBat")
  gg <- droplevels(factor(group_train))
  if (nlevels(gg) < 2L) stop("ComBat fit has only one biological group")
  mm <- stats::model.matrix(~ gg)
  batch_design <- stats::model.matrix(~ 0 + bb)
  full_design <- cbind(batch_design, mm[, -1L, drop=FALSE])
  if (qr(full_design)$rank < ncol(full_design))
    stop("Batch/group confounding within a nested training fold. ",
         "Do not silently switch to no correction. Review batch x group table.")
  cb <- tryCatch(suppressMessages(sva::ComBat(
                   dat=t(X_train), batch=bb, mod=mm,
                   par.prior=TRUE, prior.plots=FALSE,
                   mean.only=TRUE)),
                 error=function(e) stop("Training-only ComBat failed: ",
                                       conditionMessage(e)))
  corrected <- t(cb)
  dimnames(corrected) <- dimnames(X_train)
  if (!identical(dim(corrected), dim(X_train)) ||
      any(!is.finite(corrected))) stop("Invalid ComBat training output")
  delta <- X_train - corrected
  batches <- c("Run2", "Run3")
  shift <- vapply(batches, function(b) {
    colMeans(delta[bb == b, , drop=FALSE])
  }, numeric(ncol(X_train)))
  dimnames(shift) <- list(colnames(X_train), batches)
  # Under mean-only ComBat, within-batch x - corrected is the SAME protein
  # vector for all training specimens, up to floating-point rounding.
  max_deviation <- max(vapply(batches, function(b) {
    max(abs(sweep(delta[bb == b, , drop=FALSE], 2L, shift[, b], "-")))
  }, numeric(1)))
  if (!is.finite(max_deviation) || max_deviation > 1e-6)
    stop("Mean-only ComBat is not representable as fixed batch shifts; ",
         "maximum within-batch difference=", signif(max_deviation, 5),
         ". Stopping rather than leaking labels.")
  list(X=corrected, shifts=shift, max_deviation=max_deviation,
       mean_abs_shift=mean(abs(shift)), max_abs_shift=max(abs(shift)))
}

make_preprocessor <- function(X_train, group_train, batch_train, mode, seed) {
  stopifnot(mode %in% c("no_combat", "combat_mean"),
            nrow(X_train) == length(group_train),
            length(group_train) == length(batch_train))
  imp <- train_minprob(X_train, seed)
  xx <- imp$X
  shifts <- NULL
  cbdiag <- list(max_deviation=NA_real_, mean_abs_shift=NA_real_,
                 max_abs_shift=NA_real_)
  if (mode == "combat_mean") {
    cb <- fit_meanonly_combat(xx, batch_train, group_train)
    xx <- cb$X
    shifts <- cb$shifts
    cbdiag <- cb[c("max_deviation", "mean_abs_shift", "max_abs_shift")]
  }
  vars <- apply(xx, 2L, stats::var)
  good <- which(colSums(is.finite(xx)) == nrow(xx) &
                is.finite(vars) & vars > 0)
  if (length(good) < max(PANEL_GRID))
    stop("Too few proteins after zero-variance/nonfinite filter")
  xx <- xx[, good, drop=FALSE]
  mu <- colMeans(xx)
  ss <- apply(xx, 2L, stats::sd)
  scaled <- sweep(sweep(xx, 2L, mu, "-"), 2L, ss, "/")
  if (any(!is.finite(scaled))) stop("Nonfinite training standardization")
  list(X=scaled, kept=colnames(xx), mu=mu, sd=ss,
       minprob_sigma=imp$sigma, n_imputed=imp$n_imputed,
       mode=mode, shifts=shifts, combat_diagnostics=cbdiag)
}
score_preprocessor <- function(pre, X_test, batch_test, seed) {
  stopifnot(nrow(X_test) == length(batch_test))
  xx <- apply_minprob(X_test, pre$minprob_sigma, seed)
  if (pre$mode == "combat_mean") {
    batches <- as.character(batch_test)
    if (!all(batches %in% colnames(pre$shifts)))
      stop("Unseen held-out batch; cannot apply frozen ComBat")
    for (b in unique(batches)) {
      ii <- which(batches == b)
      xx[ii, ] <- sweep(xx[ii, , drop=FALSE], 2L,
                        pre$shifts[, b], "-")
    }
  }
  xx <- xx[, pre$kept, drop=FALSE]
  z <- sweep(sweep(xx, 2L, pre$mu, "-"), 2L, pre$sd, "/")
  if (any(!is.finite(z))) stop("Nonfinite held-out standardization")
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

# Outer-training-only inner grouped CV: all adaptive choices re-fitted.
inner_select <- function(X, y, family, batch, group, outer_fam, mode) {
  families <- sort(unique(as.character(family)))
  if (INNER_MODE == "lofo") {
    inner_fold <- match(as.character(family), families)
  } else {
    set.seed(family_seed(outer_fam, 10000L))
    inner_fold <- make_group_3fold(y, family)
    if (is.null(inner_fold)) stop("Unable to construct grouped three-fold inner CV")
  }
  m_grid <- if (PANEL_MODE == "fixed15") 15L else if (PANEL_MODE == "fixed30") 30L else PANEL_GRID
  predictions <- matrix(NA_real_, nrow=nrow(X), ncol=length(m_grid),
                        dimnames=list(rownames(X), as.character(m_grid)))
  audit <- vector("list", length(unique(inner_fold)))
  for (j in seq_along(sort(unique(inner_fold)))) {
    k <- sort(unique(inner_fold))[j]
    tr <- which(inner_fold != k); te <- which(inner_fold == k)
    if (sum(y[tr] == 1L) < 2L || sum(y[tr] == 0L) < 2L)
      stop("Inner split has <2 samples of a class; fold=", k)
    seed <- family_seed(paste(outer_fam, paste(sort(unique(family[te])), collapse = ","), sep = "|"), 5000000L)
    pre <- make_preprocessor(X[tr, , drop=FALSE],
                             group[tr], batch[tr], mode, seed+1L)
    rank <- stability_rank(pre$X, y[tr], family[tr], seed+2L)
    zd <- score_preprocessor(pre, X[te, , drop=FALSE], batch[te], seed+50L)
    measures <- vector("list", length(m_grid))
    for (mi in seq_along(m_grid)) {
      m <- m_grid[mi]
      ids <- head(rank$ids, m)
      fit <- fit_en_original(pre$X[, ids, drop=FALSE], y[tr], family[tr],
                             seed=seed+1000L+mi)
      predictions[te, as.character(m)] <- as.numeric(predict(
        fit$fit, newx=zd[, ids, drop=FALSE],
        s=fit$lambda, type="response"))
      measures[[mi]] <- data.frame(
        inner_fold=k, n_train=length(tr), pos_train=sum(y[tr] == 1L),
        n_held=length(te), m=m, lambda_method=fit$lambda_method,
        lambda=fit$lambda, reason=fit$fallback_reason,
        ranking_fallbacks=rank$n_fallback,
        ranking_success=rank$n_success,
        raised_pos_draws=rank$n_raised_pos,
        combat_max_delta_deviation=pre$combat_diagnostics$max_deviation,
        combat_mean_abs_shift=pre$combat_diagnostics$mean_abs_shift)
    }
    audit[[j]] <- do.call(rbind, measures)
  }
  if (any(!is.finite(predictions))) stop("Missing inner out-of-fold predictions")
  grid <- do.call(rbind, lapply(m_grid, function(m) {
    p <- predictions[, as.character(m)]
    tau <- select_threshold(p, y)
    data.frame(m=m, AUC=auc(y, p), sensitivity=tau$sens,
               specificity=tau$spec, threshold=tau$tau,
               threshold_policy=tau$policy)
  }))
  best <- order(-grid$AUC, -grid$sensitivity,
                -grid$specificity, grid$m)[1L]
  list(m=as.integer(grid$m[best]), threshold=grid$threshold[best],
       threshold_policy=grid$threshold_policy[best],
       metrics=grid, audit=do.call(rbind, audit))
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

# Two arms share identical family order, random seeds, and evaluation data.
run_outer <- function(dat, outdir, mode) {
  X <- dat$X; y <- dat$y; sm <- dat$metadata
  family <- as.character(sm$family_id)
  batch <- as.character(sm$batch)
  group <- as.character(sm$group)
  families_all <- sort(unique(family))
  families <- families_all
  if (OUTER_LIMIT > 0L) families <- head(families, OUTER_LIMIT)
  dir.create(file.path(outdir, "outer_cache"), recursive=TRUE, showWarnings=FALSE)
  results <- vector("list", length(families))
  for (i in seq_along(families)) {
    fam <- families[i]
    fn <- file.path(outdir, "outer_cache", paste0("family_", substr(sprintf("%08x", family_seed(fam, 0L)), 1, 8), ".rds"))
    if (file.exists(fn)) {
      one <- readRDS(fn)
      if (!identical(as.character(one$held_family), fam) ||
          !identical(as.character(one$mode), mode))
        stop("Incompatible cached fold: ", fn)
      results[[i]] <- one
      cat(sprintf("[%s] Outer %d/%d: family %s [cached]\n",
                  mode, i, length(families), fam))
      next
    }
    tr <- which(family != fam); te <- which(family == fam)
    cat(sprintf("[%s] Outer %d/%d: held family %s, n=%d (SM=%d)\n",
                mode, i, length(families), fam, length(te), sum(y[te] == 1L)))
    inner <- inner_select(X[tr, , drop=FALSE], y[tr], family[tr],
                          batch[tr], group[tr], outer_fam=fam, mode=mode)
    seed <- family_seed(fam, 10000000L)
    pre <- make_preprocessor(X[tr, , drop=FALSE], group[tr], batch[tr],
                             mode, seed+1L)
    rank <- stability_rank(pre$X, y[tr], family[tr], seed+2L)
    ids <- head(rank$ids, inner$m)
    fit <- fit_en_original(pre$X[, ids, drop=FALSE], y[tr], family[tr],
                           seed=seed+3L)
    zz <- score_preprocessor(pre, X[te, , drop=FALSE], batch[te], seed+4L)
    prob <- as.numeric(predict(fit$fit, newx=zz[, ids, drop=FALSE],
                               s=fit$lambda, type="response"))
    obs <- data.frame(SampleID=rownames(X)[te], family_id=family[te],
                      batch=batch[te], group=group[te], truth=y[te],
                      probability=prob, threshold=inner$threshold,
                      predicted_class=as.integer(prob >= inner$threshold),
                      selected_m=inner$m, selected_lambda=fit$lambda,
                      lambda_method=fit$lambda_method,
                      lambda_fallback_reason=fit$fallback_reason,
                      inner_threshold_policy=inner$threshold_policy)
    coefs <- extract_outer_coefficients(fit, ids, fam, rank$frequency)
    one <- list(mode=mode, held_family=fam, predictions=obs, coefficients=coefs,
                inner_grid=cbind(held_family=fam, inner$metrics),
                lambda_audit=cbind(held_family=fam, inner$audit),
                features=data.frame(held_family=fam, rank=seq_along(ids),
                                    FeatureID=ids,
                                    stability_frequency=unname(rank$frequency[ids])),
                stability=data.frame(held_family=fam,
                                     successful_refits=rank$n_success,
                                     fallback_lambda_refits=rank$n_fallback,
                                     positive_draw_adjustments=rank$n_raised_pos,
                                     failed_refits=rank$n_failed),
                combat=data.frame(held_family=fam, mode=mode,
                                  max_within_batch_shift_deviation=pre$combat_diagnostics$max_deviation,
                                  mean_abs_batch_shift=pre$combat_diagnostics$mean_abs_shift,
                                  max_abs_batch_shift=pre$combat_diagnostics$max_abs_shift))
    saveRDS(one, fn)
    results[[i]] <- one
    completed <- results[!vapply(results, is.null, logical(1))]
    write.csv(do.call(rbind, lapply(completed, `[[`, "predictions")),
              file.path(outdir, "outer_predictions_PARTIAL.csv"), row.names=FALSE)
  }
  if (length(families) < length(families_all)) {
    message("[SMOKE TEST ONLY] ", length(families), "/", length(families_all),
            " outer families; NOT a final performance estimate.")
    return(invisible(NULL))
  }
  pred <- do.call(rbind, lapply(results, `[[`, "predictions"))
  if (nrow(pred) != nrow(X) || anyDuplicated(pred$SampleID) ||
      !setequal(pred$SampleID, rownames(X)))
    stop("Outer predictions do not exactly match expected 44 samples")
  pred <- pred[match(rownames(X), pred$SampleID), , drop=FALSE]
  if (any(!is.finite(pred$probability))) stop("Non-finite outer predictions")
  write.csv(pred, file.path(outdir, "outer_predictions.csv"), row.names=FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "inner_grid")),
            file.path(outdir, "outer_inner_grid.csv"), row.names=FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "lambda_audit")),
            file.path(outdir, "outer_inner_lambda_audit.csv"), row.names=FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "features")),
            file.path(outdir, "outer_selected_features.csv"), row.names=FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "coefficients")),
            file.path(outdir, "outer_coefficients.csv"), row.names=FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "stability")),
            file.path(outdir, "outer_stability_diagnostics.csv"), row.names=FALSE)
  write.csv(do.call(rbind, lapply(results, `[[`, "combat")),
            file.path(outdir, "outer_combat_diagnostics.csv"), row.names=FALSE)
  total <- family_bootstrap_ci(pred, N_BOOT, SEED+77777L)
  write.csv(data.frame(metric=names(total$point),
                       estimate=as.numeric(total$point)),
            file.path(outdir, "outer_performance.csv"), row.names=FALSE)
  write.csv(total$ci, file.path(outdir, "outer_95CI.csv"), row.names=FALSE)
  run3_pred <- pred[pred$batch == "Run3", , drop=FALSE]
  within3 <- family_bootstrap_ci(run3_pred, N_BOOT, SEED+77778L)
  write.csv(data.frame(metric=names(within3$point),
                       estimate=as.numeric(within3$point)),
            file.path(outdir, "outer_run3_subset_performance.csv"), row.names=FALSE)
  write.csv(within3$ci, file.path(outdir, "outer_run3_subset_95CI.csv"),
            row.names=FALSE)
  cat("\n[", mode, "] OUTER RESULTS - ALL 44\n", sep="")
  print(data.frame(metric=names(total$point),
                   estimate=as.numeric(total$point)), row.names=FALSE)
  cat("\n[", mode, "] RUN3 SUBSET OF SAME 44-SAMPLE MODELS\n", sep="")
  print(data.frame(metric=names(within3$point),
                   estimate=as.numeric(within3$point)), row.names=FALSE)
  cat("\nSelected panel sizes (outer held-out family counts):\n")
  print(table(vapply(results, function(z) z$predictions$selected_m[1L], integer(1))))
  if (mode == "combat_mean") {
    dd <- do.call(rbind, lapply(results, `[[`, "combat"))
    cat("\nComBat frozen mean shift max deviation: ",
        signif(max(dd$max_within_batch_shift_deviation), 3), "\n", sep="")
  }
  invisible(list(all=total, run3=within3))
}

# Validate the reference (optional) and all paired predictions before comparing.
find_run3_prediction_file <- function() {
  fp <- Sys.getenv("CITRIN_RUN3_PREDICTIONS", "")
  if (nzchar(fp)) {
    if (dir.exists(fp)) fp <- file.path(fp, "outer_predictions.csv")
    if (!file.exists(fp)) stop("CITRIN_RUN3_PREDICTIONS missing: ", fp)
    return(fp)
  }
  roots <- unique(c(file.path(getwd(), "result", "run3_original_nested_cv"),
                    file.path(getwd(), "..", "result", "run3_original_nested_cv")))
  candidates <- unlist(lapply(roots[dir.exists(roots)], function(r) {
    list.files(r, pattern="^outer_predictions[.]csv$",
               full.names=TRUE, recursive=TRUE)
  }), use.names=FALSE)
  candidates <- unique(candidates)
  if (length(candidates) == 1L) return(candidates)
  if (length(candidates) > 1L)
    message("Multiple Run3 reference files found. Set CITRIN_RUN3_PREDICTIONS ",
            "explicitly; no automatic reference will be selected.")
  ""
}

read_predictions <- function(path) {
  if (!file.exists(path)) stop("Prediction file is missing: ", path)
  p <- read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)
  required <- c("SampleID", "family_id", "batch", "group", "truth",
                "probability", "threshold", "predicted_class")
  if (!all(required %in% names(p)))
    stop("Missing columns in: ", path, ": ",
         paste(setdiff(required, names(p)), collapse=", "))
  p$SampleID <- normalize_id(p$SampleID)
  p$family_id <- as.character(p$family_id)
  if (anyDuplicated(p$SampleID) || any(!is.finite(p$probability)) ||
      any(!p$truth %in% c(0, 1)) || any(!p$predicted_class %in% c(0, 1)))
    stop("Invalid predictions: ", path)
  p
}

match_predictions <- function(ref, alternative, expect_same_samples=TRUE) {
  if (expect_same_samples && !setequal(ref$SampleID, alternative$SampleID))
    stop("Predictions have different SampleIDs; cannot compare as paired data")
  match_idx <- match(ref$SampleID, alternative$SampleID)
  if (anyNA(match_idx)) stop("Reference has missing predictions in other analysis")
  a <- alternative[match_idx, , drop=FALSE]
  if (!identical(as.integer(ref$truth), as.integer(a$truth)) ||
      !identical(as.character(ref$family_id), as.character(a$family_id)) ||
      !identical(as.character(ref$group), as.character(a$group)) ||
      !identical(as.character(ref$batch), as.character(a$batch)))
    stop("Group/family/batch mismatches between predictions. Resolve before comparison.")
  a
}

metrics_row <- function(label, p) {
  mm <- metrics(p$truth, p$probability, p$predicted_class)
  data.frame(model=label, n=nrow(p), n_positive=sum(p$truth == 1L),
             n_negative=sum(p$truth == 0L), n_families=length(unique(p$family_id)),
             AUC=mm["AUC"], sensitivity=mm["sensitivity"],
             specificity=mm["specificity"], balanced_accuracy=mm["balanced_accuracy"],
             accuracy=mm["accuracy"], Brier=mm["Brier"],
             TP=mm["TP"], TN=mm["TN"], FP=mm["FP"], FN=mm["FN"],
             row.names=NULL)
}

# Fixed out-of-fold predictions: paired family bootstrap, not model refitting.
paired_auc_delta <- function(p1, p2, reps, seed) {
  p2 <- match_predictions(p1, p2)
  est <- auc(p1$truth, p1$probability) - auc(p2$truth, p2$probability)
  if (reps <= 0L) return(c(delta_AUC=est, lower=NA_real_, upper=NA_real_, n_valid=0))
  fams <- unique(p1$family_id)
  by_family <- split(seq_len(nrow(p1)), p1$family_id)
  diffs <- rep(NA_real_, reps)
  set.seed(seed)
  for (k in seq_len(reps)) {
    sample_fams <- sample(fams, length(fams), replace=TRUE)
    idx <- unlist(by_family[sample_fams], use.names=FALSE)
    if (length(unique(p1$truth[idx])) < 2L) next
    diffs[k] <- auc(p1$truth[idx], p1$probability[idx]) -
                auc(p2$truth[idx], p2$probability[idx])
  }
  vv <- diffs[is.finite(diffs)]
  if (length(vv) <= 20L)
    return(c(delta_AUC=est, lower=NA_real_, upper=NA_real_, n_valid=length(vv)))
  c(delta_AUC=est, lower=unname(quantile(vv, 0.025)),
    upper=unname(quantile(vv, 0.975)), n_valid=length(vv))
}

run_comparison <- function(dat, base_dir) {
  stage1 <- file.path(base_dir, "stage1_no_combat", "outer_predictions.csv")
  stage2 <- file.path(base_dir, "stage2_frozen_mean_combat", "outer_predictions.csv")
  if (!file.exists(stage1) || !file.exists(stage2)) {
    message("[STAGE3] Both completed 44-sample stages are required. ",
            "No comparison written yet.")
    return(invisible(NULL))
  }
  a <- read_predictions(stage1)
  b <- match_predictions(a, read_predictions(stage2))
  if (!setequal(a$SampleID, rownames(dat$X)) || nrow(a) != nrow(dat$X))
    stop("Predictions do not cover exactly the current 44-sample cohort")
  run3 <- a$batch == "Run3"
  a3 <- a[run3, , drop=FALSE]
  b3 <- b[run3, , drop=FALSE]
  rows <- list(metrics_row("stage1_no_combat_all44", a),
               metrics_row("stage2_meanonly_combat_all44", b),
               metrics_row("stage1_no_combat_run3_subset", a3),
               metrics_row("stage2_meanonly_combat_run3_subset", b3))
  paired_deltas <- list()
  mk_delta <- function(label, d1, d2, seed) {
    data.frame(comparison=label, as.list(paired_auc_delta(d1, d2, N_BOOT, seed)),
               row.names=NULL)
  }
  paired_deltas[[1L]] <- mk_delta("stage2 minus stage1 (all44)", b, a, SEED+9001L)
  paired_deltas[[2L]] <- mk_delta("stage2 minus stage1 (Run3)", b3, a3, SEED+9002L)
  run3_fp <- find_run3_prediction_file()
  if (nzchar(run3_fp)) {
    ref <- read_predictions(run3_fp)
    ref <- match_predictions(a3, ref)
    rows[[length(rows)+1L]] <- metrics_row("previous_run3_only_model", ref)
    paired_deltas[[length(paired_deltas)+1L]] <-
      mk_delta("stage1 Run3 subset minus Run3-only model", a3, ref, SEED+9003L)
    paired_deltas[[length(paired_deltas)+1L]] <-
      mk_delta("stage2 Run3 subset minus Run3-only model", b3, ref, SEED+9004L)
  } else {
    message("[STAGE3] Run3-only reference not specified/found; ",
            "set CITRIN_RUN3_PREDICTIONS and rerun with CITRIN_STAGE='compare'.")
  }
  comp_dir <- file.path(base_dir, "stage3_comparison")
  dir.create(comp_dir, recursive=TRUE, showWarnings=FALSE)
  summary <- do.call(rbind, rows)
  write.csv(summary, file.path(comp_dir, "model_performance_comparison.csv"),
            row.names=FALSE)
  # Cluster-bootstrap intervals for every model subset on fixed held-out scores.
  comparison_data <- list(a, b, a3, b3)
  if (nzchar(run3_fp)) comparison_data[[5L]] <- ref
  cis <- lapply(seq_along(comparison_data), function(j) {
    ci <- family_bootstrap_ci(comparison_data[[j]], N_BOOT,
                              SEED+93000L+j)$ci
    cbind(model=summary$model[j], ci)
  })
  write.csv(do.call(rbind, cis),
            file.path(comp_dir, "model_performance_family_CI.csv"),
            row.names=FALSE)
  write.csv(do.call(rbind, paired_deltas),
            file.path(comp_dir, "paired_family_bootstrap_auc_differences.csv"),
            row.names=FALSE)
  # All paired held-out predictions: inspect FN/FP, especially Run2 controls.
  paired <- data.frame(SampleID=a$SampleID, family_id=a$family_id,
                       batch=a$batch, group=a$group, truth=a$truth,
                       stage1_prob=a$probability, stage1_threshold=a$threshold,
                       stage1_pred=a$predicted_class,
                       stage1_m=a$selected_m,
                       stage2_prob=b$probability, stage2_threshold=b$threshold,
                       stage2_pred=b$predicted_class,
                       stage2_m=b$selected_m)
  write.csv(paired, file.path(comp_dir, "paired_44_predictions.csv"),
            row.names=FALSE)
  if (nzchar(run3_fp)) {
    ref <- match_predictions(a3, read_predictions(run3_fp))
    paired3 <- paired[paired$batch == "Run3", , drop=FALSE]
    paired3$run3_only_prob <- ref$probability
    paired3$run3_only_threshold <- ref$threshold
    paired3$run3_only_pred <- ref$predicted_class
    paired3$run3_only_m <- ref$selected_m
    write.csv(paired3, file.path(comp_dir, "paired_run3_predictions.csv"),
              row.names=FALSE)
  }
  run2_controls <- paired[paired$batch == "Run2", , drop=FALSE]
  write.csv(data.frame(model=c("stage1_no_combat", "stage2_meanonly_combat"),
                       Run2_n=nrow(run2_controls),
                       Run2_FP=c(sum(run2_controls$stage1_pred == 1L),
                                 sum(run2_controls$stage2_pred == 1L)),
                       Run2_specificity=c(mean(run2_controls$stage1_pred == 0L),
                                          mean(run2_controls$stage2_pred == 0L))),
            file.path(comp_dir, "run2_negative_only_errors.csv"), row.names=FALSE)
  # Report which panel sizes were selected for each held-out family.
  fampanel <- data.frame(held_family=sort(unique(a$family_id)),
                         stage1_m=NA_integer_, stage2_m=NA_integer_)
  for (j in seq_len(nrow(fampanel))) {
    f <- fampanel$held_family[j]
    aa <- unique(a$selected_m[a$family_id == f])
    bb <- unique(b$selected_m[b$family_id == f])
    if (length(aa) != 1L || length(bb) != 1L)
      stop("Same family yielded multiple selected panel sizes, ID=", f)
    fampanel$stage1_m[j] <- aa
    fampanel$stage2_m[j] <- bb
  }
  write.csv(fampanel,
            file.path(comp_dir, "family_fold_selected_panel_sizes.csv"),
            row.names=FALSE)
  for (mode in c("stage1_no_combat", "stage2_frozen_mean_combat")) {
    d <- file.path(base_dir, mode)
    lf <- read.csv(file.path(d, "outer_inner_lambda_audit.csv"))
    st <- read.csv(file.path(d, "outer_stability_diagnostics.csv"))
    pp <- read_predictions(file.path(d, "outer_predictions.csv"))
    unique_inner <- lf[!duplicated(paste(lf$held_family, lf$inner_fold)), , drop=FALSE]
    audit <- data.frame(model=mode, n_inner_panel_lambda_fits=nrow(lf),
                        n_lambda_fallbacks=sum(grepl("^fallback", lf$lambda_method)),
                        n_outer_lambda_fallbacks=sum(grepl("^fallback", pp$lambda_method)),
                        n_pooled_auc_lambda_fits=sum(lf$lambda_method == "grouped_3fold_pooled_auc"),
                        n_inner_stability_lambda_fallbacks=sum(unique_inner$ranking_fallbacks),
                        n_inner_stability_positive_draw_adjustments=sum(unique_inner$raised_pos_draws),
                        n_outer_stability_lambda_fallbacks=sum(st$fallback_lambda_refits),
                        n_outer_stability_positive_draw_adjustments=sum(st$positive_draw_adjustments),
                        n_stability_failed_refits=sum(st$failed_refits),
                        n_outer_families=nrow(st))
    write.csv(audit, file.path(d, "diagnostic_summary.csv"), row.names=FALSE)
  }
  cat("\n[STAGE 3] Performance comparison:\n")
  print(summary, row.names=FALSE, digits=4)
  cat("\n[STAGE 3] Paired AUC differences (bootstrap uses fixed out-of-fold predictions):\n")
  print(do.call(rbind, paired_deltas), row.names=FALSE, digits=4)
  cat("\n[STAGE 3] Results at: ", comp_dir, "\n", sep="")
  invisible(summary)
}

hash_config <- function(dat) {
  fhash <- unname(tools::md5sum(c(RAW_FILE, SAMPLE_FILE)))
  cfg <- paste(SCRIPT_VERSION, paste(fhash, collapse="|"),
               paste(dat$metadata$SampleID, dat$metadata$family_id,
                     dat$metadata$group, dat$metadata$batch,
                     sep=":", collapse="|"),
               SEED, B_STABILITY, PANEL_MODE, INNER_MODE, SUBSAMPLE_UNIT,
               SPEC_FLOOR, MINPROB_Q, LAMBDA_FALLBACK_RATIO,
               RUN2, RUN3, sep="|")
  temp <- tempfile(); writeLines(cfg, temp, useBytes=TRUE)
  on.exit(unlink(temp), add=TRUE)
  substr(unname(tools::md5sum(temp)), 1L, 12L)
}

main <- function() {
  dat <- read_44()
  dir.create(OUT_BASE, recursive=TRUE, showWarnings=FALSE)
  key <- hash_config(dat)
  base_dir <- file.path(OUT_BASE,
                        paste0("B", B_STABILITY, "_", PANEL_MODE,
                               "_", INNER_MODE, "_", key))
  dir.create(base_dir, recursive=TRUE, showWarnings=FALSE)
  config <- data.frame(name=c("script_version", "data_dir", "n_samples",
                              "n_positive_SM", "n_negative", "n_families",
                              "run2_n", "run3_n", "n_proteins",
                              "stability_B", "panel_mode", "inner_mode", "subsample_unit",
                              "seed", "lambda_fallback_ratio", "minprob_q",
                              "min_unique_peptides", "stage2_combat",
                              "CV_outer", "stage1_combat", "run3_reference_path"),
                       value=c(SCRIPT_VERSION, DATA_DIR, nrow(dat$X),
                               sum(dat$y), sum(dat$y == 0L),
                               length(unique(dat$metadata$family_id)),
                               sum(dat$metadata$batch == "Run2"),
                               sum(dat$metadata$batch == "Run3"),
                               ncol(dat$X), B_STABILITY, PANEL_MODE,
                               INNER_MODE, SUBSAMPLE_UNIT, SEED, LAMBDA_FALLBACK_RATIO,
                               MINPROB_Q, MIN_UNIQUE_PEPTIDES,
                               "sva::ComBat(mod=~group, mean.only=TRUE, par.prior=TRUE); frozen shifts",
                               "leave_one_family_out", "none",
                               Sys.getenv("CITRIN_RUN3_PREDICTIONS", "")))
  write.csv(config, file.path(base_dir, "analysis_parameters.csv"), row.names=FALSE)
  write.csv(dat$metadata[, c("SampleID", "group", "batch", "family_id")],
            file.path(base_dir, "cohort44_used.csv"), row.names=FALSE)
  capture.output(sessionInfo(), file=file.path(base_dir, "sessionInfo.txt"))
  cat("\n[OUTPUT] ", normalizePath(base_dir, mustWork=FALSE), "\n", sep="")
  if (STAGE %in% c("stage1", "both")) {
    d <- file.path(base_dir, "stage1_no_combat")
    dir.create(d, recursive=TRUE, showWarnings=FALSE)
    run_outer(dat, d, "no_combat")
  }
  if (STAGE %in% c("stage2", "both")) {
    if (!requireNamespace("sva", quietly=TRUE))
      stop("Stage2 requires Bioconductor sva: BiocManager::install('sva')")
    d <- file.path(base_dir, "stage2_frozen_mean_combat")
    dir.create(d, recursive=TRUE, showWarnings=FALSE)
    run_outer(dat, d, "combat_mean")
  }
  if (STAGE %in% c("both", "compare")) {
    if (OUTER_LIMIT > 0L) {
      message("[STAGE3] Comparison skipped for smoke-test outer_limit>0")
    } else {
      run_comparison(dat, base_dir)
    }
  }
  cat("\n[FINISHED] Check output files under ", base_dir, "\n", sep="")
  invisible(base_dir)
}

# Explicit source() runs main() once, matching previous Run3 script.
main()
