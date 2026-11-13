#!/usr/bin/env Rscript
# ============================================================================
# CITRIN: NESTED, FAMILY-HELD-OUT, SENSITIVITY-FIRST 30-PROTEIN ENSEMBLE
# v1.1 - 2026-11-02
#
# This is NEW, post-hoc algorithm development, not a replacement for the
# predeclared v2.2 Fixed-30 validation. Do NOT claim independent external
# validation from this same cohort. No outer-test outcome is used for feature
# selection, ensemble choice, probability threshold, or model fitting.
#
# PRIMARY: consensus of 30 proteins selected WITHIN each training partition
#          across 10 independently seeded stability rankings; an average of
#          10 elastic-net classifiers fitted on that SAME set of 30 proteins.
# SECONDARY: diverse30 (each of 10 models uses its own 30-protein subset),
#          adaptive (chooses between these two using inner CV only), and
#          an alternative >=90% inner specificity decision rule.
#
# Threshold objective: 100% INNER-OOF sample recall; among such thresholds,
# choose the largest, hence the fewest INNER-OOF false positives for a fixed
# score. It does NOT guarantee 100% OUTER-OOF sensitivity.
#
# Use v2.2 training-fitted MinProb, frozen mean-only ComBat (optional),
# ACTB normalization, stability ranking, elastic-net, and lambda fitting.
# Upstream DIA-NN quantification is not cross-validated.
#
# Parallel use: independent workers each handle disjoint outer families.
# A single post-run aggregation process gathers completed family caches.
# ============================================================================

options(stringsAsFactors = FALSE)
VERSION <- "citrin_nested_ensemble_familyLOFO_1.1"
env <- function(key, default) { z <- Sys.getenv(key, ""); if (nzchar(z)) z else default }
canon <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
parse_words <- function(s) { trimws(strsplit(s, ",", fixed = TRUE)[[1L]]) }

PROJ <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
DATA <- normalizePath(env("CITRIN_DATA_DIR", file.path(dirname(PROJ), "data")),
                      winslash = "/", mustWork = TRUE)
FAMILY_FILE <- normalizePath(env("CITRIN_FAMILY_FILE", file.path(DATA, "family_map.csv")),
                             winslash = "/", mustWork = TRUE)
HARNESS_FILE <- normalizePath(file.path(PROJ, "citrin_nested_ensemble_familyLOFO_v1_1.R"),
                              winslash = "/", mustWork = TRUE)
ENGINE_CANDIDATES <- c(file.path(PROJ, "citrin_44_nested_fixed_panel_v2_2.R"),
                       file.path(PROJ, "Claude", "citrin_44_nested_fixed_panel_v2_2.R"),
                       file.path(getwd(), "citrin_44_nested_fixed_panel_v2_2.R"))
ENGINE_CANDIDATES <- ENGINE_CANDIDATES[file.exists(ENGINE_CANDIDATES)]
ENGINE_PATH <- normalizePath(env("CITRIN_ENSEMBLE_ENGINE",
                                 if (length(ENGINE_CANDIDATES)) ENGINE_CANDIDATES[1L] else ""),
                             winslash = "/", mustWork = TRUE)
OUTPUT <- normalizePath(env("CITRIN_ENSEMBLE_OUTPUT",
                             file.path(PROJ, "results_rerun", "citrin_nested_ensemble_v1_1_lofo")),
                        winslash = "/", mustWork = FALSE)
MODE <- tolower(env("CITRIN_ENSEMBLE_COMBAT", "no_combat"))
SEEDS <- as.integer(parse_words(env("CITRIN_ENSEMBLE_SEEDS",
                                    "2025,2026,2027,2028,2029,2030,2031,2032,2033,2034")))
B_STABILITY <- as.integer(env("CITRIN_ENSEMBLE_STABILITY_B", "30"))
SUBSAMPLE <- tolower(env("CITRIN_ENSEMBLE_SUBSAMPLE", "sample"))
# Family LOFO is mandatory in v1.1; the obsolete 4-fold setting is ignored.
if (nzchar(Sys.getenv("CITRIN_ENSEMBLE_INNER_FOLDS", "")))
  message("[INFO] Ignoring obsolete CITRIN_ENSEMBLE_INNER_FOLDS: v1.1 uses inner family-LOFO only.")
INNER_VALIDATION <- "family_lofo"
OUTER_LIMIT <- as.integer(env("CITRIN_ENSEMBLE_OUTER_LIMIT", "0"))
WORKERS <- as.integer(env("CITRIN_ENSEMBLE_WORKERS", "1"))
WORKER <- as.integer(env("CITRIN_ENSEMBLE_WORKER", "1"))
AGGREGATE_ONLY <- identical(env("CITRIN_ENSEMBLE_AGGREGATE_ONLY", "false"), "true")
FINAL_MODEL <- identical(env("CITRIN_ENSEMBLE_FIT_FINAL", "true"), "true")
PANEL_SIZE <- 30L
PANEL_NAMES <- c("consensus30", "diverse30")
POLICIES <- c("full_inner_recall", "inner_spec90")
ORIGINAL15 <- c("Q06203", "Q8WUD1", "Q99829", "P01236", "Q9UHG3", "O15013",
                "O94979", "P55157", "Q9UBQ7", "Q9H115", "P78509", "P06733",
                "Q9Y5S1", "Q96BS2", "P56556")

stopifnot(MODE %in% c("no_combat", "combat_mean"),
          SUBSAMPLE %in% c("sample", "family"),
          length(SEEDS) >= 2L, !anyNA(SEEDS), !anyDuplicated(SEEDS), all(SEEDS >= 0L),
          B_STABILITY >= 1L, is.finite(B_STABILITY),
          is.finite(OUTER_LIMIT), OUTER_LIMIT >= 0L,
          is.finite(WORKERS), WORKERS >= 1L, is.finite(WORKER), WORKER >= 1L,
          WORKER <= WORKERS)
stopifnot(file.exists(file.path(DATA, "raw_data.xlsx")),
          file.exists(file.path(DATA, "sample_info.csv")))

# Source the EXACT v2.2 engine definitions, explicitly skipping its main().
# No changes to the engine, no contamination by its old results/cache.
Sys.setenv(CITRIN_DATA_DIR = DATA, CITRIN_FAMILY_FILE = FAMILY_FILE,
           CITRIN_STAGE = if (MODE == "no_combat") "stage1" else "stage2",
           CITRIN_PANEL_MODE = "fixed30", CITRIN_SUBSAMPLE_UNIT = SUBSAMPLE,
           CITRIN_B = as.character(B_STABILITY), CITRIN_SEED = "2025",
           CITRIN_N_BOOT = "0", CITRIN_OUTER_LIMIT = "0")
exprs <- parse(file = ENGINE_PATH, keep.source = FALSE)
if (!identical(deparse(exprs[[length(exprs)]]), "main()"))
  stop("Engine source no longer ends in main(): inspect before reusing functions.")
engine <- new.env(parent = globalenv())
for (j in seq_len(length(exprs) - 1L)) eval(exprs[[j]], envir = engine)
if (!identical(engine$SCRIPT_VERSION, "two_run_fixed_panel_2.2"))
  stop("Expected unmodified 44-sample engine v2.2")

# A family-held-out dataset is loaded, but only training labels are supplied
# to ALL adaptive steps below. The test's y values are read only for reporting.
dat <- engine$read_44()
X <- dat$X; y <- as.integer(dat$y)
meta <- dat$metadata
family <- as.character(meta$family_id)
batch <- as.character(meta$batch)
group <- as.character(meta$group)
fams <- sort(unique(family))
stopifnot(nrow(X) == 44L, sum(y == 1L) == 7L,
          length(fams) == 35L, length(unique(family[y == 1L])) == 5L,
          identical(rownames(X), meta$SampleID))

# The immutable plan covers all inputs and choices that affect the estimates.
dir.create(OUTPUT, recursive = TRUE, showWarnings = FALSE)
kv <- data.frame(name = c("version", "harness_md5", "engine_md5", "raw_md5", "sample_info_md5",
                          "family_md5", "engine_path", "data_dir", "family_file",
                          "mode", "seeds", "stability_B", "subsample", "inner_validation",
                          "panel_size", "outer_limit", "n_outer_families", "positive_families",
                          "optimization", "policies"),
                 value = as.character(c(VERSION, unname(tools::md5sum(HARNESS_FILE)),
                   unname(tools::md5sum(ENGINE_PATH)),
                   unname(tools::md5sum(file.path(DATA, "raw_data.xlsx"))),
                   unname(tools::md5sum(file.path(DATA, "sample_info.csv"))),
                   unname(tools::md5sum(FAMILY_FILE)), ENGINE_PATH, DATA, FAMILY_FILE,
                   MODE, paste(SEEDS, collapse = ","), B_STABILITY, SUBSAMPLE,
                   INNER_VALIDATION, PANEL_SIZE, OUTER_LIMIT, length(fams),
                   length(unique(family[y == 1L])),
                   "100pct_inner_sample_recall_min_FP_within_threshold_family_grouped",
                   paste(POLICIES, collapse = ","))), stringsAsFactors = FALSE)

# Process-owned directory locks. On the same host, a dead PID's abandoned
# lock may be recovered. Unknown/missing-owner and foreign-host locks are NOT
# deleted automatically. This avoids removing a potentially live lock.
HOST_ID <- as.character(Sys.info()[["nodename"]])
process_alive <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (length(pid) != 1L || is.na(pid) || pid <= 1L) return(NA)
  ans <- suppressWarnings(tryCatch(
    system2("/bin/kill", c("-0", as.character(pid)), stdout = FALSE, stderr = FALSE),
    error = function(e) NA_integer_))
  if (is.na(ans)) return(NA)
  identical(as.integer(ans), 0L)
}
read_lock_owner <- function(path) {
  file <- file.path(path, "owner.txt")
  zz <- tryCatch(readLines(file, warn = FALSE), error = function(e) character())
  if (length(zz) < 3L) return(NULL)
  if (!startsWith(zz[1L], "host=") || !startsWith(zz[2L], "pid=")) return(NULL)
  list(lines = zz, host = sub("^host=", "", zz[1L]),
       pid = suppressWarnings(as.integer(sub("^pid=", "", zz[2L]))))
}
write_lock_owner <- function(path) {
  fp <- file.path(path, "owner.txt")
  tmp <- tempfile(".lock_owner_", tmpdir = path)
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(paste0("host=", HOST_ID), paste0("pid=", Sys.getpid()),
               paste0("utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))),
             tmp, useBytes = TRUE)
  if (!file.rename(tmp, fp)) stop("Unable to record lock owner: ", path)
}
recover_dead_lock <- function(path) {
  ow <- read_lock_owner(path)
  if (is.null(ow) || !identical(ow$host, HOST_ID) || !identical(process_alive(ow$pid), FALSE))
    return(FALSE)
  # Atomic rename wins the claim to a stale lock. Recheck owner inside moved
  # directory before removing it; never delete a different owner's lock.
  trash <- paste0(path, ".abandoned_by_", Sys.getpid(), "_", sample.int(999999L, 1L))
  if (!file.rename(path, trash)) return(FALSE)
  after <- read_lock_owner(trash)
  if (is.null(after) || !identical(ow$lines, after$lines) ||
      !identical(process_alive(after$pid), FALSE)) {
    if (!dir.exists(path)) file.rename(trash, path)
    stop("Owner changed during stale-lock recovery; inspect: ", trash)
  }
  unlink(trash, recursive = TRUE)
  message("[RECOVERED DEAD LOCK] ", path, " (PID ", ow$pid, ")")
  TRUE
}
acquire_dir_lock <- function(path, timeout = 120) {
  t0 <- Sys.time()
  repeat {
    if (dir.create(path, showWarnings = FALSE)) {
      tryCatch(write_lock_owner(path), error = function(e) {
        unlink(path, recursive = TRUE)
        stop(e)
      })
      return(invisible(TRUE))
    }
    if (recover_dead_lock(path)) next
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > timeout) {
      owner <- read_lock_owner(path)
      msg <- if (is.null(owner)) "owner missing/unreadable" else
        paste0("host=", owner$host, " pid=", owner$pid,
               " alive=", process_alive(owner$pid))
      stop("Timed out waiting for lock: ", path, " [", msg,
           "]. If another worker is alive, do not remove the lock. ",
           "If ownership cannot be verified, inspect manually.")
    }
    Sys.sleep(0.5)
  }
}

atomic_csv <- function(z, path) {
  tmp <- tempfile(".tmp_", tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(z, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Atomic file rename failed: ", path)
}
atomic_rds <- function(z, path) {
  tmp <- tempfile(".tmp_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(z, tmp)
  if (!file.rename(tmp, path)) stop("Atomic RDS rename failed: ", path)
}
ensure_plan <- function() {
  lock <- file.path(OUTPUT, ".plan_lock")
  acquire_dir_lock(lock)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  fp <- file.path(OUTPUT, "immutable_analysis_plan.csv")
  if (file.exists(fp)) {
    old <- utils::read.csv(fp, colClasses = "character", check.names = FALSE)
    if (!isTRUE(all.equal(old, kv, check.attributes = FALSE)))
      stop("Input or analytic plan changed. Set a NEW CITRIN_ENSEMBLE_OUTPUT directory; ",
           "do not mix protocols or datasets.")
  } else {
    atomic_csv(kv, fp)
  }
}
ensure_plan()
PLAN_MD5 <- unname(tools::md5sum(file.path(OUTPUT, "immutable_analysis_plan.csv")))

# Integer-safe, deterministic seeds keyed to exact family/context/seed member.
key_seed <- function(...) engine$family_seed(paste(..., sep = "|"), 11111L)

# Strict inner leave-one-FAMILY-out: a distinct fold for each remaining family.
# For an outer-positive-family holdout, inner models have three positive
# training families; for outer-negative-family holdout, they have four.
# The training size is consequently close to the outer fit size and, in this
# cohort, stays >=30 (the cv.glmnet grouped-AUC lambda-selection branch).
make_inner_folds <- function(yy, ff, bb, context) {
  ff <- as.character(ff)
  lev <- sort(unique(ff))
  if (length(unique(ff[yy == 1L])) < 4L)
    stop("Inner family LOFO requires at least 4 positive training families: ", context)
  z <- match(ff, lev)
  if (anyNA(z) || length(unique(z)) != length(lev)) stop("Invalid family fold IDs")
  for (i in seq_along(lev)) {
    tr <- which(z != i); te <- which(z == i)
    if (!length(te) || any(ff[tr] %in% ff[te]) ||
        length(unique(ff[tr][yy[tr] == 1L])) < 3L ||
        sum(yy[tr] == 0L) < 2L)
      stop("Ineligible family-LOFO inner fold ", i, " in ", context)
    if (length(tr) < 30L)
      stop("Inner LOFO has n<30 in ", context, ": fold ", i,
           "; lambda-selection branch would change; inspect data")
    if (MODE == "combat_mean" &&
        (length(unique(bb[tr])) != 2L ||
         any(table(factor(bb[tr], levels = c("Run2", "Run3"))) < 2L)))
      stop("Insufficient training batches for frozen ComBat in fold ", i)
  }
  z
}

# Exact training-only sensitivity-first threshold, with optimal FP among all
# thresholds attaining 100% inner-OOF *sample* sensitivity at fixed scores.
full_recall_threshold <- function(s, yy) {
  if (any(!is.finite(s)) || !any(yy == 1L) || !any(yy == 0L)) stop("Invalid inner score/labels")
  tau <- min(s[yy == 1L])
  pred <- as.integer(s >= tau)
  stopifnot(all(pred[yy == 1L] == 1L))
  list(tau = tau, tp = sum(pred[yy == 1L]), fp = sum(pred[yy == 0L]),
       sens = 1, spec = mean(pred[yy == 0L] == 0L))
}
# Comparator policy, fitted using the SAME inner-OOF predictions. A floor of
# 0.9 is an internal training specification, not guaranteed outer specificity.
spec90_threshold <- function(s, yy) {
  grid <- sort(unique(c(0, s, 1 + 1e-12)))
  nneg <- sum(yy == 0L)
  candidates <- lapply(grid, function(t) {
    d <- s >= t
    c(tau = t, tp = sum(d & yy == 1L), fp = sum(d & yy == 0L))
  })
  g <- as.data.frame(do.call(rbind, candidates))
  g <- g[g$fp / nneg <= 0.10 + 1e-12, , drop = FALSE]
  if (!nrow(g)) stop("No inner threshold meeting specificity 90%")
  best <- order(-g$tp, g$fp, g$tau)[1L]
  z <- g[best, , drop = FALSE]
  list(tau = z$tau, tp = z$tp, fp = z$fp,
       sens = z$tp / sum(yy == 1L), spec = 1 - z$fp / nneg)
}

# Fit two alternatives from exactly the same K stability rankings, using ONLY
# the supplied training subset and WITHOUT test labels.
fit_ensemble <- function(xtr, ytr, ftr, btr, gtr, xnew, bnew,
                         context, keep_objects = FALSE) {
  stopifnot(nrow(xtr) == length(ytr), nrow(xnew) == length(bnew),
            !any(rownames(xtr) %in% rownames(xnew)),
            length(unique(ftr[ytr == 1L])) >= 3L)
  K <- length(SEEDS)
  bases <- vector("list", K)
  div_prob <- matrix(NA_real_, nrow(xnew), K,
                     dimnames = list(rownames(xnew), as.character(SEEDS)))
  member_diag <- vector("list", K)
  diverse_features <- vector("list", K)
  diverse_coefficients <- vector("list", K)
  for (k in seq_len(K)) {
    seed <- key_seed("fit", context, SEEDS[k])
    pre <- engine$make_preprocessor(xtr, gtr, btr, MODE, seed + 1L)
    rank <- engine$stability_rank(pre$X, ytr, ftr, seed + 2L)
    ids <- head(rank$ids, PANEL_SIZE)
    if (length(ids) != PANEL_SIZE) stop("Fewer than 30 available proteins")
    zz <- engine$score_preprocessor(pre, xnew, bnew, seed + 4L)
    fit <- engine$fit_en_original(pre$X[, ids, drop = FALSE], ytr, ftr,
                                  seed = seed + 5L)
    div_prob[, k] <- as.numeric(stats::predict(
      fit$fit, newx = zz[, ids, drop = FALSE], s = fit$lambda, type = "response"))
    bases[[k]] <- list(seed_member = SEEDS[k], seed = seed, pre = pre, zz = zz,
                       rank = rank, diverse_ids = ids,
                       diverse_fit = if (keep_objects) fit else NULL)
    diverse_features[[k]] <- data.frame(
      member_seed = SEEDS[k], rank = seq_along(ids), FeatureID = ids,
      stability_frequency = as.numeric(rank$frequency[ids]))
    cf_div <- as.matrix(stats::coef(fit$fit, s = fit$lambda))
    diverse_coefficients[[k]] <- data.frame(
      member_seed = SEEDS[k], FeatureID = rownames(cf_div),
      coefficient = as.numeric(cf_div[, 1L]),
      nonzero = as.numeric(cf_div[, 1L]) != 0,
      lambda_method = fit$lambda_method, lambda = fit$lambda)
    member_diag[[k]] <- data.frame(
      member_seed = SEEDS[k], stability_success = rank$n_success,
      stability_lambda_fallbacks = rank$n_fallback, stability_failures = rank$n_failed,
      diverse_lambda_method = fit$lambda_method, diverse_lambda = fit$lambda)
  }
  # Select ONE shared 30-protein list from the K training rankings. A protein
  # must be usable by every member's independently fitted preprocessor.
  eligible <- Reduce(intersect, lapply(bases, function(z) colnames(z$pre$X)))
  if (length(eligible) < PANEL_SIZE) stop("<30 mutually usable features")
  counts <- setNames(integer(length(eligible)), eligible)
  avg_stab <- setNames(numeric(length(eligible)), eligible)
  for (z in bases) {
    ids <- intersect(z$diverse_ids, eligible)
    counts[ids] <- counts[ids] + 1L
    avg_stab <- avg_stab + z$rank$frequency[eligible] / K
  }
  oo <- order(-counts, -avg_stab, names(counts))
  selected <- names(counts)[oo[seq_len(PANEL_SIZE)]]
  sel <- data.frame(rank = seq_along(selected), FeatureID = selected,
                    members_top30 = as.integer(counts[selected]),
                    mean_stability_frequency = as.numeric(avg_stab[selected]),
                    in_original15_exploratory = selected %in% ORIGINAL15)
  con_prob <- matrix(NA_real_, nrow(xnew), K,
                     dimnames = list(rownames(xnew), as.character(SEEDS)))
  member_coefs <- vector("list", K)
  for (k in seq_len(K)) {
    z <- bases[[k]]
    fit <- engine$fit_en_original(z$pre$X[, selected, drop = FALSE], ytr, ftr,
                                  seed = z$seed + 6L)
    con_prob[, k] <- as.numeric(stats::predict(
      fit$fit, newx = z$zz[, selected, drop = FALSE],
      s = fit$lambda, type = "response"))
    cf <- as.matrix(stats::coef(fit$fit, s = fit$lambda))
    member_coefs[[k]] <- data.frame(member_seed = SEEDS[k], FeatureID = rownames(cf),
                                    coefficient = as.numeric(cf[, 1L]),
                                    nonzero = as.numeric(cf[, 1L]) != 0,
                                    lambda_method = fit$lambda_method,
                                    lambda = fit$lambda)
    member_diag[[k]]$consensus_lambda_method <- fit$lambda_method
    member_diag[[k]]$consensus_lambda <- fit$lambda
    if (keep_objects) {
      bases[[k]]$consensus_fit <- fit
      bases[[k]]$zz <- NULL
      bases[[k]]$rank <- NULL
    }
  }
  out <- list(
    scores = cbind(consensus30 = rowMeans(con_prob), diverse30 = rowMeans(div_prob)),
    selected30 = sel, member_diagnostics = do.call(rbind, member_diag),
    member_coefficients = do.call(rbind, member_coefs),
    diverse_member_features = do.call(rbind, diverse_features),
    diverse_member_coefficients = do.call(rbind, diverse_coefficients))
  if (keep_objects) out$member_models <- bases
  if (any(!is.finite(out$scores)) || any(out$scores < 0 | out$scores > 1))
    stop("Nonfinite/out-of-range ensemble scores")
  out
}

# Use family-LOFO cross-fitting inside the outer training population.
# The outer family is never in any inner dataset.
calibrate_inner <- function(x, yy, ff, bb, gg, held_outer) {
  folds <- make_inner_folds(yy, ff, bb, paste0("outer=", held_outer))
  p <- matrix(NA_real_, nrow(x), length(PANEL_NAMES),
              dimnames = list(rownames(x), PANEL_NAMES))
  diag_list <- list()
  split_log <- list()
  for (j in sort(unique(folds))) {
    tr <- which(folds != j); te <- which(folds == j)
    stopifnot(!any(ff[tr] %in% ff[te]), length(unique(ff[tr][yy[tr] == 1L])) >= 3L)
    ctx <- paste0("outer=", held_outer, "/inner=", j)
    message("[INNER] ", ctx, " n_train=", length(tr), " n_test=", length(te),
            " n_pos_train=", sum(yy[tr] == 1L))
    z <- fit_ensemble(x[tr, , drop = FALSE], yy[tr], ff[tr], bb[tr], gg[tr],
                      x[te, , drop = FALSE], bb[te], ctx)
    p[te, ] <- z$scores
    diag_list[[length(diag_list) + 1L]] <- cbind(inner_fold = j,
      z$member_diagnostics)
    split_log[[length(split_log) + 1L]] <- data.frame(
      SampleID = rownames(x)[te], family_id = ff[te], inner_fold = j,
      n_inner_train = length(tr), n_inner_train_positive = sum(yy[tr] == 1L),
      n_inner_train_positive_families = length(unique(ff[tr][yy[tr] == 1L])),
      truth_train_side = yy[te])
  }
  stopifnot(!anyNA(p), all(is.finite(p)))
  # One inner fold per training family; each family is held out exactly once.
  stopifnot(length(unique(folds)) == length(unique(ff)),
            all(vapply(split(folds, ff), function(a) length(unique(a)) == 1L,
                       logical(1L))))
  evals <- do.call(rbind, lapply(PANEL_NAMES, function(m) {
    sf <- full_recall_threshold(p[, m], yy)
    ss <- spec90_threshold(p[, m], yy)
    data.frame(strategy = m, inner_auc = engine$auc(yy, p[, m]),
               full_threshold = sf$tau, full_inner_TP = sf$tp,
               full_inner_FP = sf$fp, full_inner_sensitivity = sf$sens,
               full_inner_specificity = sf$spec, spec90_threshold = ss$tau,
               spec90_inner_TP = ss$tp, spec90_inner_FP = ss$fp,
               spec90_inner_sensitivity = ss$sens,
               spec90_inner_specificity = ss$spec)
  }))
  # Sensitivity is fixed to 100% in inner OOF by construction. Choose method
  # with fewest false positives. AUC breaks ties; then favor consensus30.
  chosen <- evals$strategy[order(evals$full_inner_FP, -evals$inner_auc,
                                 match(evals$strategy, PANEL_NAMES))[1L]]
  list(inner_oof = p, inner_folds = folds, inner_split_log = do.call(rbind, split_log),
       inner_member_audit = do.call(rbind, diag_list), threshold_table = evals,
       train_size_audit = data.frame(
         n_inner_folds = length(unique(folds)),
         n_outer_train = nrow(x),
         min_inner_train = min(vapply(split_log, function(a) a$n_inner_train[1L], numeric(1L))),
         max_inner_train = max(vapply(split_log, function(a) a$n_inner_train[1L], numeric(1L))),
         min_inner_positive_families = min(vapply(split_log, function(a) a$n_inner_train_positive_families[1L], numeric(1L)))),
       adaptive_choice = chosen)
}

# Outer-family fit and score. Held-out y IS NOT SUPPLIED to this function.
run_one_family <- function(held) {
  tr <- which(family != held); te <- which(family == held)
  stopifnot(!any(family[tr] %in% family[te]), length(unique(family[tr][y[tr] == 1L])) >= 4L)
  message("[OUTER] Held family ", held, " n_held=", length(te),
          " n_train_pos_families=", length(unique(family[tr][y[tr] == 1L])))
  inner <- calibrate_inner(X[tr, , drop = FALSE], y[tr], family[tr],
                           batch[tr], group[tr], held)
  final <- fit_ensemble(X[tr, , drop = FALSE], y[tr], family[tr],
                        batch[tr], group[tr], X[te, , drop = FALSE],
                        batch[te], context = paste0("outer=", held, "/fit_all"))
  z <- list()
  for (strategy in c(PANEL_NAMES, "adaptive")) {
    effective <- if (strategy == "adaptive") inner$adaptive_choice else strategy
    ii <- which(inner$threshold_table$strategy == effective)
    s <- final$scores[, effective]
    for (policy in POLICIES) {
      tau <- if (policy == "full_inner_recall")
        inner$threshold_table$full_threshold[ii] else
        inner$threshold_table$spec90_threshold[ii]
      z[[length(z) + 1L]] <- data.frame(
        SampleID = rownames(X)[te], family_id = family[te],
        batch = batch[te], group = group[te], truth = y[te],
        strategy = strategy, fitted_strategy = effective, policy = policy,
        probability = as.numeric(s), threshold = as.numeric(tau),
        predicted_class = as.integer(s >= tau),
        inner_training_positive_samples = sum(y[tr] == 1L),
        inner_false_positives_full = inner$threshold_table$full_inner_FP[ii])
    }
  }
  raw_inner <- data.frame(held_outer_family = held,
    SampleID = rownames(X)[tr], family_id = family[tr], truth = y[tr],
    inner_fold = inner$inner_folds, inner_score_consensus30 = inner$inner_oof[, "consensus30"],
    inner_score_diverse30 = inner$inner_oof[, "diverse30"])
  list(held_family = held, version = VERSION, plan_md5 = PLAN_MD5,
       predictions = do.call(rbind, z),
       inner_predictions = raw_inner,
       thresholds = cbind(held_family = held, inner$threshold_table,
                          chosen_adaptive = inner$adaptive_choice),
       selected_features = cbind(held_family = held, final$selected30),
       member_coefficients = cbind(held_family = held, final$member_coefficients),
       diverse_member_features = cbind(held_family = held, final$diverse_member_features),
       diverse_member_coefficients = cbind(held_family = held, final$diverse_member_coefficients),
       member_audit = cbind(held_family = held, final$member_diagnostics),
       inner_member_audit = cbind(held_family = held, inner$inner_member_audit),
       inner_split_log = cbind(held_family = held, inner$inner_split_log),
       inner_training_size_audit = cbind(held_family = held, inner$train_size_audit))
}

# One cache per family, written atomically. Workers are assigned disjoint sets
# by position; output CSVs are written ONLY by the aggregation process.
CACHE <- file.path(OUTPUT, "outer_cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
family_file <- function(fam) file.path(CACHE,
  paste0("family_", gsub("[^a-zA-Z0-9_-]", "_", fam), ".rds"))
validate_family_cache <- function(z, fam) {
  if (!is.list(z) || !identical(as.character(z$held_family), as.character(fam)) ||
      !identical(z$version, VERSION) || !identical(z$plan_md5, PLAN_MD5))
    return(FALSE)
  p <- z$predictions
  ids <- rownames(X)[family == fam]
  if (!is.data.frame(p) || nrow(p) != length(ids) * length(POLICIES) * 3L ||
      !setequal(unique(as.character(p$SampleID)), ids) ||
      anyDuplicated(p[c("SampleID", "strategy", "policy")])) return(FALSE)
  if (!setequal(unique(p$strategy), c(PANEL_NAMES, "adaptive")) ||
      !setequal(unique(p$policy), POLICIES) ||
      !is.data.frame(z$inner_training_size_audit) ||
      !all(c("n_inner_folds", "n_outer_train", "min_inner_train",
             "max_inner_train", "min_inner_positive_families") %in%
           names(z$inner_training_size_audit)) ||
      anyNA(z$inner_training_size_audit$min_inner_train) ||
      any(z$inner_training_size_audit$min_inner_train < 30L) ||
      any(z$inner_training_size_audit$min_inner_positive_families < 3L)) return(FALSE)
  truth_expected <- y[match(as.character(p$SampleID), rownames(X))]
  if (anyNA(truth_expected) || any(p$truth != truth_expected) ||
      any(!is.finite(p$probability)) || any(!is.finite(p$threshold))) return(FALSE)
  TRUE
}
run_worker <- function() {
  # In a smoke test, prioritize a positive family, the hardest small-sample case.
  smoke_order <- c(sort(unique(family[y == 1L])), setdiff(fams, unique(family[y == 1L])))
  allf <- if (OUTER_LIMIT > 0L) head(smoke_order, OUTER_LIMIT) else fams
  mine <- which(((seq_along(allf) - 1L) %% WORKERS) == (WORKER - 1L))
  cat(sprintf("Ensemble %s: %d outer families; worker %d/%d receives %d\n",
              MODE, length(allf), WORKER, WORKERS, length(mine)))
  for (index in mine) {
    fam <- allf[index]; fp <- family_file(fam)
    if (file.exists(fp)) {
      x <- tryCatch(readRDS(fp), error = function(e) NULL)
      if (!is.null(x) && validate_family_cache(x, fam)) {
        cat("[SKIP] Family ", fam, " validated cache\n", sep = "")
        next
      }
      stop("Invalid family cache: ", fp, " ; do not overwrite without inspection")
    }
    lock <- paste0(fp, ".lock")
    hours <- local({
      acquire_dir_lock(lock, timeout = 90)
      on.exit(unlink(lock, recursive = TRUE), add = TRUE)
      t0 <- Sys.time()
      z <- run_one_family(fam)
      tmp <- tempfile(pattern = ".fold_", tmpdir = CACHE, fileext = ".rds")
      on.exit(unlink(tmp), add = TRUE)
      saveRDS(z, tmp)
      if (!file.rename(tmp, fp)) stop("Unable to atomically write cache for family ", fam)
      as.numeric(difftime(Sys.time(), t0, units = "hours"))
    })
    cat(sprintf("[FINISHED] Family %s in %.2f h\n", fam, hours))
  }
}

summarize_all <- function(z) {
  yy <- z$truth; dd <- z$predicted_class; ss <- z$probability
  pf <- unique(as.character(z$family_id[yy == 1L]))
  family_y <- tapply(yy, z$family_id, max)
  family_score <- tapply(ss, z$family_id, max)
  data.frame(n = nrow(z), n_positive = sum(yy == 1L),
    n_negative = sum(yy == 0L), n_families = length(unique(z$family_id)),
    n_positive_families = length(pf),
    AUC = engine$auc(yy, ss),
    family_auc_descriptive = engine$auc(as.integer(family_y), as.numeric(family_score)),
    sensitivity = mean(dd[yy == 1L] == 1L),
    specificity = mean(dd[yy == 0L] == 0L),
    TP = sum(dd == 1L & yy == 1L), TN = sum(dd == 0L & yy == 0L),
    FP = sum(dd == 1L & yy == 0L), FN = sum(dd == 0L & yy == 1L),
    positive_families_detected = sum(vapply(pf, function(f)
      any(dd[yy == 1L & as.character(z$family_id) == f] == 1L), logical(1L))),
    Brier = mean((ss - yy)^2))
}

# An observational ROC frontier is intentionally named EXPLORATORY because
# choosing a threshold after seeing outer-test labels would be leakage.
exploratory_frontier <- function(z) {
  tvals <- sort(unique(c(0, z$probability, 1)))
  out <- do.call(rbind, lapply(tvals, function(t) {
    d <- z$probability >= t
    data.frame(threshold_posthoc = t, tp = sum(d & z$truth == 1L),
               fn = sum(!d & z$truth == 1L), fp = sum(d & z$truth == 0L),
               tn = sum(!d & z$truth == 0L))
  }))
  out[!duplicated(out[c("tp", "fp")]), , drop = FALSE]
}

aggregate_results <- function() {
  smoke_order <- c(sort(unique(family[y == 1L])), setdiff(fams, unique(family[y == 1L])))
  allf <- if (OUTER_LIMIT > 0L) head(smoke_order, OUTER_LIMIT) else fams
  ok <- vapply(allf, function(f) file.exists(family_file(f)), logical(1L))
  cat("[AGGREGATE] Complete outer families ", sum(ok), "/", length(allf), "\n", sep = "")
  if (OUTER_LIMIT > 0L) {
    if (!all(ok)) stop("Smoke test incomplete")
    cat("Smoke test only: NO final performance estimates were produced.\n")
    return(invisible(NULL))
  }
  if (!all(ok)) stop("Incomplete outer LOFO: ", paste(allf[!ok], collapse = ", "))
  z <- lapply(allf, function(f) {
    v <- readRDS(family_file(f))
    if (!validate_family_cache(v, f)) stop("Invalid/incompatible family cache: ", family_file(f))
    v
  })
  pred <- do.call(rbind, lapply(z, `[[`, "predictions"))
  expected_groups <- expand.grid(SampleID = rownames(X), strategy = c(PANEL_NAMES, "adaptive"),
                                 policy = POLICIES, stringsAsFactors = FALSE)
  given_groups <- pred[, c("SampleID", "strategy", "policy")]
  if (anyDuplicated(given_groups) || !setequal(
      apply(given_groups, 1L, paste, collapse = "|"),
      apply(expected_groups, 1L, paste, collapse = "|")))
    stop("Outer prediction matrix incomplete/duplicated")
  stopifnot(all(vapply(split(pred, interaction(pred$strategy, pred$policy)),
       function(a) nrow(a) == 44L && sum(a$truth == 1L) == 7L, logical(1L))))
  outfiles <- c(predictions = "outer_predictions.csv", thresholds = "outer_threshold_audit.csv",
                selected_features = "outer_consensus30_features.csv",
                member_coefficients = "outer_consensus_member_coefficients.csv",
                diverse_member_features = "outer_diverse_member_features.csv",
                diverse_member_coefficients = "outer_diverse_member_coefficients.csv",
                member_audit = "outer_member_audit.csv",
                inner_member_audit = "outer_inner_member_audit.csv",
                inner_predictions = "outer_inner_predictions.csv",
                inner_split_log = "outer_inner_split_log.csv",
                inner_training_size_audit = "outer_inner_training_size_audit.csv")
  atomic_csv(pred, file.path(OUTPUT, outfiles["predictions"]))
  for (nm in setdiff(names(outfiles), "predictions"))
    atomic_csv(do.call(rbind, lapply(z, `[[`, nm)), file.path(OUTPUT, outfiles[nm]))
  stats <- do.call(rbind, lapply(split(pred, interaction(pred$strategy, pred$policy)), function(d) {
    cbind(strategy = d$strategy[1L], policy = d$policy[1L], summarize_all(d))
  }))
  rownames(stats) <- NULL
  atomic_csv(stats, file.path(OUTPUT, "outer_performance_all_policies.csv"))
  primary <- pred[pred$strategy == "consensus30" & pred$policy == "full_inner_recall", , drop = FALSE]
  primary <- primary[match(rownames(X), primary$SampleID), , drop = FALSE]
  atomic_csv(primary, file.path(OUTPUT, "outer_primary_predictions.csv"))
  atomic_csv(primary[primary$truth == 1L, , drop = FALSE],
             file.path(OUTPUT, "outer_primary_positive7.csv"))
  posfam <- unique(family[y == 1L])
  famdet <- data.frame(family_id = posfam,
                       members = vapply(posfam, function(f)
                           paste(rownames(X)[y == 1L & family == f], collapse = ";"), character(1L)),
                       detected = vapply(posfam, function(f)
                           any(primary$predicted_class[primary$family_id == f & primary$truth == 1L] == 1L), logical(1L)))
  atomic_csv(famdet, file.path(OUTPUT, "outer_positive_family_detection.csv"))
  atomic_csv(exploratory_frontier(primary), file.path(OUTPUT, "POSTHOC_DO_NOT_CLAIM_VALIDATED_threshold_frontier.csv"))
  cat("\n[PRIMARY] consensus30 + 100% inner sample recall threshold\n")
  print(stats[stats$strategy == "consensus30" & stats$policy == "full_inner_recall", , drop = FALSE], row.names = FALSE)
  cat("[NOTICE] No guarantee of seven-of-seven outer sensitivity.\n")

  if (FINAL_MODEL && !file.exists(file.path(OUTPUT, "final_training_model_NOT_VALIDATED.rds"))) {
    # A final deployment candidate may use ALL 44 training labels, AFTER
    # outer evaluation. It is NOT a new independent assessment of performance.
    final_inner <- calibrate_inner(X, y, family, batch, group, "FINAL_ALL44")
    # A scratch prediction row is used ONLY because fit_ensemble has a scoring
    # interface. It is a copy of one training vector with a distinct identifier.
    # All resulting scratch predictions are discarded and play NO role in the
    # selected panel, fitted coefficients, thresholds, or validation metrics.
    dummy <- X[1L, , drop = FALSE]
    rownames(dummy) <- "FINAL_FIT_ONLY_IGNORE_SCRATCH_PREDICTION"
    final_fit <- fit_ensemble(X, y, family, batch, group, dummy, batch[1L],
                              "FINAL_ALL44/fit", keep_objects = TRUE)
    # The deployable primary strategy is ONE consensus set of 30 proteins and
    # K separately fitted elastic-net models restricted to those 30 features.
    tf <- final_inner$threshold_table
    tf$chosen_adaptive <- final_inner$adaptive_choice
    atomic_csv(tf, file.path(OUTPUT, "final_training_inner_thresholds_NOT_VALIDATION.csv"))
    atomic_csv(final_fit$selected30,
               file.path(OUTPUT, "final_consensus30_proteins_NOT_INDEPENDENTLY_VALIDATED.csv"))
    atomic_csv(final_fit$member_coefficients,
               file.path(OUTPUT, "final_consensus_member_coefficients.csv"))
    atomic_csv(final_fit$diverse_member_features,
               file.path(OUTPUT, "final_diverse_member_features.csv"))
    deployment <- list(version = VERSION, engine_version = engine$SCRIPT_VERSION,
                       feature_set = final_fit$selected30$FeatureID,
                       member_seeds = SEEDS, fit_mode = MODE,
                       models = lapply(final_fit$member_models, function(v)
                         list(preprocessor = v$pre, consensus_fit = v$consensus_fit,
                              diverse_fit = v$diverse_fit, diverse_ids = v$diverse_ids)),
                       inner_thresholds = tf,
                       optimization_status = "TRAINING-DERIVED NOT EXTERNALLY VALIDATED")
    atomic_rds(deployment, file.path(OUTPUT, "final_training_model_NOT_VALIDATED.rds"))
    cat("[FINAL MODEL] Consensus30 feature list and fitted objects saved. ",
        "NO performance was computed from the final training fit.\n")
  }
  invisible(stats)
}

if (AGGREGATE_ONLY) aggregate_results() else {
  run_worker()
  if (WORKERS == 1L) aggregate_results() else
    cat("[WORKER DONE] Run again with CITRIN_ENSEMBLE_AGGREGATE_ONLY=true after all workers exit.\n")
}
