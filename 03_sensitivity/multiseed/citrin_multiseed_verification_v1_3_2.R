#!/usr/bin/env Rscript
# Citrin multi-seed verification harness v1.3.2 (2026-10-31)
# v1.3.2: (e) Rscript passes its --file= path with every space encoded as "~+~";
#         script_dir() now decodes it, so the harness also works from folders
#         whose path contains spaces (v1.3.1 crashed at start in "Mac (4)").
#         (f) options(warn = 1): warnings are written to the log as they occur.
# Model engines are unmodified v2.2; this file manages parallel workers,
# immutable run plans, safe resumption, output validation and aggregation.
# In parallel mode run one process per VERIFY_WORKER. After all exit, call
# VERIFY_AGGREGATE_ONLY=1 to produce the definitive summary.
#
# v1.3 changes relative to ChatGPT's v1.2 (plan, validation, aggregation unchanged):
#  (a) DYNAMIC work claiming: every worker walks the whole plan (longest runs
#      first) and claims the next unclaimed run with a PID-stamped lock, so
#      load balance no longer depends on the cost estimates. (v1.2 assumed a
#      44-sample fixed30 run costs 0.3x a select run; measured: ~0.9x, which
#      would have made two workers finish ~6 h after the others.)
#  (b) Stale locks (left by a killed worker) are detected by PID liveness and
#      removed automatically; a live worker's lock is never touched.
#  (c) A run whose RUN_COMPLETE exists but fails validation is reported as
#      "invalid" in run_status.csv instead of aborting the aggregation.
#  (d) The launcher detaches itself and runs under caffeinate so the Mac does
#      not sleep and closing the terminal does not kill the workers.

HARNESS_VERSION <- "multiseed_parallel_1.3.2"
options(stringsAsFactors = FALSE, warn = 1)

script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
  f <- gsub("~+~", " ", f, fixed = TRUE)   # Rscript encodes spaces this way
  if (length(f)) return(dirname(normalizePath(f[1], mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
ROOT <- script_dir()
getenv <- function(n, default) {
  v <- Sys.getenv(n, "")
  if (nzchar(v)) v else default
}
split_csv <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
SEEDS <- as.integer(split_csv(getenv("VERIFY_SEEDS",
                   "2025,2026,2027,2028,2029,2030,2031,2032,2033,2034")))
MODES <- split_csv(getenv("VERIFY_MODES", "fixed30,select"))
SUBS <- split_csv(getenv("VERIFY_SUBSAMPLE", "sample"))
SCOPE <- split_csv(getenv("VERIFY_SCOPE", "44,Run3"))
N_BOOT <- as.integer(getenv("VERIFY_N_BOOT", "0"))
LIMIT <- as.integer(getenv("VERIFY_OUTER_LIMIT", "0"))
WORKERS <- as.integer(getenv("VERIFY_WORKERS", "1"))
WORKER <- as.integer(getenv("VERIFY_WORKER", "1"))
AGGREGATE_ONLY <- identical(getenv("VERIFY_AGGREGATE_ONLY", "0"), "1")
REQUIRE_COMPLETE <- identical(getenv("VERIFY_REQUIRE_COMPLETE", "0"), "1")
DIR <- normalizePath(getenv("VERIFY_OUTPUT", file.path(getwd(), "result", "citrin_multiseed_v1_2")),
                     winslash = "/", mustWork = FALSE)
DATA <- normalizePath(getenv("CITRIN_DATA_DIR", file.path(getwd(), "..", "data")),
                      winslash = "/", mustWork = TRUE)
FAMILY <- normalizePath(getenv("CITRIN_FAMILY_FILE", file.path(DATA, "family_map.csv")),
                        winslash = "/", mustWork = TRUE)
SCRIPT44 <- normalizePath(getenv("VERIFY_SCRIPT44", file.path(ROOT, "citrin_44_nested_fixed_panel_v2_2.R")),
                          winslash = "/", mustWork = TRUE)
SCRIPTR3 <- normalizePath(getenv("VERIFY_SCRIPTR3", file.path(ROOT, "citrin_Run3_nested_fixed_panel_v2_2.R")),
                          winslash = "/", mustWork = TRUE)

stopifnot(length(SEEDS) > 0L, !anyNA(SEEDS), !anyDuplicated(SEEDS), all(SEEDS >= 0L),
          length(MODES) > 0L, !anyDuplicated(MODES), all(MODES %in% c("fixed30", "select", "fixed15")),
          length(SUBS) > 0L, !anyDuplicated(SUBS), all(SUBS %in% c("sample", "family")),
          length(SCOPE) > 0L, !anyDuplicated(SCOPE), all(SCOPE %in% c("44", "Run3")),
          is.finite(N_BOOT), N_BOOT >= 0L, is.finite(LIMIT), LIMIT >= 0L,
          is.finite(WORKERS), WORKERS >= 1L, is.finite(WORKER), WORKER >= 1L, WORKER <= WORKERS,
          file.exists(file.path(DATA, "raw_data.xlsx")), file.exists(file.path(DATA, "sample_info.csv")))

# Base-R directory locks avoid simultaneous plan initialization, even if all
# workers launch at once. A killed process can leave a stale lock; investigate
# running processes before manually removing a stale .plan_lock/.run_lock.
acquire_lock <- function(path, wait_seconds = 90) {
  started <- Sys.time()
  repeat {
    if (dir.create(path, showWarnings = FALSE)) return(invisible(TRUE))
    if (as.numeric(difftime(Sys.time(), started, units = "secs")) > wait_seconds)
      stop("Timed out acquiring lock: ", path,
           ". If stale, verify no active worker before removing it.")
    Sys.sleep(0.2)
  }
}
# Non-blocking PID-stamped lock for dynamic run claiming. Returns TRUE when
# this process now owns the lock. A lock whose recorded PID is no longer
# alive is treated as stale and removed; a live owner's lock is respected.
pid_alive <- function(pid) {
  if (!is.finite(pid) || pid <= 0) return(FALSE)
  if (.Platform$OS.type == "windows") return(TRUE)   # no liveness check: never steal
  identical(suppressWarnings(system2("kill", c("-0", as.character(pid)),
                                     stdout = FALSE, stderr = FALSE)), 0L)
}
try_claim <- function(path) {
  for (attempt in 1:2) {
    if (dir.create(path, showWarnings = FALSE)) {
      writeLines(c(as.character(Sys.getpid()), format(Sys.time(), tz = "UTC", usetz = TRUE),
                   sprintf("worker=%d", WORKER)), file.path(path, "owner"))
      return(TRUE)
    }
    owner <- suppressWarnings(as.integer(tryCatch(readLines(file.path(path, "owner"), n = 1L),
                                                  error = function(e) NA_character_)))
    if (is.na(owner)) { Sys.sleep(0.5); next }          # lock being written right now
    if (pid_alive(owner)) return(FALSE)                  # another live worker owns it
    cat("[STALE LOCK] removing ", path, " (pid ", owner, " not running)\n", sep = "")
    unlink(path, recursive = TRUE)
  }
  FALSE
}
atomic_lines <- function(lines, path) {
  tmp <- tempfile(pattern = ".tmp_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  writeLines(lines, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}
atomic_csv <- function(df, path) {
  tmp <- tempfile(pattern = ".tmp_", tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  write.csv(df, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Atomic CSV rename failed: ", path)
}

# All behavior-relevant environment overrides are fixed explicitly here.
# The engine version and the input-file checksums are part of the plan.
dir.create(DIR, recursive = TRUE, showWarnings = FALSE)
plan <- expand.grid(scope = SCOPE, panel = MODES, subsample = SUBS, seed = SEEDS,
                    stringsAsFactors = FALSE)
plan$run_id <- sprintf("%s_%s_%s_seed%d", plan$scope, plan$panel, plan$subsample, plan$seed)
stopifnot(!anyDuplicated(plan$run_id))
plan$script <- ifelse(plan$scope == "44", SCRIPT44, SCRIPTR3)
plan$data_dir <- DATA
plan$family_file <- FAMILY
plan$outer_limit <- LIMIT
plan$n_boot <- N_BOOT
plan$inner_mode <- "lofo"
plan$stability_B <- 30L
plan$run2_label <- "Run2"
plan$run3_label <- "Run3"
plan$expect_44 <- "true"
plan$expect_35 <- "true"
plan$run3_final_model <- "false"
plan$harness_version <- HARNESS_VERSION
plan$script_md5 <- vapply(plan$script, function(f) unname(tools::md5sum(f)), character(1L))
plan$raw_md5 <- unname(tools::md5sum(file.path(DATA, "raw_data.xlsx")))
plan$sample_info_md5 <- unname(tools::md5sum(file.path(DATA, "sample_info.csv")))
plan$family_md5 <- unname(tools::md5sum(FAMILY))
rownames(plan) <- NULL

ensure_plan <- function() {
  lock <- file.path(DIR, ".plan_lock")
  acquire_lock(lock)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  path <- file.path(DIR, "preregistered_plan.csv")
  plan_chr <- as.data.frame(lapply(plan, as.character), stringsAsFactors = FALSE)
  if (file.exists(path)) {
    old <- read.csv(path, check.names = FALSE, colClasses = "character")
    if (!isTRUE(all.equal(old, plan_chr, check.attributes = FALSE)))
      stop("Existing preregistered_plan.csv differs. Use a NEW VERIFY_OUTPUT folder; ",
           "do not mix datasets, algorithms, or verification protocols.")
  } else {
    atomic_csv(plan, path)
  }
}
ensure_plan()
# Balance workers by estimated cost rather than (plan_row %% WORKERS).
# With 44 and Run3 alternating in expand.grid, modulo-6 gives three workers
# ONLY expensive 44 runs and three ONLY short Run3 runs. Greedy assignment
# spreads expensive tasks across workers. Costs affect scheduling only.
# Measured on the user's Mac (2026-10-31): 44-sample select 5.2 h, 44-sample
# fixed30 ~4.8 h (stability selection dominates, so the panel grid matters
# little), Run3 select 0.6 h, Run3 fixed30 ~0.5 h. Used ONLY to order the
# queue (longest first); workers claim runs dynamically, see execute_run().
cost <- ifelse(plan$scope == "44",
               ifelse(plan$panel == "select", 5.2, 4.8),
               ifelse(plan$panel == "select", 0.6, 0.5))
long_first <- order(-cost, plan$run_id)
cat(sprintf("PLAN: %d runs; worker %d/%d claims runs dynamically; scope=%s; output=%s\n",
            nrow(plan), WORKER, WORKERS, paste(SCOPE, collapse = ","), DIR))
cat(sprintf("Total estimated CPU hours %.1f; ideal wall-clock with %d workers ~%.1f h\n",
            sum(cost), WORKERS, sum(cost) / WORKERS))

ENGINE_ENV <- c("CITRIN_SEED", "CITRIN_PANEL_MODE", "CITRIN_SUBSAMPLE_UNIT",
                "CITRIN_OUT_DIR", "CITRIN_STAGE", "CITRIN_N_BOOT", "CITRIN_OUTER_LIMIT",
                "CITRIN_DATA_DIR", "CITRIN_FAMILY_FILE", "CITRIN_INNER_MODE",
                "CITRIN_RUN3_PREDICTIONS", "CITRIN_B", "CITRIN_EXPECT_44",
                "CITRIN_EXPECT_35", "CITRIN_RUN2", "CITRIN_RUN3", "CITRIN_FINAL_MODEL")

run_one <- function(r, work) {
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  oldenv <- Sys.getenv(ENGINE_ENV, unset = NA_character_)
  on.exit({
    for (j in seq_along(oldenv)) {
      if (is.na(oldenv[j])) Sys.unsetenv(names(oldenv)[j])
      else do.call(Sys.setenv, setNames(list(oldenv[j]), names(oldenv)[j]))
    }
  }, add = TRUE)
  Sys.setenv(CITRIN_SEED = as.character(r$seed),
             CITRIN_PANEL_MODE = r$panel, CITRIN_SUBSAMPLE_UNIT = r$subsample,
             CITRIN_DATA_DIR = DATA, CITRIN_FAMILY_FILE = FAMILY,
             CITRIN_OUT_DIR = file.path(work, "result", "two_run_nested_cv"),
             CITRIN_STAGE = "both", CITRIN_N_BOOT = as.character(N_BOOT),
             CITRIN_OUTER_LIMIT = as.character(LIMIT), CITRIN_INNER_MODE = "lofo",
             CITRIN_B = "30", CITRIN_EXPECT_44 = "true", CITRIN_EXPECT_35 = "true",
             CITRIN_RUN2 = "Run2", CITRIN_RUN3 = "Run3", CITRIN_FINAL_MODEL = "false")
  Sys.unsetenv("CITRIN_RUN3_PREDICTIONS")
  setwd(work)
  con <- file(file.path(work, "console.log"), open = "at")
  sink(con, split = FALSE)
  sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)
  source(r$script, local = new.env(parent = globalenv()), chdir = FALSE)
  invisible(TRUE)
}

# Exactly one Run3 OOF file, or exactly one each of Stage1/Stage2 for the 44
# sample run. Smoke tests may produce only partial predictions by design.
prediction_files <- function(r, work) {
  resdir <- file.path(work, "result")
  if (!dir.exists(resdir)) stop("No result directory: ", resdir)
  fps <- list.files(resdir, pattern = "^outer_predictions[.]csv$",
                    recursive = TRUE, full.names = TRUE)
  fps_unix <- gsub("\\\\", "/", fps)
  if (identical(as.character(r$scope), "44")) {
    p1 <- fps[grepl("/stage1_no_combat/outer_predictions[.]csv$", fps_unix)]
    p2 <- fps[grepl("/stage2_frozen_mean_combat/outer_predictions[.]csv$", fps_unix)]
    if (length(p1) != 1L || length(p2) != 1L)
      stop("Expected exactly one Stage1 and Stage2 prediction CSV for ", r$run_id,
           "; found ", length(p1), " and ", length(p2))
    return(c(no_combat = p1, combat_mean = p2))
  }
  p3 <- fps[grepl("/run3_original_nested_cv/[^/]+/outer_predictions[.]csv$", fps_unix)]
  if (length(p3) != 1L)
    stop("Expected exactly one Run3 prediction CSV for ", r$run_id,
         "; found ", length(p3))
  c(run3_no_combat = p3)
}

validate_prediction <- function(fp, expected_n) {
  p <- read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("SampleID", "family_id", "batch", "group", "truth",
                "probability", "threshold", "predicted_class")
  if (!all(required %in% names(p))) stop("Missing prediction columns in ", fp)
  if (nrow(p) != expected_n) stop("Unexpected sample count in ", fp,
                                  ": ", nrow(p), " instead of ", expected_n)
  if (anyNA(p[, required]) || anyDuplicated(p$SampleID) ||
      !all(p$truth %in% c(0, 1)) || !all(p$predicted_class %in% c(0, 1)) ||
      !all(is.finite(p$probability)) || !all(is.finite(p$threshold)) ||
      any(p$probability < 0 | p$probability > 1)) stop("Invalid values in ", fp)
  if (sum(p$truth == 1) != 7L ||
      length(unique(as.character(p$family_id[p$truth == 1]))) != 5L)
    stop("Expected seven positive specimens from five families in ", fp)
  if (any(as.integer(p$predicted_class) != as.integer(p$probability >= p$threshold)))
    stop("predicted_class differs from recorded threshold in ", fp)
  p
}
validate_run <- function(r, work) {
  fps <- prediction_files(r, work)
  expected <- if (r$scope == "44") 44L else 35L
  pp <- lapply(fps, validate_prediction, expected_n = expected)
  # A valid prediction CSV alone is not sufficient: the manuscript audit
  # requires coefficients, selected features, lambda and stability diagnostics.
  required_audit <- c("outer_coefficients.csv", "outer_selected_features.csv",
                      "outer_inner_lambda_audit.csv", "outer_stability_diagnostics.csv",
                      "outer_95CI.csv")
  for (fp in fps) {
    needed <- file.path(dirname(fp), required_audit)
    if (any(!file.exists(needed)))
      stop("Missing audit output(s) for ", r$run_id, ": ",
           paste(basename(needed[!file.exists(needed)]), collapse = ", "))
  }
  if (length(pp) > 1L) {
    a <- pp[[1L]]; b <- pp[[2L]]
    ia <- order(a$SampleID); ib <- order(b$SampleID)
    for (col in c("SampleID", "family_id", "group", "batch", "truth")) {
      if (!identical(as.character(a[[col]][ia]), as.character(b[[col]][ib])))
        stop("Stage1 and Stage2 disagree on ", col, " for ", r$run_id)
    }
  }
  invisible(fps)
}

read_hours <- function(path) {
  x <- tryCatch(readLines(path), error = function(e) character())
  h <- suppressWarnings(as.numeric(sub("^hours=", "", grep("^hours=", x, value = TRUE))))
  if (length(h)) h[1L] else NA_real_
}
wauc <- function(y, s) {
  if (sum(y == 1) == 0L || sum(y == 0) == 0L || any(!is.finite(s))) return(NA_real_)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(rank(s, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

collect <- function(r) {
  work <- file.path(DIR, "runs", r$run_id)
  done <- file.path(work, "RUN_COMPLETE")
  if (!file.exists(done)) return(NULL)
  fps <- tryCatch(validate_run(r, work), error = function(e) e)
  if (inherits(fps, "error")) {
    cat("[INVALID] ", r$run_id, ": ", conditionMessage(fps), "\n", sep = "")
    return(NULL)
  }
  out <- lapply(names(fps), function(m) {
    fp <- fps[[m]]
    p <- validate_prediction(fp, if (r$scope == "44") 44L else 35L)
    y <- as.integer(p$truth)
    d <- as.integer(p$predicted_class)
    s <- as.numeric(p$probability)
    posfam <- unique(as.character(p$family_id[y == 1L]))
    detected <- vapply(posfam, function(f)
      any(d[y == 1L & p$family_id == f] == 1L), logical(1L))
    # Descriptive family AUC: maximum specimen score for each family.
    fam_y <- tapply(y, p$family_id, max)
    fam_s <- tapply(s, p$family_id, max)
    dd <- dirname(fp)
    st <- tryCatch(read.csv(file.path(dd, "outer_stability_diagnostics.csv")), error = function(e) NULL)
    au <- tryCatch(read.csv(file.path(dd, "outer_inner_lambda_audit.csv")), error = function(e) NULL)
    inner_fb <- NA_integer_
    if (!is.null(au) && all(c("held_family", "inner_fold", "ranking_fallbacks") %in% names(au))) {
      u <- au[!duplicated(au[, c("held_family", "inner_fold")]), , drop = FALSE]
      inner_fb <- sum(u$ranking_fallbacks)
    }
    data.frame(run_id = r$run_id, scope = r$scope, panel = r$panel,
               subsample = r$subsample, seed = r$seed, model = m,
               n = nrow(p), n_families = length(unique(p$family_id)),
               auc = wauc(y, s), family_auc = wauc(as.integer(fam_y), as.numeric(fam_s)),
               sensitivity = mean(d[y == 1L] == 1L),
               specificity = mean(d[y == 0L] == 0L),
               brier = mean((s - y)^2), tp = sum(d == 1L & y == 1L),
               fp = sum(d == 1L & y == 0L),
               positive_families = length(posfam),
               detected_positive_families = sum(detected),
               family_detection_rate = mean(detected),
               median_selected_m = if ("selected_m" %in% names(p)) median(p$selected_m) else NA_real_,
               n_outer_lambda_fallbacks = if ("lambda_method" %in% names(p)) sum(grepl("^fallback", p$lambda_method)) else NA_integer_,
               n_outer_stability_fallbacks = if (!is.null(st) && "fallback_lambda_refits" %in% names(st)) sum(st$fallback_lambda_refits) else NA_integer_,
               n_inner_stability_fallbacks = inner_fb,
               hours = read_hours(done), prediction_file = normalizePath(fp, mustWork = TRUE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

execute_run <- function(r) {
  work <- file.path(DIR, "runs", r$run_id)
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  done <- file.path(work, "RUN_COMPLETE")
  if (file.exists(done)) {
    if (LIMIT == 0L) validate_run(r, work)
    cat("[SKIP] Already validated: ", r$run_id, "\n", sep = "")
    return(TRUE)
  }
  lock <- file.path(work, ".run_lock")
  if (!try_claim(lock)) return(NA)          # another live worker is running it
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (file.exists(done)) {
    if (LIMIT == 0L) validate_run(r, work)
    return(TRUE)
  }
  if (file.exists(file.path(work, "RUN_FAILED.txt"))) unlink(file.path(work, "RUN_FAILED.txt"))
  t0 <- Sys.time()
  err <- tryCatch({
    run_one(r, work)
    if (LIMIT == 0L) validate_run(r, work)
    NULL
  }, error = function(e) conditionMessage(e))
  hours <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
  if (is.null(err)) {
    atomic_lines(c("OK", paste0("start=", format(t0, tz = "UTC", usetz = TRUE)),
                   paste0("end=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
                   sprintf("hours=%.5f", hours)), done)
    cat(sprintf("[OK] %s (%.2f h)\n", r$run_id, hours))
    return(TRUE)
  }
  atomic_lines(c(err, sprintf("hours=%.5f", hours)), file.path(work, "RUN_FAILED.txt"))
  cat("[FAILED] ", r$run_id, ": ", err, "\n", sep = "")
  FALSE
}

summarize <- function(metrics) {
  q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
  key <- interaction(metrics$scope, metrics$panel, metrics$subsample, metrics$model, drop = TRUE)
  out <- do.call(rbind, lapply(split(metrics, key), function(z) data.frame(
    scope = z$scope[1L], panel = z$panel[1L], subsample = z$subsample[1L], model = z$model[1L],
    n_seeds = length(unique(z$seed)),
    auc_mean = mean(z$auc), auc_median = median(z$auc), auc_min = min(z$auc),
    auc_q25 = q(z$auc, 0.25), auc_q75 = q(z$auc, 0.75), auc_max = max(z$auc),
    family_auc_median = median(z$family_auc),
    family_auc_min = min(z$family_auc), family_auc_max = max(z$family_auc),
    sensitivity_median = median(z$sensitivity), sensitivity_min = min(z$sensitivity),
    sensitivity_max = max(z$sensitivity),
    specificity_median = median(z$specificity), specificity_min = min(z$specificity),
    family_detection_median = median(z$family_detection_rate),
    family_detection_min = min(z$family_detection_rate),
    family_detection_max = max(z$family_detection_rate),
    inner_stability_fallbacks_median = median(z$n_inner_stability_fallbacks, na.rm = TRUE),
    hours_median = median(z$hours, na.rm = TRUE))))
  rownames(out) <- NULL
  out
}

aggregate_results <- function() {
  # Called once by the launcher after all workers have exited. This single
  # writer prevents stale per_seed_metrics.csv / summary races.
  status <- data.frame(run_id = plan$run_id,
                       complete = file.exists(file.path(DIR, "runs", plan$run_id, "RUN_COMPLETE")),
                       failed = file.exists(file.path(DIR, "runs", plan$run_id, "RUN_FAILED.txt")),
                       running = file.exists(file.path(DIR, "runs", plan$run_id, ".run_lock")),
                       stringsAsFactors = FALSE)
  status$invalid <- FALSE
  if (LIMIT == 0L) {
    for (i in which(status$complete)) {
      v <- tryCatch(validate_run(plan[i, , drop = FALSE], file.path(DIR, "runs", plan$run_id[i])),
                    error = function(e) e)
      if (inherits(v, "error")) { status$invalid[i] <- TRUE; status$complete[i] <- FALSE }
    }
  }
  atomic_csv(status, file.path(DIR, "run_status.csv"))
  n_ok <- sum(status$complete)
  n_fail <- sum((status$failed | status$invalid) & !status$complete)
  cat(sprintf("STATUS: %d/%d complete, %d failed/invalid, %d running, %d pending\n",
              n_ok, nrow(plan), n_fail, sum(status$running),
              nrow(plan) - n_ok - n_fail))
  if (LIMIT > 0L) {
    cat("SMOKE TEST: No performance summary. Use a NEW VERIFY_OUTPUT with outer_limit=0.\n")
    if (n_fail > 0L) stop("Smoke-test run(s) failed; inspect RUN_FAILED.txt and console.log")
    return(invisible(status))
  }
  if (n_ok == 0L) stop("No complete runs to summarize.")
  rows <- Filter(Negate(is.null), lapply(seq_len(nrow(plan)), function(i) collect(plan[i, ])))
  metrics <- do.call(rbind, rows)
  atomic_csv(metrics, file.path(DIR, "per_seed_metrics.csv"))
  summ <- summarize(metrics)
  completed <- n_ok == nrow(plan) && n_fail == 0L
  if (completed) {
    atomic_csv(summ, file.path(DIR, "multiseed_summary.csv"))
    unlink(file.path(DIR, "multiseed_summary_PARTIAL.csv"))
    cat("COMPLETE: all planned runs validated.\n")
  } else {
    atomic_csv(summ, file.path(DIR, "multiseed_summary_PARTIAL.csv"))
    unlink(file.path(DIR, "multiseed_summary.csv"))
    cat("PARTIAL: missing or failed runs are NOT counted as complete.\n")
  }
  print(summ, row.names = FALSE, digits = 4)
  if (n_fail > 0L) stop(n_fail, " run(s) failed; inspect run_status.csv and per-run logs.")
  if (REQUIRE_COMPLETE && !completed) stop("Final aggregation incomplete: ", n_ok, "/", nrow(plan), " run(s) complete.")
  invisible(summ)
}

if (AGGREGATE_ONLY) {
  aggregate_results()
} else {
  failures <- character(); n_done <- 0L
  for (i in long_first) {
    r <- plan[i, , drop = FALSE]
    ok <- tryCatch(execute_run(r), error = function(e) {
      cat("[WORKER ERROR] ", r$run_id, ": ", conditionMessage(e), "\n", sep = "")
      FALSE
    })
    if (is.na(ok)) next                     # claimed by another live worker
    if (isTRUE(ok)) n_done <- n_done + 1L else failures <- c(failures, r$run_id)
  }
  cat(sprintf("Worker %d: %d run(s) complete or verified, %d failed.\n",
              WORKER, n_done, length(failures)))
  if (WORKERS == 1L) aggregate_results()
  if (length(failures)) stop(length(failures), " worker run(s) could not complete")
  if (WORKERS > 1L) cat("Worker complete. Run VERIFY_AGGREGATE_ONLY=1 after all workers exit.\n")
}
