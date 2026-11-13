#!/usr/bin/env Rscript
# =============================================================================
# Citrin DBS proteomics: top-four protein ablation sensitivity analysis
# Based on the existing 10-seed consensus30 nested family-LOFO v1.1 pipeline.
#
# IMPORTANT: This script is a DRIVER, not a different classification method.
# It copies the ORIGINAL v1.1 harness, inserts a feature-exclusion step just
# after engine$read_44(), and otherwise preserves the original algorithm.
# The four proteins are removed BEFORE imputation, variance ranking, stability
# selection, elastic-net fitting, inner LOFO, and threshold selection.
#
# ALL NEW OUTPUTS ARE WRITTEN UNDER THE PROJECT-LOCAL result/ DIRECTORY:
#   <R project>/result/citrin_top4_ablation_v1_1/<scenario>/
# The existing Claude/ directory is NEVER used as an output destination.
# Its previous baseline results may be READ, not modified, during comparison.
#
# Terminal (run from the R project directory):
#   Rscript citrin_top4_ablation_familyLOFO.R prepare
#   Rscript citrin_top4_ablation_familyLOFO.R smoke all4
#   Rscript citrin_top4_ablation_familyLOFO.R run all4
#   Rscript citrin_top4_ablation_familyLOFO.R run no_PRL
#   Rscript citrin_top4_ablation_familyLOFO.R run no_MTTP
#   Rscript citrin_top4_ablation_familyLOFO.R run no_ARHGEF10
#   Rscript citrin_top4_ablation_familyLOFO.R run no_SEC31A
#   Rscript citrin_top4_ablation_familyLOFO.R compare
#
# RStudio Console (after source() or pasting definitions):
#   Sys.setenv(CITRIN_ABLATION_ACTION = "run",
#              CITRIN_ABLATION_SCENARIO = "all4")
#   main()
#
# Optional read-only location of the completed v1.1 baseline:
#   CITRIN_ABLATION_BASELINE_OUTPUT=/absolute/path/to/original/output
#
# These are POST-HOC, same-cohort sensitivity analyses, NOT external validation.
# Do not select the best-performing ablation as a new pre-specified model.
# =============================================================================

options(stringsAsFactors = FALSE)

SCENARIOS <- list(
  all4 = c(PRL = "P01236", MTTP = "P55157",
           ARHGEF10 = "O15013", SEC31A = "O94979"),
  no_PRL = c(PRL = "P01236"),
  no_MTTP = c(MTTP = "P55157"),
  no_ARHGEF10 = c(ARHGEF10 = "O15013"),
  no_SEC31A = c(SEC31A = "O94979"),
  no_ARHGEF10_SEC31A = c(ARHGEF10 = "O15013", SEC31A = "O94979"),
  baseline = character(0)
)
MAIN_SCENARIOS <- c("all4", "no_PRL", "no_MTTP", "no_ARHGEF10", "no_SEC31A")
MEMBER_SEEDS <- paste(2025:2034, collapse = ",")
ORIGINAL_HARNESS <- "citrin_nested_ensemble_familyLOFO_v1_1.R"
ORIGINAL_ENGINE <- "citrin_44_nested_fixed_panel_v2_2.R"
ORIGINAL_LAUNCHER <- "run_citrin_nested_ensemble_6workers_caffeinate.sh"
RESULT_FOLDER <- "citrin_top4_ablation_v1_1"

stop_here <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
is_file <- function(p) length(p) == 1L && nzchar(p) && file.exists(p) && !dir.exists(p)
file_md5 <- function(p) unname(tools::md5sum(p))
canonical <- function(p, must_exist = TRUE) {
  normalizePath(p, winslash = "/", mustWork = must_exist)
}

# Source files can be in the project root, bin/, or the original package folder.
# No manual copying to the root directory is necessary.
source_candidates <- function(project, filename) {
  unique(file.path(project,
                   c("", "bin", "citrin_nested_ensemble_v1_1_pkg"),
                   filename))
}
find_source_file <- function(project, filename, override = "") {
  paths <- c(if (nzchar(override)) override, source_candidates(project, filename))
  paths <- paths[vapply(paths, is_file, logical(1L))]
  if (!length(paths)) {
    stop_here("Missing source file: ", filename,
              "\nChecked project root, bin/, and citrin_nested_ensemble_v1_1_pkg/.")
  }
  canonical(paths[[1L]])
}

find_project <- function() {
  explicit <- Sys.getenv("CITRIN_ABLATION_PROJECT_DIR", "")
  start <- canonical(if (nzchar(explicit)) explicit else getwd())
  candidates <- if (nzchar(explicit)) start else {
    x <- start
    out <- character(0)
    for (i in seq_len(7L)) {
      out <- c(out, x)
      if (dirname(x) == x) break
      x <- dirname(x)
    }
    out
  }
  for (p in candidates) {
    project_file <- file.path(p, "20261029_rebuttal.Rproj")
    has_source <- any(vapply(source_candidates(p, ORIGINAL_HARNESS),
                             is_file, logical(1L)))
    if (is_file(project_file) || has_source) return(canonical(p))
  }
  stop_here("R project not found from: ", start,
            "\nOpen 20261029_rebuttal.Rproj or set CITRIN_ABLATION_PROJECT_DIR.")
}

# Fixed, project-local destination. Deliberately IGNORE an old environment
# variable that could otherwise silently redirect output to Claude/.
find_output_root <- function(project) {
  legacy <- Sys.getenv("CITRIN_ABLATION_ROOT", "")
  root <- canonical(file.path(project, "result", RESULT_FOLDER), FALSE)
  if (nzchar(legacy) && canonical(legacy, FALSE) != root) {
    message("[NOTE] Ignoring old CITRIN_ABLATION_ROOT: ", legacy,
            "\n       New outputs are restricted to: ", root)
  }
  root
}

find_data_dir <- function(project) {
  explicit <- Sys.getenv("CITRIN_DATA_DIR", "")
  options <- c(if (nzchar(explicit)) explicit,
               file.path(dirname(project), "data"), file.path(project, "data"))
  options <- unique(options)
  found <- options[vapply(options, function(d) {
    is_file(file.path(d, "raw_data.xlsx")) &&
      is_file(file.path(d, "sample_info.csv"))
  }, logical(1L))]
  if (!length(found)) {
    stop_here("Cannot find raw_data.xlsx and sample_info.csv. ",
              "Set CITRIN_DATA_DIR to the directory containing both files.")
  }
  canonical(found[[1L]])
}

find_family_map <- function(data_dir) {
  override <- Sys.getenv("CITRIN_FAMILY_FILE", "")
  path <- if (nzchar(override)) override else file.path(data_dir, "family_map.csv")
  if (!is_file(path)) stop_here("Family map not found: ", path)
  canonical(path)
}

# Insert the exclusion at a unique, audited source-code anchor. This is the
# only change to the v1.1 classifier; all other source lines are preserved.
patch_original_harness <- function(lines, ids, scenario) {
  if (!any(grepl('VERSION <- "citrin_nested_ensemble_familyLOFO_1.1"',
                 lines, fixed = TRUE))) {
    stop_here("Unexpected v1.1 harness version: check the source file.")
  }
  at <- which(trimws(lines) == "dat <- engine$read_44()")
  if (length(at) != 1L || at == length(lines) ||
      trimws(lines[at + 1L]) != "X <- dat$X; y <- as.integer(dat$y)") {
    stop_here("Source code differs from audited v1.1 near read_44(). ",
              "Refusing to patch an unknown algorithm.")
  }
  if (!scenario %in% names(SCENARIOS)) stop_here("Unknown scenario: ", scenario)
  if (any(!grepl("^[A-Z0-9]+$", ids))) stop_here("Unsafe accession name")
  values <- if (length(ids)) {
    paste0('c("', paste(ids, collapse = '\", \"'), '")')
  } else "character(0)"
  insertion <- c(
    "# ---- START AUDITED ABLATION INSERTION ----",
    paste0('ABLATION_SCENARIO <- "', scenario, '"'),
    paste0("ABLATION_EXCLUDED_ACCESSIONS <- ", values),
    "if (length(ABLATION_EXCLUDED_ACCESSIONS)) {",
    "  missing <- setdiff(ABLATION_EXCLUDED_ACCESSIONS, colnames(dat$X))",
    "  if (length(missing)) stop('Excluded accession absent from raw matrix: ',",
    "                           paste(missing, collapse = ', '))",
    "  dat$X <- dat$X[, !(colnames(dat$X) %in% ABLATION_EXCLUDED_ACCESSIONS),",
    "                 drop = FALSE]",
    "}",
    "if (any(ABLATION_EXCLUDED_ACCESSIONS %in% colnames(dat$X)))",
    "  stop('Excluded accession remains in candidate matrix')",
    "cat('[ABLATION] ', ABLATION_SCENARIO, ' removed=',",
    "    paste(ABLATION_EXCLUDED_ACCESSIONS, collapse = ','),",
    "    ' remaining_features=', ncol(dat$X), '\\n', sep = '')",
    "# ---- END AUDITED ABLATION INSERTION ----"
  )
  append(lines, insertion, after = at)
}

copy_if_identical_or_new <- function(src, dest) {
  if (is_file(dest)) {
    if (!identical(file_md5(src), file_md5(dest))) {
      stop_here("Existing staged file has different contents: ", dest,
                "\nDo not overwrite code associated with cached results. ",
                "Use a fresh result folder after inspecting the difference.")
    }
    return(invisible(dest))
  }
  if (file.exists(dest)) stop_here("Destination exists but is not a file: ", dest)
  if (!file.copy(src, dest, overwrite = FALSE) ||
      !identical(file_md5(src), file_md5(dest))) {
    stop_here("Unable to stage verified copy: ", src, " -> ", dest)
  }
  invisible(dest)
}

stage_scenario <- function(name, project, out_root, harness, engine, launcher) {
  if (!name %in% names(SCENARIOS)) stop_here("Unknown scenario: ", name)
  folder <- file.path(out_root, name)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(folder)) stop_here("Unable to create directory: ", folder)

  original_lines <- readLines(harness, warn = FALSE, encoding = "UTF-8")
  altered <- patch_original_harness(original_lines, unname(SCENARIOS[[name]]), name)
  staged_harness <- file.path(folder, ORIGINAL_HARNESS)
  if (is_file(staged_harness)) {
    old <- readLines(staged_harness, warn = FALSE, encoding = "UTF-8")
    if (!identical(old, altered)) {
      stop_here("Different staged v1.1 code already exists: ", staged_harness,
                "\nDo not mix earlier and current analysis caches.")
    }
  } else {
    writeLines(altered, staged_harness, useBytes = TRUE)
  }

  staged_engine <- file.path(folder, ORIGINAL_ENGINE)
  staged_launcher <- file.path(folder, ORIGINAL_LAUNCHER)
  copy_if_identical_or_new(engine, staged_engine)
  copy_if_identical_or_new(launcher, staged_launcher)

  manifest <- data.frame(
    scenario = name,
    excluded_proteins = paste(names(SCENARIOS[[name]]), collapse = ";"),
    excluded_accessions = paste(unname(SCENARIOS[[name]]), collapse = ";"),
    source_harness_md5 = file_md5(harness),
    staged_harness_md5 = file_md5(staged_harness),
    engine_md5 = file_md5(engine),
    launcher_md5 = file_md5(launcher),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(folder, "ablation_definition.csv")
  if (is_file(manifest_path)) {
    old <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    if (!identical(old, manifest)) {
      stop_here("Existing scenario manifest differs: ", manifest_path)
    }
  } else {
    utils::write.csv(manifest, manifest_path, row.names = FALSE)
  }
  message("[PREPARED] ", name, " | removed: ",
          if (length(SCENARIOS[[name]])) paste(names(SCENARIOS[[name]]), collapse = ", ")
          else "none (baseline)", " | ", folder)
  folder
}

# Validate that every selected feature and fitted coefficient excludes the
# forbidden accessions. For smoke mode, validate the one-family RDS cache.
verify_tables <- function(tables, excluded, label) {
  if (!length(excluded)) return(invisible(TRUE))
  for (nm in names(tables)) {
    tab <- tables[[nm]]
    if (!is.data.frame(tab) || !"FeatureID" %in% names(tab)) {
      stop_here(label, " lacks FeatureID table: ", nm)
    }
    leaked <- intersect(excluded, as.character(tab$FeatureID))
    if (length(leaked)) {
      stop_here(label, " contains excluded feature(s) in ", nm, ": ",
                paste(leaked, collapse = ", "))
    }
  }
  invisible(TRUE)
}

verify_full_run <- function(output, name) {
  exclude <- unname(SCENARIOS[[name]])
  csv_names <- c("outer_consensus30_features.csv",
             "outer_diverse_member_features.csv",
             "outer_consensus_member_coefficients.csv",
             "outer_diverse_member_coefficients.csv")
  tables <- lapply(csv_names, function(nm) {
    path <- file.path(output, nm)
    if (!is_file(path)) stop_here("Missing audit file: ", path)
    utils::read.csv(path, stringsAsFactors = FALSE)
  })
  names(tables) <- csv_names
  verify_tables(tables, exclude, name)
  message("[VERIFIED] ", name,
          ": excluded accessions absent from outer features and coefficients")
}

verify_smoke <- function(output, name) {
  cache <- list.files(file.path(output, "outer_cache"),
                      pattern = "^family_.*\\.rds$", full.names = TRUE)
  if (length(cache) != 1L) {
    stop_here("Smoke test must create one family cache; found ", length(cache))
  }
  value <- readRDS(cache[[1L]])
  tables <- list(consensus_features = value$selected_features,
                 diverse_features = value$diverse_member_features,
                 consensus_coefficients = value$member_coefficients,
                 diverse_coefficients = value$diverse_member_coefficients)
  verify_tables(tables, unname(SCENARIOS[[name]]), paste0(name, " smoke"))
  message("[VERIFIED] One-family smoke cache; no excluded features remain")
}

# Keep the v1.1 launcher (six workers, caffeinate on Mac, safe resume).
# We pass an explicit output path so its legacy Claude/ fallback is never used.
execute_scenario <- function(name, action, project, out_root,
                             harness, engine, launcher) {
  folder <- stage_scenario(name, project, out_root, harness, engine, launcher)
  data <- find_data_dir(project)
  family <- find_family_map(data)
  smoke <- identical(action, "smoke")
  output <- file.path(folder, if (smoke) "SMOKE_B3_one_outer_family" else "results")
  env <- c(
    CITRIN_DATA_DIR = data,
    CITRIN_FAMILY_FILE = family,
    CITRIN_ENSEMBLE_ENGINE = file.path(folder, ORIGINAL_ENGINE),
    CITRIN_ENSEMBLE_OUTPUT = output,
    CITRIN_ENSEMBLE_COMBAT = "no_combat",
    CITRIN_ENSEMBLE_SEEDS = MEMBER_SEEDS,
    CITRIN_ENSEMBLE_STABILITY_B = if (smoke) "3" else "30",
    CITRIN_ENSEMBLE_SUBSAMPLE = "sample",
    CITRIN_ENSEMBLE_OUTER_LIMIT = if (smoke) "1" else "0",
    CITRIN_ENSEMBLE_FIT_FINAL = "false",
    CITRIN_RUN2 = "Run2",
    CITRIN_RUN3 = "Run3",
    CITRIN_EXPECT_44 = "true"
  )
  do.call(Sys.setenv, as.list(env))
  Sys.unsetenv("CITRIN_ENSEMBLE_INNER_FOLDS")
  if (!nzchar(Sys.which("Rscript"))) stop_here("Rscript not found on PATH")

  message("[START] ", action, " ", name, "\n        OUTPUT: ", output)
  message("[METHOD] v1.1 inner/outer family-LOFO, consensus30, ",
          "10 seeds, B=", if (smoke) "3 (SMOKE)" else "30", ", no ComBat")
  status <- suppressWarnings(system2("/bin/bash",
                                     args = shQuote(file.path(folder, ORIGINAL_LAUNCHER)),
                                     stdout = "", stderr = "", wait = TRUE))
  if (length(status) != 1L || is.na(status) || as.integer(status) != 0L) {
    stop_here("Six-worker launcher failed (status=", status,
              "). Inspect logs in: ", file.path(output, "logs"),
              "\nRerun the SAME scenario to resume completed family folds.")
  }

  if (smoke) {
    verify_smoke(output, name)
    message("[SMOKE PASSED] Operational test only; do not report AUC or sensitivity")
  } else {
    perf_path <- file.path(output, "outer_performance_all_policies.csv")
    if (!is_file(perf_path)) stop_here("Completed performance file missing: ", perf_path)
    verify_full_run(output, name)
    perf <- utils::read.csv(perf_path, stringsAsFactors = FALSE)
    primary <- perf[perf$strategy == "consensus30" &
                      perf$policy == "full_inner_recall", , drop = FALSE]
    if (nrow(primary) != 1L || primary$n != 44L || primary$n_positive != 7L)
      stop_here("Unexpected primary performance row: ", perf_path)
    print(primary, row.names = FALSE)
  }
  invisible(output)
}

read_primary <- function(output, label) {
  required <- c("outer_performance_all_policies.csv",
                "outer_primary_predictions.csv", "immutable_analysis_plan.csv")
  for (name in required) {
    if (!is_file(file.path(output, name))) {
      stop_here(label, ": missing ", file.path(output, name))
    }
  }
  plan <- utils::read.csv(file.path(output, "immutable_analysis_plan.csv"),
                          colClasses = "character", stringsAsFactors = FALSE)
  key <- function(nm) {
    match <- plan$value[plan$name == nm]
    if (length(match) != 1L) stop_here(label, ": missing plan key ", nm)
    match
  }
  required_plan <- list(mode = "no_combat", seeds = MEMBER_SEEDS,
                        stability_B = "30", subsample = "sample",
                        inner_validation = "family_lofo", panel_size = "30",
                        outer_limit = "0", n_outer_families = "35",
                        positive_families = "5")
  for (nm in names(required_plan)) {
    if (!identical(key(nm), required_plan[[nm]])) {
      stop_here(label, ": unexpected ", nm, ": ", key(nm),
                " (expected ", required_plan[[nm]], ")")
    }
  }
  perf <- utils::read.csv(file.path(output, required[[1L]]),
                          stringsAsFactors = FALSE)
  primary <- perf[perf$strategy == "consensus30" &
                    perf$policy == "full_inner_recall", , drop = FALSE]
  if (nrow(primary) != 1L || primary$n != 44L ||
      primary$n_positive != 7L || primary$n_negative != 37L) {
    stop_here(label, ": invalid primary performance row")
  }
  pred <- utils::read.csv(file.path(output, required[[2L]]),
                          stringsAsFactors = FALSE)
  cols <- c("SampleID", "family_id", "batch", "group", "truth",
            "probability", "threshold", "predicted_class")
  if (!all(cols %in% names(pred)) || nrow(pred) != 44L ||
      anyDuplicated(pred$SampleID) || sum(pred$truth == 1L) != 7L ||
      any(!is.finite(pred$probability)) || any(!is.finite(pred$threshold))) {
    stop_here(label, ": invalid sample-level OOF prediction table")
  }
  list(label = label, stats = primary, prediction = pred, get = key)
}

# Existing v1.1 results in Claude/ are permitted as READ-ONLY baseline inputs.
# This function never creates or modifies files at any baseline path.
find_baseline_output <- function(project, root) {
  explicit <- Sys.getenv("CITRIN_ABLATION_BASELINE_OUTPUT", "")
  if (nzchar(explicit)) {
    if (!is_file(file.path(explicit, "outer_primary_predictions.csv"))) {
      stop_here("Explicit baseline output not found: ", explicit)
    }
    return(canonical(explicit))
  }
  candidates <- c(
    file.path(project, "result", "citrin_nested_ensemble_v1_1_no_combat_innerLOFO"),
    file.path(project, "result", "citrin_nested_ensemble_v1_1_no_combat_innerLOFO", "results"),
    file.path(root, "baseline", "results"),
    file.path(project, "Claude", "citrin_nested_ensemble_v1_1_no_combat_innerLOFO")
  )
  ok <- candidates[vapply(candidates, function(d)
    is_file(file.path(d, "outer_primary_predictions.csv")), logical(1L))]
  if (!length(ok)) {
    stop_here("Completed baseline not found. Set CITRIN_ABLATION_BASELINE_OUTPUT ",
              "to the original completed 10-seed v1.1 folder, or execute ",
              "'run baseline' first.")
  }
  canonical(ok[[1L]])
}

compare_scenarios <- function(project, root) {
  baseline_path <- find_baseline_output(project, root)
  message("[BASELINE / READ ONLY] ", baseline_path)
  baseline <- read_primary(baseline_path, "baseline")
  comparisons <- list()
  paired <- list()
  ordered <- c("baseline", setdiff(names(SCENARIOS), "baseline"))

  for (scenario in ordered) {
    if (identical(scenario, "baseline")) {
      x <- baseline
    } else {
      output <- file.path(root, scenario, "results")
      if (!is_file(file.path(output, "outer_primary_predictions.csv"))) {
        message("[SKIP] Scenario not completed: ", scenario)
        next
      }
      x <- read_primary(output, scenario)
      verify_full_run(output, scenario)
      for (nm in c("raw_md5", "sample_info_md5", "family_md5", "engine_md5")) {
        if (!identical(baseline$get(nm), x$get(nm))) {
          stop_here(scenario, ": cohort or engine differs from baseline: ", nm)
        }
      }
      b <- baseline$prediction[order(baseline$prediction$SampleID), , drop = FALSE]
      a <- x$prediction[match(b$SampleID, x$prediction$SampleID), , drop = FALSE]
      if (anyNA(a$SampleID) ||
          !identical(as.character(a$SampleID), as.character(b$SampleID)) ||
          !identical(as.character(a$family_id), as.character(b$family_id)) ||
          !identical(as.character(a$group), as.character(b$group)) ||
          !identical(as.integer(a$truth), as.integer(b$truth))) {
        stop_here(scenario, ": OOF patient/label/family mismatch")
      }
      paired[[scenario]] <- data.frame(
        scenario = scenario, SampleID = b$SampleID,
        family_id = b$family_id, group = b$group, truth = b$truth,
        baseline_score = b$probability, ablated_score = a$probability,
        delta_score = a$probability - b$probability,
        baseline_threshold = b$threshold, ablated_threshold = a$threshold,
        baseline_call = b$predicted_class, ablated_call = a$predicted_class,
        changed_classification = b$predicted_class != a$predicted_class
      )
    }
    row <- x$stats
    row$scenario <- scenario
    row$excluded_genes <- paste(names(SCENARIOS[[scenario]]), collapse = ";")
    row$delta_auc <- row$AUC - baseline$stats$AUC
    row$delta_sensitivity <- row$sensitivity - baseline$stats$sensitivity
    row$delta_specificity <- row$specificity - baseline$stats$specificity
    row$delta_brier <- row$Brier - baseline$stats$Brier
    comparisons[[scenario]] <- row
  }

  summary <- do.call(rbind, comparisons)
  columns <- c("scenario", "excluded_genes", "n", "n_positive", "n_negative",
               "AUC", "sensitivity", "specificity", "TP", "FN", "TN", "FP",
               "Brier", "delta_auc", "delta_sensitivity", "delta_specificity",
               "delta_brier")
  summary <- summary[, columns, drop = FALSE]
  rownames(summary) <- NULL
  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  summary_path <- file.path(root, "ABLA_AUROC_and_threshold_summary.csv")
  utils::write.csv(summary, summary_path, row.names = FALSE)
  if (length(paired)) {
    paired_path <- file.path(root, "ABLA_paired_OOF_per_sample.csv")
    utils::write.csv(do.call(rbind, paired), paired_path, row.names = FALSE)
  }
  cat("\n[COMPARISON] Out-of-family predictions; 44 samples, 35 families\n")
  print(summary[, c("scenario", "AUC", "TP", "FN", "TN", "FP",
                    "Brier", "delta_auc")], row.names = FALSE)
  cat("[CAUTION] Post-hoc sensitivity analysis within the original cohort.\n")
  cat("[OUTPUT] ", summary_path, "\n", sep = "")
  invisible(summary)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  action <- tolower(if (length(args)) args[[1L]] else
                      Sys.getenv("CITRIN_ABLATION_ACTION", "prepare"))
  scenario <- if (length(args) >= 2L) args[[2L]] else
                Sys.getenv("CITRIN_ABLATION_SCENARIO", "all4")
  project <- find_project()
  output_root <- find_output_root(project)

  if (action == "compare") {
    return(invisible(compare_scenarios(project, output_root)))
  }
  if (!action %in% c("prepare", "run", "smoke")) {
    stop_here("Usage: Rscript citrin_top4_ablation_familyLOFO.R ",
              "{prepare|smoke|run|compare} ",
              "{all4|no_PRL|no_MTTP|no_ARHGEF10|no_SEC31A|baseline}")
  }

  harness <- find_source_file(project, ORIGINAL_HARNESS,
                              Sys.getenv("CITRIN_ABLATION_HARNESS_SOURCE", ""))
  engine <- find_source_file(project, ORIGINAL_ENGINE,
                             Sys.getenv("CITRIN_ABLATION_ENGINE_SOURCE", ""))
  launcher <- find_source_file(project, ORIGINAL_LAUNCHER,
                               Sys.getenv("CITRIN_ABLATION_LAUNCHER_SOURCE", ""))
  message("[PROJECT] ", project)
  message("[RESULT ROOT] ", output_root)
  message("[SOURCE] harness=", harness)
  message("[SOURCE] engine=", engine)
  message("[SOURCE] launcher=", launcher)

  if (action == "prepare") {
    # A bare "prepare" stages all five principal ablations, independent of
    # any stale CITRIN_ABLATION_SCENARIO from a preceding smoke/run action.
    # To stage one scenario use: Rscript <driver> prepare <scenario>.
    targets <- if (length(args) >= 2L) {
      if (identical(scenario, "all")) names(SCENARIOS) else scenario
    } else {
      MAIN_SCENARIOS
    }
    for (s in targets) stage_scenario(s, project, output_root,
                                     harness, engine, launcher)
    message("\nNext commands:\n",
            "  Rscript citrin_top4_ablation_familyLOFO.R smoke all4\n",
            "  Rscript citrin_top4_ablation_familyLOFO.R run all4")
  } else {
    if (identical(scenario, "all")) {
      stop_here("Use one scenario at a time for run/smoke.")
    }
    execute_scenario(scenario, action, project, output_root,
                     harness, engine, launcher)
  }
  invisible(NULL)
}

main()
