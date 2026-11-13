#!/usr/bin/env Rscript
# ============================================================================
# Citrin deficiency | Figure 2 and supplementary plots from nested ensemble v1.1 (Observable revision v1.7: v1.6 + AUC 95% CI in panel C labels; legend wording for C/D)
# 2026-11-03
#
# RStudio: open 20261029_rebuttal.Rproj, then run:
#   source("script/citrin_ensemble_Figure2_v1_7.R", encoding = "UTF-8")
#
# All files are read from this project's citrin_nested_ensemble_v1_1/ folder.
# The source R script lives in script/; paths are PROJECT-ROOT based, not
# source-file based. No model refitting, feature search, or threshold tuning.
#
# IMPORTANT INTERPRETATION:
# - Consensus-30 was the analysis plan's PRIMARY strategy.
# - Diverse-30 was a comparison strategy selected for attention POST HOC.
# - Held-out predictions are family-held-out outer LOFO predictions.
# - "Diverse-30" uses 30 proteins PER MEMBER, not 30 distinct overall.
# - v1.8 (2026-11-08): panels B and E show the 30 proteins of the FINAL consensus-30 panel (final_consensus30_proteins_*.csv), not the 15 proteins of the superseded full-data refit.
# ============================================================================

options(stringsAsFactors = FALSE)

# --- Configuration (most users do not need to edit anything below) -----------
PROJECT_ROOT_OVERRIDE <- normalizePath(Sys.getenv("CITRIN_REPO", "."))   # repository root (public release)  # If running outside RStudio, put full project path here
RESULT_FOLDER <- file.path("01_nested_ensemble", "results")
OUTPUT_SUBFOLDER <- file.path("07_figures", "results", "Figure3_classifier_panels")
INCLUDE_ORIGINAL_PANELS_AB <- TRUE  # Uses ../data/raw_data.xlsx & sample_info.csv
READ_FINAL_MODEL_RDS <- FALSE       # TRUE only when the fitted RDS is needed
SAVE_PNG <- TRUE
DPI <- 300L

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop('Missing package "ggplot2". Run: install.packages("ggplot2")')
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop('Missing package "patchwork". Run: install.packages("patchwork")')
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

find_project_root <- function() {
  if (nzchar(PROJECT_ROOT_OVERRIDE)) {
    p <- normalizePath(PROJECT_ROOT_OVERRIDE, mustWork = TRUE)
    if (!file.exists(file.path(p, "20261029_rebuttal.Rproj"))) message("[NOTE] running in the public-repository layout: ", p)
    return(p)
  }
  candidates <- character(0)
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    rp <- tryCatch(rstudioapi::getActiveProject(), error = function(e) NULL)
    if (is.character(rp) && length(rp) == 1L && !is.na(rp) && nzchar(rp)) {
      candidates <- c(candidates, rp)
    }
  }
  current <- normalizePath(getwd(), mustWork = TRUE)
  for (i in seq_len(6L)) {
    candidates <- c(candidates, current)
    if (identical(dirname(current), current)) break
    current <- dirname(current)
  }
  candidates <- unique(candidates)
  matches <- candidates[file.exists(file.path(candidates, "20261029_rebuttal.Rproj"))]
  if (!length(matches)) {
    stop("Cannot find 20261029_rebuttal.Rproj. Open the RStudio project, ",
         "or set PROJECT_ROOT_OVERRIDE at the top of this script.")
  }
  normalizePath(matches[1L], mustWork = TRUE)
}

ROOT <- find_project_root()
RESULTS_DIR <- file.path(ROOT, RESULT_FOLDER)
OUT_DIR <- file.path(ROOT, OUTPUT_SUBFOLDER)
if (!dir.exists(RESULTS_DIR)) stop("Results directory not found: ", RESULTS_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("[ROOT]    ", ROOT, "\n", "[RESULTS] ", RESULTS_DIR,
    "\n[OUTPUT]  ", OUT_DIR, "\n", sep = "")

read_csv <- function(filename, required = FALSE) {
  fp <- file.path(RESULTS_DIR, filename)
  if (!file.exists(fp)) {
    if (required) stop("Required result not found: ", fp)
    message("[OPTIONAL MISSING] ", filename)
    return(NULL)
  }
  utils::read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE,
                  fileEncoding = "UTF-8-BOM")
}
require_fields <- function(df, fields, label) {
  miss <- setdiff(fields, names(df))
  if (length(miss)) stop(label, ": missing column(s) ", paste(miss, collapse = ", "))
  invisible(TRUE)
}
write_table <- function(df, filename) {
  utils::write.csv(df, file.path(OUT_DIR, filename), row.names = FALSE,
                   fileEncoding = "UTF-8")
}

# --- 1. Load the result files by exact filename -----------------------------
file_manifest <- c(
  plan                         = "immutable_analysis_plan.csv",
  performance                  = "outer_performance_all_policies.csv",
  predictions                  = "outer_predictions.csv",
  primary_predictions          = "outer_primary_predictions.csv",
  primary_positive7            = "outer_primary_positive7.csv",
  family_detection_primary     = "outer_positive_family_detection.csv",
  threshold_audit              = "outer_threshold_audit.csv",
  inner_training_sizes         = "outer_inner_training_size_audit.csv",
  inner_predictions            = "outer_inner_predictions.csv",
  inner_split_log              = "outer_inner_split_log.csv",
  inner_member_audit           = "outer_inner_member_audit.csv",
  outer_member_audit           = "outer_member_audit.csv",
  outer_consensus30_features   = "outer_consensus30_features.csv",
  outer_consensus_coefficients = "outer_consensus_member_coefficients.csv",
  outer_diverse_features       = "outer_diverse_member_features.csv",
  outer_diverse_coefficients   = "outer_diverse_member_coefficients.csv",
  final_consensus_features     = "final_consensus30_proteins_NOT_INDEPENDENTLY_VALIDATED.csv",
  final_consensus_coefficients = "final_consensus_member_coefficients.csv",
  final_diverse_features       = "final_diverse_member_features.csv",
  final_inner_thresholds       = "final_training_inner_thresholds_NOT_VALIDATION.csv",
  posthoc_threshold_frontier   = "POSTHOC_DO_NOT_CLAIM_VALIDATED_threshold_frontier.csv"
)
required_names <- c("plan", "performance", "predictions", "primary_predictions",
                    "inner_training_sizes", "outer_diverse_features")
citrin_ens <- setNames(lapply(seq_along(file_manifest), function(i) {
  read_csv(unname(file_manifest[i]), required = names(file_manifest)[i] %in% required_names)
}), names(file_manifest))
if (READ_FINAL_MODEL_RDS) {
  model_path <- file.path(RESULTS_DIR, "final_training_model_NOT_VALIDATED.rds")
  if (!file.exists(model_path)) stop("Final model RDS not found: ", model_path)
  citrin_ens$final_model_NOT_VALIDATED <- readRDS(model_path)
}
# citrin_ens remains available in RStudio's Environment for further analyses.
cat("[LOAD] Result CSVs read: ", sum(vapply(citrin_ens, is.data.frame, logical(1))), "\n", sep = "")

# --- 2. Structural and numerical integrity checks ---------------------------
pred <- citrin_ens$predictions
if (!"family_id" %in% names(pred)) {   # family identifiers are not distributed; supply family_map.csv via CITRIN_EXTERNAL_DIR
  fm <- utils::read.csv(file.path(Sys.getenv("CITRIN_EXTERNAL_DIR", file.path(ROOT, "data")), "family_map.csv"), stringsAsFactors = FALSE)
  pred$family_id <- fm$family_id[match(pred$SampleID, fm$SampleID)]
}
perf <- citrin_ens$performance
require_fields(pred, c("SampleID", "family_id", "group", "truth", "strategy",
                       "policy", "probability", "threshold", "predicted_class"),
               "outer_predictions.csv")
require_fields(perf, c("strategy", "policy", "AUC", "TP", "TN", "FP", "FN",
                       "sensitivity", "specificity"), "outer_performance_all_policies.csv")

pred$truth <- as.integer(pred$truth)
pred$predicted_class <- as.integer(pred$predicted_class)
pred$probability <- as.numeric(pred$probability)
pred$threshold <- as.numeric(pred$threshold)
pred$family_id <- as.character(pred$family_id)
pred$strategy <- as.character(pred$strategy)
pred$policy <- as.character(pred$policy)
STRATEGIES <- c("consensus30", "diverse30", "adaptive")
POLICIES <- c("full_inner_recall", "inner_spec90")
strategy_title <- c(consensus30 = "Consensus-30", diverse30 = "Diverse-30",
                    adaptive = "Adaptive")
strategy_cols <- c(consensus30 = "#4269D0", diverse30 = "#FF725C",
                   adaptive = "#3CA951")
if (!setequal(unique(pred$strategy), STRATEGIES) ||
    !setequal(unique(pred$policy), POLICIES)) {
  stop("Unexpected strategy/policy in predictions: verify correct v1.1 results folder")
}
if (anyNA(pred[c("truth", "predicted_class", "probability", "threshold")]) ||
    any(!is.finite(pred$probability)) || any(!is.finite(pred$threshold)) ||
    any(pred$probability < 0 | pred$probability > 1) ||
    any(pred$predicted_class != as.integer(pred$probability >= pred$threshold))) {
  stop("Missing/nonfinite scores or inconsistent prediction >= threshold classifications")
}
if (anyDuplicated(pred[c("SampleID", "strategy", "policy")])) {
  stop("Duplicate SampleID + strategy + policy rows in outer_predictions.csv")
}
for (st in STRATEGIES) for (pol in POLICIES) {
  z <- pred[pred$strategy == st & pred$policy == pol, , drop = FALSE]
  if (nrow(z) != 44L || sum(z$truth == 1L) != 7L ||
      sum(z$truth == 0L) != 37L || length(unique(z$family_id)) != 35L ||
      length(unique(z$family_id[z$truth == 1L])) != 5L) {
    stop("Incomplete/wrong cohort for strategy=", st, " policy=", pol)
  }
}
if (anyDuplicated(pred[c("SampleID", "strategy", "policy")])) stop("Duplicate prediction rows")

# Also check primary export is consensus, not diverse; prevents 6/7 vs 7/7 mixups.
pri <- citrin_ens$primary_predictions
require_fields(pri, c("SampleID", "strategy", "policy", "probability",
                      "threshold", "predicted_class"), "outer_primary_predictions.csv")
if (nrow(pri) != 44L || !all(pri$strategy == "consensus30") ||
    !all(pri$policy == "full_inner_recall")) {
  stop("The primary export is not 44 rows of consensus30/full_inner_recall")
}
prip <- pred[pred$strategy == "consensus30" &
               pred$policy == "full_inner_recall", , drop = FALSE]
ix <- match(pri$SampleID, prip$SampleID)
if (anyNA(ix) || any(abs(as.numeric(pri$probability) - prip$probability[ix]) > 1e-10) ||
    any(abs(as.numeric(pri$threshold) - prip$threshold[ix]) > 1e-10) ||
    any(as.integer(pri$predicted_class) != prip$predicted_class[ix])) {
  stop("Primary predictions disagree with all-policy outer_predictions.csv")
}

auc_rank <- function(y, score) {
  y <- as.integer(y); score <- as.numeric(score)
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L || anyNA(score)) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1L]) -
     n1 * (n1 + 1L) / 2) / (n1 * n0)
}
metric_rows <- lapply(STRATEGIES, function(st) {
  do.call(rbind, lapply(POLICIES, function(pol) {
    d <- pred[pred$strategy == st & pred$policy == pol, , drop = FALSE]
    dd <- d$predicted_class; yy <- d$truth
    positive_family <- unique(d$family_id[yy == 1L])
    data.frame(strategy = st, policy = pol, n = nrow(d),
               AUC = auc_rank(yy, d$probability),
               sensitivity = mean(dd[yy == 1L] == 1L),
               specificity = mean(dd[yy == 0L] == 0L),
               TP = sum(dd == 1L & yy == 1L), TN = sum(dd == 0L & yy == 0L),
               FP = sum(dd == 1L & yy == 0L), FN = sum(dd == 0L & yy == 1L),
               positive_families_detected = sum(vapply(positive_family, function(fid) {
                 any(dd[d$family_id == fid & yy == 1L] == 1L)
               }, logical(1L))))
  }))
})
metric_check <- do.call(rbind, metric_rows)
for (i in seq_len(nrow(metric_check))) {
  a <- metric_check[i, ]
  b <- perf[perf$strategy == a$strategy & perf$policy == a$policy, , drop = FALSE]
  if (nrow(b) != 1L) stop("Missing/duplicate summary row for ", a$strategy, "/", a$policy)
  if (any(abs(as.numeric(unlist(b[c("AUC", "sensitivity", "specificity")], use.names = FALSE)) -
              as.numeric(unlist(a[c("AUC", "sensitivity", "specificity")], use.names = FALSE))) > 1e-8) ||
      any(as.integer(unlist(b[c("TP", "TN", "FP", "FN")], use.names = FALSE)) !=
          as.integer(unlist(a[c("TP", "TN", "FP", "FN")], use.names = FALSE)))) {
    stop("Recalculated metrics disagree with outer_performance_all_policies.csv: ",
         a$strategy, "/", a$policy)
  }
}
write_table(metric_check, "AUDIT_recalculated_performance.csv")

plan <- citrin_ens$plan
if (!is.null(plan) && all(c("name", "value") %in% names(plan))) {
  plan_info <- setNames(as.character(plan$value), as.character(plan$name))
  if (!is.na(plan_info["inner_validation"]) &&
      plan_info["inner_validation"] != "family_lofo") {
    warning("Analysis plan is not inner family-LOFO. Inspect plan before publication.")
  }
  cat("[PLAN] inner=", plan_info["inner_validation"],
      "; stability B=", plan_info["stability_B"],
      "; mode=", plan_info["mode"], "\n", sep = "")
}
cat("[AUDIT] Predictions agree with performance CSV.\n")
print(metric_check[metric_check$policy == "full_inner_recall",
                   c("strategy", "TP", "FP", "FN", "TN", "AUC")], row.names = FALSE)

# --- 3. Figure 2C: ROC curves (Observable10 colors; no operating circles) ---
# Each ROC curve changes ONE global cutoff over pooled outer-OOF scores.
# Actual decisions use family-specific thresholds derived from outer training;
# their aggregate operating point is NOT a single-threshold point on this ROC.
# Therefore do not overlay these operating points; show decisions in Panel D.
roc_xy <- function(y, score) {
  th <- c(Inf, sort(unique(score), decreasing = TRUE), -Inf)
  data.frame(
    fpr = vapply(th, function(t) mean(score[y == 0L] >= t), numeric(1L)),
    tpr = vapply(th, function(t) mean(score[y == 1L] >= t), numeric(1L)),
    threshold = th
  )
}
main <- pred[pred$policy == "full_inner_recall", , drop = FALSE]
roc_dat <- do.call(rbind, lapply(STRATEGIES, function(st) {
  z <- main[main$strategy == st, , drop = FALSE]
  d <- roc_xy(z$truth, z$probability)
  d$strategy <- st
  d
}))
roc_dat$strategy <- factor(roc_dat$strategy, levels = STRATEGIES)
op <- metric_check[metric_check$policy == "full_inner_recall", , drop = FALSE]
op$fpr <- 1 - op$specificity
op$tpr <- op$sensitivity
op$strategy <- factor(op$strategy, levels = STRATEGIES)
# v1.7: family-cluster bootstrap 95% CIs of the pooled OOF AUC (2,000 replicates), from
# ensemble_v1_1_results/ensemble_v1_1_performance_with_CI.csv (rebuttal project; same values as Table 2).
ci_file <- file.path(RESULTS_DIR, "ensemble_v1_1_performance_with_CI.csv")
AUC_CI95 <- c(consensus30 = "0.710–0.964", diverse30 = "0.833–1.000", adaptive = "0.650–0.974")  # fallback (2026-11-03 run)
if (file.exists(ci_file)) {
  ci_tab <- utils::read.csv(ci_file, stringsAsFactors = FALSE)
  ci_tab <- ci_tab[ci_tab$comparison == "SM vs all controls", , drop = FALSE]
  for (st in names(AUC_CI95)) if (st %in% ci_tab$strategy) AUC_CI95[[st]] <- ci_tab$AUC_CI95_family_bootstrap[ci_tab$strategy == st][1]
  message("[CI] AUC 95% CIs read from ", ci_file)
}
legend_labels <- setNames(vapply(STRATEGIES, function(st) {
  zz <- op[as.character(op$strategy) == st, , drop = FALSE]
  sprintf("%s: AUC %.3f (95%% CI %s)", strategy_title[st], zz$AUC, AUC_CI95[[st]])
}, character(1L)), STRATEGIES)

plot_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 8.8))

# Observable10 subset: blue, coral red, green. Distinct line types aid reading.
strategy_linetypes <- c(consensus30 = "solid", diverse30 = "longdash",
                        adaptive = "dotdash")
pC <- ggplot(roc_dat, aes(x = fpr, y = tpr, color = strategy,
                           linetype = strategy, group = strategy)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4,
              linetype = "dashed", color = "#B9BDC5") +
  geom_step(direction = "hv", linewidth = 1.35, lineend = "butt") +
  scale_color_manual(values = strategy_cols, labels = legend_labels,
                     breaks = STRATEGIES, drop = FALSE) +
  scale_linetype_manual(values = strategy_linetypes, labels = legend_labels,
                        breaks = STRATEGIES, drop = FALSE) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(title = "C  Outer-family-held-out ROC",
       subtitle = "Pooled out-of-fold scores (exploratory nested ensembles)",
       x = "1 - Specificity", y = "Sensitivity") +
  guides(color = guide_legend(ncol = 1, byrow = TRUE,
                              override.aes = list(linewidth = 1.5))) +
  plot_theme +
  theme(legend.position = "bottom",
        legend.key.width = grid::unit(1.25, "cm"),
        legend.spacing.y = grid::unit(0.05, "cm"))

# --- 4. Figure 2D: seven LoF/Mis predictions (Observable10; continuous paired-family brackets) ---
# The signed margin compares each outer-OOF score with its OWN fold-specific,
# training-only threshold. A positive margin is a detected case.
# Family numeric IDs are internal serial identifiers, NOT biological family names.
# Therefore only SampleID is displayed; brackets connect within-family pairs.
pos <- main[main$truth == 1L, , drop = FALSE]
positive_ids <- paste0("SM_", seq_len(7L))
if (nrow(pos) != 7L * length(STRATEGIES) ||
    !setequal(unique(pos$SampleID), positive_ids) ||
    anyDuplicated(pos[c("SampleID", "strategy")])) {
  stop("Panel D requires exactly SM_1...SM_7 for each of the 3 strategies")
}

# Check family memberships from data, NOT from sequential family numbers.
family_by_sample <- vapply(positive_ids, function(s) {
  x <- unique(as.character(pos$family_id[pos$SampleID == s]))
  if (length(x) != 1L || is.na(x)) stop("Inconsistent family_id for ", s)
  x
}, character(1L))
names(family_by_sample) <- positive_ids
family_groups <- split(positive_ids, family_by_sample)
family_pairs <- Filter(function(x) length(x) > 1L, family_groups)
if (!identical(sort(unname(vapply(family_pairs, paste, collapse = "+", FUN.VALUE = character(1L)))),
               sort(c("SM_2+SM_3", "SM_4+SM_5")))) {
  stop("Expected two distinct paired families: SM_2/SM_3 and SM_4/SM_5. Review family_map.csv")
}

# ggplot uses numeric x and y so the brackets fit INSIDE the image.
pos$xnum <- match(as.character(pos$strategy), STRATEGIES)
pos$ynum <- length(positive_ids) + 1L - match(pos$SampleID, positive_ids)
pos$strategy <- factor(as.character(pos$strategy), levels = STRATEGIES)
pos$margin <- pos$probability - pos$threshold
pos$margin_label <- sprintf("%+.3f", pos$margin)
pos$detected <- factor(ifelse(pos$predicted_class == 1L, "Detected", "Missed"),
                       levels = c("Detected", "Missed"))

paired_brackets <- do.call(rbind, lapply(family_pairs, function(ss) {
  y <- length(positive_ids) + 1L - match(ss, positive_ids)
  if (length(y) != 2L || diff(sort(y)) != 1L) {
    stop("Relative samples must occupy adjacent rows in Panel D: ", paste(ss, collapse = ","))
  }
  data.frame(ylo = min(y), yhi = max(y))
}))
sample_labels <- data.frame(
  SampleID = positive_ids,
  ynum = length(positive_ids):1L
)

# Draw each pair as a SINGLE connected bracket (4 ordered vertices):
#   right tip -> top-left corner -> bottom-left corner -> right tip.
# Include most of the two tile rows to convey their family grouping.
# Do not map the internal numeric family IDs to displayed labels.
bracket_paths <- do.call(rbind, lapply(seq_len(nrow(paired_brackets)), function(i) {
  b <- paired_brackets[i, , drop = FALSE]
  data.frame(
    pair_id = i,
    x = c(0.495, 0.420, 0.420, 0.495),
    y = c(b$yhi + 0.37, b$yhi + 0.37,
          b$ylo - 0.37, b$ylo - 0.37)
  )
}))
stopifnot(nrow(bracket_paths) == 8L,
          !anyNA(bracket_paths),
          all(vapply(split(bracket_paths, bracket_paths$pair_id), nrow, integer(1)) == 4L))

# Pastel tints of Observable10: green = Detected, coral = Missed.
# On white backgrounds, 30% base color + 70% white ensures legible dark labels.
# Unlike Panel C, colors here code CLASSIFICATION, not algorithm identity.
status_cols <- c(Detected = "#C5E5CB", Missed = "#FFD5CE")  # 30% Observable10 over white
pD <- ggplot(pos, aes(x = xnum, y = ynum, fill = detected)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.94, height = 0.92) +
  geom_text(aes(label = margin_label), color = "#17212B", size = 3.4) +
  # Individual sample IDs (no arbitrary F13, F14, etc.)
  geom_text(data = sample_labels,
            aes(x = 0.30, y = ynum, label = SampleID),
            inherit.aes = FALSE, hjust = 1,
            color = "#273344", size = 3.5) +
  # One continuous square-bracket path PER paired family, avoiding seams
  # at the corners from 3 separately rendered geom_segment layers.
  geom_path(data = bracket_paths,
            aes(x = x, y = y, group = pair_id),
            inherit.aes = FALSE, color = "#3C4757", linewidth = 0.8,
            lineend = "round", linejoin = "round") +
  scale_fill_manual(values = status_cols, drop = FALSE, name = NULL) +
  scale_x_continuous(limits = c(0.0, 3.50), breaks = seq_along(STRATEGIES),
                     labels = unname(strategy_title[STRATEGIES]), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.49, 7.51), breaks = NULL,
                     expand = c(0, 0)) +
  labs(title = "D  Detection of 7 LoF/Mis specimens",
       subtitle = "Values = held-out score minus family-specific training-only threshold",
       caption = "Brackets indicate specimens belonging to the same family.",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "#4A5562"),
        plot.caption = element_text(size = 8.5, color = "#4A5562", hjust = 0),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 9.5, color = "#273344"),
        legend.position = "bottom",
        legend.text = element_text(size = 9.5),
        legend.key.width = grid::unit(0.65, "cm"),
        plot.margin = margin(9, 12, 6, 5))
write_table(pos, "PanelD_positive7_all_strategies.csv")

# --- 5. Figure 2E: feature recurrence from outer TRAINING fits --------------
original15 <- utils::read.csv(file.path(ROOT, "05_signature_proteins", "GO", "citrin_final30_panel.csv"), stringsAsFactors = FALSE)[, c("FeatureID", "Gene")]
if (nrow(original15) != 30L) stop("Expected the 30 proteins of the final consensus-30 panel")
features <- citrin_ens$outer_diverse_features
require_fields(features, c("held_family", "member_seed", "rank", "FeatureID"),
               "outer_diverse_member_features.csv")
features$held_family <- as.character(features$held_family)
features$member_seed <- as.character(features$member_seed)
key <- paste(features$held_family, features$member_seed, sep = "::")
member_counts <- table(key)
if (length(member_counts) != 350L || any(member_counts != 30L) ||
    anyDuplicated(features[c("held_family", "member_seed", "FeatureID")])) {
  stop("Feature selection data do not contain exactly 35 families x 10 members x 30")
}
counts_by_fold <- aggregate(FeatureID ~ held_family, data = features,
                            FUN = function(z) length(unique(z)))
names(counts_by_fold)[2L] <- "n_distinct_proteins_across_10_members"
write_table(counts_by_fold, "AUDIT_distinct_proteins_per_outer_fold.csv")

coefs <- citrin_ens$outer_diverse_coefficients
has_coefs <- !is.null(coefs)
if (has_coefs) require_fields(coefs, c("held_family", "member_seed", "FeatureID", "nonzero"),
                             "outer_diverse_member_coefficients.csv")
rec <- original15
rec$selected_members <- vapply(rec$FeatureID, function(id) sum(features$FeatureID == id), integer(1L))
rec$selection_rate <- rec$selected_members / length(member_counts)
if (has_coefs) {
  nz <- as.logical(coefs$nonzero)
  if (anyNA(nz)) stop("Unparseable nonzero values in outer_diverse_member_coefficients.csv")
  rec$nonzero_members <- vapply(rec$FeatureID, function(id) {
    sum(coefs$FeatureID == id & nz)
  }, integer(1L))
  rec$nonzero_rate <- rec$nonzero_members / length(member_counts)
} else {
  rec$nonzero_members <- NA_integer_
  rec$nonzero_rate <- NA_real_
}
rec <- rec[order(rec$selection_rate, decreasing = TRUE), , drop = FALSE]
rec$Gene <- factor(rec$Gene, levels = rev(as.character(rec$Gene)))
write_table(rec, "PanelE_final30_feature_selection_recurrence.csv")

# Figure 2E uses Observable10 colors with two distinct visual marks:
# pastel teal bars = selected among top 30 per ensemble member;
# blue circles = nonzero elastic-net coefficient.
# Percentage labels live just OUTSIDE the 100% axis boundary.
# This does NOT change feature frequencies or any model results.
OBS_TEAL_TINT <- "#B6E2D8"  # pastel based on Observable10 teal #6CC5B0
OBS_BLUE <- "#4269D0"       # Observable10 blue
rec$selection_pct <- sprintf("%.0f%%", 100 * rec$selection_rate)

pE <- ggplot(rec, aes(y = Gene)) +
  geom_vline(xintercept = c(.25, .50, .75, 1.00),
             color = "#E7EBF0", linewidth = .30) +
  geom_col(aes(x = selection_rate, fill = "Top-30 candidate selection"),
           width = .64, show.legend = TRUE) +
  # Keep percentage labels outside the plotting area, in the right margin;
  # DO NOT stretch the horizontal x axis beyond 100% for these labels.
  geom_text(aes(x = 1.045, label = selection_pct),
            hjust = 0, size = 3.05, color = "#1F2D3D", show.legend = FALSE) +
  scale_fill_manual(values = c("Top-30 candidate selection" = OBS_TEAL_TINT),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep,
                     breaks = c(0, .25, .50, .75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%"),
                     expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(0, 1), clip = "off", expand = FALSE) +
  labs(title = "E  Recurrence of the final consensus-30 panel proteins",
       subtitle = "350 dependent fits: 35 outer-family partitions x 10 ensemble members",
       x = "Proportion of ensemble member fits", y = NULL) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12, color = "#17212B"),
        plot.subtitle = element_text(size = 8.7, color = "#4A5562"),
        axis.text.y = element_text(size = 9.3, color = "#17212B"),
        axis.text.x = element_text(color = "#273344"),
        axis.title.x = element_text(color = "#273344"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 9, color = "#273344"),
        legend.key.width = grid::unit(.65, "cm"),
        legend.spacing.x = grid::unit(.3, "cm"),
        # Extra right gutter prevents the percentages from being clipped.
        plot.margin = margin(t = 6, r = 58, b = 6, l = 6))

if (has_coefs) {
  if (any(rec$nonzero_rate > rec$selection_rate + 1e-9, na.rm = TRUE)) {
    stop("Nonzero frequency exceeds candidate selection frequency; inspect coefficient file")
  }
  # White halo makes dots distinct at all positions; blue dot is the data mark.
  # Separate fill vs color scales produce an actual bar swatch and a blue dot
  # in the legend, not squares with superimposed empty circles.
  pE <- pE +
    geom_point(aes(x = nonzero_rate), shape = 16, size = 3.75,
               color = "white", show.legend = FALSE) +
    geom_point(aes(x = nonzero_rate, color = "Non-zero coefficient"),
               shape = 16, size = 2.9, show.legend = TRUE) +
    scale_color_manual(values = c("Non-zero coefficient" = OBS_BLUE),
                       name = NULL) +
    guides(fill = guide_legend(order = 1),
           color = guide_legend(order = 2))
} else {
  message("[PANEL E] Coefficient CSV missing: showing candidate-selection bars only")
}

# --- 6. Optional original panels A and B (do NOT use ensemble for these) ----
# This does not change their exploratory status or imply feature validation.
find_data_directory <- function() {
  candidates <- c(file.path(ROOT, "data"),
                  file.path(dirname(ROOT), "data"),
                  file.path(dirname(dirname(ROOT)), "data"))
  candidates <- unique(candidates)
  matches <- candidates[file.exists(file.path(candidates, "raw_data.xlsx")) &
                        file.exists(file.path(candidates, "sample_info.csv"))]
  if (length(matches)) normalizePath(matches[1L], mustWork = TRUE) else NULL
}
find_col <- function(df, candidates, required = TRUE) {
  norm <- function(s) gsub("[^a-z0-9]", "", tolower(s))
  ix <- match(norm(candidates), norm(names(df)))
  ix <- ix[!is.na(ix)]
  if (length(ix)) return(names(df)[ix[1L]])
  if (required) stop("Cannot locate column: ", paste(candidates, collapse = ", "))
  NULL
}

make_original_ab <- function(data_dir) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop('Original panels A/B need openxlsx: install.packages("openxlsx")')
  }
  meta <- utils::read.csv(file.path(data_dir, "sample_info.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
  meta$SampleID <- toupper(trimws(as.character(meta[[find_col(meta, c("SampleID", "sample_id"))]])))
  meta$group <- toupper(trimws(as.character(meta[[find_col(meta, "group")]])))
  if (nrow(meta) != 60L || anyDuplicated(meta$SampleID)) stop("Expect 60 unique sample IDs")
  tab <- table(factor(meta$group, levels = c("WW", "WS", "SS", "SM", "DD")))
  if (!identical(as.integer(tab), c(6L, 10L, 16L, 7L, 21L))) {
    stop("Unexpected sample_info.csv group composition")
  }
  raw <- openxlsx::read.xlsx(file.path(data_dir, "raw_data.xlsx"), check.names = FALSE)
  accession_col <- find_col(raw, c("Master Accession Number", "master_accession_number"))
  accessions <- toupper(trimws(as.character(raw[[accession_col]])))
  if (anyDuplicated(accessions)) {
    raw <- raw[!duplicated(accessions), , drop = FALSE]
    accessions <- accessions[!duplicated(accessions)]
  }
  actb_i <- match("P60709", accessions)
  gene_col <- find_col(raw, c("Master Gene Symbol", "Gene names", "Gene names (primary)",
                              "Gene", "Gene name"), required = FALSE)
  slc_i <- if (!is.null(gene_col)) {
    which(toupper(trimws(as.character(raw[[gene_col]]))) == "SLC25A13")
  } else integer(0)
  if (!length(slc_i)) slc_i <- which(accessions == "Q9UJS0")
  if (is.na(actb_i) || !length(slc_i)) stop("ACTB or SLC25A13 protein not found")
  col_names <- toupper(trimws(names(raw)))
  sample_cols <- match(meta$SampleID, col_names)
  if (anyNA(sample_cols)) stop("Sample column(s) missing from raw_data.xlsx")
  actb <- suppressWarnings(as.numeric(unlist(raw[actb_i, sample_cols, drop = FALSE],
                                              use.names = FALSE)))
  slc <- as.matrix(raw[slc_i, sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(slc) <- "numeric")
  if (nrow(slc) == 1L) slc <- as.numeric(slc[1L, ]) else
    slc <- apply(slc, 2L, median, na.rm = TRUE)
  if (any(!is.finite(actb) | actb <= 0)) stop("Invalid ACTB intensities")
  slc[!is.finite(slc) | slc <= 0] <- 0
  meta$slc_actb <- as.numeric(slc) / actb

  group_levels <- c("Healthy", "Carrier", "LoF/LoF", "LoF/Mis", "Disease control")
  group_cols <- c("Healthy" = "#8A9AA9", "Carrier" = "#D3A75D",
                  "LoF/LoF" = "#518D78", "LoF/Mis" = "#4EA0AE",
                  "Disease control" = "#A87083")
  group_map <- setNames(group_levels, c("WW", "WS", "SS", "SM", "DD"))
  meta$group_label <- factor(unname(group_map[meta$group]), levels = group_levels)
  summary_a <- aggregate(slc_actb ~ group_label, meta, median)
  panel_a <- ggplot(meta, aes(x = group_label, y = slc_actb, color = group_label)) +
    geom_jitter(width = .13, height = 0, size = 2, alpha = .85) +
    geom_point(data = summary_a, aes(x = group_label, y = slc_actb),
               inherit.aes = FALSE, shape = 95, size = 8, color = "#1C2429") +
    scale_color_manual(values = group_cols, guide = "none") +
    scale_x_discrete(labels = c("Healthy\nn=6", "Carrier\nn=10", "LoF/LoF\nn=16",
                                "LoF/Mis\nn=7", "Disease\nn=21")) +
    labs(title = "A  SLC25A13 abundance in DBS", x = NULL,
         y = "SLC25A13 / ACTB (ratio)") + plot_theme +
    theme(axis.text.x = element_text(size = 8.3))

  meta44 <- meta[meta$group != "SS", , drop = FALSE]
  ord <- match(meta44$SampleID, toupper(trimws(names(raw))))
  prot_idx <- match(original15$FeatureID, accessions)
  if (anyNA(prot_idx)) stop("A protein of the final panel is missing in raw_data.xlsx")
  x <- as.matrix(raw[prot_idx, ord, drop = FALSE])
  suppressWarnings(storage.mode(x) <- "numeric")
  rownames(x) <- original15$Gene
  colnames(x) <- meta44$SampleID
  actb44 <- actb[match(meta44$SampleID, meta$SampleID)]
  x[!is.finite(x) | x <= 0] <- NA_real_
  x <- sweep(log2(x), 2L, log2(actb44), "-")
  z <- t(vapply(seq_len(nrow(x)), function(i) {
    v <- x[i, ]; ok <- is.finite(v)
    if (sum(ok) < 3L || sd(v[ok]) == 0) return(rep(NA_real_, length(v)))
    (v - mean(v[ok])) / sd(v[ok])
  }, numeric(ncol(x))))
  rownames(z) <- rownames(x); colnames(z) <- colnames(x)
  display_z <- z
  display_z[!is.finite(display_z)] <- 0  # Display-order calculation only
  col_order <- colnames(display_z)[hclust(dist(t(display_z)))$order]
  # v1.8: rows ordered by selection stability (outer partitions in which the protein entered the consensus panel), then gene
  cons_feat0 <- utils::read.csv(file.path(RESULTS_DIR, "outer_consensus30_features.csv"))
  n_in0 <- vapply(original15$FeatureID, function(id) length(unique(cons_feat0$held_family[cons_feat0$FeatureID == id])), integer(1L))
  gene_order <- rev(original15$Gene[order(-n_in0, original15$Gene)])
  z <- z[gene_order, col_order, drop = FALSE]
  heat <- data.frame(x = rep(seq_along(col_order), each = nrow(z)),
                     y = rep(seq_len(nrow(z)), times = ncol(z)),
                     z = as.vector(z))
  heat$z <- pmin(3, pmax(-3, heat$z))
  anno <- meta44[match(col_order, meta44$SampleID), , drop = FALSE]
  anno$x <- seq_along(col_order); anno$y <- 31
  # v1.6: feature-selection stability tier per original protein, from the 35 outer
  # consensus-30 candidate panels (outer_consensus30_features.csv). Exploratory tier,
  # NOT a biological classification: A = in every fold; B = in 33-34 folds but absent
  # when one particular positive family is held out; C = in fewer than one third of folds.
  cons_feat <- utils::read.csv(file.path(RESULTS_DIR, "outer_consensus30_features.csv"))
  n_outer <- length(unique(cons_feat$held_family))
  pos_fams <- c("13", "14", "17", "21", "22")
  tier_tab <- do.call(rbind, lapply(seq_len(nrow(original15)), function(i) {
    hf <- unique(as.character(cons_feat$held_family[cons_feat$FeatureID == original15$FeatureID[i]]))
    miss <- pos_fams[!pos_fams %in% hf]
    n_in <- length(hf)
    tier <- if (n_in == n_outer) "A" else if (n_in >= 30L) "B" else "C"
    data.frame(Gene = original15$Gene[i], n_in_panel = n_in, n_outer = n_outer,
               absent_when_held_out = paste(miss, collapse = ";"), tier = tier)
  }))
  write_table(tier_tab, "PanelB_final30_stability_tiers.csv")
  tier_cols <- c(A = "#2A6F97", B = "#E9A23B", C = "#B8B8B8")
  tier_df <- data.frame(y = seq_along(gene_order),
                        tier = tier_tab$tier[match(gene_order, tier_tab$Gene)],
                        n_in = tier_tab$n_in_panel[match(gene_order, tier_tab$Gene)])
  tier_df$x <- length(col_order) + 1.35
  tier_df$label <- paste0(tier_df$tier, " ", tier_df$n_in, "/", n_outer)
  panel_b <- ggplot(heat, aes(x = x, y = y)) +
    geom_tile(aes(fill = z), width = .97, height = .97) +
    geom_point(data = anno, aes(x = x, y = y, color = group_label),
               inherit.aes = FALSE, shape = 15, size = 2.2) +
    geom_tile(data = tier_df, aes(x = x, y = y), fill = tier_cols[tier_df$tier],
              width = .6, height = .97, inherit.aes = FALSE) +
    geom_text(data = tier_df, aes(x = x + 0.55, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, size = 2.4, color = "#1C2429") +
    annotate("text", x = length(col_order) + 1.35, y = 31, label = "Selection\ntier",
             size = 2.4, hjust = 0.2, vjust = 0.4, color = "#1C2429", lineheight = 0.85) +
    scale_fill_gradient2(low = "#32678F", mid = "#F8FAFC", high = "#BD6059",
                         midpoint = 0, limits = c(-3, 3), na.value = "#D2D2D2",
                         name = "z-score") +
    scale_color_manual(values = group_cols, drop = FALSE, name = "Group") +
    scale_y_continuous(breaks = 1:30, labels = gene_order,
                       limits = c(.5, 31.5), expand = expansion(mult = 0)) +
    scale_x_continuous(breaks = NULL, expand = expansion(add = .55),
                       limits = c(0.45, length(col_order) + 5.2)) +
    labs(title = "B  Final consensus-30 panel (descriptive)",
         x = "44 DBS specimens (clustered)", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 12),
          panel.grid = element_blank(), legend.position = "bottom",
          axis.text.y = element_text(size = 7.2),
          legend.text = element_text(size = 7.8),
          legend.title = element_text(size = 8))
  list(A = panel_a, B = panel_b)
}

orig <- NULL
if (INCLUDE_ORIGINAL_PANELS_AB) {
  dd <- find_data_directory()
  if (is.null(dd)) {
    message("[INFO] Raw data not found in project/data or ../data: skipping A/B. ",
            "C/D/E will still be generated.")
  } else if (!requireNamespace("openxlsx", quietly = TRUE)) {
    message("[INFO] openxlsx not installed: skipping original A/B. ",
            "install.packages('openxlsx') to include them.")
  } else {
    orig <- tryCatch(make_original_ab(dd), error = function(e) {
      warning("Could not regenerate original panels A/B (C/D/E unaffected): ",
              conditionMessage(e))
      NULL
    })
  }
}

# --- 7. Supplementary plots: group-specific errors and practical constraints -
negative <- main[main$truth == 0L, , drop = FALSE]
neg_summary <- do.call(rbind, lapply(STRATEGIES, function(st) {
  z <- negative[negative$strategy == st, , drop = FALSE]
  do.call(rbind, lapply(c("DD", "WS", "WW"), function(gr) {
    q <- z[z$group == gr, , drop = FALSE]
    data.frame(strategy = st, group = gr, n = nrow(q),
               false_positives = sum(q$predicted_class == 1L),
               correctly_negative = sum(q$predicted_class == 0L))
  }))
}))
write_table(neg_summary, "Supplement_false_positive_groups.csv")
write_table(negative[negative$predicted_class == 1L, c(
  "SampleID", "family_id", "batch", "group", "strategy", "policy",
  "probability", "threshold", "predicted_class")],
  "Supplement_false_positive_sample_ids.csv")

neg_summary$strategy <- factor(neg_summary$strategy, levels = STRATEGIES)
neg_summary$group <- factor(neg_summary$group, levels = c("DD", "WS", "WW"))
neg_summary$label <- paste0(neg_summary$false_positives, "/", neg_summary$n)
pS1 <- ggplot(neg_summary, aes(x = group, y = false_positives, fill = strategy)) +
  geom_col(position = position_dodge(width = .8), width = .68) +
  geom_text(aes(label = label), position = position_dodge(width = .8),
            vjust = -.35, size = 3.1) +
  scale_fill_manual(values = strategy_cols, labels = strategy_title, name = NULL) +
  scale_x_discrete(labels = c(DD = "Disease controls", WS = "Carriers", WW = "Healthy adults")) +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) +
  labs(title = "Supplement: False positives by control subgroup",
       subtitle = "Training-derived full-inner-recall thresholds; n shown in numerator/denominator",
       x = NULL, y = "False-positive specimens") + plot_theme

pS2 <- ggplot(counts_by_fold, aes(x = n_distinct_proteins_across_10_members)) +
  geom_histogram(binwidth = 5, boundary = 0, color = "white", fill = "#668FA4") +
  labs(title = "Supplement: Number of distinct proteins per outer fold",
       subtitle = "Diverse-30 = ten separately selected 30-protein member models",
       x = "Distinct FeatureIDs across 10 models", y = "Number of outer families") +
  plot_theme

pol_plot <- metric_check
pol_plot$strategy <- factor(pol_plot$strategy, levels = STRATEGIES)
pol_plot$policy <- factor(pol_plot$policy, levels = POLICIES,
                          labels = c("Inner 100% recall target", "Inner >=90% specificity target"))
pol_long <- rbind(
  data.frame(pol_plot[, c("strategy", "policy")], metric = "Sensitivity",
             value = pol_plot$sensitivity),
  data.frame(pol_plot[, c("strategy", "policy")], metric = "Specificity",
             value = pol_plot$specificity)
)
pS3 <- ggplot(pol_long, aes(x = strategy, y = value, color = policy, group = policy)) +
  geom_point(position = position_dodge(width = .35), size = 2.6) +
  geom_line(position = position_dodge(width = .35), linewidth = .5) +
  scale_x_discrete(labels = strategy_title) +
  scale_y_continuous(limits = c(0, 1.03), breaks = seq(0, 1, by = .2)) +
  scale_color_manual(values = c("#BC7457", "#4A829B"), name = NULL) +
  facet_wrap(~ metric) +
  labs(title = "Supplement: Threshold-policy sensitivity analysis",
       subtitle = "Outer-family-held-out performance; each dot is one procedure",
       x = NULL, y = "Fraction") + plot_theme

# --- 8. Save figures and publication-supporting data ------------------------
save_plot <- function(p, base, width, height) {
  ggplot2::ggsave(file.path(OUT_DIR, paste0(base, ".pdf")), plot = p,
                  width = width, height = height, units = "in",
                  device = "pdf", useDingbats = FALSE, bg = "white")
  if (SAVE_PNG) ggplot2::ggsave(file.path(OUT_DIR, paste0(base, ".png")), plot = p,
                               width = width, height = height, units = "in",
                               dpi = DPI, bg = "white")
}
save_plot(pC, "Figure2_panel_C_ROC", 7.5, 6.5)
save_plot(pD, "Figure2_panel_D_positive7", 7.5, 6.0)
save_plot(pE, "Figure2_panel_E_final30_recurrence", 8.5, 9.0)
save_plot(pS1, "Supplement_control_subgroup_false_positives", 9.5, 5.4)
save_plot(pS2, "Supplement_diverse_unique_proteins_per_fold", 8.0, 5.2)
save_plot(pS3, "Supplement_threshold_policies", 9.2, 5.3)

if (!is.null(orig)) {
  save_plot(orig$A, "Figure2_panel_A_SLC25A13", 7.4, 5.4)
  save_plot(orig$B, "Figure2_panel_B_final30_heatmap", 8.7, 8.0)
  # Five-panel figure; Figures C-E use the NEW ensemble, not legacy 10-seed CV.
  layout <- paste("AABBB", "CCDDD", "EEEEE", sep = "\n")
  fig <- patchwork::wrap_plots(A = orig$A, B = orig$B,
                               C = pC, D = pD, E = pE,
                               design = layout,
                               heights = c(1.6, 1.0, 1.5))
  save_plot(fig, "Figure2_revised_AB_ensemble_CDE_final30", 15.5, 17.5)
} else {
  fig <- (pC | pD) / pE + patchwork::plot_layout(heights = c(1.12, 1.2))
  save_plot(fig, "Figure2_ensemble_CDE_only", 15.0, 10.3)
}

if (!is.null(citrin_ens$final_diverse_features)) {
  ff <- citrin_ens$final_diverse_features
  if ("FeatureID" %in% names(ff)) {
    unique_final <- length(unique(ff$FeatureID))
    write_table(data.frame(
      model = "Final all-44 trained diverse ensemble (NOT independently validated)",
      n_distinct_proteins = unique_final,
      note = "This is a model descriptor, NOT a performance estimate"
    ), "Supplement_final_fit_unique_proteins_NOT_VALIDATION.csv")
    cat("[FINAL MODEL] Distinct proteins in full-data diverse fit: ", unique_final,
        " (NOT independently validated)\n", sep = "")
  }
}

legend <- c(
  "Figure 2. Protein-level assessment of SLC25A13 and exploratory ensemble evaluation.",
  "(A) ACTB-normalized SLC25A13 abundance across 60 DBS specimens. Non-detects are",
  "displayed at zero under the assay conditions. Individual points and group medians are shown.",
  "(B) Heatmap of the 30 proteins of the final consensus-30 panel (Table S18).",
  "Protein-wise standardized log2 abundances are shown; gray is",
  "unobserved. The right-hand strip gives each protein's feature-selection stability",
  "tier across the 35 outer training partitions of the nested consensus ensemble:",
  "A, in the 30-protein candidate panel in every partition (35/35); B, in 30-34 partitions;",
  "C, in fewer than 30 partitions.",
  "Tiers are exploratory descriptors of selection stability in this cohort, not a",
  "biological classification. This heatmap is descriptive, not independent classifier validation.",
  "(C) Pooled out-of-fold ROC curves from family-held-out nested validation for",
  "Consensus-30, Diverse-30, and training-selected Adaptive ensemble strategies, with the",
  "AUC and its family-cluster bootstrap 95% confidence interval (2,000 replicates, families",
  "resampled from the saved out-of-fold predictions; conditional on the fitted models) for each.",
  "Each ROC curve sweeps a single global cutoff over pooled outer-OOF scores;",
  "individual family-specific training-only decision thresholds are not plotted.",
  "The corresponding threshold-based decisions are shown separately in panel D.",
  "(D) Classification of the seven LoF/Mis specimens from their family-held-out",
  "out-of-fold probabilities. Each cell contains the difference between the specimen's",
  "out-of-fold probability and the training-derived threshold of its outer fold;",
  "a nonnegative value is classified as positive. Adjacent brackets indicate",
  "specimens from the same family (SM_2/SM_3 and SM_4/SM_5).",
  "(E) Selection frequencies of the 30 final-panel proteins across 350",
  "correlated member fits (35 outer training partitions x 10 seeds). Bars denote",
  "candidate selection frequencies; if available, white circles denote nonzero",
  "elastic-net coefficient frequencies. These are not independent validation.",
  "All ensemble analyses are post hoc and exploratory; Consensus-30 was specified",
  "as the primary strategy. Diverse-30 uses 30 proteins PER ensemble member, with",
  "more than 30 distinct proteins used across its ten members. New external validation",
  "and fixed assay specifications are needed for clinical interpretation."
)
writeLines(legend, file.path(OUT_DIR, "Figure2_legend_EN.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "R_sessionInfo.txt"))

cat("\n[COMPLETE] Files saved in: ", OUT_DIR, "\n", sep = "")
cat("[NOTE] Original panel A/B not built? Check for ../data/raw_data.xlsx + sample_info.csv.\n")
cat("[NOTE] Post hoc Diverse-30 is NOT the prespecified primary Consensus-30.\n")
