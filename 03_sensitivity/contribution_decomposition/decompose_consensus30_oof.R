#!/usr/bin/env Rscript
## ---- repository paths (public release): run from the repository root, or set CITRIN_REPO ----------------------
REPO <- normalizePath(Sys.getenv("CITRIN_REPO", "."), mustWork = TRUE)
DATA <- normalizePath(Sys.getenv("CITRIN_DATA_DIR", file.path(REPO, "data")), mustWork = TRUE)   # distributed inputs
EXT  <- Sys.getenv("CITRIN_EXTERNAL_DIR", DATA)   # inputs not distributed here: family_map.csv, DBS-Tohoku.pr_matrix.tsv, search FASTA, DIA-NN report
## -------------------------------------------------------------------------------------------------------------
add_family <- function(df) { if (!"family_id" %in% names(df)) { fm <- read.csv(file.path(EXT, "family_map.csv"), stringsAsFactors = FALSE); df$family_id <- fm$family_id[match(df$SampleID, fm$SampleID)] }; df }   # family identifiers are not distributed
# Per-specimen, per-protein decomposition of the consensus30 outer-OOF scores (v1.1 production run).
# No new model, no new selection: for every outer fold (held family) and every member seed, the training-only
# preprocessor (MinProb + standardisation) is rebuilt with the harness's exact seeds, the saved member coefficients are
# applied to the held-out specimens, and the member logit is split into intercept + sum_j beta_j * z_j.
# Verification: mean over members of plogis(logit) must reproduce outer_predictions.csv (consensus30) to ~1e-10.
suppressPackageStartupMessages({library(openxlsx); library(ggplot2)})
PKG  <- file.path(REPO, "01_nested_ensemble")
PROD <- file.path(REPO, "01_nested_ensemble", "results")
# DATA: defined in the repository header above
OUT  <- file.path(REPO, "03_sensitivity", "contribution_decomposition", "results")
ENGINE <- file.path(PKG, "citrin_44_nested_fixed_panel_v2_2.R")
stopifnot(unname(tools::md5sum(ENGINE)) == read.csv(file.path(PROD, "immutable_analysis_plan.csv"))$value[
  read.csv(file.path(PROD, "immutable_analysis_plan.csv"))$name == "engine_md5"])
SEEDS <- 2025:2034; MODE <- "no_combat"

Sys.setenv(CITRIN_DATA_DIR = DATA, CITRIN_FAMILY_FILE = file.path(EXT, "family_map.csv"),
           CITRIN_STAGE = "stage1", CITRIN_PANEL_MODE = "fixed30", CITRIN_SUBSAMPLE_UNIT = "sample",
           CITRIN_B = "30", CITRIN_SEED = "2025", CITRIN_N_BOOT = "0", CITRIN_OUTER_LIMIT = "0")
exprs <- parse(file = ENGINE, keep.source = FALSE)
stopifnot(identical(deparse(exprs[[length(exprs)]]), "main()"))
engine <- new.env(parent = globalenv())
for (j in seq_len(length(exprs) - 1L)) eval(exprs[[j]], envir = engine)
stopifnot(identical(engine$SCRIPT_VERSION, "two_run_fixed_panel_2.2"))
key_seed <- function(...) engine$family_seed(paste(..., sep = "|"), 11111L)

dat <- engine$read_44(); X <- dat$X; y <- as.integer(dat$y); meta <- dat$metadata
family <- as.character(meta$family_id); batch <- as.character(meta$batch); group <- as.character(meta$group)
stopifnot(nrow(X) == 44L, identical(rownames(X), meta$SampleID))

coefs <- read.csv(file.path(PROD, "outer_consensus_member_coefficients.csv"), colClasses = c(held_family = "character"))
pred  <- add_family(read.csv(file.path(PROD, "outer_predictions.csv"))); pred$family_id <- as.character(pred$family_id)
pred  <- pred[pred$strategy == "consensus30" & pred$policy == "full_inner_recall", ]
fams  <- sort(unique(family))

rows <- list(); check <- list()
for (held in fams) {
  tr <- which(family != held); te <- which(family == held)
  cf_f <- coefs[coefs$held_family == held, ]
  logit_members <- matrix(NA_real_, length(te), length(SEEDS), dimnames = list(rownames(X)[te], SEEDS))
  contrib_sum <- NULL
  for (k in seq_along(SEEDS)) {
    seed <- key_seed("fit", paste0("outer=", held, "/fit_all"), SEEDS[k])
    pre <- engine$make_preprocessor(X[tr, , drop = FALSE], group[tr], batch[tr], MODE, seed + 1L)
    zz  <- engine$score_preprocessor(pre, X[te, , drop = FALSE], batch[te], seed + 4L)
    cf  <- cf_f[cf_f$member_seed == SEEDS[k], ]
    b0  <- cf$coefficient[cf$FeatureID == "(Intercept)"]
    feats <- cf$FeatureID[cf$FeatureID != "(Intercept)"]; stopifnot(length(feats) == 30L, all(feats %in% colnames(zz)))
    beta <- setNames(cf$coefficient[match(feats, cf$FeatureID)], feats)
    contrib <- sweep(zz[, feats, drop = FALSE], 2L, beta, "*")           # specimens x 30, logit units
    logit_members[, k] <- b0 + rowSums(contrib)
    contrib_sum <- if (is.null(contrib_sum)) contrib else contrib_sum + contrib
    rows[[length(rows) + 1L]] <- data.frame(held_family = held, member_seed = SEEDS[k],
      SampleID = rep(rownames(zz), each = 30L), FeatureID = rep(feats, times = nrow(zz)),
      z = as.vector(t(zz[, feats, drop = FALSE])), beta = rep(beta, times = nrow(zz)),
      contribution = as.vector(t(contrib)), intercept = b0)
  }
  p_rebuilt <- rowMeans(plogis(logit_members))
  p_saved <- pred$probability[match(rownames(X)[te], pred$SampleID)]
  check[[held]] <- data.frame(held_family = held, SampleID = rownames(X)[te], p_saved = p_saved,
                              p_rebuilt = p_rebuilt, abs_diff = abs(p_saved - p_rebuilt))
  message(sprintf("[fold %s] n_held=%d max|dp|=%.2e", held, length(te), max(abs(p_saved - p_rebuilt))))
}
check <- do.call(rbind, check); long <- do.call(rbind, rows)
write.csv(check, file.path(OUT, "verification_rebuilt_vs_saved_probability.csv"), row.names = FALSE)
stopifnot(max(check$abs_diff) < 1e-8)
cat(sprintf("[VERIFIED] consensus30 OOF probabilities reproduced for 44/44 specimens; max |diff| = %.2e\n", max(check$abs_diff)))

# ---- per specimen x protein: mean contribution over the 10 members (logit units) ----
agg <- aggregate(cbind(contribution, z, beta) ~ SampleID + FeatureID, data = long, FUN = mean)
agg$n_members_nonzero <- aggregate(beta != 0 ~ SampleID + FeatureID, data = long, FUN = sum)[[3L]]
ic <- aggregate(intercept ~ SampleID, data = long, FUN = mean)
agg <- merge(agg, ic, by = "SampleID")
agg$group <- group[match(agg$SampleID, rownames(X))]; agg$family <- family[match(agg$SampleID, rownames(X))]
agg$truth <- y[match(agg$SampleID, rownames(X))]
agg$probability <- pred$probability[match(agg$SampleID, pred$SampleID)]
agg$threshold <- pred$threshold[match(agg$SampleID, pred$SampleID)]
agg$called_positive <- pred$predicted_class[match(agg$SampleID, pred$SampleID)]
# gene symbols
raw <- read.xlsx(file.path(DATA, "raw_data.xlsx"), sheet = "data", sep.names = " ")
sym <- setNames(raw[["Master Gene Symbol"]], raw[["Master accession number"]])
agg$gene <- unname(sym[agg$FeatureID])
# mean logit per specimen and share of each protein in (logit - intercept)
tot <- aggregate(contribution ~ SampleID, data = agg, FUN = sum); names(tot)[2] <- "sum_contribution"
agg <- merge(agg, tot, by = "SampleID")
agg$mean_logit <- agg$intercept + agg$sum_contribution
agg <- agg[order(match(agg$group, c("SM","DD","WS","WW")), agg$SampleID, -agg$contribution), ]
write.csv(agg, file.path(OUT, "consensus30_oof_contributions_per_specimen_protein.csv"), row.names = FALSE)

# ---- summaries ----
top_n <- function(s, n = 5) { a <- agg[agg$SampleID == s, ]; a <- a[order(-a$contribution), ]
  paste(sprintf("%s %+.2f", a$gene[1:n], a$contribution[1:n]), collapse = "; ") }
bottom_n <- function(s, n = 3) { a <- agg[agg$SampleID == s, ]; a <- a[order(a$contribution), ]
  paste(sprintf("%s %+.2f", a$gene[1:n], a$contribution[1:n]), collapse = "; ") }
focus <- unique(c(agg$SampleID[agg$truth == 1L],
                  head(pred$SampleID[order(-pred$probability)][pred$truth[order(-pred$probability)] == 0], 6)))
summ <- data.frame(SampleID = focus, group = group[match(focus, rownames(X))], family = family[match(focus, rownames(X))],
  probability = round(pred$probability[match(focus, pred$SampleID)], 3),
  threshold = round(pred$threshold[match(focus, pred$SampleID)], 3),
  call = ifelse(pred$predicted_class[match(focus, pred$SampleID)] == 1L, "+", "-"),
  mean_logit = round(sapply(focus, function(s) agg$mean_logit[agg$SampleID == s][1]), 2),
  intercept = round(sapply(focus, function(s) agg$intercept[agg$SampleID == s][1]), 2),
  top5_positive_contributions = sapply(focus, top_n), top3_negative_contributions = sapply(focus, bottom_n))
write.csv(summ, file.path(OUT, "summary_positive_and_top_FP_specimens.csv"), row.names = FALSE)

four <- c(PRL = "P01236", MTTP = "P55157", ARHGEF10 = "O15013", SEC31A = "O94979")
f4 <- agg[agg$FeatureID %in% four & agg$truth == 1L, c("SampleID", "family", "gene", "z", "beta", "contribution", "n_members_nonzero", "probability", "called_positive")]
f4 <- f4[order(f4$SampleID, match(f4$gene, names(four))), ]
write.csv(f4, file.path(OUT, "four_proteins_contributions_SM.csv"), row.names = FALSE)
# share of the four in the positive part of each SM specimen's logit
share <- sapply(unique(f4$SampleID), function(s) { a <- agg[agg$SampleID == s, ]
  c(pos_total = sum(a$contribution[a$contribution > 0]), four_sum = sum(a$contribution[a$FeatureID %in% four])) })
share <- data.frame(SampleID = colnames(share), t(round(share, 2)))
share$four_share_of_positive_logit <- round(share$four_sum / share$pos_total, 2)
write.csv(share, file.path(OUT, "four_proteins_share_SM.csv"), row.names = FALSE)

# protein-level: mean |contribution| across 44, and across the 7 SM
pl <- aggregate(cbind(abs_c = abs(contribution)) ~ FeatureID + gene, data = agg, FUN = mean)
pl$mean_contrib_SM <- aggregate(contribution ~ FeatureID, data = agg[agg$truth == 1L, ], FUN = mean)$contribution[match(pl$FeatureID, aggregate(contribution ~ FeatureID, data = agg[agg$truth == 1L, ], FUN = mean)$FeatureID)]
pl$mean_contrib_ctrl <- aggregate(contribution ~ FeatureID, data = agg[agg$truth == 0L, ], FUN = mean)$contribution[match(pl$FeatureID, aggregate(contribution ~ FeatureID, data = agg[agg$truth == 0L, ], FUN = mean)$FeatureID)]
pl$n_SM_with_contrib_gt_0.5 <- sapply(pl$FeatureID, function(f) sum(agg$contribution[agg$FeatureID == f & agg$truth == 1L] > 0.5))
pl$n_folds_in_panel <- sapply(pl$FeatureID, function(f) length(unique(agg$SampleID[agg$FeatureID == f])))
pl$n_SM_in_panel <- sapply(pl$FeatureID, function(f) length(unique(agg$SampleID[agg$FeatureID == f & agg$truth == 1L])))
pl <- pl[order(-pl$n_folds_in_panel, -pl$abs_c), ]
write.csv(pl, file.path(OUT, "protein_level_contribution_summary.csv"), row.names = FALSE)

# ---- figure: heat map of mean contributions, 7 SM + 6 top FP controls x proteins with largest mean |contribution| ----
# proteins selected in (almost) every outer panel, ranked by mean |contribution|; the four discussed proteins always shown
core <- pl[pl$n_folds_in_panel >= 40L, ]; core <- core[order(-core$abs_c), ]
top_prot <- unique(c(head(core$FeatureID, 14), unname(four), "Q9H0X9"))  # + OSBPL5 (selected only when family 14 or 17 is held out)
hm <- agg[agg$SampleID %in% focus & agg$FeatureID %in% top_prot, ]
hm$SampleID <- factor(hm$SampleID, levels = focus)
hm$gene <- factor(hm$gene, levels = rev(pl$gene[match(top_prot, pl$FeatureID)]))
hm$gene <- factor(ifelse(hm$gene == "OSBPL5", "OSBPL5*", as.character(hm$gene)), levels = c(setdiff(levels(hm$gene), "OSBPL5"), "OSBPL5*")[match(c(rev(pl$gene[match(top_prot, pl$FeatureID)])), c(setdiff(levels(hm$gene), "OSBPL5"), "OSBPL5"))])
lab <- sprintf("%s\np=%.2f %s", focus, pred$probability[match(focus, pred$SampleID)],
               ifelse(pred$predicted_class[match(focus, pred$SampleID)] == 1L, "(+)", "(-)"))
lim <- max(abs(hm$contribution))
p <- ggplot(hm, aes(x = SampleID, y = gene, fill = contribution)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(abs(contribution) >= 0.3, sprintf("%+.1f", contribution), "")), size = 2.6) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "#F5F5F5", high = "#D7301F", limits = c(-lim, lim),
                       name = "Mean contribution to logit\n(10 members, outer OOF)") +
  scale_x_discrete(labels = lab) +
  geom_vline(xintercept = 7.5, linewidth = 0.6, colour = "#17212B") +
  labs(x = "LoF/Mis specimens (left of line) and the six highest-scoring controls", y = NULL,
       title = "consensus30: which proteins produced each outer-OOF score",
       caption = "Blank cells: protein not in that fold's consensus30 panel. *OSBPL5 entered the panel only when family 14 or 17 was held out.",
       subtitle = "Contribution = coefficient x standardised value, averaged over the 10 member models of the held-out fold") +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_text(size = 7.2, lineheight = 0.9), panel.grid = element_blank(),
        legend.position = "right", plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "FigureS_consensus30_contribution_heatmap.pdf"), p, width = 9.5, height = 5.4, useDingbats = FALSE)
ggsave(file.path(OUT, "FigureS_consensus30_contribution_heatmap.png"), p, width = 9.5, height = 5.4, dpi = 300, bg = "white")

cat("\n[SM + top FP summary]\n"); print(summ[, c("SampleID","family","probability","call","mean_logit","intercept","top5_positive_contributions")], row.names = FALSE)
cat("\n[four proteins, SM]\n"); print(f4, row.names = FALSE); print(share, row.names = FALSE)
cat("\n[protein level, proteins in >=40 specimens' panels]\n"); print(pl[pl$n_folds_in_panel >= 40L, ], row.names = FALSE)
cat("\n[fold-specific proteins with large single contributions]\n"); print(head(pl[pl$n_folds_in_panel < 10L, ][order(-pl$abs_c[pl$n_folds_in_panel < 10L]), ], 8), row.names = FALSE)
