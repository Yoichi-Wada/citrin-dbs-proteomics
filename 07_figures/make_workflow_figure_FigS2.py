#!/usr/bin/env python3
"""Supplementary Fig. S2 — analytical workflow (manuscript version, v4 2026-11-09: values after the WW_1–3 correction, final table numbers).
Derived from v3 (2026-11-06).
v2 incorporates the code-verified corrections of the 2026-11-06 review: specificity-constrained threshold rule,
adaptive-strategy selection in the inner loop, scope of the nested procedure, held-out wording, split of sensitivity
analyses, internal-evaluation heading, Stage 1 operational definition, final-model specification, neutral control labels.
Verified against citrin_nested_ensemble_familyLOFO_v1_1.R (spec90_threshold, adaptive choice, consensus ordering) and
citrin_44_nested_fixed_panel_v2_2.R (stability subsampling: >=3 positives per subsample)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results"); os.makedirs(OUT, exist_ok=True)
plt.rcParams["font.family"] = "Arial"
plt.rcParams["mathtext.fontset"] = "custom"; plt.rcParams["mathtext.rm"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

C = {
    "data":   dict(fc="#ECECEC", ec="#7F7F7F"),   # label-free data handling / study design
    "train":  dict(fc="#DDEEDA", ec="#3C7A3A"),   # estimated from the outer training partition only
    "held":   dict(fc="#FFF4D6", ec="#B8860B"),   # held-out family
    "report": dict(fc="#DCE6F1", ec="#376092"),   # evaluation / reported output
    "final":  dict(fc="#EFE3F3", ec="#7B4F9D"),   # frozen final model
}
FS, FS_T, LH = 8.0, 8.6, 0.19
W, H = 13.0, 12.0
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")


def box(x, y_top, w, title, lines, kind, h_pad=0.14):
    n = len(lines) + (1 if title else 0)
    h = n * LH + 2 * h_pad
    ax.add_patch(FancyBboxPatch((x, y_top - h), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                linewidth=1.1, **C[kind], zorder=2))
    yy = y_top - h_pad - LH * 0.78
    if title:
        ax.text(x + 0.12, yy, title, fontsize=FS_T, weight="bold", va="center", ha="left", zorder=3); yy -= LH
    for ln in lines:
        bold = ln.startswith("**")
        ax.text(x + 0.12, yy, ln.strip("*"), fontsize=FS, va="center", ha="left", zorder=3,
                weight="bold" if bold else "normal"); yy -= LH
    return y_top - h


def varrow(x, y1, y2, label=None):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>", mutation_scale=11, linewidth=1.0, color="#404040", zorder=1))
    if label:
        ax.text(x + 0.1, (y1 + y2) / 2, label, fontsize=FS - 0.4, va="center", ha="left", color="#404040", style="italic")


def harrow(x1, x2, y):
    ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>", mutation_scale=11, linewidth=1.0, color="#404040", zorder=1))


def frame(x, y_top, w, h, label, ec, fc="none", lw=1.4, ls="-", z=0.5):
    ax.add_patch(FancyBboxPatch((x, y_top - h), w, h, boxstyle="round,pad=0.02,rounding_size=0.1",
                                linewidth=lw, edgecolor=ec, facecolor=fc, linestyle=ls, zorder=z))
    ax.text(x + 0.12, y_top - 0.15, label, fontsize=FS_T, weight="bold", color=ec, va="center", ha="left", zorder=3)


GAP = 0.24
ax.text(0.3, H - 0.3, "Analytical workflow for SLC25A13 detection and nested family-grouped evaluation of an exploratory LoF/Mis classifier",
        fontsize=10.2, weight="bold", va="center")

# ================================================================== LEFT COLUMN
xl, wl = 0.3, 5.75
y = H - 0.72
yb = box(xl, y, wl, "Specimens and proteomic data", [
    "60 DBS: 16 LoF/LoF; 7 LoF/Mis (5 families, 2 sibling pairs); 10 carriers (parents of LoF/LoF",
    "patients); 6 healthy adults; 21 disease controls; 3 analytical runs",
    "DIA-NN 2.2.0 (1 % FDR, cross-run normalization): 7,474 protein groups → 5,199 with ≥ 3 unique",
    "peptides (technical filter, applied once, not cross-validated)",
    "Zero = not quantified; quantified values log2-transformed and ACTB-centered per specimen (Stage 2)"], "data")
varrow(xl + wl / 2, yb, yb - GAP); y = yb - GAP

yb = box(xl, y, wl, "Stage 1 — Direct SLC25A13 detection  (Fig. 2)", [
    "Rule: no SLC25A13-assigned precursor passing DIA-NN identification in the raw, non-imputed",
    "data (peptides shared with SLC25A12 not counted)",
    "Observed: not detected in 16/16 LoF/LoF; quantified in 44/44 others — a cohort-level",
    "observation, not an independently validated rule",
    "Carriers 0.35×, LoF/Mis 0.28× the median of the 27 non-carrier controls (healthy + disease)"], "data")
varrow(xl + wl / 2, yb, yb - GAP, "SLC25A13 quantified  →  Stage 2 analysis set"); y = yb - GAP

yb = box(xl, y, wl, "Stage 2 — Exploratory multiprotein classification  (Fig. 3)", [
    "44 specimens, 35 families: 7 LoF/Mis vs 37 controls (carriers, healthy, disease controls);",
    "Run 2 n = 9 (controls only), Run 3 n = 35 (all 7 LoF/Mis)",
    "Classifier preprocessing, feature selection, fitting and threshold calibration restricted to",
    "each outer training partition (right); panel size 30 fixed, inherited from same-cohort exploration"], "train")
y_s2_top, y_s2_bot = y, yb
varrow(xl + wl / 2, yb, yb - GAP); y = yb - GAP

yb = box(xl, y, wl, "Internal out-of-fold evaluation of the model-development procedure (44 specimens, 35 families)", [
    "Computed once from the 44 held-out predictions (Table 2, Fig. 3): AUC, sensitivity, specificity,",
    "balanced accuracy, Brier score, positive families detected",
    "**Primary consensus-30 ensemble: AUC 0.834 (95 % CI 0.602–0.971), sensitivity 7/7,",
    "**specificity 31/37, 5/5 families detected",
    "CI: family-cluster bootstrap of saved predictions (no refitting) or exact binomial (specimen-level)",
    "Two-stage evaluation of all 60 (Stage 1 observed, Stage 2 out-of-fold): retrospective, neither",
    "stage independently validated (Table S16)"], "report")
varrow(xl + wl / 2, yb, yb - GAP); y = yb - GAP

yb = box(xl, y, wl, "Sensitivity and reference analyses (Tables S10, S13–S15)", [
    "**Nested refitting (single models, 10 seeds):  panel size chosen in the inner loop (5–50);",
    "training-only mean-only batch adjustment; Run 3 only (n = 35); family-level subsampling;",
    "ensemble re-run after excluding PRL, MTTP, ARHGEF10 and SEC31A, individually and jointly",
    "**From saved out-of-fold predictions:  LoF/Mis vs carriers / disease controls /",
    "age-restricted controls; family-level performance; selection recurrence over 35 panels",
    "**Reference:  100 family-structured label permutations of the complete score-generation",
    "procedure of the primary ensemble (inner threshold calibration omitted); P = (b + 1)/(N + 1)"], "report")
varrow(xl + wl / 2, yb, yb - GAP); y = yb - GAP

yb = box(xl, y, wl, "Frozen research ensemble — trained once on all 44 development specimens", [
    "One 30-protein consensus panel; 10 coefficient vectors + intercepts; member-specific",
    "preprocessing constants and imputation seeds; training-derived thresholds 0.034 / 0.210",
    "Not independently validated; not evaluated on the 44 development specimens;",
    "application to the 16 LoF/LoF (parents of six among the 44) is descriptive",
    "Specification: Table S17 (machine-readable with the analysis code); LoF/LoF application: Table S18"], "final")
y_left_bottom = yb

# ================================================================== RIGHT COLUMN
xr, wr = 6.35, 6.55
yr = H - 0.72
ax.text(xr, yr + 0.02, "Nested family-held-out cross-validation  (one of 35 outer iterations)", fontsize=9.0, weight="bold", va="bottom")
outer_top = yr - 0.16
pad = 0.16
xi, wi = xr + pad, wr - 2 * pad
y = outer_top - 0.38
yb = box(xi, y, wi, "Outer training partition  (all families except the held-out family)", [
    "(i)  MinProb imputation (q = 0.01) from training specimens only",
    "(ii)  no batch correction beyond DIA-NN cross-run normalization",
    "(iii)  10 members (fixed seeds): stability selection on 30 class-stratified 80 % subsamples",
    "       (≥ 3 positives each), 1,000 highest-variance proteins per subsample, class-weighted",
    "       elastic net (α = 0.5)  →  selection-frequency ranking",
    "(iv)  panel size fixed at 30",
    "(v)  consensus panel = 30 proteins most often in members' top 30 (ties: mean frequency);",
    "       each member: standardize (training μ, σ) → elastic net, λ by family-grouped 3-fold CV",
    "(vi)  ensemble probability = mean of 10 member probabilities"], "train")
varrow(xi + wi / 2, yb, yb - GAP); y = yb - GAP

inner_top = y
y_in = inner_top - 0.38
xii, wii = xi + 0.26, wi - 0.52
yb = box(xii, y_in, wii, "Inner training set (all other families of the outer partition): steps (i)–(vi) → inner model", [], "train")
varrow(xii + wii / 2, yb, yb - GAP); y_in = yb - GAP
yb = box(xii, y_in, wii, "Inner held-out family  →  inner out-of-fold scores", []," train".strip())
varrow(xii + wii / 2, yb, yb - GAP); y_in = yb - GAP
yb = box(xii, y_in, wii, "(vii)  Calibration on the inner out-of-fold scores", [
    "primary threshold: minimum inner score among positive specimens",
    "specificity-constrained threshold: max. inner sensitivity s.t. specificity ≥ 0.90 (ties: fewer FP, lower τ)",
    "adaptive strategy: consensus or diverse panel, whichever has fewer inner FP at the",
    "   primary threshold (ties: inner AUC, then consensus)"], "train")
inner_bottom = yb - 0.2
frame(xi, inner_top, wi, inner_top - inner_bottom,
      "Inner family-LOFO calibration (34 folds; within the outer training partition): thresholds and adaptive strategy",
      ec="#2F6DA8", fc="#E6F0FA", lw=2.0, ls=(0, (5, 2.5)), z=0.8)
varrow(xi + wi / 2, inner_bottom, inner_bottom - GAP); y = inner_bottom - GAP

yb = box(xi, y, wi, "Held-out family  (excluded from all fold-specific fitting and calibration)", [
    "impute from the specimen's own 1 % quantile with training σ → standardize (training μ, σ)",
    "→ mean of 10 member probabilities → call at the fold's threshold = one out-of-fold prediction"], "held")
outer_bottom = yb - pad
frame(xr, outer_top, wr, outer_top - outer_bottom, "Outer loop  (35 iterations, one per family)", ec="#3C7A3A", fc="#FAFCF9")

yb2 = box(xr, outer_bottom - GAP, wr, "Ensemble strategies  (fixed before the final run)", [
    "primary: consensus panel;  alternative \"diverse\": each member keeps its own top 30;",
    "alternative \"adaptive\": chosen per fold in the inner loop (vii)"], "data")

harrow(xl + wl, xr - 0.05, (y_s2_top + y_s2_bot) / 2)

# ================================================================== legend
ly = min(y_left_bottom, yb2) - 0.42
items = [("data", "label-free data handling / study design"), ("train", "estimated from the outer training partition only"),
         ("held", "held-out family"), ("report", "evaluation / reported output"), ("final", "frozen final model")]
lx = 0.3
for kind, lab in items:
    ax.add_patch(FancyBboxPatch((lx, ly - 0.09), 0.32, 0.18, boxstyle="round,pad=0.01,rounding_size=0.03", linewidth=1.0, **C[kind]))
    ax.text(lx + 0.42, ly, lab, fontsize=FS, va="center", ha="left"); lx += 2.6 if kind != "train" else 3.1
ax.text(0.3, ly - 0.3, "LoF, loss-of-function; Mis, missense; LOFO, leave-one-family-out; CV, cross-validation; FP, false positives; MinProb, left-censored imputation at the 1 % quantile.",
        fontsize=FS - 0.4, va="center", color="#404040")

y_floor = ly - 0.6
ax.set_ylim(y_floor, H); fig.set_size_inches(W, H - y_floor)
for ext in ("pdf", "png", "svg"):
    fig.savefig(os.path.join(OUT, f"FigureS2_analytical_workflow.{ext}"), dpi=220 if ext == "png" else None, bbox_inches="tight", pad_inches=0.12)
print("saved; left bottom", round(y_left_bottom, 2), "right bottom", round(yb2, 2))
