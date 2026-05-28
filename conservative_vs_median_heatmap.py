"""
plot_performance_heatmap.py
----------------------------
Generates a two-panel heatmap comparing performance levels of single biomarkers
under two evaluation schemes:
  - Left panel  : Worst-case (conservative) — from summary.csv
  - Right panel : Median-based              — from summary_median.csv

Acceptance rule: a biomarker PASSES only if ALL 8 metrics are at least
"Minimally acceptable". Any single "Not acceptable" = FAIL.

Usage
-----
  python plot_performance_heatmap.py \
      --summary     summary.csv \
      --summary_med summary_median.csv \
      --out         heatmap_comparison.pdf

Dependencies: pandas, matplotlib, numpy
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm

# ──────────────────────────────────────────────────────────────────────────────
# 1.  Threshold definitions
# ──────────────────────────────────────────────────────────────────────────────
# Performance levels: 0=Not acceptable, 1=Minimally acceptable, 2=Good, 3=Excellent

def classify_higher_better(value, thresholds=(0.7, 0.8, 0.9)):
    """AUC1, AUC2, Pairwise — higher is better.
    EPS guards against floating-point boundary issues.
    """
    if pd.isna(value):
        return np.nan
    EPS = 1e-4
    t1, t2, t3 = thresholds
    if value < t1 - EPS:   return 0
    elif value < t2 - EPS: return 1
    elif value < t3 - EPS: return 2
    else:                  return 3

def classify_lr_pos(value, thresholds=(2.0, 5.0, 10.0)):
    """LR_pos_th1, LR_pos_th2 — higher is better.
    EPS guards against floating-point boundary issues.
    """
    if pd.isna(value):
        return np.nan
    EPS = 1e-4
    t1, t2, t3 = thresholds
    if value < t1 - EPS:   return 0
    elif value < t2 - EPS: return 1
    elif value < t3 - EPS: return 2
    else:                  return 3

def classify_lr_neg(value, thresholds=(0.5, 0.2, 0.1)):
    """LR_neg_th1, LR_neg_th2 — lower is better.
    EPS guards against floating-point values marginally above a boundary
    (e.g. 0.5000005 from bootstrap averaging) being misclassified.
    """
    if pd.isna(value):
        return np.nan
    EPS = 1e-4
    t1, t2, t3 = thresholds          # t1 > t2 > t3
    if value > t1 + EPS:   return 0
    elif value > t2 + EPS: return 1
    elif value > t3 + EPS: return 2
    else:                  return 3

def classify_cl_error(value, thresholds=(1.0, 0.5, 0.3)):
    """Mean_classification_error — lower is better.
    EPS guards against floating-point boundary issues.
    """
    if pd.isna(value):
        return np.nan
    EPS = 1e-4
    t1, t2, t3 = thresholds          # t1 > t2 > t3
    if value > t1 + EPS:   return 0
    elif value > t2 + EPS: return 1
    elif value > t3 + EPS: return 2
    else:                  return 3

# Map column name → classifier function
METRIC_CLASSIFIERS = {
    "AUC1":                 classify_higher_better,
    "AUC2":                 classify_higher_better,
    "pairwise":             classify_higher_better,
    "LR_pos_th1":           classify_lr_pos,
    "LR_pos_th2":           classify_lr_pos,
    "LR_neg_th1":           classify_lr_neg,
    "LR_neg_th2":           classify_lr_neg,
    "Mean_classification_error": classify_cl_error,
}

# Display labels for columns (x-axis)
METRIC_LABELS = {
    "AUC1":                 "AUC₁",
    "AUC2":                 "AUC₂",
    "pairwise":             "Pairwise",
    "LR_pos_th1":           "LR⁺ th₁",
    "LR_pos_th2":           "LR⁺ th₂",
    "LR_neg_th1":           "LR⁻ th₁",
    "LR_neg_th2":           "LR⁻ th₂",
    "Mean_classification_error": "Cl. error",
}

METRIC_COLS   = list(METRIC_CLASSIFIERS.keys())
METRIC_XLABELS = [METRIC_LABELS[m] for m in METRIC_COLS]

# ──────────────────────────────────────────────────────────────────────────────
# 2.  Colour scheme
# ──────────────────────────────────────────────────────────────────────────────
LEVEL_COLORS = [
    "#d62728",   # 0 = Not acceptable       (red)
    "#ff7f0e",   # 1 = Minimally acceptable (orange)
    "#2ca02c",   # 2 = Good                 (green)
    "#1f77b4",   # 3 = Excellent            (blue)
]
CMAP = ListedColormap(LEVEL_COLORS)
NORM = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], CMAP.N)

LEVEL_NAMES = ["Not acceptable", "Minimally acceptable", "Good", "Excellent"]

# ──────────────────────────────────────────────────────────────────────────────
# 3.  Core processing
# ──────────────────────────────────────────────────────────────────────────────
def detect_biomarker_col(df):
    """Return the name of the biomarker identifier column."""
    for c in ["feature_1", "feature_2", "Feature_Pair"]:
        if c in df.columns:
            return c
    # Fallback: first string column
    for c in df.columns:
        if df[c].dtype == object:
            return c
    return None


def build_level_matrix(df):
    """
    Returns
    -------
    level_mat : np.ndarray  shape (n_biomarkers, n_metrics)   float, NaN possible
    pass_mask : np.ndarray  shape (n_biomarkers,)              bool
    biomarkers : list of str
    """
    bio_col    = detect_biomarker_col(df)
    biomarkers = df[bio_col].tolist() if bio_col else [f"Biomarker_{i+1}" for i in range(len(df))]
    n = len(biomarkers)
    m = len(METRIC_COLS)
    level_mat = np.full((n, m), np.nan)

    for i, col in enumerate(METRIC_COLS):
        if col not in df.columns:
            continue
        for j, val in enumerate(df[col]):
            level_mat[j, i] = METRIC_CLASSIFIERS[col](val)

    # A biomarker passes iff every non-NaN metric is ≥ 1 (Minimally acceptable)
    pass_mask = np.array([
        np.all(level_mat[j, ~np.isnan(level_mat[j, :])] >= 1)
        for j in range(n)
    ])
    return level_mat, pass_mask, biomarkers


def annotate_values(ax, df, biomarkers):
    """Write numeric value in each cell."""
    bio_col = detect_biomarker_col(df)
    for i, col in enumerate(METRIC_COLS):
        if col not in df.columns:
            continue
        for j, bio in enumerate(biomarkers):
            val = df.loc[df[bio_col] == bio, col].values if bio_col else df[[col]].iloc[[j]].values.flatten()
            if len(val) == 0 or pd.isna(val[0]):
                txt = "NA"
            else:
                txt = f"{val[0]:.2f}"
            ax.text(i, j, txt, ha="center", va="center",
                    fontsize=6.5, color="white", fontweight="bold")


# ──────────────────────────────────────────────────────────────────────────────
# 4.  Plotting
# ──────────────────────────────────────────────────────────────────────────────
def draw_panel(ax, level_mat, pass_mask, biomarkers, title, show_yticklabels=True):
    n_bio, n_met = level_mat.shape

    im = ax.imshow(level_mat, cmap=CMAP, norm=NORM,
                   aspect="auto", interpolation="none")

    # Grid lines
    ax.set_xticks(np.arange(-0.5, n_met, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_bio, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)

    # x-axis
    ax.set_xticks(range(n_met))
    ax.set_xticklabels(METRIC_XLABELS, rotation=40, ha="right", fontsize=9)

    # y-axis (biomarker names; color by pass/fail)
    ax.set_yticks(range(n_bio))
    if show_yticklabels:
        labels = ax.set_yticklabels(biomarkers, fontsize=9)
        for label, passed in zip(labels, pass_mask):
            label.set_color("black" if passed else "#d62728")
            label.set_fontweight("bold" if passed else "normal")
    else:
        ax.set_yticklabels([])

    # Highlight rows that PASS with a left bracket
    if show_yticklabels:
        for j, passed in enumerate(pass_mask):
            if passed:
                ax.annotate("", xy=(-0.52, j), xycoords=("axes fraction", "data"),
                            xytext=(-0.52, j),
                            annotation_clip=False)

    # Title
    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)

    return im


def make_heatmap(df_wc, df_med, out_path):
    mat_wc,  pass_wc,  bios_wc  = build_level_matrix(df_wc)
    mat_med, pass_med, bios_med = build_level_matrix(df_med)

    # Align biomarkers:
    # If summary_median has no identifier column, inherit names from summary.csv
    if any(b.startswith("Biomarker_") for b in bios_med):
        bios_med = bios_wc
    assert len(bios_wc) == len(bios_med), (
        f"Row count mismatch: summary has {len(bios_wc)}, "
        f"summary_median has {len(bios_med)} rows."
    )
    n_bio = len(bios_wc)

    # ── Reorder rows: PASS (median) first, then FAIL ──────────────────────
    # Sorting key: (0=PASS, 1=FAIL) so PASS biomarkers float to top.
    # Within each group, original CSV order is preserved.
    sort_idx   = sorted(range(n_bio), key=lambda i: (0 if pass_med[i] else 1))
    mat_wc     = mat_wc[sort_idx, :]
    mat_med    = mat_med[sort_idx, :]
    pass_wc    = pass_wc[sort_idx]
    pass_med   = pass_med[sort_idx]
    bios_wc    = [bios_wc[i]  for i in sort_idx]
    bios_med   = [bios_med[i] for i in sort_idx]
    # Reorder df rows so annotate_values stays in sync
    bio_col_wc  = detect_biomarker_col(df_wc)
    bio_col_med = detect_biomarker_col(df_med)
    df_wc  = df_wc.iloc[sort_idx].reset_index(drop=True)
    df_med = df_med.iloc[sort_idx].reset_index(drop=True)

    # ── figure layout ──────────────────────────────────────────────────────
    fig_w = 15.0
    fig_h = max(5.0, 0.45 * n_bio + 2.5)
    fig, axes = plt.subplots(
        1, 2,
        figsize=(fig_w, fig_h),
        gridspec_kw={"wspace": 0.2, "width_ratios": [1, 1]}
    )

    im = draw_panel(axes[0], mat_wc,  pass_wc,  bios_wc,
                    "Conservative-based", show_yticklabels=True)
    draw_panel(axes[1], mat_med, pass_med, bios_med,
               "Median-based", show_yticklabels=True)

    # Numeric annotations
    annotate_values(axes[0], df_wc,  bios_wc)
    annotate_values(axes[1], df_med, bios_med)

    # ── pass / fail summary text ───────────────────────────────────────────
    n_pass_wc  = int(pass_wc.sum())
    n_pass_med = int(pass_med.sum())
    for ax, n_pass, scheme in zip(axes,
                                   [n_pass_wc, n_pass_med],
                                   ["Worst-case", "Median"]):
        ax.set_xlabel(
            f"Passed: {n_pass}/{n_bio} biomarkers",
            fontsize=9, labelpad=6,
            color="#2ca02c" if n_pass > 0 else "#d62728"
        )

    # ── shared legend ──────────────────────────────────────────────────────
    patches = [
        mpatches.Patch(color=LEVEL_COLORS[i], label=LEVEL_NAMES[i])
        for i in range(4)
    ]
    # "Pass" row indicator
    pass_patch = mpatches.Patch(
        facecolor="none", edgecolor="black", linewidth=0,
        label="Bold black name = biomarker passes all criteria"
    )
    fail_patch = mpatches.Patch(
        facecolor="none", edgecolor="#d62728", linewidth=0,
        label="Red name = at least one metric rejected"
    )

    fig.legend(
        handles=patches + [pass_patch, fail_patch],
        loc="lower center",
        ncol=3,
        fontsize=9,
        frameon=True,
        bbox_to_anchor=(0.5, -0.06),
        title="Performance level",
        title_fontsize=10
    )

    fig.suptitle(
        "Single-biomarker TdP risk prediction: performance level comparison\n"
        "Conservative-based vs. Median-based evaluation scheme",
        fontsize=11, fontweight="bold", y=0.97
    )

    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved → {out_path}")

    # ── console summary ───────────────────────────────────────────────────
    print("\n── Acceptance summary ──────────────────────────────────────")
    print(f"{'Biomarker':<22} {'Worst-case':>12} {'Median':>8}")
    print("-" * 45)
    for bio, pw, pm in zip(bios_wc, pass_wc, pass_med):
        wc_str  = "PASS" if pw  else "FAIL"
        med_str = "PASS" if pm  else "FAIL"
        flag = " ←" if (not pw and pm) else ""  # gained under median
        print(f"{bio:<22} {wc_str:>12} {med_str:>8}{flag}")
    print("-" * 45)
    print(f"{'TOTAL':<22} {n_pass_wc:>12} {n_pass_med:>8}")
    gained = int((~pass_wc & pass_med).sum())
    lost   = int((pass_wc & ~pass_med).sum())
    print(f"\nGained under median scheme : {gained} biomarker(s)")
    print(f"Lost   under median scheme : {lost} biomarker(s)")


# ──────────────────────────────────────────────────────────────────────────────
# 5.  Entry point
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Performance heatmap: worst-case vs median")
    parser.add_argument("--summary",     required=True, help="Path to summary.csv (worst-case)")
    parser.add_argument("--summary_med", required=True, help="Path to summary_median.csv (median)")
    parser.add_argument("--out",         default="heatmap_comparison.pdf",
                        help="Output file path (.pdf / .png / .svg)")
    args = parser.parse_args()

    df_wc  = pd.read_csv(args.summary)
    df_med = pd.read_csv(args.summary_med)

    make_heatmap(df_wc, df_med, args.out)