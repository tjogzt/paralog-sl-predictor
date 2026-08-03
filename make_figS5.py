#!/usr/bin/env python3
"""make_figS5.py — regenerate Supplementary Fig. S5 (Mutation Type Stratification)
from the CURRENT artifact output/muttype_all_results.csv.

Replaces the stale figure (panel-d legend read All=0.735 / Missense=0.609 /
Truncating=0.726) with the audited values computed from the signed DD columns
with NaN->0 (roc_auc_score(is_known_sl, col.fillna(0))):
  Breast: All=0.465, Missense=0.437, Truncating=0.460.

Layout/style matches the existing figure (R_figS5_muttype.R): 2x2 composite,
180x180 mm, Arial 7 pt base, white background, bold panel letters a-d.

Outputs (overwritten):
  output/figures/FigS5_MutationType.pdf / .svg / .tiff (300 dpi) / .png (300 dpi)
  output/figures/FigS5_panel_{a,b,c,d}.pdf
  ../figure_review/FigS5_MutationType.png (copy)

Usage: python3 make_figS5.py   (run from repo root)
"""
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42  # TrueType (journals reject Type 3)
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.size"] = 7
matplotlib.rcParams["axes.labelsize"] = 7
matplotlib.rcParams["axes.linewidth"] = 0.4
matplotlib.rcParams["xtick.labelsize"] = 7
matplotlib.rcParams["ytick.labelsize"] = 7
matplotlib.rcParams["xtick.major.width"] = 0.3
matplotlib.rcParams["ytick.major.width"] = 0.3
matplotlib.rcParams["legend.fontsize"] = 7
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

ROOT = Path(__file__).resolve().parent
CSV = ROOT / "output" / "muttype_all_results.csv"
OUT = ROOT / "output" / "figures"
REVIEW = ROOT.parent / "figure_review"

RED = "#CB181D"; BLUE = "#2171B5"; GRAY = "#636363"
DARK = "#252525"; PURPLE = "#8E44AD"

MM = 1 / 25.4
FIGSIZE = (180 * MM, 180 * MM)      # 180 x 180 mm composite (~7.09 in)
PANELSIZE = (90 * MM, 90 * MM)      # 90 x 90 mm single panels

CANCERS = ["Ovarian", "Endometrial", "Colorectal", "Breast"]
# Audited expectations (verified from the current artifact) — compared, never
# hardcoded into the figure.
EXPECTED_PANEL_A = {
    "Ovarian":    {"dd_all": 0.567, "dd_trunc": 0.633, "dd_miss": 0.471},
    "Endometrial": {"dd_all": 0.463, "dd_trunc": 0.600, "dd_miss": 0.365},
    "Colorectal": {"dd_all": 0.793, "dd_trunc": 0.857, "dd_miss": 0.686},
    "Breast":     {"dd_all": 0.465, "dd_trunc": 0.460, "dd_miss": 0.437},
}
EXPECTED_PANEL_D_LABELS = ["All (AUC=0.465)", "Missense (AUC=0.437)",
                           "Truncating (AUC=0.460)"]

# ── Data ────────────────────────────────────────────────────────────────────
if not CSV.exists():
    sys.exit(f"ERROR: {CSV} not found — run mutation_type_analysis.py first")
mut = pd.read_csv(CSV)
mut["is_known_sl"] = mut["is_known_sl"].astype(bool)

# Panel a: per-cancer AUROCs on the signed DD columns (NaN -> 0)
panel_a = {c: {} for c in CANCERS}
print("=== Panel a: per-cancer DD AUROC (signed, NaN->0) ===")
for canc in CANCERS:
    sub = mut[mut["cancer"] == canc]
    y = sub["is_known_sl"].astype(int)
    for col in ["dd_all", "dd_trunc", "dd_miss"]:
        panel_a[canc][col] = roc_auc_score(y, sub[col].fillna(0))
    exp = EXPECTED_PANEL_A[canc]
    got = {k: round(v, 3) for k, v in panel_a[canc].items()}
    ok = all(abs(panel_a[canc][k] - exp[k]) < 5e-4 for k in exp)
    print(f"  {canc:11s} All={panel_a[canc]['dd_all']:.3f} "
          f"Trunc={panel_a[canc]['dd_trunc']:.3f} "
          f"Miss={panel_a[canc]['dd_miss']:.3f}  (expected match: {ok})")

# Panel b: |dd_trunc| - |dd_miss| over rows with both non-NaN
both = mut[mut["dd_trunc"].notna() & mut["dd_miss"].notna()]
diff = (both["dd_trunc"].abs() - both["dd_miss"].abs()).values
mean_diff = float(np.mean(diff))
mean_label = f"Mean = {mean_diff:+.3f}"
print(f"=== Panel b: n={len(diff)}, mean={mean_diff:.6f} -> label '{mean_label}'")

# Panel c: known SL pairs, max |DD| per pair across cancers, current-figure
# ordering (max|dd_trunc| descending, top 8)
known = (mut[mut["is_known_sl"]]
         .groupby(["driver", "paralog"], as_index=False)
         .agg(dd_trunc=("dd_trunc", lambda s: float(s.abs().max())),
              dd_miss=("dd_miss", lambda s: float(s.abs().max()))))
known = (known.sort_values("dd_trunc", ascending=False, kind="stable")
         .head(8).reset_index(drop=True))
known["label"] = known["driver"] + "\u2192" + known["paralog"]
print("=== Panel c: known SL pairs (max|dd_trunc| desc, top 8) ===")
for _, r in known.iterrows():
    print(f"  {r['label']:22s} trunc={r['dd_trunc']:.3f} miss={r['dd_miss']:.3f}")

# Panel d: Breast ROCs from the signed columns (NaN -> 0)
br = mut[mut["cancer"] == "Breast"]
y_br = br["is_known_sl"].astype(int).values
roc_series = []  # (name, color, fpr, tpr, legend label)
for col, name, color in [("dd_all", "All", GRAY),
                         ("dd_miss", "Missense", BLUE),
                         ("dd_trunc", "Truncating", RED)]:
    scores = br[col].fillna(0).values
    auc = roc_auc_score(y_br, scores)
    fpr, tpr, _ = roc_curve(y_br, scores)
    roc_series.append((name, color, fpr, tpr, f"{name} (AUC={auc:.3f})"))
d_labels = [s[4] for s in roc_series]
print(f"=== Panel d: Breast ROC labels: {d_labels}")
if d_labels != EXPECTED_PANEL_D_LABELS:
    sys.exit("DISCREPANCY: panel-d legend labels computed from the artifact "
             f"are {d_labels}, expected {EXPECTED_PANEL_D_LABELS}. "
             "Stopping without writing files.")


# ── Panel builders (draw onto a given Axes) ──────────────────────────────────
def _despine(ax, keep_left=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if not keep_left:
        ax.spines["left"].set_visible(False)


def draw_panel_a(ax):
    x = np.arange(len(CANCERS))
    series = [("dd_all", "All", GRAY), ("dd_trunc", "Truncating", RED),
              ("dd_miss", "Missense", BLUE)]
    w = 0.55 / 3  # R: geom_col(width=0.55, position_dodge(0.7)) -> touching bars
    for i, (col, name, color) in enumerate(series):
        vals = [panel_a[c][col] for c in CANCERS]
        ax.bar(x + (i - 1) * w, vals, width=w, color=color, label=name)
    ax.axhline(0.5, linewidth=0.3, color=GRAY, linestyle=(0, (4, 3)), alpha=0.3)
    ax.set_xticks(x)
    ax.set_xticklabels(CANCERS)
    ax.set_ylim(0, 0.9)
    ax.set_yticks(np.arange(0, 0.81, 0.2))
    ax.set_ylabel("DD AUROC")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
              frameon=False, handlelength=1.2, handleheight=1.0,
              columnspacing=1.2)
    _despine(ax)


def draw_panel_b(ax):
    ax.hist(diff, bins=25, color=PURPLE, alpha=0.7, edgecolor="white",
            linewidth=0.2)
    ax.axvline(0, linewidth=0.5, color=DARK)
    ax.axvline(mean_diff, linewidth=0.5, color=RED, linestyle="dashed")
    lo, hi = diff.min(), diff.max()
    span = hi - lo
    ax.set_xlim(lo - 0.03 * span, hi + 0.14 * span)  # R: expand mult c(0.03, .14)
    ymax = ax.get_ylim()[1]
    ax.set_ylim(0, ymax * 1.14)
    ax.text(0.97, 0.94, mean_label, transform=ax.transAxes, ha="right",
            va="top", color=RED, fontsize=7)
    ax.set_xlabel("|DD_trunc| \u2212 |DD_miss|")
    ax.set_ylabel("Frequency")
    _despine(ax)


def draw_panel_c(ax):
    n = len(known)
    ys = {i: n - 1 - i for i in range(n)}  # first pair (EP300) on top
    h = 0.27
    for i, r in known.iterrows():
        y = ys[i]
        ax.barh(y + 0.175, r["dd_miss"], height=h, color=BLUE)    # Missense
        ax.barh(y - 0.175, r["dd_trunc"], height=h, color=RED)    # Truncating
        ax.text(0, y + 0.45, r["label"], fontsize=7, color=DARK,
                ha="left", va="center")
    ax.set_xlim(0, 0.62)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6])
    ax.set_ylim(-0.6, n - 1 + 0.95)
    ax.set_yticks([])
    ax.set_xlabel("|DD|")
    ax.legend(handles=[Patch(color=RED, label="Truncating"),
                       Patch(color=BLUE, label="Missense")],
              loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2,
              frameon=False, handlelength=1.2, handleheight=1.0,
              columnspacing=1.2)
    _despine(ax, keep_left=False)


def draw_panel_d(ax):
    ax.plot([0, 1], [0, 1], linewidth=0.3, color=GRAY, linestyle="dashed",
            alpha=0.5)
    for name, color, fpr, tpr, label in roc_series:
        ax.plot(fpr, tpr, linewidth=0.6, color=color, label=label)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Breast Cancer", fontsize=7, fontweight="bold", loc="left")
    ax.legend(loc="lower right", frameon=False, handlelength=1.6,
              borderpad=0.2, labelspacing=0.3)
    _despine(ax)


# ── Single-panel PDFs (90 x 90 mm, no panel letters) ────────────────────────
OUT.mkdir(parents=True, exist_ok=True)
for letter, draw in zip("abcd", [draw_panel_a, draw_panel_b,
                                 draw_panel_c, draw_panel_d]):
    fig, ax = plt.subplots(figsize=PANELSIZE)
    draw(ax)
    fig.savefig(OUT / f"FigS5_panel_{letter}.pdf",
                bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"  wrote output/figures/FigS5_panel_{letter}.pdf")

# ── Composite 2x2 (180 x 180 mm) with panel letters ─────────────────────────
fig, axes = plt.subplots(2, 2, figsize=FIGSIZE)
fig.subplots_adjust(left=0.085, right=0.99, top=0.955, bottom=0.085,
                    wspace=0.45, hspace=0.65)
for ax, letter, draw in zip(axes.flat, "abcd", [draw_panel_a, draw_panel_b,
                                                draw_panel_c, draw_panel_d]):
    draw(ax)
    pos = ax.get_position()
    fig.text(pos.x0 - 0.026, pos.y1 + 0.012, letter, fontsize=9,
             fontweight="bold", ha="left", va="bottom")

fig.savefig(OUT / "FigS5_MutationType.pdf")
fig.savefig(OUT / "FigS5_MutationType.svg")
fig.savefig(OUT / "FigS5_MutationType.tiff", dpi=300,
            pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(OUT / "FigS5_MutationType.png", dpi=300)
plt.close(fig)
for ext in ["pdf", "svg", "tiff", "png"]:
    print(f"  wrote output/figures/FigS5_MutationType.{ext}")

REVIEW.mkdir(parents=True, exist_ok=True)
shutil.copyfile(OUT / "FigS5_MutationType.png", REVIEW / "FigS5_MutationType.png")
print(f"  copied PNG -> {REVIEW / 'FigS5_MutationType.png'}")
print("FigS5 regeneration complete.")
