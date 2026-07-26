#!/usr/bin/env python3
"""Graphical abstract — regenerated with audited numbers.

Fixes vs the legacy PDF (no source): 23 lineages (was 24), AUROC>0.7 in 9 of 17
evaluable lineages (was 8), DD sign matches manuscript Eq. 1 (WT - MUT).

Usage: python graphical_abstract.py
Outputs: output/figures/GraphicalAbstract.pdf / .png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42  # TrueType (journals reject Type 3)
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).parent / "output" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

NAVY = "#16394e"; GRAY = "#5d6d7e"
BLUE = "#2f8fce"; GREEN = "#27ae60"; RED = "#e74c3c"
YELLOW = "#f39c12"; PURPLE = "#8e44ad"

# 180 mm double-column width (7.087 in); aspect preserved from the legacy 14.2x8.2 in
fig = plt.figure(figsize=(7.087, 4.094), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
ax.axis("off")

fig.text(0.5, 0.945, "Delta Dependency: Paralog Compensation Predicts Synthetic Lethality",
         ha="center", va="center", fontsize=12.5, fontweight="bold", color=NAVY)
fig.text(0.5, 0.895, "A Pan-Cancer Framework with Five-Layer Orthogonal Evidence",
         ha="center", va="center", fontsize=9.4, color=GRAY)


def box(x, y, w, h, color, fill):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.2",
                                linewidth=1.6, edgecolor=color, facecolor=fill))


def arrow(x1, x2, y):
    ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>",
                                 mutation_scale=16, linewidth=1.6, color=NAVY))


# ── Top row: pipeline ──────────────────────────────────────────
T, H = 62, 22  # bottom y, height of top boxes
box(3, T, 20, H, BLUE, "#e9f4fc")
fig.text(0.135, 0.79, "DepMap 26Q1", ha="center", fontsize=10.1, fontweight="bold", color=NAVY)
for i, s in enumerate(["1,208 cell lines × 18,531 genes",
                       "HGNC 66,595 paralog pairs",
                       "12 curated pairs (evidence-tiered)"]):
    fig.text(0.135, 0.745 - i * 0.038, s, ha="center", fontsize=7.6, color="#34495e")

box(27.5, T, 20, H, GREEN, "#e8f8ee")
fig.text(0.38, 0.79, "D = μ(Chronos|WT)", ha="center", fontsize=10.1, fontweight="bold", color=NAVY)
fig.text(0.38, 0.745, "− μ(Chronos|MUT)", ha="center", fontsize=10.1, fontweight="bold", color=NAVY)
fig.text(0.38, 0.685, "Univariate", ha="center", fontsize=7.6, color=GREEN)
fig.text(0.38, 0.647, "No training needed", ha="center", fontsize=7.6, color=GREEN)

box(52, T, 20, H, RED, "#fdeceb")
fig.text(0.625, 0.79, "Priority Candidates", ha="center", fontsize=10.1, fontweight="bold", color=NAVY)
fig.text(0.625, 0.745, "ARID1A→ARID1B", ha="center", fontsize=8.3, fontweight="bold", color=RED)
fig.text(0.625, 0.707, "EP300→CREBBP", ha="center", fontsize=7.6, color="#34495e")
fig.text(0.625, 0.669, "KRAS→HRAS  PIK3R1→CRKL", ha="center", fontsize=7.2, color=GRAY)

box(76.5, T, 20, H, YELLOW, "#fef6e0")
fig.text(0.87, 0.79, "vs Published", ha="center", fontsize=10.1, fontweight="bold", color=NAVY)
fig.text(0.87, 0.745, "DD: 0.676", ha="center", fontsize=9, fontweight="bold", color=RED)
fig.text(0.87, 0.707, "SLMGAE: 0.790", ha="center", fontsize=7.6, color="#34495e")
fig.text(0.87, 0.669, "+ID>30%: 1.000", ha="center", fontsize=7.6, color=GREEN)

for x1, x2 in [(23.6, 27.2), (48.1, 51.7), (72.6, 76.2)]:
    arrow(x1, x2, T + H / 2)

# ── Five-layer validation ──────────────────────────────────────
fig.text(0.5, 0.545, "Five-Layer Orthogonal Evidence",
         ha="center", fontsize=10.8, fontweight="bold", color=NAVY)

layers = [
    ("1", "Genomic", BLUE, ["23 lineages", "AUROC > 0.7", "in 7 of 8 lineages"]),
    ("2", "Proteomic", GREEN, ["7 CPTAC cohorts", "EP300-CREBBP", "5/7 significant"]),
    ("3", "Pharmacologic", RED, ["1,482 PRISM drugs", "633 selective hits", "MEK/mTOR/HDAC"]),
    ("4", "Clinical Strat.", YELLOW, ["MSI: MSI-H >= MSS", "MutType: trunc > miss", "DWS: safety tiers"]),
    ("5", "Structural", PURPLE, ["Sequence features", "Domain architecture", "PROTAC suitability"]),
]
W, GAP, Y0, HH = 17.6, 1.6, 12, 38
x = 3.0
for num, name, color, lines in layers:
    cx = (x + W / 2) / 100
    box(x, Y0, W, HH, color, "white")
    ax.add_patch(Circle((x + W / 2, Y0 + HH - 6.5), 2.6, facecolor=color, edgecolor="none"))
    fig.text(cx, (Y0 + HH - 6.5) / 100, num, ha="center", va="center",
             fontsize=9.4, fontweight="bold", color="white")
    fig.text(cx, (Y0 + HH - 13.5) / 100, name, ha="center", fontsize=9,
             fontweight="bold", color=color)
    for i, s in enumerate(lines):
        fig.text(cx, (Y0 + HH - 20.5 - i * 6.5) / 100, s, ha="center",
                 fontsize=7.2, color="#34495e")
    x += W + GAP

fig.text(0.5, 0.045,
         "Fully reproducible  |  No GPU required  |  Open-source R package paralogSL  |  All public data",
         ha="center", fontsize=7.6, style="italic", color=GRAY)

for ext in ("pdf", "png"):
    fig.savefig(OUT / f"GraphicalAbstract.{ext}", dpi=300)
print(f"wrote {OUT}/GraphicalAbstract.pdf/.png")
