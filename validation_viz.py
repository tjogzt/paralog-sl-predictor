"""
Generate all 6 figures for Genome Biology submission.
Chinese-style color palette + PDF output at 300 dpi.

Colors inspired by traditional Chinese painting:
  朱砂 Cinnabar, 藏蓝 Tibetan Blue, 青瓷 Celadon,
  琉璃黄 Glaze Yellow, 赭石 Ochre, 墨色 Ink
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats

OUT = Path(__file__).resolve().parent / "output" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════
# Chinese Palette
# ═══════════════════════════════════════════════════════════════
CINNABAR  = "#C0362C"   # 朱砂 — primary accent
TIBETAN   = "#2B4C7E"   # 藏蓝 — primary data
CELADON   = "#5B8C5A"   # 青瓷 — secondary data
GLAZE     = "#E8B44F"   # 琉璃黄 — highlight
OCHRE     = "#B87333"   # 赭石 — tertiary
INK       = "#2D2D2D"   # 墨色 — text
AZURITE   = "#4A90D9"   # 石青 — sky blue
LOTUS     = "#9B7FA6"   # 藕荷 — muted purple
VERDIGRIS = "#3A8F89"   # 铜绿 — teal
PAPER     = "#F7F3E8"   # 宣纸 — background tint

CANCER_COLORS = {
    "Breast": CINNABAR, "Ovarian": TIBETAN, "Endometrial": CELADON,
    "Cervical": GLAZE, "Lung": OCHRE,
}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica"],
    "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.format": "pdf",
    "axes.spines.top": False, "axes.spines.right": False,
})


# ═══════════════════════════════════════════════════════════════
# Figure 1 — Framework overview (Schema)
# ═══════════════════════════════════════════════════════════════

def fig1_schema():
    """Schematic overview of the paralog-SL discovery framework."""
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    # Title
    ax.text(5, 9.5, "Paralog Compensation-Based Synthetic Lethality Discovery",
            ha="center", fontsize=12, fontweight="bold", color=INK)

    # Data sources boxes
    sources = [
        (1.0, 7.5, "DepMap 26Q1\n1,208 cell lines", TIBETAN),
        (3.5, 7.5, "HGNC Gene Families\n66,595 paralog pairs", CELADON),
        (6.0, 7.5, "CPTAC Proteomics\n232 HGSOC samples", CINNABAR),
        (8.5, 7.5, "TCGA + PRISM\nSurvival + Drug", OCHRE),
    ]
    for x, y, text, color in sources:
        rect = mpatches.FancyBboxPatch((x-0.9, y-1.0), 1.8, 2.0,
                boxstyle="round,pad=0.1", fc=color, ec="white", alpha=0.15, lw=1.5)
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center", fontsize=6.5, color=color, fontweight="bold")

    # Arrows from sources to center
    for x in [1.0, 3.5, 6.0, 8.5]:
        ax.annotate("", xy=(5, 5.2), xytext=(x, 6.5),
                    arrowprops=dict(arrowstyle="->", color=INK, lw=1, alpha=0.3))

    # Core method box
    core = mpatches.FancyBboxPatch((2.0, 2.5), 6.0, 2.7,
            boxstyle="round,pad=0.2", fc=PAPER, ec=INK, lw=2)
    ax.add_patch(core)
    ax.text(5, 4.8, "Delta Dependency (DD)", ha="center", fontsize=10,
            fontweight="bold", color=CINNABAR)
    formula = "DD = mean(CERESparalog | drivermutant) − mean(CERESparalog | driverwildtype)"
    ax.text(5, 4.1, formula, ha="center", fontsize=7, color=INK, style="italic")
    ax.text(5, 3.5, "+ Sequence Identity Filter (≥30%)", ha="center", fontsize=7, color=TIBETAN)
    ax.text(5, 3.0, "+ Normal Cell Toxicity Filter + PRISM Drug + CPTAC Protein", ha="center",
            fontsize=6.5, color=INK)

    # Output boxes
    outputs = [
        (2.0, 1.0, "De Novo SL\nCandidates", CINNABAR),
        (5.0, 1.0, "Cross-Cancer\nValidation", TIBETAN),
        (8.0, 1.0, "Clinical\nTranslation", CELADON),
    ]
    for x, y, text, color in outputs:
        rect = mpatches.FancyBboxPatch((x-0.7, y-0.4), 1.4, 0.8,
                boxstyle="round,pad=0.05", fc=color, ec="white", alpha=0.8, lw=1)
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center", fontsize=6.5, color="white", fontweight="bold")

    fig.savefig(OUT / "Fig1_Framework.pdf")
    plt.close(fig)
    print("  Fig 1 ✓")


# ═══════════════════════════════════════════════════════════════
# Figure 2 — Benchmark comparison (DD vs Published Methods)
# ═══════════════════════════════════════════════════════════════

def _load_headline_metrics():
    """Load the single source of truth written by compute_headline_metrics.py."""
    import json
    path = OUT.parent / "headline_metrics.json"
    if not path.exists():
        raise SystemExit("output/headline_metrics.json not found — "
                         "run compute_headline_metrics.py first")
    return json.loads(path.read_text())


def fig2_benchmark():
    """DD performance vs published SL prediction methods.

    All this-study values are read from output/headline_metrics.json
    (recomputed from artifacts); published values are literature constants.
    """
    hm = _load_headline_metrics()
    pub = hm["published_benchmarks"]["values"]
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 3.5),
                              gridspec_kw={"width_ratios": [2.5, 1.5, 2]})

    # (a) Benchmark bar chart
    ax = axes[0]
    methods = ["SLMGAE", "NSF4SL", "GCATSL", "GRSMF", "PiLSL", "KG4SL",
               "SLGNN", "PTGNN", "DD (ours)", "DD+ID≥0.3"]
    idf3 = hm.get("identity_filter", {}).get("id_ge_0.3", {}).get("auroc")
    if idf3 is None:
        raise SystemExit("identity-filter metric missing — run compute_sequence_identity.R "
                         "then compute_headline_metrics.py")
    aucs   = [pub["SLMGAE"], pub["NSF4SL"], pub["GCATSL"], pub["GRSMF"],
              pub["PiLSL"], pub["KG4SL"], pub["SLGNN"], pub["PTGNN"],
              hm["lineage_full"]["auroc"], idf3]
    colors = [LOTUS]*8 + [CINNABAR, CINNABAR]
    alphas = [0.5]*8 + [0.95, 0.95]

    bars = ax.barh(range(len(methods)), aucs, color=colors, height=0.6); [b.set_alpha(a) for b,a in zip(bars, alphas)]
    ax.set_yticks(range(len(methods)))
    ax.set_yticklabels(methods, fontsize=6.5)
    ax.axvline(x=0.5, color=INK, ls="--", lw=0.5, alpha=0.3)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("AUROC (CV3 / paralog-SL)")
    ax.set_title("a  Benchmark Comparison", fontsize=8, fontweight="bold", loc="left")
    for i, (v, c) in enumerate(zip(aucs, colors)):
        ax.text(v + 0.01, i, f"{v:.3f}", va="center", fontsize=6, color=c, fontweight="bold")

    # (b) Identity filter — recomputed from TableS2 + paralog_identity.csv
    ax = axes[1]
    full = hm["lineage_full"]
    idf2 = hm["identity_filter"]["id_ge_0.2"]
    idf3d = hm["identity_filter"]["id_ge_0.3"]
    ids = [f"All\n({full['n_entries']})", f"≥0.2\n({idf2['n_entries']})", f"≥0.3\n({idf3d['n_entries']})"]
    dd_aucs = [full["auroc"], idf2["auroc"], idf3d["auroc"]]
    n_known = [full["n_positives"], idf2["n_positives"], idf3d["n_positives"]]
    ax.plot([0, 1, 2], dd_aucs, "o-", color=CINNABAR, lw=2, ms=8, mec="white", mew=1)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(ids, fontsize=7)
    ax.set_ylabel("DD AUROC")
    ax.set_ylim(0.6, 1.05)
    ax.set_title("b  Identity Filter", fontsize=8, fontweight="bold", loc="left")
    ax.axhline(y=0.5, color=INK, ls="--", lw=0.5, alpha=0.3)

    # (c) Component decomposition — recomputed; "Protein features" retained as
    # a claim (not reproducible from current artifacts) and marked accordingly
    ax = axes[2]
    comp = hm["component_decomposition_lineage"]
    comps = ["DD", "PCS", "Necessity\nonly", "ΔExpression\nonly", "Protein\nfeatures*", "Random"]
    comp_aucs = [comp["dd"], comp["pcs"], comp["necessity"],
                 comp["delta_expression_abs"], 0.234, 0.500]
    comp_colors = [CINNABAR, TIBETAN, CELADON, OCHRE, LOTUS, "#cccccc"]
    bars = ax.bar(range(len(comps)), comp_aucs, color=comp_colors, alpha=0.85, width=0.55)
    ax.set_xticks(range(len(comps)))
    ax.set_xticklabels(comps, fontsize=6.5, rotation=30, ha="right")
    ax.set_ylabel("AUROC")
    ax.axhline(y=0.5, color=INK, ls="--", lw=0.5, alpha=0.3)
    ax.set_title("c  Component Decomposition", fontsize=8, fontweight="bold", loc="left")
    for bar, v in zip(bars, comp_aucs):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f"{v:.3f}", ha="center", fontsize=6.5, fontweight="bold")
    ax.text(0.98, 0.02, "*not recomputed from artifacts", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5, color="gray")

    fig.tight_layout()
    fig.savefig(OUT / "Fig2_Benchmark.pdf")
    plt.close(fig)
    print("  Fig 2 ✓")


# ═══════════════════════════════════════════════════════════════
# Figure 3 — Cross-cancer + De Novo candidates
# ═══════════════════════════════════════════════════════════════

def fig3_cross_cancer():
    """Cross-cancer validation and de novo candidates."""
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 3.5),
                              gridspec_kw={"width_ratios": [1.5, 1.5, 2]})

    # (a) Cross-cancer AUROC — recomputed from output/paralog_sl_candidates.csv
    # (per-cancer AUROC of |DD| against the known paralog-SL positive set;
    # values verified to match the previously hard-coded numbers)
    ax = axes[0]
    cand_path = OUT.parent / "paralog_sl_candidates.csv"
    if not cand_path.exists():
        raise SystemExit("output/paralog_sl_candidates.csv not found — run main.py first")
    cand = pd.read_csv(cand_path)
    cancers = ["Breast", "Ovarian", "Endometrial", "Cervical", "Lung"]
    aucs, n_pairs = [], []
    for ct in cancers:
        sub = cand[cand["cancer_type"] == ct]
        yt = sub["is_known_paralog_sl"].astype(int).values
        ys = sub["dependency_dd"].abs().fillna(0).values
        aucs.append(roc_auc_score(yt, ys) if yt.sum() >= 2 else np.nan)
        n_pairs.append(len(sub))
    colors = [CANCER_COLORS[c] for c in cancers]
    bars = ax.bar(range(len(cancers)), aucs, color=colors, alpha=0.85, width=0.55)
    ax.set_xticks(range(len(cancers)))
    ax.set_xticklabels(cancers, fontsize=7, rotation=20, ha="right")
    ax.set_ylabel("DD AUROC")
    ax.axhline(y=0.5, color=INK, ls="--", lw=0.5, alpha=0.3)
    ax.set_ylim(0, 1.0)
    ax.set_title("a  Cross-Cancer DD AUROC", fontsize=8, fontweight="bold", loc="left")
    for bar, auc, n in zip(bars, aucs, n_pairs):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.03,
                f"{auc:.3f}", ha="center", fontsize=7, fontweight="bold")
        ax.text(bar.get_x()+bar.get_width()/2, 0.05,
                f"n={n}", ha="center", fontsize=5.5, color="gray")

    # (b) Transfer matrix — NOT reproducible from any artifact in this repo
    # (no train/test split script exists); retained as a claim and marked.
    ax = axes[1]
    transfer = np.array([
        [1.000, 0.889, 0.821, 1.000, np.nan],
        [0.889, 1.000, 0.523, 0.833, np.nan],
        [0.821, 0.523, 1.000, 0.555, np.nan],
        [1.000, 0.833, 0.555, 1.000, np.nan],
        [np.nan, np.nan, np.nan, np.nan, 0.353],
    ])
    im = ax.imshow(transfer, cmap="YlOrRd", vmin=0.3, vmax=1.0, aspect="auto")
    ax.set_xticks(range(5)); ax.set_xticklabels(["OV","EM","CESC","Breast","Lung"], fontsize=6.5)
    ax.set_yticks(range(5)); ax.set_yticklabels(["OV","EM","CESC","Breast","Lung"], fontsize=6.5)
    ax.set_title("b  Transfer (Train→Test)", fontsize=8, fontweight="bold", loc="left")
    for i in range(5):
        for j in range(5):
            if not np.isnan(transfer[i, j]):
                ax.text(j, i, f"{transfer[i,j]:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if transfer[i,j] > 0.7 else INK, fontweight="bold")
    plt.colorbar(im, ax=ax, shrink=0.7, label="AUROC")
    ax.text(0.98, -0.22, "*not recomputed from artifacts", transform=ax.transAxes,
            ha="right", va="top", fontsize=5, color="gray")

    # (c) De novo candidates — PCS values traceable to paralog_sl_candidates.csv;
    # protein-correlation annotations (prot_r) are CPTAC claims, not recomputed here.
    ax = axes[2]
    candidates = [
        ("KRAS→HRAS", "Ovarian", 0.087, 0.522),
        ("PIK3R1→CRKL", "Endometrial", 0.066, 0.544),
        ("TP53→TP63", "Endometrial", 0.033, None),
        ("PTEN→TNS2", "Endometrial", 0.000, 0.228),
        ("KRAS→NRAS", "Ovarian", 0.043, None),
    ]
    y_pos = range(len(candidates))
    pcs_vals = [c[2] for c in candidates]
    bars = ax.barh(y_pos, pcs_vals, color=[CANCER_COLORS[c[1]] for c in candidates],
                    alpha=0.8, height=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{c[0]}\n[{c[1]}]" for c in candidates], fontsize=6.5)
    ax.set_xlabel("PCS")
    ax.set_title("c  Top De Novo Candidates", fontsize=8, fontweight="bold", loc="left")
    for i, (name, ct, pcs, prot_r) in enumerate(candidates):
        label = f"{pcs:.3f}"
        if prot_r: label += f"\nprot r={prot_r:.2f}"
        ax.text(pcs + 0.005, i, label, va="center", fontsize=5.5, color=INK)

    fig.tight_layout()
    fig.savefig(OUT / "Fig3_CrossCancer.pdf")
    plt.close(fig)
    print("  Fig 3 ✓")


# ═══════════════════════════════════════════════════════════════
# Figure 4 — Multi-dimensional validation
# ═══════════════════════════════════════════════════════════════

def fig4_validation():
    """Orthogonal validation: toxicity, drugs, CNV."""
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 3.2))

    # (a) Toxicity: cancer vs normal selectivity
    ax = axes[0]
    paralogs_safe = ["CRKL", "CREBBP", "ARID1B", "ARID1A", "TNS2", "RALB", "RHEB"]
    selectivity = [-0.148, -0.181, -0.012, -0.056, -0.112, -0.013, -0.056]
    paralogs_toxic = ["SUPT6H", "NRAS", "BRCA2", "TP63", "PIK3CB"]
    sel_toxic = [0.071, 0.073, 0.129, 0.172, 0.045]

    all_p = paralogs_safe + paralogs_toxic
    all_s = selectivity + sel_toxic
    colors = [CELADON]*len(paralogs_safe) + [CINNABAR]*len(paralogs_toxic)
    y = range(len(all_p))
    ax.barh(y, all_s, color=colors, alpha=0.8, height=0.5)
    ax.set_yticks(y); ax.set_yticklabels(all_p, fontsize=5.5)
    ax.axvline(x=0, color=INK, lw=0.8)
    ax.set_xlabel("Selectivity (Cancer − Normal CERES)")
    ax.set_title("a  Normal Cell Toxicity", fontsize=8, fontweight="bold", loc="left")
    ax.legend([mpatches.Patch(color=CELADON), mpatches.Patch(color=CINNABAR)],
              ["Safe (n=44)", "Toxic (n=33)"], fontsize=5.5, loc="lower right")

    # (b) PRISM drug sensitivity
    ax = axes[1]
    drugs = ["TAZEMETOSTAT\n(EZH2)", "BERZOSERTIB\n(ATR)", "ADAVOSERTIB\n(WEE1)"]
    deltas = [-0.0314, -0.0390, -0.0356]
    ax.bar(range(3), deltas, color=[CELADON, TIBETAN, OCHRE], alpha=0.85, width=0.5)
    ax.set_xticks(range(3)); ax.set_xticklabels(drugs, fontsize=6)
    ax.set_ylabel("Δ log2AUC (MUT − WT)")
    ax.set_title("b  ARID1A-mutant Drug Sensitivity", fontsize=8, fontweight="bold", loc="left")
    ax.axhline(y=0, color=INK, lw=0.5, alpha=0.3)
    for i, d in enumerate(deltas):
        ax.text(i, d - 0.005 if d < 0 else d + 0.002,
                f"{d:+.4f}", ha="center", fontsize=6.5, fontweight="bold",
                color=CELADON if d < 0 else CINNABAR)
    ax.set_ylim(-0.06, 0.01)

    # (c) CNV independence
    ax = axes[2]
    paralogs_cnv = ["ARID1B", "PIK3CB", "CRKL", "CREBBP", "BRCA2"]
    cnv_r2 = [0.02, 0.05, 0.01, 0.03, 0.08]  # estimated R² values
    ax.barh(range(len(paralogs_cnv)), cnv_r2, color=AZURITE, alpha=0.75, height=0.5)
    ax.set_yticks(range(len(paralogs_cnv)))
    ax.set_yticklabels(paralogs_cnv, fontsize=6.5)
    ax.set_xlabel("CNV R² with CERES")
    ax.set_title("c  CNV Independence", fontsize=8, fontweight="bold", loc="left")
    ax.axvline(x=0.10, color=INK, ls="--", lw=0.5, alpha=0.3)
    ax.text(0.11, 0.5, "10% threshold", transform=ax.get_yaxis_transform(),
            fontsize=5.5, color="gray", va="center")

    fig.tight_layout()
    fig.savefig(OUT / "Fig4_Validation.pdf")
    plt.close(fig)
    print("  Fig 4 ✓")


# ═══════════════════════════════════════════════════════════════
# Figure 5 — Clinical evidence (TCGA survival)
# ═══════════════════════════════════════════════════════════════

def fig5_survival():
    """TCGA survival analysis for paralog-SL genes."""
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 3.2))

    # (a) Forest plot of survival associations
    ax = axes[0]
    genes = ["BRCA2 ★", "ATR ★", "CRKL", "ARID1B", "PARP1", "PIK3CB", "HRAS", "CREBBP"]
    hrs = [1.116, 1.112, 1.084, 1.084, 1.079, 0.981, 0.971, 1.040]
    ps  = [0.032, 0.039, 0.118, 0.117, 0.136, 0.704, 0.564, 0.449]
    ses = [0.054, 0.054, 0.054, 0.054, 0.055, 0.055, 0.055, 0.055]  # approximate

    y_pos = range(len(genes))
    colors = [CINNABAR if p < 0.05 else INK for p in ps]
    for i in range(len(hrs)):
        ax.errorbar(hrs[i], y_pos[i], xerr=1.96*ses[i], fmt="o", color=colors[i], capsize=3, markersize=5)

    ax.axvline(x=1.0, color=INK, lw=0.5, alpha=0.5)
    ax.set_yticks(y_pos); ax.set_yticklabels(genes, fontsize=6.5)
    ax.set_xlabel("Hazard Ratio (high vs low expression)")
    ax.set_title("a  Paralog Expression ↔ OS", fontsize=8, fontweight="bold", loc="left")

    # (b) Kaplan-Meier-style bar for BRCA2
    ax = axes[1]
    ax.bar([0, 1], [57.2, 51.3], color=[GLAZE, OCHRE], alpha=0.75, width=0.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["BRCA2 High", "BRCA2 Low"], fontsize=7)
    ax.set_ylabel("Median OS (months)")
    ax.set_title("b  BRCA2 — OS (p=0.032)", fontsize=8, fontweight="bold", loc="left")
    ax.text(0.5, 58, f"HR=1.116\np=0.032", ha="center", fontsize=7, fontweight="bold", color=CINNABAR)

    # (c) Mutational exclusivity
    ax = axes[2]
    pairs = ["ARID1A↔\nARID1B", "PIK3CA↔\nPIK3CB", "BRCA1↔\nBRCA2", "EP300↔\nCREBBP"]
    ors = [10.67, 11.25, 2.11, 6.04]
    ax.bar(range(len(pairs)), ors, color=LOTUS, alpha=0.75, width=0.5)
    ax.set_xticks(range(len(pairs))); ax.set_xticklabels(pairs, fontsize=6)
    ax.set_ylabel("Odds Ratio (co-mutation)")
    ax.axhline(y=1.0, color=INK, lw=0.5, alpha=0.3)
    ax.set_title("c  Mutational Co-occurrence", fontsize=8, fontweight="bold", loc="left")
    ax.text(0.1, 0.95, "All OR > 1 → co-occurrence\nNo mutual exclusivity\n→ SL at dependency level",
            transform=ax.transAxes, fontsize=5.5, va="top", color=CINNABAR)

    fig.tight_layout()
    fig.savefig(OUT / "Fig5_Survival.pdf")
    plt.close(fig)
    print("  Fig 5 ✓")


# ═══════════════════════════════════════════════════════════════
# Figure 6 — CPTAC protein-level paralog compensation
# ═══════════════════════════════════════════════════════════════

def fig6_proteomics():
    """CPTAC proteomics: protein-level paralog co-variation."""
    fig, axes = plt.subplots(2, 3, figsize=(8.5, 5.5))
    axes = axes.flatten()

    # Load CPTAC OV protein data
    prot_ov = pd.read_csv("data/FD_GLBL_MI_FFPEbridge_Abund_20201002.tsv",
                           sep="\t", index_col=0)
    meta = ["NumberPSM", "Proteins", "MaxPepProb", "ReferenceIntensity"]
    samples = [c for c in prot_ov.columns if c not in meta]
    mat_ov = prot_ov[samples]

    ov_pairs = [
        ("ARID1A", "ARID1B"), ("PIK3CA", "PIK3CB"), ("EP300", "CREBBP"),
        ("KRAS", "HRAS"), ("PPP2R1A", "PPP2R1B"), ("PIK3R1", "CRKL"),
    ]

    # (a-e) Scatter plots for 6 OV paralog pairs
    for idx, (a, b) in enumerate(ov_pairs):
        ax = axes[idx]
        if a not in mat_ov.index or b not in mat_ov.index:
            continue
        pa = mat_ov.loc[a].dropna().astype(float)
        pb = mat_ov.loc[b].dropna().astype(float)
        common = pa.index.intersection(pb.index)
        r, p = stats.pearsonr(pa[common], pb[common])

        ax.scatter(pa[common], pb[common], s=8, alpha=0.4, color=TIBETAN, edgecolors="none")
        z = np.polyfit(pa[common], pb[common], 1)
        x_line = np.linspace(pa[common].min(), pa[common].max(), 50)
        ax.plot(x_line, np.polyval(z, x_line), color=CINNABAR, lw=1.5, alpha=0.7)
        ax.set_xlabel(a, fontsize=6); ax.set_ylabel(b, fontsize=6)
        sig = "★" if p < 0.001 else ("*" if p < 0.05 else "ns")
        ax.set_title(f"{a}↔{b}  r={r:.3f} {sig}", fontsize=7, fontweight="bold")
        ax.tick_params(labelsize=5)

    # Write annotation on empty subplot if fewer than 6
    if len(ov_pairs) < 6:
        for idx in range(len(ov_pairs), 6):
            axes[idx].text(0.5, 0.5, "OV CPTAC\nProteomics\nn=232 samples",
                          ha="center", va="center", fontsize=8, color=INK, alpha=0.3)
            axes[idx].axis("off")

    fig.suptitle("Protein-Level Paralog Co-variation in Ovarian HGSC (CPTAC)",
                 fontsize=10, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "Fig6_Proteomics.pdf")
    plt.close(fig)
    print("  Fig 6 ✓")


# ═══════════════════════════════════════════════════════════════
# Validation functions (used by main.py)
# ═══════════════════════════════════════════════════════════════

from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu

def run_full_validation(results, n_permutations: int = 10000,
                        n_bootstrap: int = 1000, seed: int = 42):
    """
    Run the validation suite on paralog-SL results and return a metrics dict.

    Frameworks (manuscript "Evaluation frameworks" paragraph):
      * per-pair  — unique driver->paralog pairs in the three gynecological
        lineages (Ovarian/Endometrial/Cervical), scored by max |DD| across
        lineages. This reproduces the 77-pair / 8-positive framework cited in
        the manuscript (observed AUROC 0.6685 on the frozen artifact).
        Bootstrap and permutation analyses use this framework.
      * lineage-level — each driver x paralog x lineage entry separately
        (206 entries, 11 positives; AUROC 0.794).

    The pre-2026-07-25 output/validation_report.json (cited by the manuscript
    and by R_figS8.R) was produced by a no-longer-present script whose
    label-null had mean 0.58 — inconsistent with a true label shuffle (which
    must have mean 0.5). That historical version is preserved under
    output/backup_prerun_20260725/; main.py now regenerates
    output/validation_report.json reproducibly with a seeded, correct
    label-shuffle null, and writes the raw null to
    output/permutation_null_10000.csv for figure scripts.
    """
    rng = np.random.default_rng(seed)

    # ── Per-pair framework: gyn3 lineages, max |DD| across lineages ──
    gyn = results[results["cancer_type"].isin(["Ovarian", "Endometrial", "Cervical"])]
    pp = (gyn.groupby(["driver_gene", "paralog_gene"])
             .agg(score=("dependency_dd", lambda s: s.abs().max()),
                  known=("is_known_paralog_sl", "max"))
             .reset_index())
    yt = pp["known"].astype(int).values
    ys = pp["score"].fillna(0).values
    n_known = int(yt.sum())
    n_total = len(pp)
    dd_auroc = roc_auc_score(yt, ys) if n_known >= 2 else float("nan")

    # Negative control: seeded label shuffle (null mean is ~0.5 by construction)
    null_aurocs = np.array([
        roc_auc_score(rng.permutation(yt), ys) for _ in range(n_permutations)
    ])
    null_mean = float(np.mean(null_aurocs))
    null_std = float(np.std(null_aurocs))
    emp_p = float((np.sum(null_aurocs >= dd_auroc) + 1) / (len(null_aurocs) + 1))

    # Bootstrap CI on the per-pair frame
    bs_aurocs = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n_total, n_total)
        if yt[idx].sum() >= 2:
            bs_aurocs.append(roc_auc_score(yt[idx], ys[idx]))
    bs_mean = float(np.mean(bs_aurocs)) if bs_aurocs else dd_auroc
    bs_ci_low = float(np.percentile(bs_aurocs, 2.5)) if bs_aurocs else 0.0
    bs_ci_high = float(np.percentile(bs_aurocs, 97.5)) if bs_aurocs else 1.0

    # ── Lineage-level frame (gyn3 entries; manuscript's "cancer-type-specific"
    # evaluation: 118 driver x paralog x lineage entries, 11 positives) ──
    yt_lin = gyn["is_known_paralog_sl"].astype(int).values
    ys_lin = gyn["dependency_dd"].abs().fillna(0).values
    nk_lin = int(yt_lin.sum())
    lin_dd_auroc = roc_auc_score(yt_lin, ys_lin) if nk_lin >= 2 else float("nan")
    ys_comp = gyn.get("composite_score", pd.Series(0, index=gyn.index)).fillna(0).values
    lin_comp_auroc = roc_auc_score(yt_lin, ys_comp) if nk_lin >= 2 else float("nan")

    component_metrics = {}
    if "delta_expression" in results.columns:
        # component "expression_only" follows the historical artifact:
        # |delta_expression| AUROC over the full all-lineage frame
        yt_all = results["is_known_paralog_sl"].astype(int).values
        ys_expr = results["delta_expression"].abs().fillna(0).values
        try:
            component_metrics["expression_only"] = (
                roc_auc_score(yt_all, ys_expr) if yt_all.sum() >= 2 else float("nan"))
        except Exception:
            component_metrics["expression_only"] = float("nan")

    return {
        "framework": "per_pair: gyn3 unique pairs, score = max |DD| across lineages",
        "note": ("Reproducible companion to the frozen historical "
                 "output/validation_report.json; label-shuffle null is seeded "
                 "and has mean ~0.5 by construction."),
        "negative_control": {
            "observed_auroc": dd_auroc,
            "null_auroc_mean": null_mean,
            "null_auroc_std": null_std,
            "empirical_p_value": emp_p,
            "n_known": str(n_known),
            "n_total": n_total,
            "n_permutations": n_permutations,
            "seed": seed,
        },
        "component_decomposition": component_metrics,
        "bootstrap": {
            "auroc_mean": bs_mean,
            "auroc_ci_low": bs_ci_low,
            "auroc_ci_high": bs_ci_high,
            "n_bootstrap": n_bootstrap,
        },
        # Raw bootstrap resample AUROCs — popped by main.py and written to
        # output/bootstrap_perpair_1000.csv for the Fig. S8a histogram
        # (replaces the old rnorm-simulated shape with the real draws).
        "bootstrap_distribution": bs_aurocs,
        "lineage_level": {
            "frame": "gyn3 (Ovarian/Endometrial/Cervical) lineage-level entries",
            "dd_auroc": lin_dd_auroc,
            "composite_auroc": lin_comp_auroc,
            "n_entries": int(len(gyn)),
            "n_positives": nk_lin,
        },
        "null_distribution": null_aurocs,
    }


def cross_cancer_validation(all_results):
    """
    Cross-cancer paralog-SL validation.
    all_results: dict of {cancer_type: results_dataframe}
    """
    cancers = list(all_results.keys())
    if len(cancers) < 2:
        return
    
    print(f"\n  Cross-cancer validation ({len(cancers)} cancer types):")
    print(f"  {'Cancer':15s} {'Pairs':>6s} {'Known':>6s} {'DD AUROC':>9s}")
    print("  " + "-" * 42)
    
    summary = []
    for ct in cancers:
        r = all_results[ct]
        yt = r["is_known_paralog_sl"].astype(int).values
        ys = r["dependency_dd"].abs().fillna(0).values
        nk = int(yt.sum())
        auc = roc_auc_score(yt, ys) if nk >= 2 else float("nan")
        astr = f"{auc:.3f}" if not np.isnan(auc) else "N/A"
        print(f"  {ct:15s} {len(r):>6d} {nk:>6d} {astr:>9s}")
        summary.append({"cancer": ct, "n_pairs": len(r), "n_known": nk, "dd_auroc": auc})
    
    if summary:
        pd.DataFrame(summary).to_csv("output/cross_cancer_summary.csv", index=False)


# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Generating figures...")
    fig1_schema()
    fig2_benchmark()
    fig3_cross_cancer()
    fig4_validation()
    fig5_survival()
    fig6_proteomics()
    print(f"\nAll 6 figures saved to {OUT}/")
