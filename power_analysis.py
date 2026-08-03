#!/usr/bin/env python3
"""
power_analysis.py
=================
Analytic (no fabricated data) power calculation for the AUROC benchmark,
using the binormal AUROC model with the Hanley-McNeil (1982) variance.

Formulas (recorded in the output JSON metadata)
-----------------------------------------------
Hanley & McNeil (1982), "The meaning and use of the area under a receiver
operating characteristic (ROC) curve", Radiology 143:29-36:

    Q1 = A / (2 - A)            Q2 = 2 A^2 / (1 + A)
    Var(A) = [ A(1-A) + (n1-1)(Q1 - A^2) + (n0-1)(Q2 - A^2) ] / (n1 n0)

with n1 = # positives, n0 = # controls. Under the null A = 0.5 this
reduces to  Var0 = 0.25 * (1/n1 + 1/n0).

One-sided test of A vs 0.5 at alpha = 0.05:
    reject when Ahat > T,  T = 0.5 + z(1-alpha) * SE0
    power  = P(Ahat > T | A) = Phi( (A - T) / SE1 )
with SE0 = sqrt(Var0), SE1 = sqrt(Var at the alternative A), and
z(1-0.05) = 1.644854 (scipy.stats.norm.ppf(0.95)).

Required sample size: smallest integer n1 with power >= 0.80, found by
integer search (bisection-free linear scan; monotone in n1 in the
configurations used here). The equal-variance closed form
    n_per_group = 0.5 * ((z_{1-alpha} + z_{power}) / (A - 0.5))^2   (1:1)
is reported alongside as the closed-form cross-check that reproduces the
manuscript's "an estimated >= 25 validated positive pairs would be needed
to achieve 80% power ... for detecting AUROC = 0.85 against a null of
0.5" (Evaluation frameworks, Methods).

Questions answered (all from observed universe sizes; no fabricated data):
  (i)   power of the CURRENT positive-set sizes — 8 positive entries vs
        102 controls (lineage-level frame) and 6 positive pairs vs 66
        (per-pair frame) — to detect the observed AUROC 0.629 vs 0.5.
  (ii)  positives needed for 80% power at AUROC 0.629: (a) with the
        control count fixed at 102, (b) at a 1:8 positive:negative ratio.
  (iii) positives needed at AUROC 0.85 for reference (manuscript: ~25).

Output: output/power_analysis.json
Usage:  python3 power_analysis.py   (run from repo root)
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent
JSON_OUT = ROOT / "output" / "power_analysis.json"

ALPHA = 0.05          # one-sided
TARGET_POWER = 0.80
A_OBS = 0.629         # primary lineage-level AUROC (TableS2, 110 entries)
A_REF = 0.85          # manuscript reference effect
Z_A = float(norm.ppf(1 - ALPHA))      # 1.644854
Z_B = float(norm.ppf(TARGET_POWER))   # 0.841621


def hm_var(a, n1, n0):
    """Hanley-McNeil variance of the AUROC under the binormal model."""
    q1 = a / (2.0 - a)
    q2 = 2.0 * a * a / (1.0 + a)
    return (a * (1.0 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0)


def power(a, n1, n0):
    """One-sided power of AUROC A vs 0.5 at alpha (Hanley-McNeil, two variances)."""
    se0 = np.sqrt(hm_var(0.5, n1, n0))
    se1 = np.sqrt(hm_var(a, n1, n0))
    threshold = 0.5 + Z_A * se0
    return float(norm.cdf((a - threshold) / se1))


def n_pos_required(a, ratio=None, n0_fixed=None, cap=100_000):
    """Smallest n1 with power >= TARGET_POWER.

    ratio:     n0 = round(ratio * n1)  (e.g. 8 for a 1:8 pos:neg ratio)
    n0_fixed:  hold the control count constant and grow only positives
    """
    for n1 in range(2, cap):
        n0 = n0_fixed if n0_fixed is not None else int(round(ratio * n1))
        if power(a, n1, n0) >= TARGET_POWER:
            return n1, n0, power(a, n1, n0)
    return None, None, None


def main():
    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "power_analysis.py",
        "method": "binormal AUROC model, Hanley-McNeil (1982) variance; one-sided alpha=0.05",
        "formulas": {
            "variance": "Var(A) = [A(1-A) + (n1-1)(Q1-A^2) + (n0-1)(Q2-A^2)] / (n1*n0), "
                        "Q1 = A/(2-A), Q2 = 2A^2/(1+A); under H0 (A=0.5): Q1=Q2=1/3, so "
                        "Var0 = [0.25 + (n1+n0-2)/12] / (n1*n0)",
            "rejection_threshold": "T = 0.5 + z(1-alpha)*SE0, z(0.95) = 1.644854 (one-sided)",
            "power": "power = Phi((A - T)/SE1), SE1 = sqrt(Var(A)) at the alternative",
            "required_n": "smallest integer n1 with power >= 0.80 (integer scan)",
            "equal_variance_closed_form": "n_per_group = 0.5*((z_(1-alpha)+z_power)/(A-0.5))^2 "
                                          "(1:1 ratio; cross-check for the manuscript's ~25)",
        },
        "inputs": {
            "alpha_one_sided": ALPHA,
            "target_power": TARGET_POWER,
            "auroc_observed": A_OBS,
            "auroc_reference": A_REF,
            "current_sizes": {
                "lineage_level_frame": {"n_pos": 8, "n_neg": 102,
                                        "note": "110 TableS2 entries = 8 positive + 102 control"},
                "per_pair_frame": {"n_pos": 6, "n_neg": 66,
                                   "note": "72 unique pairs = 6 positive + 66 control"},
            },
        },
    }

    # (i) Power of the current positive-set sizes at the observed AUROC
    cur = {}
    for name, n1, n0 in (("lineage_level_8v102", 8, 102), ("per_pair_6v66", 6, 66)):
        cur[name] = {
            "n_pos": n1, "n_neg": n0,
            "auroc_alternative": A_OBS,
            "power": power(A_OBS, n1, n0),
            "se0_under_null": float(np.sqrt(hm_var(0.5, n1, n0))),
            "se1_under_alternative": float(np.sqrt(hm_var(A_OBS, n1, n0))),
        }
    out["i_power_current_sizes"] = cur
    print(f"(i) power @ AUROC {A_OBS}: 8v102 -> {cur['lineage_level_8v102']['power']:.4f}, "
          f"6v66 -> {cur['per_pair_6v66']['power']:.4f}")

    # (ii) Positives needed for 80% power at AUROC 0.629
    n_fix, n0_fix, p_fix = n_pos_required(A_OBS, n0_fixed=102)
    n_rat, n0_rat, p_rat = n_pos_required(A_OBS, ratio=8)
    out["ii_required_for_80pct_at_observed"] = {
        "same_control_count_102": {"n_pos_required": n_fix, "n_neg": n0_fix,
                                   "power_achieved": p_fix},
        "ratio_1_to_8": {"n_pos_required": n_rat, "n_neg": n0_rat,
                         "power_achieved": p_rat,
                         "note": "n_neg = 8 x n_pos"},
    }
    print(f"(ii) n_pos for 80% @ {A_OBS}: n0=102 fixed -> {n_fix}; "
          f"1:8 ratio -> {n_rat} (vs {n0_rat} controls)")

    # (iii) Reference: AUROC 0.85 (manuscript: ~25 positives)
    n_hm, n0_hm, p_hm = n_pos_required(A_REF, ratio=1)   # full H-M, 1:1
    n_cf = 0.5 * ((Z_A + Z_B) / (A_REF - 0.5)) ** 2      # equal-variance closed form, 1:1
    out["iii_reference_auroc_0.85"] = {
        "hanley_mcneil_two_variance_1to1": {
            "n_pos_required": n_hm, "n_neg": n0_hm, "power_achieved": p_hm,
            "note": "integer scan with distinct null/alternative variances; "
                    "more efficient because SE1 < SE0 at high A"},
        "equal_variance_closed_form_1to1": {
            "n_per_group_continuous": float(n_cf),
            "n_per_group_ceiling": int(np.ceil(n_cf)),
            "formula": "n = 0.5*((1.644854+0.841621)/(0.85-0.5))^2",
        },
        "manuscript_quote": "an estimated >= 25 validated positive pairs would be needed to "
                            "achieve 80% power ... for detecting AUROC = 0.85 against a null "
                            "of 0.5 (Evaluation frameworks, Methods)",
        "match_note": "the closed-form continuous estimate (25.2) reproduces the "
                      "manuscript's '~25'; the ceiling integer is 26",
    }
    print(f"(iii) n @ {A_REF}: H-M 1:1 -> {n_hm}; closed form -> {n_cf:.2f} "
          f"(ceil {int(np.ceil(n_cf))}) [manuscript ~25]")

    JSON_OUT.write_text(json.dumps(out, indent=2, allow_nan=False, default=str))
    print(f"\nWrote {JSON_OUT}")


if __name__ == "__main__":
    main()
