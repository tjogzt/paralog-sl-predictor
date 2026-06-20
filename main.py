#!/usr/bin/env python3
"""Paralog SL Predictor v2 — with negative controls, bootstrap, cross-cancer validation."""

import argparse, sys, json
import pandas as pd
from config import RESULTS_FILE, SUMMARY_FILE
from data_loader import load_all_data
from pcs import run_full_analysis
from validation_viz import run_full_validation, cross_cancer_validation


def rank_and_export(results, results_file, summary_file):
    results_file = results_file or RESULTS_FILE
    summary_file = summary_file or SUMMARY_FILE
    results["novelty"] = results["is_known_paralog_sl"].map({True: "Known", False: "Novel"})
    display = ["driver_gene","paralog_gene","cancer_type","pcs","delta_expression",
               "necessity","dependency_dd","composite_score","novelty",
               "q_value","mutation_frequency","n_mut","n_wt","is_known_paralog_sl"]
    available = [c for c in display if c in results.columns]
    ranked = results.sort_values(["is_known_paralog_sl","composite_score"],
                                  ascending=[False,False])[available]
    ranked.to_csv(results_file, index=False)
    with open(summary_file, "w") as f:
        f.write("="*60+"\nParalog SL Predictor — Analysis Summary\n"+"="*60+"\n\n")
        f.write(f"Total pairs analyzed: {len(ranked)}\n")
        f.write(f"Known paralog-SL recovered: {ranked['is_known_paralog_sl'].sum()}\n\n")
        f.write("--- Top 20 ---\n\n")
        for i, (_, r) in enumerate(ranked.head(20).iterrows(), 1):
            flag = "★" if r["is_known_paralog_sl"] else "·"
            f.write(f"{i:2d}. {flag} {r['driver_gene']:10s} → {r['paralog_gene']:10s}  "
                    f"PCS={r['pcs']:.4f}  Score={r['composite_score']:.3f}  [{r['cancer_type']}]\n")
    return ranked


def main():
    p = argparse.ArgumentParser(description="Paralog SL Predictor v2")
    p.add_argument("--cancer", type=str, default=None, choices=["Ovarian","Endometrial","Cervical"])
    p.add_argument("--no-viz", action="store_true")
    args = p.parse_args()

    print("\n"+"█"*60+"\n  Paralog SL Predictor v2\n  + Negative Controls + Bootstrap + Cross-Cancer\n"+"█"*60+"\n")

    try:
        dep, expr, mod, mut, para, sldb = load_all_data()
    except FileNotFoundError as e:
        print(f"[ERROR] {e}"); sys.exit(1)

    cancer_types = [args.cancer] if args.cancer else None
    results = run_full_analysis(dep, expr, mod, mut, para, cancer_types=cancer_types)
    if results.empty:
        print("No results."); sys.exit(1)

    # ── Validation ──
    vr = run_full_validation(results)
    with open("output/validation_report.json","w") as f:
        json.dump({k:v for k,v in vr.items() if not isinstance(v,pd.DataFrame)}, f, indent=2, default=str)

    # ── Cross-cancer validation ──
    all_res = {}
    for ct in ["Ovarian","Endometrial","Cervical"]:
        r = run_full_analysis(dep, expr, mod, mut, para, cancer_types=[ct])
        if not r.empty: all_res[ct] = r
    cross_cancer_validation(all_res)

    # ── Export ──
    ranked = rank_and_export(results, RESULTS_FILE, SUMMARY_FILE)

    print(f"\n{'='*60}\nOutputs:\n  Results:  {RESULTS_FILE}\n  Summary:  {SUMMARY_FILE}\n  Validation: output/validation_report.json\n{'='*60}\n")

    top_novel = ranked[ranked["novelty"]=="Novel"].head(10)
    if not top_novel.empty:
        print("Top 10 Novel Candidate SL Pairs:")
        for i,(_,r) in enumerate(top_novel.iterrows(),1):
            print(f"  {i:2d}. {r['driver_gene']:10s} → {r['paralog_gene']:10s}  "
                  f"PCS={r['pcs']:.4f}  Score={r['composite_score']:.3f}  [{r['cancer_type']}]")

if __name__ == "__main__":
    main()
