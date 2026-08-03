#!/usr/bin/env python3
"""
cptac_rna_purity.py — Reviewer-response analyses for CPTAC paralog co-variation
================================================================================

Two analyses on top of the existing CPTAC protein-level results:

1. Matched mRNA co-variation: for each of the 7 CPTAC cohorts, fetch mRNA
   expression for the union of genes in the 26 analyzed paralog pairs and
   compute per-pair Pearson correlations of log-scale mRNA on samples matched
   to the proteomics samples. BH-correct within each pair across cohorts
   (same convention as the protein analysis).

2. Tumor-purity partial correlation: for each cohort fetch ESTIMATE_TUMORPURITY
   (fallback CIBERSORT_ABSOLUTE_SCORE) and recompute each pair's protein
   correlation controlling for purity (residualize each protein on purity,
   correlate residuals), alongside the unadjusted r recomputed on the SAME
   purity-matched subset.

Inputs (existing artifacts, read-only):
  - output/cptac_pair_matrix.csv          (the 26 analyzed pairs)
  - output/cptac_{cohort}_correlations.csv (full-sample protein r/p/n)
  - data/cptac_cache/{COHORT}_protein_data.json (cached protein matrices)

Outputs (NEW files):
  - output/cptac_rna_covariation.{json,csv}
  - output/cptac_purity_partial.{json,csv}
  - data/cptac_cache/{COHORT}_mrna_data.json
  - data/cptac_cache/{COHORT}_purity.json

All numbers come from live cBioPortal API responses or existing local
artifacts. No simulated data. Seed 42 fixed (nothing stochastic here, but set
for reproducibility).
"""

import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats
from scipy.stats import false_discovery_control

sys.path.insert(0, str(Path(__file__).resolve().parent))
from download_proteomics import CPTAC_STUDIES, GENE_ENTREZ, BASE_URL  # noqa: E402
from config import DATA_DIR, OUTPUT_DIR  # noqa: E402

np.random.seed(42)

SLEEP = 0.3          # polite delay between API calls
MIN_RNA_SAMPLES = 10     # same convention as protein analysis
MIN_PURITY_SAMPLES = 15
CACHE_DIR = DATA_DIR / "cptac_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

PURITY_PRIMARY = "ESTIMATE_TUMORPURITY"
PURITY_FALLBACK = "CIBERSORT_ABSOLUTE_SCORE"
# Task-specified chain (1–2) plus documented equivalents for cohorts whose
# studies name purity differently (3–5). Normalization ignores case, spaces,
# and underscores so e.g. "CIBERSORT_ABSOLUTE _SCORE" matches fallback 2.
PURITY_CHAIN = [
    "ESTIMATE_TUMORPURITY",
    "CIBERSORT_ABSOLUTE_SCORE",
    "TUMOR_PURITY_BYESTIMATE_RNASEQ",
    "TUMOR_PURITY",
    "TSNET_PURITY",
]


def _norm(s):
    return "".join(ch for ch in s.upper() if ch not in " _-")

# ── API helpers ────────────────────────────────────────────────


def api(method, url, params=None, json_body=None, timeout=90, retries=5):
    """GET/POST with retry + 0.3 s politeness sleep. Raises on hard failure."""
    last_err = None
    for attempt in range(retries):
        try:
            r = requests.request(
                method, url, params=params, json=json_body, timeout=timeout,
                headers={"Content-Type": "application/json"},
            )
            time.sleep(SLEEP)
            if r.ok:
                return r.json()
            last_err = f"HTTP {r.status_code}: {r.text[:200]}"
            if r.status_code in (400, 404):
                break
        except Exception as e:  # SSL hiccups etc.
            last_err = f"{type(e).__name__}: {str(e)[:150]}"
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"{method} {url} failed after {retries} tries: {last_err}")


# ── Gene set ───────────────────────────────────────────────────


def load_pairs():
    """The 26 analyzed pairs from the existing pair matrix artifact."""
    df = pd.read_csv(OUTPUT_DIR / "cptac_pair_matrix.csv")
    pairs = []
    for p in df["pair"]:
        a, b = p.split("↔")
        pairs.append((a.strip(), b.strip()))
    return pairs


# ── Per-cohort discovery ───────────────────────────────────────


def discover_rna_profile(study):
    """Return (profile_id, datatype, note).

    Preference order (per task spec, extended for heterogeneous CPTAC naming):
      1. {study}_rna_seq_v2_mrna            (continuous; task-specified preferred)
      2. {study}_rna_seq_mrna               (continuous)
      3. {study}_mrna                       (continuous)
      4. any other continuous MRNA_EXPRESSION (excluding microRNA '*_mirna')
      5. {study}_mrna_median_Zscores        (task-specified fallback)
      6. any other Z-SCORE MRNA_EXPRESSION  (excluding microRNA)
    """
    profiles = api("GET", f"{BASE_URL}/studies/{study}/molecular-profiles")
    mrna = [p for p in profiles if p["molecularAlterationType"] == "MRNA_EXPRESSION"]
    by_id = {p["molecularProfileId"]: p for p in mrna}
    notes = {
        f"{study}_rna_seq_v2_mrna": "preferred rna_seq_v2_mrna (continuous)",
        f"{study}_rna_seq_mrna": "rna_seq_mrna (continuous; naming variant)",
        f"{study}_mrna": "study_mrna (continuous; naming variant)",
    }
    for cand in (f"{study}_rna_seq_v2_mrna", f"{study}_rna_seq_mrna", f"{study}_mrna"):
        p = by_id.get(cand)
        if p and p["datatype"] != "Z-SCORE":
            return cand, p["datatype"], notes[cand]
    for p in mrna:  # any other continuous, no microRNA
        if p["datatype"] != "Z-SCORE" and "_mirna" not in p["molecularProfileId"]:
            return (p["molecularProfileId"], p["datatype"],
                    "other continuous MRNA_EXPRESSION profile")
    zscore = f"{study}_mrna_median_Zscores"
    if zscore in by_id:
        return zscore, by_id[zscore]["datatype"], "fallback mrna_median_Zscores"
    for p in mrna:  # any other z-score, no microRNA
        if p["datatype"] == "Z-SCORE" and "_mirna" not in p["molecularProfileId"]:
            return (p["molecularProfileId"], p["datatype"],
                    "other Z-SCORE MRNA_EXPRESSION profile")
    return None, None, "rna_unavailable"


def choose_mrna_sample_list(study, profile_id):
    """Pick the sample list covering mRNA-measured samples."""
    lists = api("GET", f"{BASE_URL}/studies/{study}/sample-lists")
    by_id = {l["sampleListId"]: l for l in lists}
    for cat in ("all_cases_with_mrna_rnaseq_data", "all_cases_with_mrna_data"):
        for l in lists:
            if l["category"] == cat:
                return l["sampleListId"]
    for cand in (f"{study}_rna_seq_mrna", f"{study}_mrna", profile_id, f"{study}_all"):
        if cand in by_id:
            return cand
    return f"{study}_all"


def fetch_mrna(profile_id, sample_list, gene_entrez):
    """One POST for all genes; returns {symbol: {sampleId: value}}."""
    ids = sorted(set(gene_entrez.values()))
    data = api("POST",
               f"{BASE_URL}/molecular-profiles/{profile_id}/molecular-data/fetch",
               json_body={"sampleListId": sample_list, "entrezGeneIds": ids})
    e2s = {v: k for k, v in gene_entrez.items()}
    out = {}
    for d in data:
        sym = e2s.get(d.get("entrezGeneId"))
        if sym is None:
            continue
        val = d.get("value", "")
        try:
            fval = float(val)
        except (TypeError, ValueError):
            continue
        if np.isnan(fval):
            continue
        out.setdefault(sym, {})[d["sampleId"]] = fval
    return out


def fetch_sample_patient_map(study):
    samples = api("GET", f"{BASE_URL}/studies/{study}/samples",
                  params={"projection": "SUMMARY", "pageSize": 10000000})
    return {s["sampleId"]: s["patientId"] for s in samples}


def discover_purity_attribute(study):
    """Return (attribute_id, level) or (None, None). level in {SAMPLE, PATIENT}.

    Walks PURITY_CHAIN with normalized matching; returns the raw attribute id
    (which may contain spaces, e.g. UCEC's "CIBERSORT_ABSOLUTE _SCORE")."""
    attrs = api("GET", f"{BASE_URL}/studies/{study}/clinical-attributes",
                params={"projection": "SUMMARY"})
    table = {a["clinicalAttributeId"]: a["patientAttribute"] for a in attrs}
    norm_map = {_norm(k): k for k in table}
    for cand in PURITY_CHAIN:
        raw = norm_map.get(_norm(cand))
        if raw is not None:
            return raw, ("PATIENT" if table[raw] else "SAMPLE")
    return None, None


def fetch_purity_values(study, attribute, level):
    """Fetch all values for one clinical attribute."""
    data = api("GET", f"{BASE_URL}/studies/{study}/clinical-data",
               params={"clinicalDataType": level, "attributeId": attribute,
                       "projection": "SUMMARY", "pageSize": 10000000})
    out = {}
    key = "patientId" if level == "PATIENT" else "sampleId"
    for d in data:
        try:
            out[d[key]] = float(d["value"])
        except (TypeError, ValueError, KeyError):
            continue
    return out


# ── Statistics ─────────────────────────────────────────────────


def pearson(x, y):
    r, p = stats.pearsonr(x, y)
    return float(r), float(p)


def partial_corr_purity(a, b, purity):
    """Partial Pearson r of a,b controlling for purity via residualization."""
    X = np.column_stack([np.ones(len(purity)), purity])
    ra = a - X @ np.linalg.lstsq(X, a, rcond=None)[0]
    rb = b - X @ np.linalg.lstsq(X, b, rcond=None)[0]
    return pearson(ra, rb)


def bh_within_pair(records, p_key, q_key):
    """BH-adjust p across cohorts within each pair (protein-analysis convention).
    records: list of dicts with 'pair' and p_key. Adds q_key in place."""
    by_pair = {}
    for i, rec in enumerate(records):
        p = rec.get(p_key)
        if p is not None and not (isinstance(p, float) and np.isnan(p)):
            by_pair.setdefault(rec["pair"], []).append(i)
    for _, idxs in by_pair.items():
        ps = np.array([records[i][p_key] for i in idxs], dtype=float)
        qs = false_discovery_control(ps, method="bh")
        for i, q in zip(idxs, qs):
            records[i][q_key] = float(q)


# ── Main ───────────────────────────────────────────────────────


def main():
    t0 = time.time()
    pairs = load_pairs()
    genes = sorted({g for pr in pairs for g in pr})
    gene_entrez = {}
    for g in genes:
        if g not in GENE_ENTREZ:
            raise KeyError(f"{g} missing from GENE_ENTREZ mapping")
        gene_entrez[g] = GENE_ENTREZ[g]
    print(f"{len(pairs)} pairs, {len(genes)} unique genes")

    # Full-sample protein results (existing artifacts) + cached protein matrices
    protein_csv = {}   # cohort -> {(a,b): row dict}
    protein_mat = {}   # cohort -> {gene: pd.Series}
    sanity = []
    for cohort in CPTAC_STUDIES:
        csv_path = OUTPUT_DIR / f"cptac_{cohort.lower()}_correlations.csv"
        df = pd.read_csv(csv_path)
        protein_csv[cohort] = {(r["gene_a"], r["gene_b"]): r for _, r in df.iterrows()}
        raw = json.load(open(CACHE_DIR / f"{cohort}_protein_data.json"))
        protein_mat[cohort] = {g: pd.Series(v, dtype=float) for g, v in raw.items() if v}
        # sanity: recompute full-sample protein r from cache, compare to CSV
        for (a, b), row in protein_csv[cohort].items():
            if row["status"] == "ok" and a in protein_mat[cohort] and b in protein_mat[cohort]:
                common = protein_mat[cohort][a].index.intersection(protein_mat[cohort][b].index)
                rr, _ = pearson(protein_mat[cohort][a][common], protein_mat[cohort][b][common])
                sanity.append(abs(rr - row["r"]))
    print(f"protein cache-vs-CSV max |Δr| sanity check: {max(sanity):.2e}")

    availability = {}
    rna_records = []      # pair×cohort mRNA tests
    purity_records = []   # pair×cohort purity partial tests

    for cohort, cfg in CPTAC_STUDIES.items():
        study = cfg["study"]
        print(f"\n{'─'*60}\n{cohort} ({study})")
        info = {"study": study}

        # ── mRNA profile + sample list ──
        profile_id, datatype, note = discover_rna_profile(study)
        info["rna_profile"] = profile_id
        info["rna_datatype"] = datatype
        info["rna_profile_note"] = note
        print(f"  mRNA profile: {profile_id or 'UNAVAILABLE'} ({note})")

        mrna = {}
        log_note = None
        if profile_id:
            sample_list = choose_mrna_sample_list(study, profile_id)
            info["rna_sample_list"] = sample_list
            mrna = fetch_mrna(profile_id, sample_list, gene_entrez)
            n_vals = sum(len(v) for v in mrna.values())
            all_vals = np.array([x for v in mrna.values() for x in v.values()]) if n_vals else np.array([0.0])
            if datatype == "Z-SCORE":
                log_note = "z-scores used as-is"
            elif np.nanmax(np.abs(all_vals)) > 30:
                mrna = {g: {s: float(np.log2(x + 1.0)) for s, x in v.items()}
                        for g, v in mrna.items()}
                log_note = "log2(x+1) applied to raw continuous values"
            else:
                log_note = "values already log-scale; used as-is"
            info["mrna_transform"] = log_note
            info["mrna_genes_fetched"] = len(mrna)
            info["mrna_samples_per_gene_median"] = (
                float(np.median([len(v) for v in mrna.values()])) if mrna else 0)
            print(f"  mRNA genes: {len(mrna)}/{len(genes)}  [{log_note}]")
            # cache raw matrix
            with open(CACHE_DIR / f"{cohort}_mrna_data.json", "w") as f:
                json.dump({
                    "cohort": cohort, "study": study, "profile": profile_id,
                    "datatype": datatype, "sample_list": sample_list,
                    "transform": log_note,
                    "fetched_utc": datetime.now(timezone.utc).isoformat(),
                    "values": mrna,
                }, f)
        else:
            info["rna_sample_list"] = None

        # ── purity ──
        purity_attr, purity_level = discover_purity_attribute(study)
        info["purity_attribute"] = purity_attr
        info["purity_level"] = purity_level
        sample_purity = {}
        if purity_attr:
            sp_map = fetch_sample_patient_map(study)
            vals = fetch_purity_values(study, purity_attr, purity_level)
            if purity_level == "SAMPLE":
                sample_purity = vals
            else:  # PATIENT-level: map to each of the patient's samples
                for sid, pid in sp_map.items():
                    if pid in vals:
                        sample_purity[sid] = vals[pid]
            info["purity_samples_with_value"] = len(sample_purity)
            with open(CACHE_DIR / f"{cohort}_purity.json", "w") as f:
                json.dump({
                    "cohort": cohort, "study": study,
                    "attribute": purity_attr, "level": purity_level,
                    "fetched_utc": datetime.now(timezone.utc).isoformat(),
                    "raw_values": vals,
                    "sample_values": sample_purity,
                }, f)
            print(f"  purity: {purity_attr} ({purity_level}), "
                  f"{len(sample_purity)} samples with values")
        else:
            info["purity_samples_with_value"] = 0
            with open(CACHE_DIR / f"{cohort}_purity.json", "w") as f:
                json.dump({"cohort": cohort, "study": study, "attribute": None,
                           "fetched_utc": datetime.now(timezone.utc).isoformat(),
                           "raw_values": {}, "sample_values": {}}, f)
            print("  purity: UNAVAILABLE")

        availability[cohort] = info
        prot = protein_mat[cohort]

        # ── per-pair tests ──
        for a, b in pairs:
            pair_key = f"{a}↔{b}"
            csv_row = protein_csv[cohort].get((a, b))
            prot_full_r = float(csv_row["r"]) if csv_row is not None and csv_row["status"] == "ok" else None
            prot_full_p = float(csv_row["p"]) if csv_row is not None and csv_row["status"] == "ok" else None
            prot_full_n = int(csv_row["n"]) if csv_row is not None and csv_row["status"] == "ok" else None

            # (1) mRNA co-variation on protein-matched samples
            if mrna and a in mrna and b in mrna:
                common_m = set(mrna[a]) & set(mrna[b])
                if a in prot and b in prot:
                    matched = sorted(common_m & set(prot[a].index) & set(prot[b].index))
                    mode = "protein_matched"
                else:
                    matched = sorted(common_m)
                    mode = "mrna_only_pair_protein_missing"
                if len(matched) >= MIN_RNA_SAMPLES:
                    xa = np.array([mrna[a][s] for s in matched])
                    xb = np.array([mrna[b][s] for s in matched])
                    mr, mp = pearson(xa, xb)
                    rec = {"cohort": cohort, "pair": pair_key, "gene_a": a, "gene_b": b,
                           "rna_profile": profile_id, "rna_datatype": datatype,
                           "transform": log_note, "match_mode": mode,
                           "n_mrna": len(matched), "mrna_r": mr, "mrna_p": mp,
                           "protein_n_full": prot_full_n,
                           "protein_r_full": prot_full_r, "protein_p_full": prot_full_p}
                    if mode == "protein_matched":
                        pr, pp = pearson(np.array([prot[a][s] for s in matched]),
                                         np.array([prot[b][s] for s in matched]))
                        rec["protein_r_matched"] = pr
                        rec["protein_p_matched"] = pp
                    else:
                        rec["protein_r_matched"] = None
                        rec["protein_p_matched"] = None
                    rna_records.append(rec)

            # (2) purity partial correlation (protein level)
            if sample_purity and a in prot and b in prot:
                matched = sorted(set(prot[a].index) & set(prot[b].index) & set(sample_purity))
                if len(matched) >= MIN_PURITY_SAMPLES:
                    xa = np.array([prot[a][s] for s in matched])
                    xb = np.array([prot[b][s] for s in matched])
                    pur = np.array([sample_purity[s] for s in matched])
                    ur, up = pearson(xa, xb)
                    pr_, pp_ = partial_corr_purity(xa, xb, pur)
                    purity_records.append({
                        "cohort": cohort, "pair": pair_key, "gene_a": a, "gene_b": b,
                        "purity_attribute": purity_attr, "purity_level": purity_level,
                        "n_matched": len(matched),
                        "protein_r_full": prot_full_r, "protein_p_full": prot_full_p,
                        "protein_n_full": prot_full_n,
                        "protein_r_matched": ur, "protein_p_matched": up,
                        "partial_r": pr_, "partial_p": pp_,
                        "delta_r": pr_ - ur,
                        "abs_delta_r": abs(pr_ - ur),
                    })

    print(f"\nmRNA tests: {len(rna_records)}; purity tests: {len(purity_records)}")

    # ── BH within pair across cohorts ──
    bh_within_pair(rna_records, "mrna_p", "mrna_q_bh")
    bh_within_pair(rna_records, "protein_p_full", "protein_q_bh_full")
    bh_within_pair(rna_records, "protein_p_matched", "protein_q_bh_matched")
    bh_within_pair(purity_records, "protein_p_full", "protein_q_bh_full")
    bh_within_pair(purity_records, "protein_p_matched", "protein_q_bh_matched")
    bh_within_pair(purity_records, "partial_p", "partial_q_bh")

    # ── Summary metrics ──
    summary = {"per_cohort": {}, "overall": {}}
    for cohort in CPTAC_STUDIES:
        recs = [r for r in rna_records if r["cohort"] == cohort]
        prec = [r for r in purity_records if r["cohort"] == cohort]
        paired = [r for r in recs if r.get("protein_r_matched") is not None]
        entry = {
            "n_pairs_mrna_tested": len(recs),
            "median_mrna_r": float(np.median([r["mrna_r"] for r in recs])) if recs else None,
            "median_protein_r_full_sample": (
                float(np.median([r["protein_r_full"] for r in recs if r["protein_r_full"] is not None]))
                if any(r["protein_r_full"] is not None for r in recs) else None),
            "median_protein_r_matched": (
                float(np.median([r["protein_r_matched"] for r in paired])) if paired else None),
            "n_pairs_purity_tested": len(prec),
            "median_abs_delta_r_purity": (
                float(np.median([r["abs_delta_r"] for r in prec])) if prec else None),
        }
        summary["per_cohort"][cohort] = entry

    all_paired = [r for r in rna_records if r.get("protein_r_matched") is not None]
    summary["overall"] = {
        "n_pair_cohort_mrna_tests": len(rna_records),
        "n_pair_cohort_with_matched_protein": len(all_paired),
        "median_protein_r_matched": (
            float(np.median([r["protein_r_matched"] for r in all_paired])) if all_paired else None),
        "median_mrna_r": float(np.median([r["mrna_r"] for r in all_paired])) if all_paired else None,
        "median_abs_protein_r_matched": (
            float(np.median([abs(r["protein_r_matched"]) for r in all_paired])) if all_paired else None),
        "median_abs_mrna_r": (
            float(np.median([abs(r["mrna_r"]) for r in all_paired])) if all_paired else None),
        "median_abs_delta_r_purity": (
            float(np.median([r["abs_delta_r"] for r in purity_records])) if purity_records else None),
        "n_pair_cohort_purity_tests": len(purity_records),
    }

    # ── EP300–CREBBP per cohort ──
    ep = {}
    for cohort in CPTAC_STUDIES:
        r_rec = next((r for r in rna_records if r["cohort"] == cohort and r["pair"] == "EP300↔CREBBP"), None)
        p_rec = next((r for r in purity_records if r["cohort"] == cohort and r["pair"] == "EP300↔CREBBP"), None)
        ep[cohort] = {
            "protein_r_full": (r_rec or p_rec or {}).get("protein_r_full"),
            "mrna_r": r_rec["mrna_r"] if r_rec else None,
            "mrna_q_bh": r_rec.get("mrna_q_bh") if r_rec else None,
            "n_mrna": r_rec["n_mrna"] if r_rec else None,
            "purity_adjusted_protein_r": p_rec["partial_r"] if p_rec else None,
            "purity_adjusted_q_bh": p_rec.get("partial_q_bh") if p_rec else None,
            "n_purity_matched": p_rec["n_matched"] if p_rec else None,
        }

    # ── Sign test: protein vs mRNA magnitude (matched samples) ──
    diffs_mag, diffs_signed = [], []
    for r in all_paired:
        diffs_mag.append(abs(r["protein_r_matched"]) - abs(r["mrna_r"]))
        diffs_signed.append(r["protein_r_matched"] - r["mrna_r"])

    def sign_test(diffs):
        pos = sum(1 for d in diffs if d > 0)
        neg = sum(1 for d in diffs if d < 0)
        n = pos + neg
        p = float(stats.binomtest(pos, n, 0.5).pvalue) if n else None
        return {"n": n, "n_positive": pos, "n_negative": neg, "p_two_sided": p}

    st_mag = sign_test(diffs_mag)
    st_signed = sign_test(diffs_signed)

    # ── Conclusion crossings under purity adjustment ──
    def crossings(before_key, after_key="partial_q_bh"):
        lost, gained = [], []
        for r in purity_records:
            qb, qa = r.get(before_key), r.get(after_key)
            if qb is None or qa is None:
                continue
            base = {"cohort": r["cohort"], "pair": r["pair"],
                    "q_before": qb, "q_after": qa,
                    "r_full_sample": r["protein_r_full"],
                    "n_full_sample": r["protein_n_full"],
                    "r_matched_unadjusted": r["protein_r_matched"],
                    "r_partial": r["partial_r"],
                    "n_matched": r["n_matched"]}
            if qb < 0.05 <= qa:
                lost.append(base)
            elif qa < 0.05 <= qb:
                gained.append(base)
        return lost, gained

    lost_same, gained_same = crossings("protein_q_bh_matched")
    lost_full, gained_full = crossings("protein_q_bh_full")

    conclusion_flags = {
        "protein_stronger_than_mrna": {
            "sign_test_abs_r": st_mag,
            "sign_test_signed_r": st_signed,
            "flag": bool(st_mag["p_two_sided"] is not None
                         and st_mag["p_two_sided"] < 0.05
                         and st_mag["n_positive"] > st_mag["n_negative"]),
        },
        "purity_changes_conclusions": {
            "same_sample_comparison": {
                "lost_significance": lost_same, "gained_significance": gained_same,
                "flag": bool(lost_same or gained_same),
            },
            "full_sample_comparison": {
                "lost_significance": lost_full, "gained_significance": gained_full,
                "flag": bool(lost_full or gained_full),
            },
        },
    }

    meta = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "api_base": BASE_URL,
        "seed": 42,
        "api_sleep_seconds": SLEEP,
        "min_rna_samples": MIN_RNA_SAMPLES,
        "min_purity_samples": MIN_PURITY_SAMPLES,
        "n_pairs": len(pairs),
        "n_genes": len(genes),
        "genes": genes,
        "bh_convention": "BH within each pair across the 7 cohorts (same as protein analysis)",
        "purity_attribute_preference": PURITY_CHAIN,
        "deviations": [
            "mRNA profile discovery extended beyond the two task-specified ids: "
            "continuous {study}_rna_seq_mrna / {study}_mrna and other non-microRNA "
            "MRNA_EXPRESSION profiles are used when {study}_rna_seq_v2_mrna is absent; "
            "Z-score profiles are last-resort. Profile used is recorded per cohort.",
            "Purity attribute chain extended (normalized matching) to "
            "TUMOR_PURITY_BYESTIMATE_RNASEQ / TUMOR_PURITY / TSNET_PURITY for studies "
            "lacking ESTIMATE_TUMORPURITY and CIBERSORT_ABSOLUTE_SCORE; attribute used "
            "is recorded per cohort.",
        ],
        "protein_cache_vs_csv_max_abs_dr": float(max(sanity)),
        "runtime_seconds": round(time.time() - t0, 1),
    }

    # ── Write outputs ──
    rna_json = {
        "meta": meta, "availability": availability,
        "summary": summary, "ep300_crebbp_per_cohort": ep,
        "conclusion_flags": conclusion_flags,
        "tests": rna_records,
    }
    with open(OUTPUT_DIR / "cptac_rna_covariation.json", "w") as f:
        json.dump(rna_json, f, indent=2)
    pd.DataFrame(rna_records).to_csv(OUTPUT_DIR / "cptac_rna_covariation.csv", index=False)

    purity_json = {
        "meta": meta, "availability": availability,
        "summary": summary, "ep300_crebbp_per_cohort": ep,
        "conclusion_flags": conclusion_flags,
        "tests": purity_records,
    }
    with open(OUTPUT_DIR / "cptac_purity_partial.json", "w") as f:
        json.dump(purity_json, f, indent=2)
    pd.DataFrame(purity_records).to_csv(OUTPUT_DIR / "cptac_purity_partial.csv", index=False)

    # ── Console summary ──
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for c in CPTAC_STUDIES:
        a = availability[c]
        s = summary["per_cohort"][c]
        print(f"  {c:5s} rna={a['rna_profile'] or 'NONE':36s} "
              f"purity={a['purity_attribute'] or 'NONE':26s} "
              f"medProteinR={s['median_protein_r_matched']} medMrnaR={s['median_mrna_r']}")
    ov = summary["overall"]
    print(f"\n  overall median protein r (matched) = {ov['median_protein_r_matched']:.4f}")
    print(f"  overall median mRNA r              = {ov['median_mrna_r']:.4f}")
    print(f"  sign test |r|: {st_mag['n_positive']}/{st_mag['n']} protein>mRNA, p={st_mag['p_two_sided']:.3e}")
    print(f"  median |Δr| purity adjustment      = {ov['median_abs_delta_r_purity']:.4f}")
    print(f"  purity crossings (same-sample): lost={len(lost_same)}, gained={len(gained_same)}")
    print(f"  purity crossings (full-sample):   lost={len(lost_full)}, gained={len(gained_full)}")
    print("\n  EP300↔CREBBP per cohort:")
    for c, v in ep.items():
        print(f"    {c:5s} protein_r={v['protein_r_full']} mrna_r={v['mrna_r']} "
              f"purity_adj_r={v['purity_adjusted_protein_r']}")
    print(f"\nDone in {meta['runtime_seconds']} s")


if __name__ == "__main__":
    main()
