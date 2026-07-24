#!/usr/bin/env python3
"""
Fetch pre-computed AlphaFold2 structures from AlphaFold Protein Structure DB.
URL: https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-model_v4.pdb

This uses the official DeepMind AF2 predictions — no local model needed.
Reference: Varadi et al., NAR 2024 (AlphaFold Protein Structure Database)
"""

import os, json, csv, time
from pathlib import Path
import requests
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
AF_DIR = BASE / "outputs" / "af2_structures"
AF_DIR.mkdir(parents=True, exist_ok=True)

# UniProt IDs for paralog proteins
UNIPROT_IDS = {
    "ARID1A": "O14497", "ARID1B": "Q8NFD5",
    "EP300": "Q09472", "CREBBP": "Q92793",
    "PIK3CA": "P42336", "PIK3CB": "P42338",
    "PPP2R1A": "P30153", "PPP2R1B": "P30154",
    "KRAS": "P01116", "HRAS": "P01112",
    "BRCA1": "P38398", "BRCA2": "P51587",
    "SMARCA4": "P51532", "SMARCA2": "P51531",
    "NF1": "P21359", "RASA2": "Q15283",
    "KMT2D": "O14686", "KMT2C": "Q8NEZ4",
    "STK11": "Q15831", "SIK1": "P57059",
    "FBXW7": "Q969H0", "FBXW2": "Q9UKT8",
    "PTEN": "P60484", "PIK3R1": "P27986",
}

KEY_PAIRS = [
    ("KRAS", "HRAS"),
    ("PPP2R1A", "PPP2R1B"),
    ("PIK3CA", "PIK3CB"),
    ("ARID1A", "ARID1B"),
    ("EP300", "CREBBP"),
]

AF_URL = "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb"

def fetch_af2_structure(uniprot_id, gene_name):
    """Download pre-computed AF2 structure."""
    url = AF_URL.format(uid=uniprot_id)
    out_path = AF_DIR / f"{gene_name}.pdb"
    
    if out_path.exists():
        print(f"  {gene_name}: already downloaded")
        return str(out_path)
    
    try:
        r = requests.get(url, timeout=120, stream=True)
        if r.status_code == 404:
            # Try v3 (older format)
            url_v3 = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v3.pdb"
            r = requests.get(url_v3, timeout=120, stream=True)
        
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        
        size_mb = out_path.stat().st_size / 1e6
        print(f"  {gene_name}: downloaded ({size_mb:.1f} MB)")
        return str(out_path)
    except Exception as e:
        print(f"  {gene_name}: FAILED - {e}")
        return None

def parse_pdb_metrics(pdb_path):
    """Extract pLDDT and other metrics from AF2 PDB file."""
    plddt_values = []
    ca_coords = []
    atom_lines = 0
    
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and " CA " in line:
                atom_lines += 1
                bfactor = float(line[60:66])
                plddt_values.append(bfactor)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ca_coords.append([x, y, z])
    
    if not plddt_values:
        return None
    
    plddt = np.array(plddt_values)
    return {
        "n_residues": len(plddt),
        "mean_plddt": float(plddt.mean()),
        "median_plddt": float(np.median(plddt)),
        "plddt_gt_90": float((plddt > 90).mean()),
        "plddt_gt_70": float((plddt > 70).mean()),
        "plddt_lt_50": float((plddt < 50).mean()),
        "coords": np.array(ca_coords),
    }

def compute_pairwise_tm(coords_a, coords_b):
    """TM-score approximation between two structures."""
    n = min(len(coords_a), len(coords_b))
    d0 = 1.24 * max(n - 15, 1) ** (1/3) - 1.8
    d0 = max(d0, 0.5)
    diff = coords_a[:n] - coords_b[:n]
    rmsd = np.sqrt((diff ** 2).sum(axis=1).mean())
    tm = 1.0 / (1.0 + (rmsd / d0) ** 2)
    return float(tm), float(rmsd)

def compute_contact_map(coords, threshold=8.0):
    """Compute binary contact map (Cα-Cα distance < threshold)."""
    n = len(coords)
    # Use a subset for very large proteins
    if n > 1500:
        idx = np.linspace(0, n-1, 1500, dtype=int)
        coords = coords[idx]
        n = 1500
    
    diff = coords[None, :, :] - coords[:, None, :]
    dist = np.sqrt((diff ** 2).sum(axis=2))
    contacts = (dist < threshold).astype(int)
    # Remove diagonal and i,i+1 (sequential)
    for i in range(n):
        contacts[i, max(0,i-1):min(n,i+2)] = 0
    return contacts

def compute_contact_similarity(cm_a, cm_b):
    """Jaccard similarity of contact maps."""
    n = min(cm_a.shape[0], cm_b.shape[0])
    intersection = (cm_a[:n,:n] & cm_b[:n,:n]).sum()
    union = (cm_a[:n,:n] | cm_b[:n,:n]).sum()
    return float(intersection / union) if union > 0 else 0.0

# ============================================================================
# Main
# ============================================================================
print("=" * 60)
print("Fetching AlphaFold2 Pre-computed Structures")
print("=" * 60)

# Phase 1: Download structures
results = {}
for gene, uid in UNIPROT_IDS.items():
    path = fetch_af2_structure(uid, gene)
    if path:
        metrics = parse_pdb_metrics(path)
        if metrics:
            results[gene] = {**metrics, "pdb_path": path}
            print(f"    pLDDT={metrics['mean_plddt']:.1f}  n={metrics['n_residues']}")

print(f"\nDownloaded: {len(results)}/{len(UNIPROT_IDS)} structures")

# Save summary
df = pd.DataFrame([
    {"gene": g, "uniprot": UNIPROT_IDS.get(g, ""), **{k:v for k,v in r.items() if k != "coords"}}
    for g, r in results.items()
])
df.to_csv(AF_DIR / "af2_summary.csv", index=False)
print(f"Summary: {AF_DIR}/af2_summary.csv")

# Phase 2: Pairwise structural similarity
print("\n" + "=" * 60)
print("Pairwise Structural Similarity (AF2)")
print("=" * 60)

pairs = []
for gene_a, gene_b in KEY_PAIRS:
    if gene_a not in results or gene_b not in results:
        print(f"  {gene_a}-{gene_b}: SKIP (structure missing)")
        continue
    
    ra, rb = results[gene_a], results[gene_b]
    tm, rmsd = compute_pairwise_tm(ra["coords"], rb["coords"])
    cm_a = compute_contact_map(ra["coords"])
    cm_b = compute_contact_map(rb["coords"])
    contact_sim = compute_contact_similarity(cm_a, cm_b)
    
    # Domain-level pLDDT comparison (first 200 aa vs first 200 aa as rough domain proxy)
    n_domain = min(200, ra["n_residues"], rb["n_residues"])
    domain_tm, domain_rmsd = compute_pairwise_tm(
        ra["coords"][:n_domain], rb["coords"][:n_domain]
    )
    
    pairs.append({
        "gene_a": gene_a, "gene_b": gene_b,
        "len_a": ra["n_residues"], "len_b": rb["n_residues"],
        "plddt_a": ra["mean_plddt"], "plddt_b": rb["mean_plddt"],
        "tm_score": tm, "rmsd_ca": rmsd,
        "domain_tm_score": domain_tm, "domain_rmsd": domain_rmsd,
        "contact_similarity": contact_sim,
        "plddt_ratio": ra["mean_plddt"] / max(rb["mean_plddt"], 0.01),
    })
    
    print(f"  {gene_a:12s} - {gene_b:12s}  "
          f"TM={tm:.3f}  RMSD={rmsd:.1f}Å  "
          f"domain_TM={domain_tm:.3f}  contact_sim={contact_sim:.3f}  "
          f"pLDDT_a={ra['mean_plddt']:.1f}  pLDDT_b={rb['mean_plddt']:.1f}")

df_pairs = pd.DataFrame(pairs)
df_pairs.to_csv(AF_DIR / "af2_pairwise.csv", index=False)
print(f"\nPairwise metrics: {AF_DIR}/af2_pairwise.csv")

print("\n" + "=" * 60)
print("Complete. All structures from AlphaFold DB (Varadi et al., 2024).")
print("=" * 60)
