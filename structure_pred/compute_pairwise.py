"""Compute pairwise structural similarity from ESMFold PDB files."""
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parent
PDB_DIR = BASE / "outputs" / "monomers"

def load_coords(pdb_path):
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and " CA " in line:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    return np.array(coords)

def compute_tm(coords_a, coords_b):
    n = min(len(coords_a), len(coords_b))
    d0 = 1.24 * max(n - 15, 1) ** (1/3) - 1.8
    d0 = max(d0, 0.5)
    diff = coords_a[:n] - coords_b[:n]
    rmsd = np.sqrt((diff**2).sum(axis=1).mean())
    tm = 1.0 / (1.0 + (rmsd / d0)**2)
    return float(tm), float(rmsd)

KEY_PAIRS = [
    ("KRAS", "HRAS"), ("PPP2R1A", "PPP2R1B"),
    ("PIK3CA", "PIK3CB"), ("STK11", "SIK1"),
    ("FBXW7", "FBXW2"), ("SMARCA4", "SMARCA2"),
    ("NF1", "RASA2"),
]

results = pd.read_csv(BASE / "outputs" / "esmfold_results.csv")

pairs = []
for ga, gb in KEY_PAIRS:
    pa = PDB_DIR / f"{ga}.pdb"
    pb = PDB_DIR / f"{gb}.pdb"
    if not pa.exists() or not pb.exists():
        continue
    
    coords_a = load_coords(pa)
    coords_b = load_coords(pb)
    tm, rmsd = compute_tm(coords_a, coords_b)
    
    ra = results[results.gene == ga].iloc[0]
    rb = results[results.gene == gb].iloc[0]
    
    pairs.append({
        "gene_a": ga, "gene_b": gb,
        "len_a_pred": ra["pred_length"], "len_b_pred": rb["pred_length"],
        "plddt_a": ra["mean_plddt"], "plddt_b": rb["mean_plddt"],
        "tm_score": round(tm, 4), "rmsd_ca": round(rmsd, 2),
        "plddt_mean": round((ra["mean_plddt"] + rb["mean_plddt"]) / 2, 3),
    })
    print(f"  {ga:12s}-{gb:12s}  TM={tm:.4f}  RMSD={rmsd:.1f}Å  pLDDT_mean={pairs[-1]['plddt_mean']:.3f}")

df = pd.DataFrame(pairs)
df.to_csv(BASE / "outputs" / "esmfold_pairwise.csv", index=False)
print(f"\nSaved: esmfold_pairwise.csv ({len(pairs)} pairs)")
