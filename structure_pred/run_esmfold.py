#!/usr/bin/env python3
"""
ESMFold Structure Prediction Pipeline for Paralog-SL Proteins
=============================================================
Predicts monomer structures via ESMFold (facebook/esmfold_v1) for all
22 paralog proteins, and linker-based dimer predictions for 5 key
paralog-SL pairs.

Outputs:
  - monomers/        PDB files for individual proteins
  - dimers/          PDB files for complex predictions (linker-based)
  - esmfold_results.csv  Summary metrics (pLDDT, tm_score, etc.)
  - esmfold_pairwise.csv Paralog pair structural similarity metrics

Requirements: transformers, torch, biopython, numpy, pandas
Model: facebook/esmfold_v1 (~3 GB download, ~6 GB RAM at inference)
Estimated time: 30-60 min for 22 monomers + 5 dimers on M4 Max
"""

import os, sys, json, time, csv
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

# ============================================================================
# Configuration
# ============================================================================
BASE = Path(__file__).resolve().parent
INPUT_FASTA = BASE / "inputs" / "all_paralog_proteins.fasta"
INPUT_CSV = BASE / "inputs" / "paralog_pairs.csv"
OUTPUT_MONOMERS = BASE / "outputs" / "monomers"
OUTPUT_DIMERS = BASE / "outputs" / "dimers"
RESULTS_CSV = BASE / "outputs" / "esmfold_results.csv"
PAIRWISE_CSV = BASE / "outputs" / "esmfold_pairwise.csv"

for d in [OUTPUT_MONOMERS, OUTPUT_DIMERS]:
    d.mkdir(parents=True, exist_ok=True)

# Key pairs for dimer prediction
KEY_PAIRS = [
    ("KRAS", "HRAS"),
    ("PPP2R1A", "PPP2R1B"),
    ("PIK3CA", "PIK3CB"),
    ("ARID1A", "ARID1B"),
    ("EP300", "CREBBP"),
]

# Poly-Gly linker for dimer prediction
LINKER = "G" * 25

# ============================================================================
# FASTA loading
# ============================================================================
def load_fasta(path):
    """Load FASTA file, return dict {gene_name: sequence}."""
    sequences = {}
    current_gene = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_gene:
                    sequences[current_gene] = "".join(current_seq)
                header = line[1:]
                current_gene = header.split("|")[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_gene:
            sequences[current_gene] = "".join(current_seq)
    return sequences

# ============================================================================
# ESMFold inference
# ============================================================================
def predict_structure(model, tokenizer, sequence, gene_name, output_dir, device):
    """
    Run ESMFold inference on a single protein sequence.
    Returns: dict with metrics
    """
    t0 = time.time()
    
    with torch.no_grad():
        # Tokenize
        tokenized = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.to(device) for k, v in tokenized.items()}
        
        # Forward pass
        output = model(**tokenized)
    
    # Extract pLDDT (per-residue confidence)
    plddt = output.plddt.cpu().numpy().flatten()
    mean_plddt = float(plddt.mean())
    
    # Extract predicted positions
    positions = output.positions[-1].cpu().numpy().flatten().reshape(-1, 3)
    
    # Save as PDB
    pdb_path = os.path.join(output_dir, f"{gene_name}.pdb")
    _write_pdb(positions, plddt, sequence, pdb_path)
    
    elapsed = time.time() - t0
    print(f"  {gene_name}: {len(sequence)} aa, pLDDT={mean_plddt:.3f}, {elapsed:.0f}s")
    
    return {
        "gene": gene_name,
        "length": len(sequence),
        "mean_plddt": mean_plddt,
        "time_seconds": elapsed,
        "pdb_path": str(pdb_path),
    }

def _write_pdb(positions, plddt, sequence, path):
    """Write predicted structure as PDB file (minimal format)."""
    with open(path, "w") as f:
        f.write(f"REMARK  ESMFold predicted structure\n")
        f.write(f"REMARK  Mean pLDDT: {plddt.mean():.3f}\n")
        atom_num = 0
        for i, (pos, conf, aa) in enumerate(zip(positions, plddt, sequence)):
            atom_num += 1
            res_num = i + 1
            x, y, z = pos[0], pos[1], pos[2]
            f.write(
                f"ATOM  {atom_num:5d}  CA  {aa:3s} A{res_num:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{conf:6.2f}           C  \n"
            )
        f.write("END\n")

# ============================================================================
# Structural similarity metrics
# ============================================================================
def compute_tm_score(coords_a, coords_b):
    """Simple TM-score approximation (length-normalized RMSD)."""
    n = min(len(coords_a), len(coords_b))
    d0 = 1.24 * (n - 15) ** (1/3) - 1.8
    diff = coords_a[:n] - coords_b[:n]
    rmsd = np.sqrt((diff ** 2).sum(axis=1).mean())
    tm = 1.0 / (1.0 + (rmsd / d0) ** 2) if d0 > 0 else 0.0
    return float(tm), float(rmsd)

def compute_pairwise_metrics(results, sequences):
    """Compute structural similarity between paralog pairs."""
    pairs = []
    for gene_a, gene_b in KEY_PAIRS:
        res_a = results.get(gene_a)
        res_b = results.get(gene_b)
        if not res_a or not res_b:
            continue
        
        # Load coordinates
        coords_a = _load_coords(res_a["pdb_path"])
        coords_b = _load_coords(res_b["pdb_path"])
        
        tm, rmsd = compute_tm_score(coords_a, coords_b)
        plddt_ratio = res_a["mean_plddt"] / max(res_b["mean_plddt"], 0.01)
        
        pairs.append({
            "gene_a": gene_a, "gene_b": gene_b,
            "len_a": res_a["length"], "len_b": res_b["length"],
            "plddt_a": res_a["mean_plddt"], "plddt_b": res_b["mean_plddt"],
            "tm_score": tm, "rmsd_ca": rmsd,
            "plddt_ratio": plddt_ratio,
        })
    return pd.DataFrame(pairs)

def _load_coords(pdb_path):
    """Extract CA atom coordinates from PDB file."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and " CA " in line:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    return np.array(coords)

# ============================================================================
# Main pipeline
# ============================================================================
def main():
    print("=" * 60)
    print("ESMFold Paralog-SL Structure Prediction Pipeline")
    print("=" * 60)
    
    # Load sequences
    sequences = load_fasta(INPUT_FASTA)
    print(f"\nLoaded {len(sequences)} protein sequences from FASTA")
    
    # Select proteins for monomer prediction (top paralog proteins, exclude very large)
    # Focus on all 22 proteins we have
    target_genes = list(sequences.keys())
    
    # Device selection: prefer MPS (Apple Silicon GPU)
    if torch.backends.mps.is_available():
        device = torch.device("mps")
        print(f"\nDevice: MPS (Apple Silicon GPU)")
    else:
        device = torch.device("cpu")
        print(f"\nDevice: CPU")
    
    # Load model
    print("Loading ESMFold model (facebook/esmfold_v1)...")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1",
        low_cpu_mem_usage=True,
    ).to(device)
    model.eval()
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    print("Model loaded.\n")
    
    # ================================================================
    # Phase 1: Monomer predictions
    # ================================================================
    print("=" * 60)
    print("Phase 1: Monomer Structure Predictions")
    print("=" * 60)
    
    results = {}
    for gene in target_genes:
        seq = sequences[gene]
        if len(seq) > 3000:
            print(f"  {gene}: SKIPPED ({len(seq)} aa, >3000 threshold)")
            continue
        
        try:
            # Use TRUNCATION for very long proteins (>2000aa) - focus on functional domains
            if len(seq) > 2000:
                print(f"  {gene}: TRUNCATING to first 2000 aa (total {len(seq)} aa)")
                seq = seq[:2000]  # Truncate for computational feasibility
            
            res = predict_structure(model, tokenizer, seq, gene, 
                                    OUTPUT_MONOMERS, device)
            results[gene] = res
        except Exception as e:
            print(f"  {gene}: FAILED - {e}")
    
    # Save monomer results
    df_results = pd.DataFrame(results.values())
    df_results.to_csv(RESULTS_CSV, index=False)
    print(f"\nMonomer results saved: {RESULTS_CSV}")
    print(f"  Successful: {len(results)}/{len(target_genes)} proteins")
    print(f"  Mean pLDDT: {df_results['mean_plddt'].mean():.3f}" if len(df_results) > 0 else "")
    
    # ================================================================
    # Phase 2: Dimer predictions (linker-based)
    # ================================================================
    print("\n" + "=" * 60)
    print("Phase 2: Dimer Structure Predictions (linker-based)")
    print("=" * 60)
    
    dimer_results = {}
    for gene_a, gene_b in KEY_PAIRS:
        if gene_a not in sequences or gene_b not in sequences:
            print(f"  {gene_a}-{gene_b}: SKIPPED (sequence missing)")
            continue
        
        seq_a = sequences[gene_a]
        seq_b = sequences[gene_b]
        total_len = len(seq_a) + len(seq_b) + len(LINKER)
        
        # Skip if too large
        if total_len > 3000:
            print(f"  {gene_a}-{gene_b}: SKIPPED ({total_len} aa total, >3000 threshold)")
            continue
        
        # Truncate large proteins
        if len(seq_a) > 1500:
            seq_a = seq_a[:1500]
        if len(seq_b) > 1500:
            seq_b = seq_b[:1500]
        
        combined_seq = seq_a + LINKER + seq_b
        pair_name = f"{gene_a}_{gene_b}"
        
        try:
            res = predict_structure(model, tokenizer, combined_seq, pair_name,
                                    OUTPUT_DIMERS, device)
            res["gene_a"] = gene_a
            res["gene_b"] = gene_b
            dimer_results[pair_name] = res
        except Exception as e:
            print(f"  {gene_a}-{gene_b}: FAILED - {e}")
    
    print(f"\nDimer predictions: {len(dimer_results)}/{len(KEY_PAIRS)} successful")
    
    # ================================================================
    # Phase 3: Pairwise metrics
    # ================================================================
    print("\n" + "=" * 60)
    print("Phase 3: Pairwise Structural Similarity")
    print("=" * 60)
    
    df_pairwise = compute_pairwise_metrics(results, sequences)
    df_pairwise.to_csv(PAIRWISE_CSV, index=False)
    print(f"Pairwise metrics saved: {PAIRWISE_CSV}")
    
    if len(df_pairwise) > 0:
        print("\nTop structurally similar pairs:")
        for _, row in df_pairwise.sort_values("tm_score", ascending=False).iterrows():
            print(f"  {row['gene_a']:12s} - {row['gene_b']:12s}  "
                  f"TM={row['tm_score']:.3f}  RMSD={row['rmsd_ca']:.1f}Å  "
                  f"pLDDT_a={row['plddt_a']:.3f}  pLDDT_b={row['plddt_b']:.3f}")
    
    print("\n" + "=" * 60)
    print("Pipeline complete.")
    print(f"  Monomers: {OUTPUT_MONOMERS}/ ({len(results)} PDB files)")
    print(f"  Dimers:   {OUTPUT_DIMERS}/ ({len(dimer_results)} PDB files)")
    print(f"  Results:  {RESULTS_CSV}")
    print(f"  Pairwise: {PAIRWISE_CSV}")
    print("=" * 60)

if __name__ == "__main__":
    main()
