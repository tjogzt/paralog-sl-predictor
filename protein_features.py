"""
Direction C: Protein Feature Computation from FASTA + Half-life Data
=====================================================================
Computes protein-level features from the human proteome FASTA file
(no API calls needed after initial download).

Features computed:
  Sequence-based (from ProtParam):
    - length, molecular_weight, isoelectric_point
    - gravy (hydrophobicity), instability_index
    - lysine_count (ubiquitination potential)
    - cysteine_count (disulfide bond potential)
    - aromaticity, aliphatic_index

  Half-life (from Mathieson et al. 2018):
    - HeLa cell half-life (hours)

Usage: python3.10 protein_features.py
"""

import gzip
import re
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

DATA_DIR = Path(__file__).resolve().parent / "data"
FASTA_FILE = "/tmp/human_proteome.fasta.gz"
HALFLIFE_FILE = "/tmp/halflife_data.xlsx"
OUTPUT = DATA_DIR / "protein_features.csv"

# ── 1. Parse human proteome FASTA ─────────────────────────────

def parse_fasta(fasta_path: str) -> dict:
    """Parse UniProt reference proteome FASTA.
    Returns: {gene_symbol: sequence}
    """
    print("Parsing human proteome FASTA...")
    
    gene_to_seq = {}
    current_gene = None
    current_seq = []
    
    with gzip.open(fasta_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                # Save previous
                if current_gene:
                    gene_to_seq[current_gene] = "".join(current_seq)
                
                # Parse gene symbol from header
                # Format: >sp|P12345|GENE_HUMAN ... GN=GENE ...
                header = line[1:].strip()
                gn_match = re.search(r"GN=(\S+)", header)
                if gn_match:
                    current_gene = gn_match.group(1).upper()
                else:
                    # Try to extract from ID
                    parts = header.split("|")
                    if len(parts) >= 3:
                        current_gene = parts[2].split("_")[0].upper()
                    else:
                        current_gene = None
                current_seq = []
            elif current_gene:
                current_seq.append(line.strip())
        
        if current_gene:
            gene_to_seq[current_gene] = "".join(current_seq)
    
    # Take longest isoform per gene
    gene_seqs = {}
    for gene, seq in gene_to_seq.items():
        if gene not in gene_seqs or len(seq) > len(gene_seqs[gene]):
            gene_seqs[gene] = seq
    
    print(f"  Parsed {len(gene_seqs)} unique genes")
    return gene_seqs


# ── 2. Compute sequence-based features ─────────────────────────

def compute_sequence_features(sequence: str) -> dict:
    """Compute protein features from amino acid sequence."""
    seq = sequence.upper()
    length = len(seq)
    
    if length == 0:
        return {}
    
    # Amino acid composition
    aa_counts = {}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        aa_counts[aa] = seq.count(aa)
    
    # Lysine count (ubiquitination potential)
    lysine_count = seq.count("K")
    lysine_density = lysine_count / length
    
    # Cysteine count (structural stability / disulfide bonds)
    cysteine_count = seq.count("C")
    cysteine_density = cysteine_count / length
    
    # Molecular weight (approximate)
    aa_weights = {
        "A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13, "F": 165.19,
        "G": 75.07, "H": 155.16, "I": 131.17, "K": 146.19, "L": 131.17,
        "M": 149.21, "N": 132.12, "P": 115.13, "Q": 146.15, "R": 174.20,
        "S": 105.09, "T": 119.12, "V": 117.15, "W": 204.23, "Y": 181.19,
    }
    mol_weight = sum(aa_weights.get(aa, 0) * count for aa, count in aa_counts.items())
    mol_weight += 18.02  # water
    
    # Isoelectric point (approximate pI)
    # pKa values for charged residues
    pKa = {"D": 3.9, "E": 4.3, "H": 6.0, "C": 8.3, "Y": 10.1, "K": 10.5, "R": 12.5}
    # Approximate pI using Henderson-Hasselbalch iteration
    # Simplified: use net charge at different pH
    def net_charge(ph):
        charge = 0
        charge += aa_counts.get("D", 0) * (-1 / (1 + 10**(pKa["D"] - ph)))
        charge += aa_counts.get("E", 0) * (-1 / (1 + 10**(pKa["E"] - ph)))
        charge += aa_counts.get("C", 0) * (-1 / (1 + 10**(pKa["C"] - ph)))
        charge += aa_counts.get("Y", 0) * (-1 / (1 + 10**(pKa["Y"] - ph)))
        charge += aa_counts.get("H", 0) * (1 / (1 + 10**(ph - pKa["H"])))
        charge += aa_counts.get("K", 0) * (1 / (1 + 10**(ph - pKa["K"])))
        charge += aa_counts.get("R", 0) * (1 / (1 + 10**(ph - pKa["R"])))
        charge += 1  # N-terminus
        charge -= 1  # C-terminus
        return charge
    
    # Binary search for pI (where net_charge ≈ 0)
    lo, hi = 0.0, 14.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if net_charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    pI = (lo + hi) / 2
    
    # GRAVY (hydrophobicity)
    hydropathy = {
        "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
        "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
        "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
        "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
    }
    gravy = sum(hydropathy.get(aa, 0) * count for aa, count in aa_counts.items()) / length
    
    # Aromaticity
    aromatic = (aa_counts.get("F", 0) + aa_counts.get("Y", 0) + aa_counts.get("W", 0)) / length
    
    # Aliphatic index
    aliphatic = (aa_counts.get("A", 0) + 2.9 * aa_counts.get("V", 0) + 
                 3.9 * aa_counts.get("I", 0) + 3.9 * aa_counts.get("L", 0)) / length * 100
    
    # Disorder propensity (simple metric: fraction of disorder-promoting residues)
    disorder_promoting = {"P", "E", "S", "Q", "K", "A", "G"}
    disorder_ratio = sum(aa_counts.get(aa, 0) for aa in disorder_promoting) / length
    
    # Instability index (simplified - counts certain dipeptides)
    # Full calculation requires dipeptide weight table; simplified version:
    instable_dp = {"DP", "GT", "NN", "PG", "PP", "QQ", "RS", "TV"}
    dipeptide_count = 0
    for i in range(length - 1):
        if seq[i:i+2] in instable_dp:
            dipeptide_count += 1
    instability = (dipeptide_count / max(length - 1, 1)) * 100  # simplified
    
    return {
        "length": length,
        "molecular_weight": round(mol_weight, 1),
        "isoelectric_point": round(pI, 2),
        "gravy": round(gravy, 3),
        "aromaticity": round(aromatic, 4),
        "aliphatic_index": round(aliphatic, 1),
        "instability_index": round(instability, 1),
        "lysine_count": lysine_count,
        "lysine_density": round(lysine_density, 5),
        "cysteine_count": cysteine_count,
        "cysteine_density": round(cysteine_density, 5),
        "disorder_propensity": round(disorder_ratio, 4),
    }


# ── 3. Load half-life data ─────────────────────────────────────

def load_halflife_data(path: str) -> dict:
    """Load protein half-life data from Mathieson et al. 2018."""
    if not Path(path).exists():
        print("  Half-life file not found. Skipping half-life feature.")
        return {}
    
    try:
        df = pd.read_excel(path)
        print(f"  Half-life data: {len(df)} rows, cols: {list(df.columns)[:5]}")
        
        # Try to find gene symbol and half-life columns
        hl_map = {}
        
        # Look for columns with gene/name info
        gene_col = None
        for c in df.columns:
            if "gene" in str(c).lower() or "symbol" in str(c).lower() or "uniprot" in str(c).lower():
                gene_col = c
                break
        
        # Look for half-life column
        hl_col = None
        for c in df.columns:
            if "half" in str(c).lower() or "t1/2" in str(c).lower() or "turnover" in str(c).lower():
                hl_col = c
                break
        
        if gene_col and hl_col:
            for _, row in df.iterrows():
                gene = str(row[gene_col]).upper().strip()
                try:
                    hl = float(row[hl_col])
                    hl_map[gene] = hl
                except (ValueError, TypeError):
                    continue
        
        print(f"  Mapped {len(hl_map)} genes with half-life data")
        return hl_map
    except Exception as e:
        print(f"  Error loading half-life: {e}")
        return {}


# ── 4. Main ────────────────────────────────────────────────────

def main():
    # Parse proteome
    gene_seqs = parse_fasta(FASTA_FILE)
    
    # Collect genes of interest
    from config import DRIVER_GENES, KNOWN_PARALOG_SL
    
    all_genes = set()
    for genes in DRIVER_GENES.values():
        all_genes.update(g.upper() for g in genes)
    for a, b in KNOWN_PARALOG_SL:
        all_genes.add(a.upper())
        all_genes.add(b.upper())
    
    # Also add all paralog genes from our 66K table
    para_df = pd.read_csv(DATA_DIR / "ensembl_paralogs.csv", nrows=0)
    # (just use the driver set for now - enough for integration)
    
    print(f"Computing features for {len(all_genes)} genes...")
    
    # Load half-life
    hl_map = load_halflife_data(HALFLIFE_FILE)
    
    # Compute features
    records = []
    found = 0
    for gene in sorted(all_genes):
        if gene in gene_seqs:
            feats = compute_sequence_features(gene_seqs[gene])
            feats["gene"] = gene
            feats["half_life_hours"] = hl_map.get(gene)
            records.append(feats)
            found += 1
        else:
            records.append({"gene": gene, "length": None})
    
    df = pd.DataFrame(records)
    df.to_csv(OUTPUT, index=False)
    
    print(f"  Genes with sequences: {found}/{len(all_genes)}")
    print(f"  With half-life: {df['half_life_hours'].notna().sum()}")
    print(f"  Mean lysine density: {df['lysine_density'].mean():.4f}")
    print(f"  Mean disorder propensity: {df['disorder_propensity'].mean():.3f}")
    print(f"Saved to {OUTPUT}")


if __name__ == "__main__":
    main()
