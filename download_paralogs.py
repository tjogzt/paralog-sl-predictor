"""
BioMart Paralogs Downloader
============================
Downloads the complete set of human paralogs from Ensembl BioMart.
Replaces the built-in 180-pair table with thousands of real pairs.

Uses Ensembl BioMart REST API with XML query.
"""

import time
import requests
import pandas as pd
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"
OUTPUT_FILE = DATA_DIR / "ensembl_paralogs.csv"

BIOMART_URL = "https://www.ensembl.org/biomart/martservice"

def fetch_all_human_paralogs():
    """
    Fetch all human paralog pairs via BioMart.
    
    Strategy: Query hsapiens_gene_ensembl joined with ensembl_compara
    to get paralog relationships with homology type and percent identity.
    """
    print("Fetching complete human paralog set from Ensembl BioMart...")
    
    # BioMart XML query for human paralogs
    # This queries the compara database for paralog relationships
    xml_query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="CSV" header="1" uniqueRows="1" 
       datasetConfigVersion="0.6" completionStamp="1">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
  </Dataset>
</Query>"""
    
    # Step 1: Get all human genes with gene names
    print("  Step 1: Fetching all human gene names...")
    r = requests.post(BIOMART_URL, data={"query": xml_query}, timeout=120)
    r.raise_for_status()
    lines = r.text.strip().split("\n")
    print(f"    Got {len(lines)-1} genes")
    
    # Parse to dict: ensembl_id → gene_symbol
    gene_map = {}
    for line in lines[1:]:  # skip header
        parts = line.split(",")
        if len(parts) >= 2:
            eid = parts[0].strip()
            symbol = parts[1].strip().upper()
            if symbol:
                gene_map[eid] = symbol
    
    # Step 2: Get paralog relationships from compara
    # We use the homology table to get paralogs
    print("  Step 2: Fetching paralog relationships from Ensembl Compara...")
    
    # Build a query for paralogs
    # We filter for human-human paralogs (homology type includes 'paralog')
    xml_query2 = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="CSV" header="1" uniqueRows="1"
       datasetConfigVersion="0.6" completionStamp="1">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
    <Filter name="homolog_paralog_ortholog" value="paralog"/>
    <Attribute name="hsapiens_homolog_ensembl_gene"/>
    <Attribute name="hsapiens_homolog_associated_gene_name"/>
    <Attribute name="hsapiens_homolog_orthology_type"/>
    <Attribute name="hsapiens_homolog_perc_id"/>
    <Attribute name="hsapiens_homolog_perc_id_r1"/>
  </Dataset>
</Query>"""
    
    # Try with compara homology type filter
    try:
        r2 = requests.post(BIOMART_URL, data={"query": xml_query2}, timeout=300)
        r2.raise_for_status()
        results_text = r2.text
    except Exception as e:
        print(f"    BioMart query failed: {e}")
        print("    Trying alternative approach via Ensembl REST API...")
        return fetch_via_rest_api(gene_map)
    
    lines2 = results_text.strip().split("\n")
    print(f"    Got {len(lines2)-1} paralog relationships")
    
    if len(lines2) <= 1:
        print("    BioMart returned empty. Trying REST API fallback...")
        return fetch_via_rest_api(gene_map)
    
    # Parse results
    records = []
    for line in lines2[1:]:
        parts = line.split(",")
        if len(parts) >= 6:
            gene_a = parts[1].strip() if parts[1].strip() else gene_map.get(parts[0].strip(), parts[0].strip())
            gene_b = parts[3].strip() if parts[3].strip() else parts[2].strip()
            orthology_type = parts[4].strip()
            perc_id = parts[5].strip() if parts[5].strip() else ""
            
            gene_a = gene_a.upper()
            gene_b = gene_b.upper()
            
            if gene_a and gene_b and gene_a != gene_b:
                try:
                    perc_id_val = float(perc_id) if perc_id else None
                except ValueError:
                    perc_id_val = None
                
                records.append({
                    "gene_A": gene_a,
                    "gene_B": gene_b,
                    "homology_type": orthology_type,
                    "identity_pct": perc_id_val,
                })
    
    if len(records) < 100:
        print(f"    Only {len(records)} records from BioMart. Trying REST API fallback...")
        return fetch_via_rest_api(gene_map)
    
    df = pd.DataFrame(records).drop_duplicates()
    print(f"    Final: {len(df)} paralog pairs among {df['gene_A'].nunique()} genes")
    
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"    Saved to {OUTPUT_FILE}")
    return df


def fetch_via_rest_api(gene_map=None):
    """
    Fallback: Use Ensembl REST API to fetch paralogs for a curated
    set of cancer-relevant genes + all driver genes from config.
    
    This is slower than BioMart but more reliable for targeted queries.
    """
    print("    Using REST API fallback (targeted gene set)...")
    
    from config import DRIVER_GENES
    
    # Build comprehensive query set
    query_genes = set()
    for genes in DRIVER_GENES.values():
        query_genes.update(g.upper() for g in genes)
    
    # Add common cancer genes with known paralogs
    extra_genes = [
        "AKT1", "AKT2", "AKT3", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2",
        "RAF1", "BRAF", "ARAF", "NRAS", "HRAS", "CDK2", "CDK4", "CDK6",
        "CCND1", "CCND2", "CCND3", "BCL2", "BCL2L1", "MCL1",
        "MYC", "MYCN", "MYCL", "E2F1", "E2F2", "E2F3",
        "SMAD2", "SMAD3", "SMAD4", "STAT1", "STAT2", "STAT3",
        "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
        "FGFR1", "FGFR2", "FGFR3", "FGFR4",
        "ERBB2", "ERBB3", "ERBB4",
        "PARP1", "PARP2", "PARP3",
        "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC6",
        "DNMT1", "DNMT3A", "DNMT3B",
        "CASP3", "CASP7", "CASP8", "CASP9",
        "BAX", "BAK1", "BOK", "BID", "BIM",
        "RAD51", "RAD51B", "RAD51C", "RAD51D", "XRCC2", "XRCC3",
        "BRIP1", "PALB2", "BARD1",
        "MLH1", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2",
        "RPA1", "RPA2", "RPA3",
        "TOP2A", "TOP2B",
        "HIF1A", "HIF2A", "ARNT",
        "NFKB1", "NFKB2", "RELA", "RELB", "REL",
        "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2",
        "WNT1", "WNT2", "WNT3", "WNT5A", "WNT5B",
        "LATS1", "LATS2", "YAP1", "TAZ",
        "RB1", "RBL1", "RBL2",
        "CHEK1", "CHEK2", "ATM", "ATR",
        "WEE1", "PKMYT1",
        "EZH1", "EZH2",
        "SMARCA2", "SMARCA4", "ARID1A", "ARID1B",
        "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCD2", "SMARCD3",
        "KMT2A", "KMT2B", "KMT2C", "KMT2D",
        "SETD2", "SETDB1",
        "CUL1", "CUL2", "CUL3", "CUL4A", "CUL4B", "CUL5",
        "FBXW7", "FBXW2", "FBXW4", "FBXW5", "FBXW8",
        "SPOP", "SPOPL",
        "KEAP1", "KLHL20",
        "VHL",
        "MDM2", "MDM4",
        "BIRC2", "BIRC3", "XIAP",
        "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG",
        "PIK3R1", "PIK3R2", "PIK3R3",
        "MTOR", "RPTOR", "RICTOR",
        "TSC1", "TSC2",
        "AXIN1", "AXIN2",
        "APC", "APC2",
        "TCF7", "TCF7L1", "TCF7L2",
        "LEF1",
        "PPP2CA", "PPP2CB", "PPP2R1A", "PPP2R1B",
        "PPP2R2A", "PPP2R2B", "PPP2R2C", "PPP2R2D",
        "PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D", "PPP2R5E",
        "STK11", "SIK1", "SIK2", "SIK3",
        "AMPK", "NUAK1", "NUAK2",
        "EP300", "CREBBP",
        "BRD2", "BRD3", "BRD4", "BRDT",
        "TET1", "TET2", "TET3",
        "IDH1", "IDH2",
        "GATA1", "GATA2", "GATA3", "GATA4", "GATA6",
        "RUNX1", "RUNX2", "RUNX3",
        "CEBPA", "CEBPB", "CEBPD", "CEBPE",
        "AR", "ESR1", "ESR2", "PGR", "NR3C1",
        "KIT", "PDGFRA", "PDGFRB", "FLT3", "CSF1R",
        "MET", "MST1R",
        "IGF1R", "INSR", "INSRR",
        "TGFBR1", "TGFBR2", "ACVR1", "ACVR2A", "ACVR2B",
        "BMPR1A", "BMPR1B", "BMPR2",
        "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6",
    ]
    query_genes.update(g.upper() for g in extra_genes)
    
    server = "https://rest.ensembl.org"
    all_paralogs = []
    
    print(f"    Querying {len(query_genes)} genes via REST API...")
    
    for i, gene in enumerate(sorted(query_genes)):
        if i % 10 == 0:
            print(f"      Progress: {i}/{len(query_genes)}")
        try:
            ext = f"/homology/symbol/homo_sapiens/{gene}?type=paralogues;format=json"
            r = requests.get(f"{server}{ext}",
                            headers={"Content-Type": "application/json"},
                            timeout=15)
            if not r.ok:
                continue
            
            data = r.json()
            homologues = data.get("data", [])
            if not homologues:
                continue
            
            for h in homologues[:5]:  # cap at 5 paralogs per gene
                h_type = h.get("homology_type", "")
                target = h.get("target", {})
                target_species = target.get("species", "")
                
                # Only keep human-human paralogs
                if "homo_sapiens" not in target_species.lower():
                    continue
                
                target_gene = target.get("id", "")
                perc_id = h.get("perc_id", 0) or 0
                
                # Fetch target symbol
                try:
                    time.sleep(0.05)
                    t_r = requests.get(
                        f"{server}/lookup/id/{target_gene}?expand=0",
                        headers={"Content-Type": "application/json"},
                        timeout=10)
                    if t_r.ok:
                        t_data = t_r.json()
                        target_symbol = t_data.get("display_name", "").upper()
                        if target_symbol and target_symbol != gene:
                            all_paralogs.append({
                                "gene_A": gene,
                                "gene_B": target_symbol,
                                "homology_type": h_type,
                                "identity_pct": float(perc_id) if perc_id else None,
                            })
                except Exception:
                    pass
        except Exception:
            continue
    
    if not all_paralogs:
        print("    REST API returned no results. Using built-in paralog table.")
        from data_loader import _build_builtin_paralogs
        df = _build_builtin_paralogs()
    else:
        df = pd.DataFrame(all_paralogs).drop_duplicates()
    
    print(f"    Final: {len(df)} paralog pairs among {df['gene_A'].nunique()} genes")
    df.to_csv(OUTPUT_FILE, index=False)
    return df


if __name__ == "__main__":
    df = fetch_all_human_paralogs()
    print(f"\nDone. {len(df)} paralog pairs saved to {OUTPUT_FILE}")
