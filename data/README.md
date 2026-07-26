# Data Sources

Large data files are excluded from this repository via `.gitignore`. Download them manually before running the pipeline.

## Required Files

| File | Source | Size | URL |
|------|--------|------|-----|
| `CRISPRGeneEffect.csv` | DepMap 26Q1 | 420 MB | https://depmap.org/portal/download/ |
| `OmicsSomaticMutations.csv` | DepMap 26Q1 | 554 MB | https://depmap.org/portal/download/ |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | DepMap 26Q1 | 483 MB | https://depmap.org/portal/download/ |
| `OmicsCNGene.csv` | DepMap 26Q1 | ~400 MB | https://depmap.org/portal/download/ |
| `Model.csv` | DepMap 26Q1 | 681 KB | https://depmap.org/portal/download/ |
| `PRISM_log2AUC.csv` | DepMap PRISM Repurposing | ~50 MB | https://depmap.org/portal/download/ |

## Included Files (small, tracked by Git)

| File | Description |
|------|-------------|
| `ensembl_paralogs.csv` | HGNC paralog pair definitions |
| `protein_features.csv` | UniProt-derived protein features |
| `uniprot_sequences.rds` | Cached UniProt protein sequences |
| `cptac_cache/*.json` | CPTAC protein abundance per cohort |
| `OmicsGlobalSignatures.csv` | DepMap 26Q1 official genomic signatures incl. MSIsensor2 MSIscore (MSI annotation; MSI-H = score > 20) |

## Download Instructions

```bash
# 1. Visit https://depmap.org/portal/download/
# 2. Select release "DepMap Public 26Q1"
# 3. Download the files listed above
# 4. Place them in this data/ directory
```
