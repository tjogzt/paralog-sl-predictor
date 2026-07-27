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

## File Integrity (SHA256) and Download Dates

SHA256 checksums of the exact files used in this study, for verification after download (`shasum -a 256 <file>`):

| File | SHA256 | Downloaded |
|------|--------|-----------|
| `CRISPRGeneEffect.csv` | `e610a4cefb13a82b5b256b47eb08b63ff14843f8dbd0fb164bc0a32688e5b89e` | 2026-05-24 |
| `OmicsSomaticMutations.csv` | `4d50634373578621bfb9f8ca69b38cfffc3f9b7abcdefd653044a675b856c520` | 2026-05-24 |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | `2a71dc94110efcc0221eae821bb93a9f03b54bea16f005818911a09d33383d56` | 2026-05-24 |
| `OmicsCNGene.csv` | `4851d3e939d48837a39a0f01294deb90fa507a85703586a927b77474f999134c` | 2026-05-24 |
| `Model.csv` | `ea4e0b2a3bc806f81df62689a5ae75f1a100135727a3d7b8a4c7ccc8815183f8` | 2026-05-24 |
| `PRISM_log2AUC.csv` | `19d3b51f214d8bfad627503804feee07a50afeb9eaf46df14667d9e173040f29` | 2026-05-24 |
| `OmicsGlobalSignatures.csv` | `d210b664edf42598129a3bb982af4a3b4aa22100411bb61ecb9baa7f81cf6011` | 2026-07-26 |
| `ensembl_paralogs.csv` | `6288928be7a49c494e548a9fa8d84e8f6a1b098f891956e14e8ef0b46dfa4308` | 2026-05-24 |
| `FD_GLBL_MI_FFPEbridge_Abund_20201002.tsv` (CPTAC UCEC) | `577eee78c55a6cd025725d320fb703a071410e5ad3f83a34efc01e713684e4e9` | 2026-05-25 |
