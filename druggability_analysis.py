import os, json
from collections import Counter

# Single-letter → three-letter amino acid mapping
AA_MAP = {
    'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS',
    'I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN',
    'R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'
}
HYDROPHOBIC = {'ALA','VAL','LEU','ILE','MET','PHE','TRP','PRO'}

def parse_pdb(filepath):
    """Parse PDB handling both standard (3-letter) and ESMFold (1-letter) formats."""
    residues = {}
    for line in open(filepath):
        if not line.startswith('ATOM'):
            continue
        atom_name = line[12:16].strip()
        if atom_name != 'CA':
            continue
        raw_resname = line[17:20].strip()
        # Handle single-letter codes (ESMFold format)
        if len(raw_resname) == 1 and raw_resname in AA_MAP:
            resname = AA_MAP[raw_resname]
        else:
            resname = raw_resname
        resnum = int(line[22:26])
        plddt = float(line[60:66])
        residues[resnum] = {'name': resname, 'plddt': plddt}
    return residues

monomer_dir = 'paralog_sl_predictor/structure_pred/outputs/monomers'

targets = [
    ('SMARCA2', 'Paralog of SMARCA4'),
    ('SMARCA4', 'Driver'),
    ('PIK3CB', 'Paralog of PIK3CA'),
    ('PIK3CA', 'Driver'),
    ('PPP2R1B', 'Paralog of PPP2R1A'),
    ('PPP2R1A', 'Driver'),
    ('FBXW2', 'Paralog of FBXW7'),
    ('FBXW7', 'Driver'),
    ('KRAS', 'Driver'),
    ('HRAS', 'Paralog of KRAS'),
    ('NF1', 'Driver'),
    ('RASA2', 'Paralog of NF1'),
    ('STK11', 'Driver'),
    ('SIK1', 'Paralog of STK11'),
    ('BRCA1', 'Driver (functional analog)'),
]

header = "{:<10} {:>5} {:>6} {:>6} {:>4} {:>4} {:>7} {:>7} {:>7}".format(
    "Protein", "Res", "pLDDT", ">70%", "Lys", "Cys", "Hydro%", "Pockets", "Drug.")
print(header)
print("-" * 68)

results = []
for gene, role in targets:
    path = os.path.join(monomer_dir, gene + '.pdb')
    if not os.path.exists(path):
        continue
    res = parse_pdb(path)
    if len(res) < 10:
        continue
    n = len(res)
    plddts = [r['plddt'] for r in res.values()]
    mean_plddt = sum(plddts) / len(plddts)
    high_conf = sum(1 for p in plddts if p > 0.7) / n * 100  # ESMFold pLDDT is 0-1
    if max(plddts) > 1:  # AlphaFold uses 0-100 scale
        high_conf = sum(1 for p in plddts if p > 70) / n * 100

    aa_counts = Counter(r['name'] for r in res.values())
    lys = aa_counts.get('LYS', 0)
    cys = aa_counts.get('CYS', 0)
    hydro_n = sum(aa_counts.get(a, 0) for a in HYDROPHOBIC)
    hydro_pct = hydro_n / n * 100

    # Potential pocket regions: structured + mixed hydrophobic/polar
    sorted_res = sorted(res.items())
    pocket_windows = 0
    for i in range(len(sorted_res) - 10):
        chunk = sorted_res[i:i+10]
        cp = sum(r[1]['plddt'] for r in chunk) / 10
        ch = sum(1 for r in chunk if r[1]['name'] in HYDROPHOBIC) / 10
        threshold = 0.7 if max(plddts) <= 1 else 70
        if cp > threshold and 0.2 <= ch <= 0.7:
            pocket_windows += 1
    pocket_regions = max(1, pocket_windows // 12)

    # Composite druggability score
    conf_norm = high_conf / 100
    drug_score = round(
        conf_norm * 0.30 +
        min(hydro_pct / 100, 0.5) * 0.50 +  # hydrophobic = binding potential
        min(lys / n * 6, 0.25) * 0.80 +      # Lys = PROTAC conjugation
        min(pocket_regions / 6, 0.20) * 1.0,  # pocket count
        3
    )

    row = "{:<10} {:>5} {:>6.2f} {:>5.1f}% {:>4} {:>4} {:>6.1f}% {:>7} {:>7.3f}".format(
        gene, n, mean_plddt, high_conf, lys, cys, hydro_pct, pocket_regions, drug_score)
    print(row)
    results.append({
        'gene': gene, 'role': role, 'residues': n,
        'mean_plddt': round(mean_plddt, 3), 'high_conf_pct': round(high_conf, 1),
        'lys_count': lys, 'cys_count': cys, 'hydrophobic_pct': round(hydro_pct, 1),
        'pocket_regions': pocket_regions, 'druggability_score': drug_score
    })

with open('paralog_sl_predictor/output/druggability_analysis.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved. Key finding:")
arid1b = [r for r in results if r['gene'] in ('SMARCA2','PIK3CB','FBXW2','HRAS','RASA2','SIK1','PPP2R1B')]
if arid1b:
    top = max(arid1b, key=lambda x: x['druggability_score'])
    print(f"  Top druggable paralog: {top['gene']} (score={top['druggability_score']}, Lys={top['lys_count']}, pockets={top['pocket_regions']})")
