"""Complete independent verification of ALL 35 manuscript references against PubMed."""
import subprocess, xml.etree.ElementTree as ET, time, re

# All 35 references from manuscript.tex, with key search terms
refs = [
    # SL reviews
    ("ONeil2017",   "O'Neil+NJ", "2017", "Nat+Rev+Genet", "synthetic+lethality+cancer"),
    ("Huang2020",   "Huang+A", "2020", "Nat+Rev+Drug+Discov", "synthetic+lethality+drug"),
    ("Lord2017",    "Lord+CJ", "2017", "Science", "PARP+inhibitors+clinic"),
    # Landmark PARP papers
    ("Bryant2005",  "Bryant+HE", "2005", "Nature", "BRCA2+polymerase"),
    ("Farmer2005",  "Farmer+H", "2005", "Nature", "BRCA+mutant+repair"),
    # SL databases & methods
    ("Zheng2022",   "Zheng+J", "2022", "Database", "SynLethDB"),
    ("Feng2024",    "Feng+Y", "2024", "Nat+Commun", "benchmarking+synthetic+lethality"),
    ("Chen2023",    "Chen+L", "2023", "Bioinformatics", "SLGNN+synthetic+lethality"),
    ("Das2024",     "Das+S", "2024", "Bioinformatics", "DDSL+deep+double+strand"),
    ("Liu2024",     "Liu+H", "2024", "Brief+Bioinform", "MVGCN+multi+view+synthetic"),
    ("Li2023",      "Li+J", "2023", "Brief+Bioinform", "NSF4SL+negative+sampling"),
    ("Zhang2024",   "Zhang+Y", "2024", "Genome+Biol", "SLAST+sequence+learning+synthetic"),
    # Evolutionary biology
    ("Koonin2005",  "Koonin+EV", "2005", "Annu+Rev+Genet", "orthologs+paralogs+evolutionary"),
    # Paralog dependency
    ("DAntonio2023","D'Antonio+M", "2023", "Cell+Syst", "paralog+dependency+cancer"),
    # SWI/SNF and chromatin
    ("Hoffman2014", "Hoffman+GR", "2014", "Proc+Natl+Acad+Sci+USA", "BRM+SMARCA2+synthetic"),
    ("Helming2014", "Helming+KC", "2014", "Nat+Med", "ARID1B+ARID1A"),
    ("Bitler2015",  "Bitler+BG", "2015", "Nat+Med", "EZH2+ARID1A"),
    # PI3K signaling
    ("Jia2008",     "Jia+S", "2008", "Nature", "PI3K+p110beta+tumorigenesis"),
    ("Wee2008",     "Wee+S", "2008", "PNAS", "PTEN+PIK3CB"),
    # PPP2R1A
    ("Taylor2022",  "Taylor+SE", "2022", "Cancer+Res", "PPP2R1A+RNR+ribonucleotide"),
    # Data resources
    ("DepMap2026",  None, None, None, None),  # URL, skip
    ("Meyers2017",  "Meyers+RM", "2017", "Nat+Genet", "CERES+CRISPR+copy+number"),
    ("Dempster2021","Dempster+JM", "2021", "Genome+Biol", "Chronos+CRISPR"),
    ("Cerami2012",  "Cerami+E", "2012", "Cancer+Discov", "cBioPortal+cancer+genomics"),
    ("Gao2013",     "Gao+J", "2013", "Sci+Signal", "integrative+cancer+cBioPortal"),
    ("Corsello2020","Corsello+SM", "2020", "Nat+Cancer", "PRISM+viability+drug+non-oncology"),
    ("Liu2018",     "Liu+J", "2018", "Cell", "TCGA+pan-cancer+clinical"),
    ("Bailey2018",  "Bailey+MH", "2018", "Cell", "cancer+driver+genes+mutations"),
    ("ICGC2020",    None, "2020", "Nature", "pan-cancer+whole+genomes"),  # consortium
    ("Benjamini1995","Benjamini+Y", "1995", None, "false+discovery+rate"),  # stats journal
    ("Pedregosa2011","Pedregosa+F", "2011", None, "scikit-learn+machine+learning"),  # CS journal
    ("Edwards2015", "Edwards+NJ", "2015", "J+Proteome+Res", "CPTAC+data+portal+proteomics"),
    ("Virtanen2020","Virtanen+P", "2020", "Nat+Methods", "SciPy+computing"),
    ("Harris2020",  "Harris+CR", "2020", "Nature", "NumPy+array+programming"),
    ("Pacini2021",  "Pacini+C", "2021", "Nat+Commun", "cross-study+genetic+dependencies"),
]

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
results = []

for key, author, year, journal, keyword in refs:
    if key == "DepMap2026":
        results.append((key, "✅", "N/A", "URL reference (depmap.org)", "", "", ""))
        continue
    
    # Build query
    parts = []
    if author: parts.append(f"{author}%5BAuthor%5D")
    if year: parts.append(f"{year}%5Bdp%5D")
    if keyword: parts.append(f"{keyword}%5BTitle%5D")
    query = "+AND+".join(parts)
    
    try:
        url = f"{BASE}/esearch.fcgi?db=pubmed&retmax=3&term={query}"
        r1 = subprocess.run(['curl', '-s', url], capture_output=True, text=True, timeout=15)
        root = ET.fromstring(r1.stdout)
        count = int(root.findtext('.//Count', '0'))
        ids = [e.text for e in root.findall('.//Id')]
        
        if not ids:
            # Fallback: try without journal filter
            parts2 = [p for p in parts if not p.endswith('%5BJournal%5D')]
            if keyword not in [p for p in parts if '%5BTitle%5D' in p]:
                parts2.append(f"{keyword}%5BAll%5D")
            query2 = "+AND+".join(parts2)
            url2 = f"{BASE}/esearch.fcgi?db=pubmed&retmax=3&term={query2}"
            r2 = subprocess.run(['curl', '-s', url2], capture_output=True, text=True, timeout=15)
            root2 = ET.fromstring(r2.stdout)
            ids2 = [e.text for e in root2.findall('.//Id')]
            count2 = int(root2.findtext('.//Count', '0'))
            
            if ids2:
                ids = ids2
                count = count2
            else:
                results.append((key, "❌", "0", f"no PMID (count={count2})", "", "", ""))
                time.sleep(0.35)
                continue
        
        pmid = ids[0]
        
        # Fetch summary
        url3 = f"{BASE}/esummary.fcgi?db=pubmed&id={pmid}"
        r3 = subprocess.run(['curl', '-s', url3], capture_output=True, text=True, timeout=15)
        root3 = ET.fromstring(r3.stdout)
        
        title = root3.findtext('.//Item[@Name="Title"]', '')[:90]
        source = root3.findtext('.//Item[@Name="Source"]', '')
        pubdate = root3.findtext('.//Item[@Name="PubDate"]', '')[:4]
        volume = root3.findtext('.//Item[@Name="Volume"]', '')
        pages = root3.findtext('.//Item[@Name="Pages"]', '')
        doi = root3.findtext('.//Item[@Name="DOI"]', '')
        
        # First author
        auth_list = root3.findall('.//Item[@Name="AuthorList"]//Item[@Name="Author"]')
        first_author = ""
        if auth_list:
            l = auth_list[0].findtext('Item[@Name="LastName"]', '')
            f = auth_list[0].findtext('Item[@Name="ForeName"]', '')
            first_author = f"{l} {f}".strip()
        
        # Match checks
        expected_author = author.replace("+", " ").replace("'","") if author else ""
        author_ok = expected_author.lower().split()[0] in first_author.lower() if expected_author and first_author else True
        year_ok = (year == pubdate) if year else True
        journal_ok = (journal and journal.replace("+"," ").lower()[:15] in source.lower()[:30])
        
        status = "✅" if (author_ok and year_ok and journal_ok) else "⚠️"
        detail = f"{first_author[:20]} | {source[:25]} | {volume}:{pages} | {pubdate}"
        
        if not author_ok: detail += " [author mismatch]"
        if not year_ok: detail += f" [year: expected={year}, got={pubdate}]"
        if not journal_ok: detail += " [journal mismatch]"
        
        results.append((key, status, pmid, title, volume, pages, detail))
        
        time.sleep(0.35)
        
    except Exception as e:
        results.append((key, "💥", "error", str(e)[:60], "", "", ""))

# Print table
print(f"\n{'Key':20s} {'St':3s} {'PMID':11s} {'Title':50s} {'Vol:Pages':15s} {'Details'}")
print("-" * 140)
for r in results:
    key, status, pmid, title, vol, pages, detail = r
    print(f"{key:20s} {status:3s} {pmid:11s} {title:50s} {vol}:{pages:8s}  {detail[:80]}")

# Summary
n_ok = sum(1 for r in results if r[1] == "✅")
n_warn = sum(1 for r in results if r[1] == "⚠️")
n_fail = sum(1 for r in results if r[1] == "❌")
n_err = sum(1 for r in results if r[1] == "💥")
print(f"\n✅={n_ok}  ⚠️={n_warn}  ❌={n_fail}  💥={n_err}  (total={len(results)})")
