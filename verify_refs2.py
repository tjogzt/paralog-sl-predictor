"""Verify remaining references against PubMed."""
import subprocess, xml.etree.ElementTree as ET, time

refs = {
    "Chen2023":   ("Chen L", "Bioinformatics", "2023", "SLGNN"),
    "Das2024":    ("Das S", "Bioinformatics", "2024", "DDSL"),
    "Liu2024":    ("Liu H", "Brief Bioinform", "2024", "MVGCN"),
    "Li2023":     ("Li J", "Brief Bioinform", "2023", "NSF4SL"),
    "Zhang2024":  ("Zhang Y", "Genome Biol", "2024", "SLAST"),
    "Taylor2022": ("Taylor SE", "Cancer Res", "2022", "PPP2R1A"),
    "D'Antonio2023": ("D'Antonio M", "Cell Syst", "2023", "paralog"),
    "ICGC2020":   ("ICGC", "Nature", "2020", "pan-cancer"),
    "Benjamini1995":("Benjamini Y", "J R Stat Soc B", "1995", "FDR"),
    "Pedregosa2011":("Pedregosa F", "J Mach Learn Res", "2011", "scikit"),
    "Hoffman2014": ("Hoffman GR", "Proc Natl Acad Sci USA", "2014", "SMARCA2"),
}

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

for key, (author, journal, year, keyword) in refs.items():
    author_enc = author.replace("'", "").replace(" ", "+").replace("/", "")
    journal_enc = journal.replace(" ", "+")
    kw_enc = keyword.replace(" ", "+")
    query = f"{author_enc}%5BAuthor%5D+AND+{year}%5Bdp%5D+AND+{kw_enc}%5BTitle%5D"
    
    try:
        url = f"{BASE}/esearch.fcgi?db=pubmed&retmax=3&term={query}"
        r1 = subprocess.run(['curl', '-s', url], capture_output=True, text=True, timeout=15)
        root = ET.fromstring(r1.stdout)
        ids = [e.text for e in root.findall('.//Id')]
        count = root.findtext('.//Count', '0')
        
        if not ids:
            # Try without journal filter
            query2 = f"{author_enc}%5BAuthor%5D+AND+{year}%5Bdp%5D+AND+{kw_enc}%5BAll%5D"
            url2 = f"{BASE}/esearch.fcgi?db=pubmed&retmax=3&term={query2}"
            r2 = subprocess.run(['curl', '-s', url2], capture_output=True, text=True, timeout=15)
            root2 = ET.fromstring(r2.stdout)
            ids2 = [e.text for e in root2.findall('.//Id')]
            count2 = root2.findtext('.//Count', '0')
            
            if ids2:
                pmid = ids2[0]
                url3 = f"{BASE}/esummary.fcgi?db=pubmed&id={pmid}"
                r3 = subprocess.run(['curl', '-s', url3], capture_output=True, text=True, timeout=15)
                root3 = ET.fromstring(r3.stdout)
                title = root3.findtext('.//Item[@Name="Title"]', '')
                source = root3.findtext('.//Item[@Name="Source"]', '')
                yr = root3.findtext('.//Item[@Name="PubDate"]', '')[:4]
                print(f"  ⚠️  {key:20s} PMID:{pmid}: {title[:80]} [{source} {yr}] (broad search)")
            else:
                print(f"  ❌ {key:20s} not found in PubMed (count={count2})")
        else:
            # Found
            for pmid in ids[:1]:
                url3 = f"{BASE}/esummary.fcgi?db=pubmed&id={pmid}"
                r3 = subprocess.run(['curl', '-s', url3], capture_output=True, text=True, timeout=15)
                root3 = ET.fromstring(r3.stdout)
                title = root3.findtext('.//Item[@Name="Title"]', '')
                source = root3.findtext('.//Item[@Name="Source"]', '')
                yr = root3.findtext('.//Item[@Name="PubDate"]', '')[:4]
                print(f"  ✅ {key:20s} PMID:{pmid}: {title[:80]} [{source} {yr}]")
        
        time.sleep(0.4)
        
    except Exception as e:
        print(f"  💥 {key:20s} Error: {str(e)[:50]}")
