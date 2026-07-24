"""Verify manuscript references against PubMed."""
import subprocess, json, time, re, xml.etree.ElementTree as ET

refs = {
    "ONeil2017":   ("O'Neil NJ", "Nat Rev Genet", "2017"),
    "Huang2020":   ("Huang A", "Nat Rev Drug Discov", "2020"),
    "Lord2017":    ("Lord CJ", "Science", "2017"),
    "Bryant2005":  ("Bryant HE", "Nature", "2005"),
    "Farmer2005":  ("Farmer H", "Nature", "2005"),
    "Zheng2022":   ("Wang J", "Database", "2022"),
    "Feng2024":    ("Feng Y", "Nat Commun", "2024"),
    "Hao2021":     ("Hao Z", "IEEE J Biomed Health Inform", "2021"),
    "Long2021":    ("Long Y", "Bioinformatics", "2021"),
    "Cai2020":     ("Cai R", "Bioinformatics", "2020"),
    "Huang2019":   ("Huang J", "BMC Bioinformatics", "2019"),
    "Wang2021":    ("Wang S", "Bioinformatics", "2021"),
    "Wang2022":    ("Wang S", "Bioinformatics", "2022"),
    "Zhu2023":     ("Zhu Y", "Bioinformatics", "2023"),
    "Liu2022":     ("Liu X", "Bioinformatics", "2022"),
    "Long2022":    ("Long Y", "Bioinformatics", "2022"),
    "Koonin2005":  ("Koonin EV", "Annu Rev Genet", "2005"),
    "DAntonio2013":("D'Antonio M", "Cell Rep", "2013"),
    "Hoffman2014": ("Hoffman GR", "Proc Natl Acad Sci", "2014"),
    "Helming2014": ("Helming KC", "Nat Med", "2014"),
    "Jia2008":     ("Jia S", "Nature", "2008"),
    "Wee2008":     ("Wee S", "Proc Natl Acad Sci", "2008"),
    "Corsello2020":("Corsello SM", "Nat Cancer", "2020"),
    "Meyers2017":  ("Meyers RM", "Nat Genet", "2017"),
    "Dempster2021":("Dempster JM", "Genome Biol", "2021"),
    "Cerami2012":  ("Cerami E", "Cancer Discov", "2012"),
    "Gao2013":     ("Gao J", "Sci Signal", "2013"),
    "Liu2018":     ("Liu J", "Cell", "2018"),
    "Bailey2018":  ("Bailey MH", "Cell", "2018"),
    "ICGC2020":    ("ICGC/TCGA", "Nature", "2020"),
    "Benjamini1995":("Benjamini Y", "J R Stat Soc", "1995"),
    "Edwards2015": ("Edwards NJ", "J Proteome Res", "2015"),
    "Bitler2015":  ("Bitler BG", "Nat Med", "2015"),
    "Pedregosa2011":("Pedregosa F", "J Mach Learn Res", "2011"),
    "Virtanen2020":("Virtanen P", "Nat Methods", "2020"),
    "Harris2020":  ("Harris CR", "Nature", "2020"),
    "Pacini2021":  ("Pacini C", "Nat Commun", "2021"),
}

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
verified = []; mismatch = []; failed = []

for key, (author, journal, year) in refs.items():
    author_enc = author.replace("'", "").replace(" ", "+")
    journal_enc = journal.replace(" ", "+")
    query = f"{author_enc}%5BAuthor%5D+AND+{year}%5Bdp%5D+AND+{journal_enc}%5BJournal%5D"
    
    try:
        # ESearch (XML response)
        url = f"{BASE}/esearch.fcgi?db=pubmed&retmax=1&term={query}"
        r1 = subprocess.run(['curl', '-s', url], capture_output=True, text=True, timeout=15)
        root = ET.fromstring(r1.stdout)
        ids = [e.text for e in root.findall('.//Id')]
        
        if not ids:
            failed.append((key, "no PMID"))
            continue
        
        pmid = ids[0]
        
        # ESummary (XML response)
        url2 = f"{BASE}/esummary.fcgi?db=pubmed&id={pmid}"
        r2 = subprocess.run(['curl', '-s', url2], capture_output=True, text=True, timeout=15)
        root2 = ET.fromstring(r2.stdout)
        
        # Extract fields
        art_title = root2.findtext('.//Item[@Name="Title"]', '')
        art_source = root2.findtext('.//Item[@Name="Source"]', '')
        art_pubdate = root2.findtext('.//Item[@Name="PubDate"]', '')[:4]
        art_volume = root2.findtext('.//Item[@Name="Volume"]', '')
        art_pages = root2.findtext('.//Item[@Name="Pages"]', '')
        
        # First author
        auth_list = root2.findall('.//Item[@Name="AuthorList"]//Item[@Name="Author"]')
        first_author = ''
        if auth_list:
            fa_elem = auth_list[0]
            last = fa_elem.findtext('Item[@Name="LastName"]', '')
            first = fa_elem.findtext('Item[@Name="ForeName"]', '')
            first_author = f"{last} {first}" if last and first else ''
        
        # Match check
        author_surname = author.split()[0].lower().replace("'","").replace("O'Neil","oneil")
        fa_surname = first_author.split()[0].lower() if first_author else ''
        am = author_surname in fa_surname or fa_surname in author_surname
        ym = year in art_pubdate
        jm = journal.lower()[:12] in art_source.lower()[:30]
        
        if am and ym and jm:
            verified.append((key, pmid, art_source[:30], art_volume, art_pages, art_title[:80]))
        else:
            mismatch.append((key, pmid, first_author[:25], art_pubdate, art_source[:25], art_volume, art_pages, art_title[:60]))
        
        time.sleep(0.35)
        
    except Exception as e:
        failed.append((key, str(e)[:60]))

print(f"\n✅ VERIFIED: {len(verified)} / {len(refs)}")
for v in verified:
    print(f"  {v[0]:20s} PMID:{v[1]:10s} {v[2]:28s} {v[3]}:{v[4]}  [{v[5]}]")

if mismatch:
    print(f"\n⚠️  MISMATCH: {len(mismatch)}")
    for m in mismatch:
        print(f"  {m[0]:20s} PMID:{m[1]:10s} Author:{m[2]:25s} Yr:{m[3]:5s} Jrnl:{m[4]:25s} Vol:{m[5]}:{m[6]}")

if failed:
    print(f"\n❌ FAILED: {len(failed)}")
    for f in failed:
        print(f"  {f[0]:20s} {f[1]}")
