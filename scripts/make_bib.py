#!/usr/bin/env python3
"""Build a BibTeX file from PubMed metadata (NCBI E-utilities efetch)."""
import sys, time, subprocess, urllib.parse
import xml.etree.ElementTree as ET

PMIDS = {
 "bray2024global": "38572751",
 "siegel2023colorectal": "36856579",
 "sudlow2015ukbiobank": "25826379",
 "wu2018integrative": "29500431",
 "vosa2021eqtlgen": "34475573",
 "gtex2020": "32913098",
 "pietzner2021mapping": "34648354",
 "hemani2018mrbase": "29846171",
 "verbanck2018presso": "29686387",
 "bowden2015egger": "26050253",
 "bowden2016median": "27061298",
 "hartwig2017mode": "29040600",
 "burgess2013ivw": "24114802",
 "skrivankova2021strobe": "34698778",
 "zheng2020phenome": "32895551",
 "schmidt2020drugtarget": "32591531",
 "ochoa2021opentargets": "33196847",
 "kelleher2023pharos": "36624666",
 "finan2017druggable": "28356508",
 "kurki2023finngen": "36653562",
 "hao2021seurat": "34062119",
 "wolf2018scanpy": "29409532",
 "kleshchevnikov2022cell2location": "35027729",
 "stahl2016spatial": "27365449",
 "weinstein2013tcga": "24071849",
 "buniello2019gwas": "30445434",
 "denny2010phewas": "20335276",
 "hong2025integrative": "40496862",
 "tian2025dbi": "40426943",
 "wang2025discovoncol": "40381082",
 "xia2025bmccancer": "40855282",
}

URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def text(node, path):
    el = node.find(path)
    if el is None:
        return None
    return "".join(el.itertext()).strip()

def fetch(ids):
    q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids), "retmode": "xml"})
    out = subprocess.run(["curl", "-s", "--max-time", "60", URL + "?" + q],
                         capture_output=True, check=True)
    return out.stdout

records = {}
items = list(PMIDS.items())
for i in range(0, len(items), 15):
    chunk = items[i:i + 15]
    root = ET.fromstring(fetch([p for _, p in chunk]))
    for art in root.findall(".//PubmedArticle"):
        cit = art.find("MedlineCitation")
        pmid = text(cit, "PMID")
        artel = cit.find("Article")
        authors = []
        for a in artel.findall(".//AuthorList/Author"):
            last = text(a, "LastName")
            init = text(a, "Initials")
            coll = text(a, "CollectiveName")
            if last:
                authors.append("%s, %s" % (last, init) if init else last)
            elif coll:
                authors.append(coll)
        journal = text(artel, "Journal/Title") or text(artel, "Journal/ISOAbbreviation")
        year = text(artel, "Journal/JournalIssue/PubDate/Year") or (text(artel, "Journal/JournalIssue/PubDate/MedlineDate") or "")[:4]
        vol = text(artel, "Journal/JournalIssue/Volume")
        issue = text(artel, "Journal/JournalIssue/Issue")
        pages = text(artel, "Pagination/MedlinePgn") or text(artel, "Pagination/StartPage")
        doi = None
        for aid in artel.findall(".//ArticleIdList/ArticleId"):
            if aid.get("IdType") == "doi":
                doi = (aid.text or "").strip()
        pmcid = None
        for aid in art.findall(".//ArticleIdList/ArticleId"):
            if aid.get("IdType") == "pmc":
                pmcid = (aid.text or "").strip()
        records[pmid] = dict(title=text(artel, "ArticleTitle"), authors=authors, journal=journal,
                             year=year, volume=vol, issue=issue, pages=pages, doi=doi, pmcid=pmcid)
    time.sleep(0.4)

missing = [k for k, p in PMIDS.items() if p not in records]
out = []
for key, pmid in PMIDS.items():
    r = records.get(pmid)
    if not r:
        out.append("%% MISSING: %s (PMID %s)\n" % (key, pmid))
        continue
    fields = ["  title={{%s}}" % (r["title"] or "")]
    fields.append("  author={%s}" % " and ".join(r["authors"]) if r["authors"] else "  author={}")
    fields.append("  journal={%s}" % (r["journal"] or ""))
    if r["volume"]: fields.append("  volume={%s}" % r["volume"])
    if r["issue"]:  fields.append("  number={%s}" % r["issue"])
    if r["pages"]:  fields.append("  pages={%s}" % r["pages"])
    fields.append("  year={%s}" % (r["year"] or ""))
    if r["doi"]:   fields.append("  doi={%s}" % r["doi"])
    fields.append("  pmid={%s}" % pmid)
    out.append("@article{%s,\n%s\n}\n" % (key, ",\n".join(fields)))

open("manuscript/crc_smr_atlas_new.bib", "w", encoding="utf-8").write("\n".join(out))
print("wrote %d entries; missing: %s" % (len(PMIDS) - len(missing), missing))