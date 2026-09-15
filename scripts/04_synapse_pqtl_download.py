#!/usr/bin/env python3
"""
Download UKB-PPP pQTL files for the 75 significant eQTL genes.
Matches gene symbols to Synapse file names (format: GENE_UNIPROT_...).
"""
import synapseclient
import os
import sys
import re

# Config
TOKEN = "eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsib3BlbmlkIiwiZW1haWwiLCJwcm9maWxlIiwiZ2E0Z2hfcGFzc3BvcnRfdjEiLCJ2aWV3IiwiZG93bmxvYWQiLCJtb2RpZnkiLCJvZmZsaW5lX2FjY2VzcyJdLCJvaWRjX2NsYWltcyI6e319LCJ0b2tlbl90eXBlIjoiUEVSU09OQUxfQUNDRVNTX1RPS0VOIiwiaXNzIjoiaHR0cHM6Ly9yZXBvLXByb2QucHJvZC5zYWdlYmFzZS5vcmcvYXV0aC92MSIsImF1ZCI6IjAiLCJuYmYiOjE3ODU0MTQ0OTQsImlhdCI6MTc4NTQxNDQ5NCwianRpIjoiNDM0ODciLCJzdWIiOiIzNjA1ODEzIn0.uEV_Tjl_GgplNIyjiIDQp4z5aXbxZfQNISysPK6pHON07N2--KzlQm5QBT3GHLxDX_k11UCQsw-22pmf7F8g8jZmTeKUA44tJXGvMIhMjrneBbXi5BfHhShyxy7HOxpwo5ntQbrUcKoTV-x4--xjHlLFBDWgDXgmrIGh9mbFGZVaWJ5_X_o7yNguRMG1PGmAOih3JCPyvTBz6uaDsKsuyb5TeojrxKWJrAf0PypRY2XbFCXvUSyBAfK0qJ5sslO3ReDm1syYboYRAs_8A4glU5MENu0TNQYsOW9K4OYteCPqjboyt1b_vd6hRymCFVe-YB-cP6taFNMhcBBx8J57Pg"

PARENT_ID = "syn51365303"  # UKB-PPP European discovery pQTL folder
OUTDIR = "/ifs1/User/zhouman/project9-v5-crc-atlas/data/ukbppp"

# 75 significant genes (from Phase 1 SMR)
GENES = [
    "AAMP", "ABHD12B", "ACTR1B", "ARPC2", "ARPC5", "ASIC1", "ATF1",
    "B3GNTL1", "B9D2", "BMP2", "CABLES2", "CASC8", "CCDC97", "CCHCR1",
    "CCM2", "CDKN1A", "CERS5", "CLDN4", "COX14", "COX15", "CXCR1", "CXCR2",
    "EIF3H", "GDPGP1", "GPBAR1", "HCG27", "HLA-DRB5", "HLA-E", "HLA-K",
    "LAMC1", "LEMD3", "LIMA1", "LINC01270", "LOC128966610", "LTA",
    "MAPKAPK5-AS1", "METTL27", "MMP11", "MYO9A", "MZF1", "NCF2", "PELATON",
    "PGAP3", "PNKD", "POLD3", "RGL1", "SF3A3", "SLC25A28", "SMAD6", "SMAD7",
    "SMARCD1", "STAB1", "STAT6", "STIMATE", "TET2", "TMBIM1", "TMEM258",
    "TMT1A", "TNF", "TNNC1", "TRIM28", "TRIM4", "UBE2M", "UQCC5", "UTP11",
    "UTP23", "ZFP57", "ZFP90", "ZMIZ1"
]

os.makedirs(OUTDIR, exist_ok=True)


def main():
    print(f"Connecting to Synapse...")
    syn = synapseclient.Synapse()
    syn.login(authToken=TOKEN)
    print(f"Logged in as: {syn.username}")

    # List children
    print(f"\nListing files in {PARENT_ID}...")
    children = list(syn.getChildren(PARENT_ID))
    print(f"Found {len(children)} files")

    # Build mapping: gene_symbol → [(file_id, file_name)]
    file_map = {}
    for child in children:
        name = child['name']
        # Files named like: A1BG_P04217_A1BG:OID1:2008 chr1:228761401_A_G_hg38
        # or: ALB_P02768_ALB:OID30788 chr4:73404114_A_G_hg38
        # Match gene prefix before first underscore
        match = re.match(r'^([A-Z0-9]+(?:\.[0-9]+)?)_', name)
        if match:
            gene = match.group(1)
            if gene not in file_map:
                file_map[gene] = []
            file_map[gene].append((child['id'], name))

    # Match our genes
    matched = []
    unmatched = []
    for gene in GENES:
        if gene in file_map:
            matched.append((gene, file_map[gene]))
        else:
            # Try case-insensitive
            found = False
            for fgene in file_map:
                if fgene.upper() == gene.upper():
                    matched.append((gene, file_map[fgene]))
                    found = True
                    break
            if not found:
                unmatched.append(gene)

    print(f"\nMatched: {len(matched)}/{len(GENES)} genes")
    print(f"Unmatched: {len(unmatched)}/{len(GENES)} genes")
    if unmatched:
        print(f"Unmatched genes: {', '.join(unmatched)}")

    # Download matched files
    print(f"\n=== Downloading {len(matched)} genes ===")
    for gene, files in matched:
        for fid, fname in files:
            outpath = os.path.join(OUTDIR, fname)
            if os.path.exists(outpath):
                print(f"  {gene}: {fname} (already exists, skip)")
                continue
            print(f"  {gene}: downloading {fname} ({fid})...")
            try:
                syn.get(fid, downloadLocation=OUTDIR)
                print(f"    ✅ done")
            except Exception as e:
                print(f"    ❌ failed: {e}")

    # Write gene-to-file mapping
    map_path = os.path.join(OUTDIR, "gene_pqtl_mapping.tsv")
    with open(map_path, 'w') as f:
        f.write("gene\tfile_id\tfile_name\n")
        for gene, files in matched:
            for fid, fname in files:
                f.write(f"{gene}\t{fid}\t{fname}\n")
    print(f"\nMapping saved to {map_path}")

    print(f"\n=== Summary ===")
    print(f"Total genes: {len(GENES)}")
    print(f"Matched in UKB-PPP: {len(matched)}")
    print(f"Unmatched: {len(unmatched)}")
    print(f"Output dir: {OUTDIR}")


if __name__ == "__main__":
    main()
