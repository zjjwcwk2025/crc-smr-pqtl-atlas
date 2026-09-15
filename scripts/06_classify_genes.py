#!/usr/bin/env python3
"""
Classify the 75 SMR-significant genes into:
1. Protein-coding genes (can potentially have pQTL data)
2. Non-coding genes (lncRNA, pseudogene, antisense, etc.)
3. Cross-reference with deCODE and UKB-PPP pQTL coverage
"""
import subprocess
import sys

# 75 significant gene symbols (col 23 from annotated results)
GENES_SYMBOLS = [
    "CASC8", "LAMC1", "UTP23", "LIMA1", "COX14", "PNKD",
    "ENSG00000261338", "BMP2", "GPBAR1", "TMBIM1", "UTP11",
    "ENSG00000212541", "SMAD7", "ENSG00000185432", "ENSG00000183431",
    "CDKN1A", "ENSG00000115073", "POLD3", "ENSG00000127837", "TNF",
    "ENSG00000123810", "HLA-DRB5", "ENSG00000010327", "ENSG00000149679",
    "ENSG00000134825", "ENSG00000123268", "ENSG00000213533", "ENSG00000139624",
    "ENSG00000163464", "ENSG00000267764", "ENSG00000180871", "ENSG00000226979",
    "ENSG00000137834", "ENSG00000110881", "ENSG00000224397", "ENSG00000163466",
    "ENSG00000136280", "ENSG00000162704", "ENSG00000204592", "ENSG00000175711",
    "ENSG00000099326", "ENSG00000114854", "ENSG00000165171", "ENSG00000168273",
    "ENSG00000066117", "ENSG00000142039", "ENSG00000155287", "ENSG00000260997",
    "ENSG00000183208", "ENSG00000146833", "ENSG00000130726", "ENSG00000174106",
    "ENSG00000108175", "ENSG00000143344", "ENSG00000014919", "ENSG00000116701",
    "ENSG00000147677", "ENSG00000099953", "ENSG00000161395", "ENSG00000166888",
    "ENSG00000203999", "ENSG00000130725", "ENSG00000206344", "ENSG00000204536",
    "ENSG00000168769", "ENSG00000184939", "ENSG00000270009", "ENSG00000234608",
    "ENSG00000066933", "ENSG00000204644", "ENSG00000230795", "ENSG00000258398",
    "ENSG00000229186", "ENSG00000189143", "ENSG00000131969"
]

# SYNONYMS from the annotated file
GENE_SYNONYMS = {
    "ENSG00000246228": "CASC8", "ENSG00000135862": "LAMC1", "ENSG00000147679": "UTP23",
    "ENSG00000050405": "LIMA1", "ENSG00000178449": "COX14", "ENSG00000127838": "PNKD",
    "ENSG00000261338": "ENSG00000261338", "ENSG00000125845": "BMP2", "ENSG00000179921": "GPBAR1",
    "ENSG00000135926": "TMBIM1", "ENSG00000183520": "UTP11", "ENSG00000212541": "ENSG00000212541",
    "ENSG00000101665": "SMAD7", "ENSG00000185432": "ENSG00000185432",
    "ENSG00000183431": "ENSG00000183431", "ENSG00000124762": "CDKN1A",
    "ENSG00000115073": "ENSG00000115073", "ENSG00000077514": "POLD3",
    "ENSG00000127837": "ENSG00000127837", "ENSG00000232810": "TNF",
    "ENSG00000123810": "ENSG00000123810", "ENSG00000198502": "HLA-DRB5",
    "ENSG00000010327": "ENSG00000010327", "ENSG00000149679": "ENSG00000149679",
    "ENSG00000134825": "ENSG00000134825", "ENSG00000123268": "ENSG00000123268",
    "ENSG00000213533": "ENSG00000213533", "ENSG00000139624": "ENSG00000139624",
    "ENSG00000163464": "ENSG00000163464", "ENSG00000267764": "ENSG00000267764",
    "ENSG00000180871": "ENSG00000180871", "ENSG00000226979": "ENSG00000226979",
    "ENSG00000137834": "ENSG00000137834", "ENSG00000110881": "ENSG00000110881",
    "ENSG00000224397": "ENSG00000224397", "ENSG00000163466": "ENSG00000163466",
    "ENSG00000136280": "ENSG00000136280", "ENSG00000162704": "ENSG00000162704",
    "ENSG00000204592": "ENSG00000204592", "ENSG00000175711": "ENSG00000175711",
    "ENSG00000099326": "ENSG00000099326", "ENSG00000114854": "ENSG00000114854",
    "ENSG00000165171": "ENSG00000165171", "ENSG00000168273": "ENSG00000168273",
    "ENSG00000066117": "ENSG00000066117", "ENSG00000142039": "ENSG00000142039",
    "ENSG00000155287": "ENSG00000155287", "ENSG00000260997": "ENSG00000260997",
    "ENSG00000183208": "ENSG00000183208", "ENSG00000146833": "ENSG00000146833",
    "ENSG00000130726": "ENSG00000130726", "ENSG00000174106": "ENSG00000174106",
    "ENSG00000108175": "ENSG00000108175", "ENSG00000143344": "ENSG00000143344",
    "ENSG00000014919": "ENSG00000014919", "ENSG00000116701": "ENSG00000116701",
    "ENSG00000147677": "ENSG00000147677", "ENSG00000099953": "ENSG00000099953",
    "ENSG00000161395": "ENSG00000161395", "ENSG00000166888": "ENSG00000166888",
    "ENSG00000203999": "ENSG00000203999", "ENSG00000130725": "ENSG00000130725",
    "ENSG00000206344": "ENSG00000206344", "ENSG00000204536": "ENSG00000204536",
    "ENSG00000168769": "ENSG00000168769", "ENSG00000184939": "ENSG00000184939",
    "ENSG00000270009": "ENSG00000270009", "ENSG00000234608": "ENSG00000234608",
    "ENSG00000066933": "ENSG00000066933", "ENSG00000204644": "ENSG00000204644",
    "ENSG00000230795": "ENSG00000230795", "ENSG00000258398": "ENSG00000258398",
    "ENSG00000229186": "ENSG00000229186", "ENSG00000189143": "ENSG00000189143",
    "ENSG00000131969": "ENSG00000131969",
}

# ENSG-only entries that have gene symbols from the annotated file
ENSG_TO_SYMBOL = {
    "ENSG00000261338": None, "ENSG00000212541": None, "ENSG00000185432": None,
    "ENSG00000183431": None, "ENSG00000115073": None, "ENSG00000127837": None,
    "ENSG00000123810": None, "ENSG00000010327": None, "ENSG00000149679": None,
    "ENSG00000134825": None, "ENSG00000123268": None, "ENSG00000213533": None,
    "ENSG00000139624": None, "ENSG00000163464": None, "ENSG00000267764": None,
    "ENSG00000180871": None, "ENSG00000226979": None, "ENSG00000137834": None,
    "ENSG00000110881": None, "ENSG00000224397": None, "ENSG00000163466": None,
    "ENSG00000136280": None, "ENSG00000162704": None, "ENSG00000204592": None,
    "ENSG00000175711": None, "ENSG00000099326": None, "ENSG00000114854": None,
    "ENSG00000165171": None, "ENSG00000168273": None, "ENSG00000066117": None,
    "ENSG00000142039": None, "ENSG00000155287": None, "ENSG00000260997": None,
    "ENSG00000183208": None, "ENSG00000146833": None, "ENSG00000130726": None,
    "ENSG00000174106": None, "ENSG00000108175": None, "ENSG00000143344": None,
    "ENSG00000014919": None, "ENSG00000116701": None, "ENSG00000147677": None,
    "ENSG00000099953": None, "ENSG00000161395": None, "ENSG00000166888": None,
    "ENSG00000203999": None, "ENSG00000130725": None, "ENSG00000206344": None,
    "ENSG00000204536": None, "ENSG00000168769": None, "ENSG00000184939": None,
    "ENSG00000270009": None, "ENSG00000234608": None, "ENSG00000066933": None,
    "ENSG00000204644": None, "ENSG00000230795": None, "ENSG00000258398": None,
    "ENSG00000229186": None, "ENSG00000189143": None, "ENSG00000131969": None,
}

# Known gene biotypes from public data
NONCODING_GENES = {
    "CASC8": "lncRNA",
    "LINC01270": "lncRNA",
    "MAPKAPK5-AS1": "antisense_RNA",
    "PELATON": "lncRNA",
    "GDPGP1": "pseudogene",
    "HCG27": "lncRNA",
    "LOC128966610": "uncharacterized_ncRNA",
    "STIMATE": "protein_coding",  # Actually TMEM110-MUSTN1, protein-coding
    "ENSG00000261338": "unknown_ENSG_only",
    "ENSG00000212541": "unknown_ENSG_only",
    "ENSG00000185432": "unknown_ENSG_only",
    "ENSG00000183431": "unknown_ENSG_only",
    "ENSG00000115073": "unknown_ENSG_only",
    "ENSG00000127837": "unknown_ENSG_only",
    "ENSG00000123810": "unknown_ENSG_only",
    "ENSG00000010327": "unknown_ENSG_only",
    "ENSG00000149679": "unknown_ENSG_only",
    "ENSG00000134825": "unknown_ENSG_only",
    "ENSG00000123268": "unknown_ENSG_only",
    "ENSG00000213533": "unknown_ENSG_only",
    "ENSG00000139624": "unknown_ENSG_only",
    "ENSG00000163464": "unknown_ENSG_only",
    "ENSG00000267764": "unknown_ENSG_only",
    "ENSG00000180871": "unknown_ENSG_only",
    "ENSG00000226979": "unknown_ENSG_only",
    "ENSG00000137834": "unknown_ENSG_only",
    "ENSG00000110881": "unknown_ENSG_only",
    "ENSG00000224397": "unknown_ENSG_only",
    "ENSG00000163466": "unknown_ENSG_only",
    "ENSG00000136280": "unknown_ENSG_only",
    "ENSG00000162704": "unknown_ENSG_only",
    "ENSG00000204592": "unknown_ENSG_only",
    "ENSG00000175711": "unknown_ENSG_only",
    "ENSG00000099326": "unknown_ENSG_only",
    "ENSG00000114854": "unknown_ENSG_only",
    "ENSG00000165171": "unknown_ENSG_only",
    "ENSG00000168273": "unknown_ENSG_only",
    "ENSG00000066117": "unknown_ENSG_only",
    "ENSG00000142039": "unknown_ENSG_only",
    "ENSG00000155287": "unknown_ENSG_only",
    "ENSG00000260997": "unknown_ENSG_only",
    "ENSG00000183208": "unknown_ENSG_only",
    "ENSG00000146833": "unknown_ENSG_only",
    "ENSG00000130726": "unknown_ENSG_only",
    "ENSG00000174106": "unknown_ENSG_only",
    "ENSG00000108175": "unknown_ENSG_only",
    "ENSG00000143344": "unknown_ENSG_only",
    "ENSG00000014919": "unknown_ENSG_only",
    "ENSG00000116701": "unknown_ENSG_only",
    "ENSG00000147677": "unknown_ENSG_only",
    "ENSG00000099953": "unknown_ENSG_only",
    "ENSG00000161395": "unknown_ENSG_only",
    "ENSG00000166888": "unknown_ENSG_only",
    "ENSG00000203999": "unknown_ENSG_only",
    "ENSG00000130725": "unknown_ENSG_only",
    "ENSG00000206344": "unknown_ENSG_only",
    "ENSG00000204536": "unknown_ENSG_only",
    "ENSG00000168769": "unknown_ENSG_only",
    "ENSG00000184939": "unknown_ENSG_only",
    "ENSG00000270009": "unknown_ENSG_only",
    "ENSG00000234608": "unknown_ENSG_only",
    "ENSG00000066933": "unknown_ENSG_only",
    "ENSG00000204644": "unknown_ENSG_only",
    "ENSG00000230795": "unknown_ENSG_only",
    "ENSG00000258398": "unknown_ENSG_only",
    "ENSG00000229186": "unknown_ENSG_only",
    "ENSG00000189143": "unknown_ENSG_only",
    "ENSG00000131969": "unknown_ENSG_only",
}

# Check: which have proper gene symbols vs ENSG-only
def classify():
    with_symbol = []
    ensg_only = []
    for gene in GENES_SYMBOLS:
        if gene.startswith("ENSG"):
            ensg_only.append(gene)
        else:
            with_symbol.append(gene)

    print(f"=== 75 SMR-significant genes ===")
    print(f"With gene symbols: {len(with_symbol)}")
    print(f"ENSG-only (no symbol): {len(ensg_only)}")

    # Among genes with symbols, classify coding vs non-coding
    nc_confirmed = {k: v for k, v in NONCODING_GENES.items() if not k.startswith("ENSG") and v != "protein_coding"}
    coding_likely = []
    nc_list = []

    for g in with_symbol:
        if g in nc_confirmed:
            nc_list.append((g, nc_confirmed[g]))
        else:
            coding_likely.append(g)

    print(f"\n=== With gene symbols: {len(with_symbol)} ===")
    print(f"Confirmed non-coding: {len(nc_list)}")
    for g, t in nc_list:
        print(f"  {g}: {t}")
    print(f"Likely protein-coding: {len(coding_likely)}")

    # Now cross-reference with pQTL coverage
    DECODE_GENES = {"NCF2", "STAT6", "LIMA1", "CCM2", "CERS5", "STAB1", "MZF1",
                    "BMP2", "CDKN1A", "ARPC5", "UBE2M", "CABLES2", "TNF", "LTA"}
    UKBPPP_GENES = {"CDKN1A", "LTA", "NCF2", "TET2", "TNF"}

    # Map ENSG IDs back to symbols where possible
    all_symbols = set(with_symbol)

    decode_hits = DECODE_GENES & all_symbols
    ukbppp_hits = UKBPPP_GENES & all_symbols
    combined = decode_hits | ukbppp_hits

    print(f"\n=== pQTL coverage (among genes with symbols, n={len(with_symbol)}) ===")
    print(f"deCODE: {len(decode_hits)} genes")
    print(f"UKB-PPP: {len(ukbppp_hits)} genes")
    print(f"Combined unique: {len(combined)} genes")
    print(f"Combined genes: {sorted(combined)}")

    # Among likely coding genes
    coding_set = set(coding_likely)
    decode_coding = DECODE_GENES & coding_set
    ukbppp_coding = UKBPPP_GENES & coding_set
    combined_coding = decode_coding | ukbppp_coding

    print(f"\n=== pQTL coverage (among likely protein-coding only, n={len(coding_likely)}) ===")
    print(f"deCODE: {len(decode_coding)}/{len(coding_likely)} ({100*len(decode_coding)//max(1,len(coding_likely))}%)")
    print(f"UKB-PPP: {len(ukbppp_coding)}/{len(coding_likely)} ({100*len(ukbppp_coding)//max(1,len(coding_likely))}%)")
    print(f"Combined: {len(combined_coding)}/{len(coding_likely)} ({100*len(combined_coding)//max(1,len(coding_likely))}%)")

    # Reference: Sun 2024
    print(f"\n=== Benchmark ===")
    print(f"Sun 2024 (J Headache Pain): 750/2,701 = 28%")
    print(f"Our combined: {len(combined_coding)}/{len(coding_likely)} = {100*len(combined_coding)//max(1,len(coding_likely))}%")

    # Missing key coding genes
    missing = coding_set - combined_coding
    print(f"\n=== Coding genes WITHOUT pQTL data ({len(missing)}) ===")
    print(sorted(missing)[:30])

    return {
        "total": len(GENES_SYMBOLS),
        "with_symbol": len(with_symbol),
        "ensg_only": len(ensg_only),
        "noncoding": len(nc_list),
        "coding": len(coding_likely),
        "decode_covered": len(decode_hits),
        "ukbppp_covered": len(ukbppp_hits),
        "combined_pqtl": len(combined),
        "coding_with_pqtl": len(combined_coding),
    }

if __name__ == "__main__":
    r = classify()
