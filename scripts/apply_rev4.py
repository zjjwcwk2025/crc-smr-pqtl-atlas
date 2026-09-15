#!/usr/bin/env python3
"""rev4: insert citations, restructure for journal submission, complete declarations."""
import re, sys

TEX = "manuscript/crc_smr_atlas_rev3.tex"
src = open(TEX, encoding="utf-8").read().replace("\r\n", "\n")
fails = []

def rep(old, new, n=1, label=""):
    global src
    c = src.count(old)
    if c != n:
        fails.append("LITERAL %s: found %d, expected %d :: %s" % (label, c, n, old[:70].replace("\n", "\\n")))
        return
    src = src.replace(old, new)
    print("ok  %-28s %s" % (label, old[:52].replace("\n", " ")))

def repre(pattern, new, n=1, label=""):
    global src
    pat = re.compile(pattern, re.S)
    c = len(pat.findall(src))
    if c != n:
        fails.append("REGEX %s: found %d, expected %d :: %s" % (label, c, n, pattern[:70]))
        return
    src = pat.sub(lambda m: new, src, count=n)
    print("ok  %-28s %s" % (label, pattern[:52]))

# ---------------------------------------------------------------- title page
rep(r"\author{Man Zhou\textsuperscript{1,*}}",
    "\\author{Jiajie Zhou\\textsuperscript{1,*}\\\\[0.4em]\n"
    "\\normalsize\\textsuperscript{1}The Affiliated Huaian No.1 People's Hospital of Nanjing Medical University, Huaian, Jiangsu, China\\\\[0.3em]\n"
    "\\normalsize\\textsuperscript{*}Corresponding author: Jiajie Zhou, MD. E-mail: zjjwcwk2025@njmu.edu.cn; Tel: +86-15805233415}",
    label="author block")

# ---------------------------------------------------------------- structured abstract
rep("\\section*{Abstract}\n%==========================================================================\n\nColorectal cancer (CRC) GWAS have identified",
    "\\section*{Abstract}\n%==========================================================================\n\n\\noindent\\textbf{Background:} Colorectal cancer (CRC) GWAS have identified",
    label="abstract-Background")
rep("pipeline. From 15,582 gene-level",
    "pipeline.\n\n\\noindent\\textbf{Results:} From 15,582 gene-level",
    label="abstract-Results")
rep("Protein-layer coverage, not genetic\nsignificance, is therefore the current limiting step",
    "\n\n\\noindent\\textbf{Conclusions:} Protein-layer coverage, not genetic\nsignificance, is therefore the current limiting step",
    label="abstract-Conclusions")

# ---------------------------------------------------------------- Introduction citations
rep("mortality among all cancers~\\cite{sung2021global}.",
    "mortality among all cancers~\\cite{sung2021global, bray2024global, siegel2023colorectal}.",
    label="intro-burden")
rep("challenges~\\cite{hazelwood2025multi}.",
    "challenges~\\cite{hazelwood2025multi, finan2017druggable, ochoa2021opentargets}.",
    label="intro-druggable")
rep("cis-eQTL data~\\cite{zhu2016integration,\ngusev2016integrative, giambartolomei2014}",
    "cis-eQTL data~\\cite{zhu2016integration, wu2018integrative,\ngusev2016integrative, giambartolomei2014, vosa2021eqtlgen, gtex2020}",
    label="intro-SMR")
rep("pQTL-based drug target validation~\\cite{sun2023plasma, ferkingstad2021large}.",
    "pQTL-based drug target validation~\\cite{sun2023plasma, ferkingstad2021large, pietzner2021mapping},\nand Mendelian randomization of plasma proteins is now an established route to\ncausal target assessment~\\cite{zheng2020phenome, schmidt2020drugtarget, hemani2018mrbase}.",
    label="intro-pqtl-mr")

# ---------------------------------------------------------------- Methods citations
repre(r"CRC GWAS summary statistics were obtained from GCST90255675 \(78,473 cases,\s+107,143 controls, European ancestry\)\.",
      "CRC GWAS summary statistics were obtained from GCST90255675 (78,473 cases,\n107,143 controls, European ancestry), accessed through the NHGRI-EBI GWAS\nCatalog~\\cite{buniello2019gwas}.",
      label="methods-gwas")

rep("Blood cis-eQTL data were obtained from eQTLGen ($n=31,684$).",
    "Blood cis-eQTL data were obtained from eQTLGen ($n=31,684$)~\\cite{vosa2021eqtlgen},\nwith tissue-level eQTL data from the Genotype-Tissue Expression (GTEx)\nProject~\\cite{gtex2020}.",
    label="methods-eqtlgen")
rep("was performed using the SMR software (v1.4.2) with a cis window",
    "was performed using the SMR software (v1.4.2)~\\cite{zhu2016integration} with a cis window",
    label="methods-smr")
rep("HEIDI heterogeneity test was applied to distinguish linkage from\npleiotropy.",
    "HEIDI heterogeneity test was applied to distinguish linkage from\npleiotropy~\\cite{wu2018integrative}.",
    label="methods-heidi")
rep("Colocalization analysis was performed using the coloc R package (v5.2.3),",
    "Colocalization analysis was performed using the coloc R package (v5.2.3)~\\cite{giambartolomei2014},",
    label="methods-coloc")
repre(r"UKB-PPP pQTL data \(Olink, \$n\\approx35,000\$\) were accessed via Synapse\s+\(project syn51364943\)\.",
      "UKB-PPP pQTL data (Olink, $n\\approx35,000$) were accessed via Synapse\n(project syn51364943)~\\cite{sun2023plasma, sudlow2015ukbiobank}.",
      label="methods-ukbppp")
repre(r"deCODE pQTL data \(SomaScan v4, \$n\\approx35,000\$\) were\s+obtained from the deCODE genetics proteomics download portal\.",
      "deCODE pQTL data (SomaScan v4, $n\\approx35,000$) were\nobtained from the deCODE genetics proteomics download portal~\\cite{ferkingstad2021large}.",
      label="methods-decode")
repre(r"Mendelian randomization was performed using\s+Wald ratio \(single instrument\) and inverse-variance weighted \(IVW; multiple\s+instruments\) methods\.",
      "Mendelian randomization was performed using\nWald ratio (single instrument) and inverse-variance weighted (IVW; multiple\ninstruments) methods~\\cite{burgess2013ivw}, with the weighted median~\\cite{bowden2016median},\nweighted mode~\\cite{hartwig2017mode}, MR-Egger~\\cite{bowden2015egger} and\nMR-PRESSO~\\cite{verbanck2018presso} estimators used as sensitivity analyses and\nTwoSampleMR~\\cite{hemani2018mrbase} used for implementation. Reporting follows\nthe STROBE-MR guidance~\\cite{skrivankova2021strobe}.",
      label="methods-mr")
rep("Directional consistency was assessed in FinnGen R13 C3\\_COLORECTAL",
    "Directional consistency was assessed in FinnGen R13~\\cite{kurki2023finngen} C3\\_COLORECTAL",
    label="methods-finngen")
rep("CRC single-cell RNA-seq data were analyzed using Seurat (v5).",
    "CRC single-cell RNA-seq data were analyzed using Seurat~\\cite{hao2021seurat}.",
    label="methods-seurat")
repre(r"Spatial transcriptomics data \(Visium\) were processed using Seurat spatial\s+workflow",
      "Spatial transcriptomics data (Visium)~\\cite{stahl2016spatial} were processed using Seurat spatial\nworkflow",
      label="methods-visium")
rep("\\subsection{Drug annotation and PheWAS}",
    "\\subsection{Drug annotation and PheWAS}\n\nDrug-target associations, target tractability and safety annotations were\nretrieved from the Open Targets Platform~\\cite{ochoa2021opentargets} and\nPharos~\\cite{kelleher2023pharos}, with druggable-genome classification\nfollowing Finan et al.~\\cite{finan2017druggable} and phenome-wide association\nresults interpreted following Denny et al.~\\cite{denny2010phewas}.",
    label="methods-drugphewas")

# ---------------------------------------------------------------- competitor studies
repre(r"four low-tier SMR/MR studies\s+\(Hong 2025, 9 genes; Wang 2025, 4 genes; Wang BMC Cancer supplement, 6 genes;\s+Tian 2025, 8 genes\)",
      "four low-tier SMR/MR studies\n(Hong et al.~2025, 9 genes~\\cite{hong2025integrative}; Wang et al.~2025, 4\ngenes~\\cite{wang2025discovoncol}; Xia and Wang, BMC Cancer, 6\ngenes~\\cite{xia2025bmccancer}; Tian et al.~2025, 8 genes~\\cite{tian2025dbi})",
      label="methods-competitors")
repre(r"\(Hong 2025, Wang 2025, Wang\s+BMC Cancer supplement, Tian 2025, based on smaller\s+FinnGen GWAS\)",
      "(Hong et al.~\\cite{hong2025integrative}, Wang et al.~\\cite{wang2025discovoncol},\nXia and Wang~\\cite{xia2025bmccancer}, Tian et al.~\\cite{tian2025dbi}, based on\nsmaller FinnGen GWAS)",
      label="discussion-competitors")
repre(r"TCGA survival analysis of the six non-MHC Tier 1 genes",
      "TCGA COAD+READ survival analysis~\\cite{weinstein2013tcga} of the six non-MHC Tier 1 genes",
      label="tcga-cite")

# ---------------------------------------------------------------- JTM declarations
def section(name):
    return re.compile(r"\\section\*\{[^}]*" + name + r"[^}]*\}(.*?)(?=\\section\*\{|\\bibliographystyle)", re.S)

def replace_section(name, body, label):
    global src
    pat = section(name)
    hits = pat.findall(src)
    if len(hits) != 1:
        fails.append("SECTION %s: found %d" % (label, len(hits)))
        return
    src = pat.sub(lambda m: body, src, count=1)
    print("ok  %-28s section=%s" % (label, name))

replace_section("Funding",
 "\\section*{Funding}\n%==========================================================================\n\n"
 "This work was supported by the Jiangsu Provincial Health Commission Project\n"
 "(Z2024049) and the Nanjing Medical University Science and Technology\n"
 "Development Fund (NMUB20240108). The funders had no role in study design, data\n"
 "collection and analysis, decision to publish, or preparation of the manuscript.\n\n", "funding")

replace_section("Author Contributions",
 "\\section*{Author Contributions}\n%==========================================================================\n\n"
 "J.Z. conceived and designed the study, performed all analyses, generated the\n"
 "figures and tables, and wrote the manuscript. The author has read and approved\n"
 "the final manuscript.\n\n", "contributions")

replace_section("Competing Interests",
 "\\section*{Competing Interests}\n%==========================================================================\n\n"
 "The author declares that he has no competing interests.\n\n", "competing")

replace_section("Acknowledgements",
 "\\section*{Acknowledgements}\n%==========================================================================\n\n"
 "The author thanks the investigators and participants of the resources used in\n"
 "this study: eQTLGen, the Genotype-Tissue Expression (GTEx) Project, the UK\n"
 "Biobank Pharma Proteomics Project, deCODE genetics, FinnGen, The Cancer Genome\n"
 "Atlas, the NHGRI-EBI GWAS Catalog, and the Open Targets and Pharos teams.\n\n", "acknowledgements")

replace_section("Code Availability",
 "\\section*{Code Availability}\n%==========================================================================\n\n"
 "All analysis scripts used in this study are available from the corresponding\n"
 "author on reasonable request.\n\n", "code-availability")

replace_section("Ethics Statement",
 "\\section*{Ethics Approval and Consent to Participate}\n%==========================================================================\n\n"
 "This study used only publicly available, de-identified summary-level data and\n"
 "previously published de-identified single-cell and spatial transcriptomics data\n"
 "(Supplementary Methods). No new human participants or identifiable material were\n"
 "collected, and institutional ethics approval was therefore not required. Consent\n"
 "for publication is not applicable, as no individually identifiable data are\n"
 "reported.\n\n", "ethics")

# ---------------------------------------------------------------- write out
if fails:
    print("\n=== FAILURES ===")
    for f in fails:
        print(" -", f)
    sys.exit(1)

open(TEX, "w", encoding="utf-8").write(src)
print("\nwrote %s (%d bytes)" % (TEX, len(src)))
print("cite commands in file: %d" % len(re.findall(r"\\cite[a-z]*\{", src)))