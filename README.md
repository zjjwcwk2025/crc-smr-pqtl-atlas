# Limited plasma proteomic coverage constrains genetic drug-target prioritization in colorectal cancer

Analysis code for the manuscript by Jiajie Zhou (The Affiliated Huaian No.1 People's
Hospital of Nanjing Medical University). The repository contains every script needed to
reproduce the two evidence layers described in the paper: blood cis-eQTL SMR discovery
and plasma pQTL validation.

## Data sources (all public)

| Layer | Resource | Access |
|---|---|---|
| CRC GWAS | GCST90255675 (78,473 cases / 107,143 controls, EUR) | NHGRI-EBI GWAS Catalog |
| Blood cis-eQTL | eQTLGen (n = 31,684) | https://www.eqtlgen.org |
| Tissue eQTL | GTEx v8 (colon) | https://gtexportal.org |
| Plasma pQTL | UKB-PPP (Olink, n ~ 35,000) | Synapse project syn51364943 |
| Plasma pQTL | deCODE (SomaScan v4, n ~ 35,000) | https://www.decode.com/summarydata |
| Replication | FinnGen R13, endpoint C3_COLORECTAL | https://r13.finngen.fi |
| Single-cell | CRC scRNA-seq atlas (n = 50,214 cells, 10 cell types) | local processed object |
| Spatial | Visium CRC, GSE285505 (3,493 spots) | GEO |
| Survival / expression | TCGA COAD + READ | GDC |

## Software

SMR v1.4.2; coloc v5.2.3; susieR v0.14.2; TwoSampleMR; MR-PRESSO; Seurat v5;
R 4.4.3; Python 3.11. Genomic coordinates are hg38.

## Script index

- `scripts/02_parse_smr_results.R` - parse SMR output, Manhattan/QQ (Fig S1)
- `scripts/11_coloc_eqtl.R`, `scripts/11b_coloc_eqtl_all75.R` - Bayesian colocalization (Fig 1b, Fig S4/S5)
- `scripts/15_susie_finemap.R`, `scripts/15b_susie_extended.R`, `scripts/31_susie_sensitivity.R` - SuSiE fine-mapping and sensitivity (Fig S3)
- `scripts/07_ukbppp_pqtl_mr.R`, `scripts/42_decode_pqtl_mr.R`, `scripts/18_pqtl_coloc.R`, `scripts/44_decode_supplement_mr.R` - pQTL MR and colocalization (Fig 3a, Fig S20-S27)
- `scripts/24_mr_sensitivity.R`, `scripts/45_compute_fstatistics.R`, `scripts/48_mr_methods_pleiotropy.R` - instrument strength, MR-PRESSO, MR-Egger
- `scripts/12_finngen_replication.R` - FinnGen directional consistency (Fig S22)
- `scripts/40_gtex_colon_compare.R` - GTEx colon sensitivity (Fig 3b/3c)
- `scripts/16_scrna_celltype_localization.R` - single-cell cell-type specificity (Fig S7)
- `scripts/17_visium_spatial_validation.R`, `scripts/51_spatial_microenvironment.R`, `scripts/55_spatial_gene_crosscorrelation.R` - spatial layers (Fig S8, S9, S14-S16)
- `scripts/13_drug_annotation.R`, `scripts/21_deep_drug_phewas.R`, `scripts/22_phewas_safety.py`, `scripts/33_tier1_idg_druggability.R` - drug annotation, IDG grading, PheWAS (Fig 4, Fig S18/S19)
- `scripts/09_competitor_overlap.R` - overlap with prior CRC target-discovery studies (Fig 5)
- `scripts/rev3_main_figures.R`, `scripts/generate_manuscript_figures.R`, `scripts/theme_pub.R` - figure generation and the unified publication theme
- `scripts/make_bib.py` - reference entries retrieved and verified through NCBI E-utilities

## Reproducing the build

1. Supplementary first, then main (the main text cross-references the supplement):
   `pdflatex crc_smr_supplementary.tex` (twice)
2. `pdflatex crc_smr_atlas_rev3.tex` ; `bibtex crc_smr_atlas_rev3` ; `pdflatex` twice

## Licence

MIT (see LICENSE). Data remain subject to the terms of the original resources.

## Manuscript sources and submission figures

- `manuscript/crc_smr_atlas_rev3.tex` - main text (TeX source)
- `manuscript/crc_smr_supplementary.tex` - supplementary material (19 figures, 9 tables)
- `manuscript/tables/` - Table 1-3 and Supplementary Tables S1-S9
- `manuscript/cover_letter.tex` - cover letter
- `scripts/91_build_submission_figures.py` - builds one combined PDF (and a 600-dpi TIFF)
  per figure - `Figure 1`-`Figure 5` and `Figure S1`-`Figure S19` - into `submission/figures/`,
  for journals that require multi-panel figures as a single file. Panel letters (a), (b), (c)
  are preserved; captions are not embedded.

Compile the manuscript with `pdflatex manuscript/crc_smr_atlas_rev3.tex` (twice, plus
`bibtex crc_smr_atlas_rev3` for the main text). Figures are resolved through
`\graphicspath{{results/figures/}{results/phase6_scrna/}}`.
