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
R 4.4.3; Python 3.11; Poppler (pdftoppm, pdftotext). Genomic coordinates are hg38.

## Analysis scripts

- `scripts/01_convert_gwas_to_ma.R` - harmonise the CRC GWAS summary statistics for SMR
- `scripts/02_parse_smr_results.R` - parse SMR output; Manhattan and QQ plots (Fig S1)
- `scripts/03_convert_vcf_to_plink_eur.sh`, `scripts/merge_plink_eur.sh` - reference-panel preparation for HEIDI
- `scripts/04_synapse_pqtl_download.py`, `scripts/08_download_decode.sh` - pQTL panel download
- `scripts/05_heidi_from_vcf.py`, `scripts/06_run_heidi_75genes.sh` - HEIDI test for every significant probe
- `scripts/06_classify_genes.py` - probe classification and gene-level annotation
- `scripts/11_coloc_eqtl.R`, `scripts/11b_coloc_eqtl_all75.R` - Bayesian colocalization at the eQTL layer (Fig 1c)
- `scripts/15_susie_finemap.R`, `scripts/15b_susie_extended.R`, `scripts/31_susie_sensitivity.R` - SuSiE fine-mapping and sensitivity (Fig S3)
- `scripts/07_ukbppp_pqtl_mr.R`, `scripts/42_decode_pqtl_mr.R`, `scripts/18_pqtl_coloc.R`, `scripts/43_decode_coloc.R`, `scripts/44_decode_supplement_mr.R` - pQTL MR and colocalization (Fig 3a; Fig S6-S13; Tables S4, S5)
- `scripts/24_mr_sensitivity.R`, `scripts/45_compute_fstatistics.R`, `scripts/48_mr_methods_pleiotropy.R` - instrument strength, MR-Egger, leave-one-out (Fig S5, S10-S13)
- `scripts/49_mr_presso_table3.R` - MR-PRESSO for the nine genes of Table 3 at `NbDistribution = 5000`; at the package default of 1000 the smallest resolvable outlier-test p-value is `n/1000` (0.061-0.228 for these genes), above 0.05, so the default cannot flag outliers at these instrument counts. TNF (n = 2288) and LTA (n = 10213) are not feasible and are recorded as skipped in Table 3.
- `scripts/50_fast_presso_reference.R` - vectorised re-implementation of the MR-PRESSO core algorithm, used for the MR-PRESSO columns of Table 3. Verified against `MRPRESSO::mr_presso` at seed 20260921 and NbDistribution = 5000: RSSobs agrees to <= 2.3e-13 and the flagged instrument sets are identical for CCM2, LIMA1, STAT6, CDKN1A, TNFRSF1A and NCF2. The package call is about 1.5 h per 60 instruments on a single core, so this script is what produced the reported values for all nine genes.
- `scripts/12_finngen_replication.R` - FinnGen directional consistency (Fig S17)
- `scripts/40_gtex_colon_compare.R` - GTEx colon sensitivity (Fig 3b, 3c)
- `scripts/30_locus_count.R`, `scripts/32_tier1_smr_or.R` - locus counting and Tier-1 SMR odds ratios (Fig S2, S4)
- `scripts/16_scrna_celltype_localization.R` - single-cell cell-type specificity (Fig S14)
- `scripts/17_visium_spatial_validation.R`, `scripts/51_spatial_microenvironment.R`, `scripts/55_spatial_gene_crosscorrelation.R`, `scripts/58_spatial_tier1_panel.R` - spatial layers (Fig S15, S16; Tables S6, S7)
- `scripts/13_drug_annotation.R`, `scripts/21_deep_drug_phewas.R`, `scripts/22_phewas_safety.py`, `scripts/33_tier1_idg_druggability.R` - drug annotation, IDG tractability grading, PheWAS safety (Fig 4; Fig S18, S19; Table S8)
- `scripts/09_competitor_overlap.R` - overlap with prior CRC target-discovery studies (Fig 5; Table S9)
- `scripts/19_tcga_survival.R`, `scripts/20_tcga_paired_de.R`, `scripts/20_tcga_paired_tumor_normal.R`, `scripts/26_tcga_pancancer.R`, `scripts/25_pathway_enrichment.R`, `scripts/27_string_ppi.R` - supporting expression, survival, pathway and network analyses
- `scripts/make_bib.py` - reference entries retrieved and verified through NCBI E-utilities

## Figure assembly

- `scripts/theme_pub.R` - the single publication theme (fonts, sizes, palettes) used by every figure
- `scripts/rev3_main_figures.R` - draws every main figure (Fig 1-5) and every supplementary figure (Fig S1-S19) into `results/figures_rev3/`
- `scripts/generate_manuscript_figures.R` - upstream versions of several panels, kept for provenance
- `scripts/redraw_fig1_figs1.R` - Manhattan and QQ panel (Fig S1)
- `scripts/63_fig5b_upset.R` - three-set overlap of this study, Chen 2024 and Hazelwood 2025 (Fig 5b)
- `scripts/fig5a_regen.R` - competitor-overlap bar (Fig 5a)
- `scripts/figS26S27_regen.R` - MR method comparison and funnel/pleiotropy panels (Fig S12, S13); reads MR-PRESSO at NbDistribution = 5000, which is required to resolve n/NbDistribution <= 0.046 at these instrument counts
- `scripts/71_regen_small_text.R` - re-renders Fig S1 and Fig S16 at their final printed width so that no label falls below 7 pt
- `scripts/73_export_fig1a_panel.py` - exports the TikZ conceptual overview as a standalone vector panel (Fig 1a)
- `scripts/trim_figs.py`, `scripts/split_docs.py`, `scripts/fix_cap.py`, `scripts/fix_supp.py` - manuscript and caption housekeeping

## Submission files

- `scripts/91_build_submission_figures.py` - builds one combined PDF per figure
  (`Figure1`-`Figure5`, `FigureS1`-`FigureS19`) into `submission/figures/`, for journals
  that require multi-panel figures as a single file. Panel letters (a), (b), (c) are
  preserved; captions are not embedded. Set `ONLY=Figure3,Figure5` to rebuild a subset.
- `scripts/74_export_tifs.py` - exports 600-dpi LZW TIFFs of every combined figure
  (`submission/figures/tif_figures/`) and of every individual panel
  (`submission/figures/tif_panels/`, 35 panels including Fig 1a). Resolution steps down
  through 500/400/300 dpi when a file would exceed the 10 MB upload limit.
- `scripts/75_build_upload_pack.py` - writes `submission/upload_pack/main_text/` and
  `submission/upload_pack/supplementary/`. Each directory is flat (no subdirectories)
  and carries its own copies of the .tex, .bib, .bbl, tables and figure panels, so it
  compiles on its own inside the single directory that a journal's LaTeX service uses.
  `manuscript/...` input paths are rewritten to plain file names for this.

The generated PDFs and TIFFs are not tracked here; both scripts recreate them from the
figure sources.

## Reproducing the build

Figures, then the two documents. They are independent: the main text quotes the
supplementary item numbers literally (Supplementary Fig. S1-S19, Table S1-S9) instead of
reading the supplement's .aux file, so either document can be compiled on its own.

1. `Rscript scripts/rev3_main_figures.R`
2. `python3 scripts/91_build_submission_figures.py`
3. from the repository root:
   `pdflatex -interaction=nonstopmode manuscript/crc_smr_supplementary.tex` (twice)
4. `pdflatex -interaction=nonstopmode manuscript/crc_smr_atlas_rev3.tex` ;
   `bibtex crc_smr_atlas_rev3` ; `pdflatex` twice
5. `python3 scripts/75_build_upload_pack.py` - assembles the flat upload pack

Figures are resolved through
`\graphicspath{{results/figures_rev3/}{results/figures/}{results/phase6_scrna/}}`.

## Manuscript sources

- `manuscript/crc_smr_atlas_rev3.tex` - main text
- `manuscript/crc_smr_supplementary.tex` - supplementary material (19 figures, 9 tables)
- `manuscript/tables/` - Table 1-3 and Supplementary Tables S1-S9
- `manuscript/cover_letter.tex` - cover letter
- `manuscript/ai_disclosure.md` - AI-use disclosure

## Figure and table map

| Item | Source file | Script |
|---|---|---|
| Fig 1a | TikZ block in the main text | `73_export_fig1a_panel.py` |
| Fig 1b | Fig1b_study_design.pdf | `rev3_main_figures.R` |
| Fig 1c | Fig1c_coloc_pph4.pdf | `rev3_main_figures.R` |
| Fig 2a | Fig2a_coverage_funnel.pdf | `rev3_main_figures.R` |
| Fig 2b | Fig2b_eqtl_pqtl_classification.pdf | `rev3_main_figures.R` |
| Fig 2c | Fig2c_eqtl_pqtl_scatter.pdf | `rev3_main_figures.R` |
| Fig 3a | Fig3a_ccm2_convergence.pdf | `rev3_main_figures.R` |
| Fig 3b | Fig3b_gtex_colon_scatter.pdf | `rev3_main_figures.R` |
| Fig 3c | Fig3c_gtex_colon_forest.pdf | `rev3_main_figures.R` |
| Fig 4a | Fig4a_phewas_safety.pdf | `rev3_main_figures.R` |
| Fig 4b | Fig4b_drug_repurposing.pdf | `rev3_main_figures.R` |
| Fig 5a | Fig5a_competitor_overlap.pdf | `fig5a_regen.R` |
| Fig 5b | Fig5b_gene_overlap_venn.pdf | `63_fig5b_upset.R` |
| Fig S1 | FigS1a_manhattan.pdf, FigS1b_qq_plot.pdf | `71_regen_small_text.R` |
| Fig S2 | FigS2_tier1_regional_plots.pdf | `rev3_main_figures.R` |
| Fig S3 | FigS3_susie_sensitivity.pdf | `rev3_main_figures.R` |
| Fig S4 | FigS4_smr_or_forest.pdf | `rev3_main_figures.R` |
| Fig S5 | FigS5_mr_sensitivity.pdf | `rev3_main_figures.R` |
| Fig S6 | FigS6_ukbppp_pqtl_coloc.pdf | `rev3_main_figures.R` |
| Fig S7 | FigS7_decode_pqtl_coloc.pdf | `rev3_main_figures.R` |
| Fig S8 | FigS8_pqtl_mr_decode_forest.pdf | `rev3_main_figures.R` |
| Fig S9 | FigS9_pqtl_mr_ukbppp_forest.pdf | `rev3_main_figures.R` |
| Fig S10 | FigS10_pqtl_mr_scatter.pdf | `rev3_main_figures.R` |
| Fig S11 | FigS11_pqtl_mr_loo.pdf | `rev3_main_figures.R` |
| Fig S12 | FigS12_pqtl_mr_methods_forest.pdf | `figS26S27_regen.R` |
| Fig S13 | FigS13_pqtl_mr_funnel_pleiotropy.pdf | `figS26S27_regen.R` |
| Fig S14 | FigS14_celltype_dotplot.pdf | `rev3_main_figures.R` |
| Fig S15 | FigS15_spatial_tier1_genes.pdf | `rev3_main_figures.R` |
| Fig S16 | FigS16a/b/c (elbow, ME-zone bar, ME Tier-1 heatmap) | `71_regen_small_text.R` |
| Fig S17 | FigS17_finngen_replication_forest.pdf | `rev3_main_figures.R` |
| Fig S18 | FigS18_idg_druggability.pdf | `rev3_main_figures.R` |
| Fig S19 | FigS19_phewas_safety_scores.pdf | `rev3_main_figures.R` |
| Table 1 | Table1_tier1_genes.tex | `rev3_tables.R` |
| Table 2 | Table2_pqtl_coverage.tex | `rev3_tables.R` |
| Table 3 | Table3_mr_sensitivity.tex | `rev3_tables.R` |
| Table S1 | TableS1_smr_top10.tex | `rev3_tables.R` |
| Table S2 | TableS2_protein_coding.tex | `rev3_tables.R` |
| Table S3 | TableS3_susie_sensitivity.tex | `rev3_tables.R` |
| Table S4 | TableS4_ukbppp_pqtl_coloc.tex | `rev3_tables.R` |
| Table S5 | TableS5_decode_pqtl_coloc.tex | `rev3_tables.R` |
| Table S6 | TableS6_me_zone_characterization.tex | `rev3_tables.R` |
| Table S7 | TableS7_spatial_crosscorrelation.tex | `rev3_tables.R` |
| Table S8 | TableS8_phewas_safety.tex | `rev3_tables.R` |
| Table S9 | TableS9_chen2024_gene_list.tex | `rev3_tables.R` |

## Licence

MIT (see LICENSE). Data remain subject to the terms of the original resources.
