# AI Usage Disclosure Statement

## For: Nature Communications

---

## Statement of AI-Assisted Tools Used in This Research

In accordance with the journal's policy on the use of artificial intelligence (AI)-assisted technologies, the authors disclose the following:

### 1. Writing and Language Assistance
Claude (Anthropic) was used for language editing assistance, including improving sentence clarity, grammar correction, and ensuring consistent scientific terminology throughout the manuscript. All AI-suggested text was critically reviewed, verified, and approved by the authors. No text was generated de novo by AI without author review.

### 2. Code Generation and Debugging
Claude (Anthropic) was used to assist in writing, debugging, and optimizing R and Python scripts for data analysis, including:
- SMR analysis pipeline scripts (shell/R)
- Colocalization and SuSiE fine-mapping scripts (R)
- Mendelian randomization analysis scripts (R)
- Drug annotation and PheWAS screening scripts (Python)
- Single-cell and spatial transcriptomics analysis scripts (R)
- Figure generation scripts (R/ggplot2)
- Pipeline auditing and quality control scripts

All AI-generated code was reviewed, tested, and validated by the authors. All analytical decisions (statistical thresholds, parameter choices, interpretation of results) were made by the authors.

### 3. Manuscript Review
Claude (Anthropic) was used to simulate a peer review process for manuscript quality improvement prior to submission. This included:
- Methodological consistency checking
- Literature coverage assessment
- Statistical reporting standards verification
- Claims-narrative consistency review

All reviewer-simulated feedback was evaluated by the authors, and all revisions were authored by the human authors. No AI-generated review comments were accepted without critical evaluation.

### 4. Data Analysis
All primary data analyses were performed by the authors using standard open-source bioinformatics tools (SMR v1.4.2, PLINK2, coloc v5.2.3, susieR v0.14.2, Seurat v5, TwoSampleMR v0.7.6). AI assistance was limited to code generation and debugging as described in Section 2.

### 5. Author Responsibility
The authors take full responsibility for the accuracy, integrity, and originality of the work. All data, analyses, and conclusions have been independently verified by the authors. The use of AI tools did not replace human scientific judgment, interpretation, or ethical decision-making at any stage of the research.

### 6. Tools and Versions
- **Claude** (Anthropic, various versions): Writing assistance, code generation, manuscript review simulation
- **SMR** v1.4.2: Summary-based Mendelian randomization
- **PLINK2** v2.0.0a7: Genotype data processing
- **coloc** R package v5.2.3: Bayesian colocalization
- **susieR** R package v0.14.2: Sum of Single Effects fine-mapping
- **Seurat** R package v5: Single-cell and spatial transcriptomics analysis
- **TwoSampleMR** R package v0.7.6: Two-sample Mendelian randomization

---

*This disclosure is provided in accordance with Nature Portfolio's editorial policies on AI-assisted technologies. The authors confirm that no generative AI tools were used to fabricate or manipulate data, create fraudulent references, or circumvent peer review.*
