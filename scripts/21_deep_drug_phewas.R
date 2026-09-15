#!/usr/bin/env Rscript
# Phase 5b: Deep Drug Annotation + Drug Repurposing Score
# Data source: Open Targets Platform v4 GraphQL API (accessible from server)
# DGIdb blocked (HTML response), OT Genetics DNS fails — skipped gracefully
# Integrates: SMR, HEIDI, coloc, SuSiE, TCGA, scRNA/Visium for repurposing score
# Output: results/phase5b_deep_drug_annotation.csv

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

PROJ <- "/ifs1/User/zhouman/project9-v5-crc-atlas"
setwd(PROJ)

# ── 1. Load all existing results ─────────────────────────────────────────────

cat("[1/4] Loading existing results...\n")

smr <- read_tsv("results/phase1_bonferroni_significant_annotated.tsv", show_col_types = FALSE) %>%
  select(gene = Gene, ensg = ENSG_clean, chr = ProbeChr, pos = Probe_bp,
         p_SMR, b_SMR, p_HEIDI_raw = p_HEIDI) %>%
  group_by(gene) %>% slice_min(p_SMR, n = 1) %>% ungroup()

coloc <- read_csv("results/phase3_coloc_top10.csv", show_col_types = FALSE) %>%
  select(ensg, PPH4, PPH3)

susie <- read_csv("results/phase3_susie_finemap.csv", show_col_types = FALSE) %>%
  select(ensg, n_cs, top_pip_snp, converged, tier)

tcga <- read_csv("results/phase6_scrna/TableS_tcga_tumor_vs_normal.csv", show_col_types = FALSE) %>%
  select(gene, log2FC, p_value_tcga = p_value, tcga_sig = significant, tcga_dir = direction)

scrna <- read_csv("results/phase6_scrna/phase6_celltype_specificity.csv", show_col_types = FALSE) %>%
  select(gene = gene_symbol, top_celltype, specificity, expr_fibro, expr_epithelial, expr_immune)

dat <- smr %>%
  left_join(coloc, by = "ensg") %>%
  left_join(susie, by = "ensg") %>%
  left_join(tcga, by = "gene") %>%
  left_join(scrna, by = "gene") %>%
  mutate(heidi_pass = ifelse(!is.na(p_HEIDI_raw) & p_HEIDI_raw >= 0.01, TRUE, FALSE))

cat(sprintf("  %d genes loaded (coloc-passing: %d, SuSiE Tier 1: %d)\n",
            nrow(dat), sum(!is.na(dat$PPH4) & dat$PPH4 > 0.8),
            sum(!is.na(dat$tier) & dat$tier == 1)))

# ── 2. Query Open Targets Platform for ALL genes ─────────────────────────────

cat("\n[2/4] Querying Open Targets Platform...\n")

OT_URL <- "https://api.platform.opentargets.org/api/v4/graphql"

query_ot <- function(ensg, gene_sym) {
  # Rich query: all known drugs, tractability, related targets
  query <- sprintf('
    query {
      target(ensemblId: "%s") {
        approvedSymbol
        tractability {
          label modalities { modality }
        }
        knownDrugs {
          count uniqueDrugs uniqueTargets
          rows {
            phase
            status
            mechanismOfAction
            drugType
            drug { id name maximumClinicalTrialPhase }
            disease { id name }
            approvedIndication
            targetClass
          }
        }
        geneticConstraint {
          constraintType
          oe
          oeLower
          oeUpper
          expLOF
          expMissense
        }
      }
    }', ensg)

  Sys.sleep(0.25)
  tryCatch({
    r <- POST(OT_URL, body = list(query = query), encode = "json", timeout(30))
    if (status_code(r) != 200) return(NULL)
    data <- fromJSON(content(r, "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
    target <- data$data$target
    if (is.null(target) || length(target) == 0) return(NULL)

    known <- target$knownDrugs
    if (is.null(known) || length(known$rows) == 0) {
      tract <- target$tractability
      tract_str <- if (!is.null(tract) && nrow(tract) > 0) {
        paste(tract$label, collapse = "; ")
      } else ""
      return(list(
        n_drugs = 0L, n_drugs_unique = 0L,
        top_drug = NA_character_, top_phase = NA_character_,
        top_indication = NA_character_, top_moa = NA_character_,
        top_drug_type = NA_character_,
        all_drug_names = NA_character_, all_phases = NA_character_,
        n_approved = 0L, n_phase3 = 0L, n_phase2 = 0L, n_phase1 = 0L,
        has_crc_drug = FALSE, crc_drug_text = NA_character_,
        tractability = tract_str
      ))
    }

    rows <- known$rows

    # Phase mapping
    phase_order <- c("4" = 4, "3" = 3, "2" = 2, "1" = 1, "0" = 0)
    phases <- as.character(rows$phase)
    phases_num <- phase_order[phases]
    phases_num[is.na(phases_num)] <- 0

    # Sort by phase
    ord <- order(phases_num, decreasing = TRUE)
    rows <- rows[ord, , drop = FALSE]
    phases <- phases[ord]
    phases_num <- phases_num[ord]

    # Count by phase
    p4 <- sum(phases_num == 4, na.rm = TRUE)
    p3 <- sum(phases_num == 3, na.rm = TRUE)
    p2 <- sum(phases_num == 2, na.rm = TRUE)
    p1 <- sum(phases_num <= 1 & phases_num > 0, na.rm = TRUE)

    # Top drug
    top <- rows[1, ]
    drug_name <- top$drug$name %||% top$drug$id
    top_phase <- top$phase %||% "unknown"
    top_ind <- if (!is.null(top$disease$name)) top$disease$name[1] else NA_character_
    top_moa <- if (!is.null(top$mechanismOfAction) && length(top$mechanismOfAction) > 0) {
      top$mechanismOfAction[1]
    } else NA_character_
    top_type <- top$drugType %||% NA_character_

    all_names <- paste(unique(sapply(rows$drug, function(d) d$name %||% d$id)), collapse = " | ")
    all_phases <- paste(phases, collapse = ", ")

    # CRC-related drugs
    crc_keywords <- c("colorectal", "colon", "rectal", "bowel", "crc")
    has_crc <- any(sapply(rows$disease, function(d) {
      if (is.null(d$name)) return(FALSE)
      any(sapply(crc_keywords, function(kw) grepl(kw, tolower(d$name))))
    }))

    crc_text <- if (has_crc) {
      crc_rows <- rows[sapply(rows$disease, function(d) {
        if (is.null(d$name)) return(FALSE)
        any(sapply(crc_keywords, function(kw) grepl(kw, tolower(d$name))))
      }), ]
      paste(sapply(crc_rows$drug, function(d) d$name %||% d$id), collapse = "; ")
    } else NA_character_

    # Tractability
    tract <- target$tractability
    tract_str <- if (!is.null(tract) && nrow(tract) > 0) {
      paste(unique(tract$label), collapse = "; ")
    } else ""

    list(
      n_drugs = known$count %||% nrow(rows),
      n_drugs_unique = known$uniqueDrugs %||% length(unique(sapply(rows$drug, `[[`, "name"))),
      top_drug = drug_name,
      top_phase = as.character(top_phase),
      top_indication = top_ind,
      top_moa = top_moa,
      top_drug_type = top_type,
      all_drug_names = all_names,
      all_phases = all_phases,
      n_approved = p4, n_phase3 = p3, n_phase2 = p2, n_phase1 = p1,
      has_crc_drug = has_crc,
      crc_drug_text = crc_text,
      tractability = tract_str
    )
  }, error = function(e) {
    message(sprintf("  %s: OT error - %s", gene_sym, e$message))
    NULL
  })
}

# Query all genes
genes_to_query <- dat %>% select(gene, ensg) %>% distinct()
n_total <- nrow(genes_to_query)
ot_results <- list()

for (i in seq_len(n_total)) {
  g <- genes_to_query$gene[i]
  e <- genes_to_query$ensg[i]
  r <- query_ot(e, g)
  if (!is.null(r)) {
    r$gene <- g
    ot_results[[g]] <- r
    cat(sprintf("  [%d/%d] %s: %d drugs (%d approved)\n", i, n_total, g,
                r$n_drugs, r$n_approved))
  } else {
    cat(sprintf("  [%d/%d] %s: no data\n", i, n_total, g))
  }
}

# Flatten
ot_df <- bind_rows(lapply(ot_results, as_tibble))
dat <- left_join(dat, ot_df, by = "gene")

n_druggable <- sum(!is.na(dat$n_drugs) & dat$n_drugs > 0)
n_approved <- sum(!is.na(dat$n_approved) & dat$n_approved > 0)
n_crc <- sum(!is.na(dat$has_crc_drug) & dat$has_crc_drug)
cat(sprintf("\n  Druggable: %d/%d (%.0f%%)\n", n_druggable, nrow(dat), 100*n_druggable/nrow(dat)))
cat(sprintf("  Approved drugs: %d genes\n", n_approved))
cat(sprintf("  CRC-specific: %d genes\n", n_crc))

# ── 3. Drug Repurposing Potential Score ──────────────────────────────────────

cat("\n[3/4] Computing drug repurposing potential scores...\n")

# Data-driven normalization: max values from actual data
max_lfc <- max(abs(dat$log2FC), na.rm=TRUE)
if (max_lfc == 0 || is.na(max_lfc)) max_lfc <- 1.0
max_spec <- max(dat$specificity, na.rm=TRUE)
if (max_spec == 0 || is.na(max_spec)) max_spec <- 1.0

dat <- dat %>% mutate(
  # Normalized component scores (0-1)
  score_smr      = pmin(-log10(p_SMR) / -log10(min(p_SMR)), 1),
  score_coloc    = ifelse(!is.na(PPH4), PPH4, 0),
  score_susie    = case_when(
    !is.na(tier) & tier == 1 ~ 1.0,
    TRUE ~ 0
  ),
  score_heidi    = ifelse(heidi_pass, 1.0, 0.3),
  score_tcga     = ifelse(!is.na(log2FC), pmin(abs(log2FC) / max_lfc, 1), 0),
  score_celltype = ifelse(!is.na(specificity), pmin(specificity / max_spec, 1), 0),
  score_tract    = case_when(
    !is.na(tractability) & grepl("SMALL_MOLECULE|ANTIBODY", tractability) ~ 1.0,
    !is.na(tractability) & tractability != "" ~ 0.5,
    TRUE ~ 0
  ),
  # Drug evidence scoring rubric (definitional: approved > Phase 3 > other > none)
  score_drugs = case_when(
    !is.na(n_approved) & n_approved > 0 ~ 1.0,
    !is.na(n_phase3) & n_phase3 > 0 ~ 0.7,
    !is.na(n_drugs) & n_drugs > 0 ~ 0.4,
    TRUE ~ 0
  )
)

# Weighted composite
dat <- dat %>% mutate(
  score_genetic  = score_smr * 0.20 + score_coloc * 0.35 + score_susie * 0.30 + score_heidi * 0.15,
  score_biology  = score_tcga * 0.50 + score_celltype * 0.30 + score_tract * 0.20,
  score_drug_ev  = score_drugs,
  repurposing    = score_genetic * 0.40 + score_biology * 0.25 + score_drug_ev * 0.35,

  drug_tier = case_when(
    repurposing >= 0.70 ~ "A (High — multi-evidence + approved drugs)",
    repurposing >= 0.50 ~ "B (Strong — coloc/SuSiE + pipeline drugs)",
    repurposing >= 0.30 ~ "C (Moderate — coloc/SuSiE, no drugs yet)",
    TRUE ~ "D (Exploratory — limited evidence)"
  )
)

cat(sprintf("  Tier A: %d, B: %d, C: %d, D: %d\n",
            sum(grepl("^A", dat$drug_tier), na.rm = TRUE),
            sum(grepl("^B", dat$drug_tier), na.rm = TRUE),
            sum(grepl("^C", dat$drug_tier), na.rm = TRUE),
            sum(grepl("^D", dat$drug_tier), na.rm = TRUE)))

# ── 4. Save output ───────────────────────────────────────────────────────────

cat("\n[4/4] Saving results...\n")

out <- dat %>%
  select(
    gene, ensg, chr, pos,
    p_SMR, b_SMR, p_HEIDI = p_HEIDI_raw, heidi_pass,
    PPH4, PPH3, n_cs, top_pip_snp, tier,
    log2FC, tcga_pval = p_value_tcga, tcga_sig, tcga_dir,
    top_celltype, specificity, expr_fibro, expr_epithelial,
    # OT Platform drug data
    n_drugs, n_drugs_unique, n_approved, n_phase3, n_phase2, n_phase1,
    top_drug, top_phase, top_indication, top_moa, top_drug_type,
    all_drug_names, has_crc_drug, crc_drug_text,
    tractability,
    # Scores
    score_smr, score_coloc, score_susie, score_heidi,
    score_tcga, score_celltype, score_tract, score_drugs,
    score_genetic, score_biology, score_drug_ev,
    repurposing, drug_tier
  ) %>%
  arrange(desc(repurposing))

write_csv(out, "results/phase5b_deep_drug_annotation.csv")

# ── 5. Summary report ────────────────────────────────────────────────────────

cat("\n", paste(rep("=", 90), collapse = ""), "\n")
cat("  DEEP DRUG ANNOTATION + REPURPOSING SCORE\n")
cat(paste(rep("=", 90), collapse = ""), "\n\n")

# Top druggable targets
cat("── Highest Repurposing Potential ──\n\n")
top <- out %>% filter(repurposing >= 0.50) %>% arrange(desc(repurposing))
if (nrow(top) > 0) {
  for (i in seq_len(min(nrow(top), 20))) {
    g <- top[i, ]
    col_v <- ifelse(!is.na(g$PPH4) & g$PPH4 > 0.8, "coloc✓", "coloc✗")
    hei_v <- ifelse(g$heidi_pass, "HEIDI✓", "HEIDI✗")
    drug_info <- ifelse(!is.na(g$n_approved) & g$n_approved > 0,
                        sprintf("APPROVED(%d)", g$n_approved),
                        ifelse(!is.na(g$n_drugs) & g$n_drugs > 0,
                               sprintf("%d drugs", g$n_drugs), "none"))
    crc_flag <- ifelse(!is.na(g$has_crc_drug) & g$has_crc_drug, " [CRC]", "")

    cat(sprintf("  %-10s score=%.3f  SMR=%.1e  |  %s %s  |  %s%s\n",
                str_pad(g$gene, 10, "right"), g$repurposing,
                g$p_SMR, col_v, hei_v, drug_info, crc_flag))
    if (!is.na(g$top_drug)) {
      cat(sprintf("    → %s (%s) — %s\n", g$top_drug, g$top_phase, g$top_indication))
    }
  }
}

cat("\n── Drug Pipeline Summary ──\n\n")
cat(sprintf("  Approved drugs:         %d genes\n", sum(!is.na(out$n_approved) & out$n_approved > 0)))
cat(sprintf("  Phase 3:                %d genes\n", sum(!is.na(out$n_phase3) & out$n_phase3 > 0)))
cat(sprintf("  Phase 2:                %d genes\n", sum(!is.na(out$n_phase2) & out$n_phase2 > 0)))
cat(sprintf("  Any pipeline drug:      %d genes\n", sum(!is.na(out$n_drugs) & out$n_drugs > 0)))
cat(sprintf("  CRC-specific drugs:     %d genes\n", sum(!is.na(out$has_crc_drug) & out$has_crc_drug)))
cat(sprintf("  Druggable (tractable):  %d genes\n",
            sum(!is.na(out$tractability) & out$tractability != "")))

cat("\n── Genes with Approved Drugs ──\n\n")
approved <- out %>% filter(!is.na(n_approved) & n_approved > 0) %>% arrange(desc(n_approved))
for (i in seq_len(nrow(approved))) {
  g <- approved[i, ]
  cat(sprintf("  %s: %d approved drugs, top=%s → %s\n", g$gene, g$n_approved, g$top_drug, g$top_indication))
}

cat(sprintf("\nFull results: results/phase5b_deep_drug_annotation.csv (%d genes)\n", nrow(out)))
cat("Done.\n")
