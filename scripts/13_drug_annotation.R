#!/usr/bin/env Rscript
# Phase 5: Drug-target annotation via Open Targets Platform v4 API
# Query: all 75 SMR-significant genes → known clinical drugs
# Output: results/phase5_drug_annotation.csv

library(data.table)
library(dplyr)
library(jsonlite)

# === 1. Load SMR genes with Ensembl IDs ===
cat("Loading SMR genes...\n")
smr <- fread("results/phase1_bonferroni_significant_annotated.tsv") %>%
  filter(SYMBOL != "" & !is.na(SYMBOL)) %>%
  select(ensg = ENSG_clean, gene = SYMBOL, gene_chr = ProbeChr, gene_pos = Probe_bp, p_SMR, b_SMR)

# Add coloc results if available
coloc_file <- "results/phase3_coloc_top10.csv"
if(file.exists(coloc_file)) {
  coloc <- fread(coloc_file) %>% select(ensg, PPH4)
  smr <- smr %>% left_join(coloc, by="ensg")
} else {
  smr$PPH4 <- NA_real_
}

# Keep unique genes
genes <- smr %>% filter(!duplicated(gene))
cat(sprintf("Gene set: %d unique genes\n", nrow(genes)))

# === 2. Query Open Targets API (batch) ===
# OT v4: target(ensemblId:) → drugAndClinicalCandidates
cat("Querying Open Targets API...\n")

query_one_gene <- function(ensg, gene) {
  cat(sprintf("  Querying %s (%s)...", gene, ensg))

  query <- sprintf('query{target(ensemblId:"%s"){id approvedSymbol drugAndClinicalCandidates{rows{drug{name maximumClinicalStage drugType mechanismsOfAction{rows{mechanismOfAction}} indications{rows{disease{name}}}}}} associatedDiseases{count rows{disease{name} datasourceScores{componentId score}}}}}}', ensg)

  # Use curl for HTTP
  cmd <- sprintf('curl -s --max-time 30 "https://api.platform.opentargets.org/api/v4/graphql" -H "Content-Type: application/json" -d \'%s\'', query)
  resp <- tryCatch(system(cmd, intern=TRUE, ignore.stderr=TRUE), error=function(e) NULL)

  if(is.null(resp) || length(resp) == 0) {
    cat(" FAILED\n")
    return(NULL)
  }

  json_str <- paste(resp, collapse="")
  parsed <- tryCatch(fromJSON(json_str), error=function(e) NULL)

  if(is.null(parsed) || !is.null(parsed$errors)) {
    cat(" ERROR\n")
    return(NULL)
  }

  target <- parsed$data$target
  drugs <- target$drugAndClinicalCandidates$rows

  if(length(drugs) == 0) {
    cat(" 0 drugs\n")
    return(data.frame(
      gene = gene, ensg = ensg,
      n_drugs = 0, top_drug = NA_character_, top_stage = NA_character_,
      top_type = NA_character_, top_moa = NA_character_,
      top_indication = NA_character_, all_drugs = NA_character_,
      n_approved = 0, n_phase4 = 0, n_phase3 = 0, n_phase2 = 0,
      n_phase1 = 0, n_phase0 = 0,
      stringsAsFactors = FALSE
    ))
  }

  n_drugs <- nrow(drugs)
  cat(sprintf(" %d drugs\n", n_drugs))

  # Extract top drug info
  drugs <- drugs %>%
    mutate(
      name = drug$name,
      stage = drug$maximumClinicalStage,
      dtype = drug$drugType,
      moa = sapply(drug$mechanismsOfAction$rows, function(x) {
        if(length(x) > 0 && nrow(x) > 0) paste(x$mechanismOfAction[1:min(2,nrow(x))], collapse="; ") else NA_character_
      }),
      top_indication = sapply(drug$indications$rows, function(x) {
        if(length(x) > 0 && nrow(x) > 0) x$disease$name[1] else NA_character_
      })
    )

  # Stage counts
  stage_counts <- table(factor(drugs$stage, levels=c(0,1,2,3,4)))

  data.frame(
    gene = gene, ensg = ensg,
    n_drugs = n_drugs,
    top_drug = drugs$name[1],
    top_stage = as.character(drugs$stage[1]),
    top_type = as.character(drugs$drugType[1]),
    top_moa = as.character(drugs$moa[1]),
    top_indication = as.character(drugs$top_indication[1]),
    all_drugs = paste(drugs$name, collapse=", "),
    all_stages = paste(drugs$stage, collapse=","),
    n_approved = as.integer(stage_counts["4"]),
    n_phase4   = as.integer(stage_counts["4"]),
    n_phase3   = as.integer(stage_counts["3"]),
    n_phase2   = as.integer(stage_counts["2"]),
    n_phase1   = as.integer(stage_counts["1"]),
    n_phase0   = as.integer(stage_counts["0"]),
    stringsAsFactors = FALSE
  )
}

# Query all genes (limit to first 30 for speed, expand later)
# Focus on coloc-significant genes first
genes_q <- genes %>% arrange(p_SMR) %>% head(30)

results <- list()
for(i in 1:nrow(genes_q)) {
  res <- query_one_gene(genes_q$ensg[i], genes_q$gene[i])
  if(!is.null(res)) {
    results[[length(results)+1]] <- res
  }
  Sys.sleep(0.3)  # rate limiting
}

cat(sprintf("\nQueried %d genes successfully\n", length(results)))

# === 3. Merge with SMR results ===
if(length(results) > 0) {
  drug_table <- rbindlist(results)

  # Merge with SMR stats
  final <- genes_q %>%
    left_join(drug_table, by=c("gene","ensg")) %>%
    mutate(
      n_drugs = ifelse(is.na(n_drugs), 0, n_drugs),
      # Drug repurposing potential score
      drug_score = case_when(
        n_drugs > 0 & n_approved > 0 ~ "Known drug target",
        n_drugs > 0 & n_phase3 > 0 ~ "Phase 3+",
        n_drugs > 0 ~ "Pipeline drug target",
        TRUE ~ "Novel biology"
      ),
      # CRC-specific check
      crc_drug = grepl("colorect|colon|rectal|bowel", top_indication, ignore.case=TRUE)
    ) %>%
    select(gene, ensg, p_SMR, b_SMR, PPH4,
           n_drugs, top_drug, top_stage, top_type, top_moa, top_indication,
           crc_drug, drug_score, all_drugs, all_stages,
           n_approved, n_phase3, n_phase2, n_phase1, n_phase0)

  # === 4. Output ===
  final <- final %>% arrange(p_SMR)
  fwrite(final, "results/phase5_drug_annotation.csv")

  cat("\n\n=== Drug Annotation Summary ===\n")
  cat(sprintf("Total genes queried: %d\n", nrow(final)))
  cat(sprintf("With known drugs: %d (%.0f%%)\n",
              sum(final$n_drugs > 0), 100*sum(final$n_drugs > 0)/nrow(final)))
  cat(sprintf("Approved drugs: %d genes\n", sum(final$n_approved > 0)))
  cat(sprintf("CRC-specific drugs: %d genes\n", sum(final$crc_drug, na.rm=TRUE)))

  cat("\n\n=== Top 15 Drug-Target Results ===\n")
  top_druggable <- final %>%
    filter(n_drugs > 0) %>%
    arrange(p_SMR) %>%
    head(15)

  cat(sprintf("%-10s %12s %12s %6s %-20s %-10s %-6s %s\n",
              "Gene", "p_SMR", "b_SMR", "Drugs", "Top Drug", "Stage", "CRC?", "Drug Score"))
  for(i in 1:nrow(top_druggable)) {
    r <- top_druggable[i, ]
    cat(sprintf("%-10s %12.2e %12.4f %6d %-20s %-10s %-6s %s\n",
                r$gene, r$p_SMR, r$b_SMR, r$n_drugs,
                substr(r$top_drug, 1, 20), r$top_stage,
                ifelse(!is.na(r$crc_drug) && r$crc_drug, "YES", "no"),
                r$drug_score))
  }

  cat("\n✅ Saved: results/phase5_drug_annotation.csv\n")
} else {
  cat("\n⚠ No drug data retrieved\n")
}
