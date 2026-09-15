#!/usr/bin/env Rscript
# Phase 2: eQTL vs pQTL 对比矩阵 (UKB-PPP 5 + deCODE 12 genes) + pQTL coloc
library(data.table)
library(coloc)

# ===== 1. Load eQTL SMR (annotated with SYMBOL) =====
cat("=== 1. eQTL SMR ===\n")
smr <- fread("results/phase1_bonferroni_significant_annotated.tsv")
smr <- smr[SYMBOL != "" & !is.na(SYMBOL)]
smr_dt <- smr[, .(SYMBOL, b_SMR, se_SMR, p_SMR, p_HEIDI, nsnp_HEIDI)]

# Genes with pQTL data
ukbppp_genes <- c("CDKN1A","LTA","NCF2","TET2","TNF")
decode_genes <- c("NCF2","STAT6","LIMA1","CCM2","CERS5","STAB1","MZF1",
                  "BMP2","ARPC5","UBE2M","CABLES2","CDKN1A","TNF","LTA",
                  "TNFRSF1B","TNFRSF1A")

smr_sub <- smr_dt[SYMBOL %in% c(ukbppp_genes, decode_genes)]
cat(sprintf("SMR entries for pQTL genes: %d\n", nrow(smr_sub)))

# ===== 2. Load pQTL MR results =====
cat("\n=== 2. pQTL MR ===\n")

ukbppp_mr <- data.table(
  SYMBOL = c("CDKN1A","LTA","NCF2","TET2","TNF"),
  pqtl_source = "UKB-PPP",
  pqtl_wald_b = c(0.130, -0.014, -0.071, NA, -0.295),
  pqtl_wald_p = c(0.020, 0.506, 0.307, NA, 0.0005),
  pqtl_significant = c(TRUE, FALSE, FALSE, FALSE, TRUE)
)

decode_mr <- fread("results/phase2_decode/phase2_decode_pqtl_mr_combined.csv")
# Add genes with 0 instruments that were skipped
decode_no_instr <- data.table(
  SYMBOL = c("CERS5","STAB1","MZF1","BMP2","ARPC5","UBE2M","CABLES2"),
  gene = c("CERS5","STAB1","MZF1","BMP2","ARPC5","UBE2M","CABLES2"),
  n_instruments = 0,
  wald_b = NA, wald_se = NA, wald_pval = NA,
  wald_significant = FALSE, ivw_b = NA, ivw_se = NA, ivw_pval = NA,
  ivw_significant = FALSE, cochran_q_pval = NA
)
decode_mr <- rbind(decode_mr, decode_no_instr, fill=TRUE)

decode_mr_sub <- decode_mr[, .(SYMBOL = gene,
                                pqtl_source = "deCODE",
                                pqtl_wald_b = wald_b,
                                pqtl_wald_p = wald_pval,
                                pqtl_significant = wald_pval < 0.05,
                                n_instruments)]

# ===== 3. Build contrast matrix =====
cat("\n=== 3. eQTL vs pQTL Contrast Matrix ===\n")

# Merge eQTL + UKB-PPP
ukbppp_merged <- merge(smr_sub, ukbppp_mr, by="SYMBOL", all.y=TRUE)

# Merge eQTL + deCODE (only genes with pQTL in smr)
decode_merged <- merge(smr_sub[SYMBOL %in% decode_mr_sub$SYMBOL], decode_mr_sub, by="SYMBOL", all.y=TRUE)

# Combine
all_merged <- rbind(ukbppp_merged, decode_merged, fill=TRUE)

# Direction
all_merged[, eqtl_direction := fifelse(b_SMR > 0, "+", "−")]
all_merged[, pqtl_direction := fifelse(pqtl_wald_b > 0, "+", "−")]
all_merged[, direction_match := fifelse(is.na(eqtl_direction)|is.na(pqtl_direction),
                                        NA, eqtl_direction == pqtl_direction)]

# Classification
all_merged[, classification := "C: pQTL not significant"]
all_merged[pqtl_significant == TRUE & direction_match == TRUE,
           classification := "A: pQTL concordant"]
all_merged[pqtl_significant == TRUE & direction_match == FALSE,
           classification := "B: pQTL discordant"]

setorder(all_merged, -pqtl_significant, pqtl_source, SYMBOL)

cat(sprintf("%-12s %8s %9s %9s %10s %10s %8s %6s %6s %8s %s\n",
            "Gene","Source","b_SMR","p_SMR","b_pQTL","p_pQTL","pQTLsig","eQTL","pQTL","Concord","Class"))
cat(strrep("─", 100), "\n")

for (i in seq_len(nrow(all_merged))) {
  r <- all_merged[i]
  cat(sprintf("%-12s %8s %9.3f %9.2e %9.3f %9.2e %5s %6s %6s %6s %s\n",
              r$SYMBOL, r$pqtl_source,
              r$b_SMR, r$p_SMR,
              r$pqtl_wald_b, r$pqtl_wald_p,
              ifelse(isTRUE(r$pqtl_significant), "✓", "—"),
              ifelse(is.na(r$eqtl_direction), "NA", r$eqtl_direction),
              ifelse(is.na(r$pqtl_direction), "NA", r$pqtl_direction),
              ifelse(is.na(r$direction_match), "—",
                     ifelse(r$direction_match, "✓", "⚠")),
              gsub("..*: ","", r$classification)))
}

cat(sprintf("\n=== Summary ===\n"))
cat(sprintf("Total: %d | pQTL sig: %d | Concordant: %d | Discordant: %d\n",
            nrow(all_merged),
            sum(all_merged$pqtl_significant, na.rm=TRUE),
            sum(all_merged$direction_match & all_merged$pqtl_significant, na.rm=TRUE),
            sum(!all_merged$direction_match & all_merged$pqtl_significant, na.rm=TRUE)))
cat(sprintf("No cis-pQTL (deCODE): %d genes\n",
            sum(all_merged$pqtl_source=="deCODE" & is.na(all_merged$pqtl_wald_b))))
cat(sprintf("No instruments (UKB-PPP): %d genes\n",
            sum(all_merged$pqtl_source=="UKB-PPP" & is.na(all_merged$pqtl_wald_b))))

fwrite(all_merged, "results/phase2_decode/eqtl_pqtl_contrast_matrix.csv")

# ===== 4. deCODE pQTL coloc =====
cat("\n\n=== 4. deCODE pQTL coloc ===\n")

GWAS_FILE <- "data/gwas/GCST90255675.h.tsv.gz"
DECODE_DIR <- "data/decode"
CIS_WINDOW <- 1000000L

gene_list <- list(
  STAT6    = list(chr=12, pos=57491338, file="10372_18_STAT6_STAT6.txt.gz"),
  LIMA1    = list(chr=12, pos=50623450, file="11543_84_LIMA1_LIMA1.txt.gz"),
  CCM2     = list(chr=7,  pos=44976160, file="12347_29_CCM2_CCM2.txt.gz"),
  TNFRSF1A = list(chr=12, pos=6291050,  file="2654_19_TNFRSF1A_TNF_sR_I.txt.gz"),
  TNFRSF1B = list(chr=1,  pos=12155565, file="3152_57_TNFRSF1B_TNF_sR_II.txt.gz")
)

gwas <- fread(GWAS_FILE)
setnames(gwas, old=c("chromosome","base_pair_location","effect_allele","other_allele",
                      "beta","standard_error","effect_allele_frequency","p_value"),
          new=c("chr","bp","ea","oa","beta","se","eaf","pval"))
gwas[, varbeta := se^2]
comp <- function(x) chartr("ACGTacgt","TGCAtgca",x)

coloc_results <- list()

for (gname in names(gene_list)) {
  g <- gene_list[[gname]]
  cat(sprintf("\n--- %s ---\n", gname))

  filter_cmd <- sprintf(
    "zcat %s | awk -F'\\t' 'NR==1 || ($1==\"chr%d\" && $2>=%d && $2<=%d)'",
    file.path(DECODE_DIR, g$file), g$chr, g$pos-CIS_WINDOW, g$pos+CIS_WINDOW)
  pqtl <- fread(cmd=filter_cmd)
  if (nrow(pqtl) <= 1) { cat("No cis SNPs\n"); next }

  setnames(pqtl, old=c("Chrom","Pos","effectAllele","otherAllele","Beta","SE","Pval","N","ImpMAF"),
           new=c("chr_raw","bp","ea","oa","beta","se","pval","n","maf"))
  pqtl[, chr := as.integer(gsub("chr","",chr_raw))]
  pqtl[, varbeta := se^2]
  pqtl[, MAF := pmax(maf, 0.001, na.rm=TRUE)]

  setkey(pqtl, chr, bp)
  gw_chr <- gwas[chr == g$chr]
  setkey(gw_chr, bp)

  matched <- gw_chr[pqtl, on=.(chr, bp), nomatch=0L,
    .(chr, bp.pqtl=i.bp, bp.gwas=bp,
      pqtl.beta=i.beta, pqtl.se=i.se, pqtl.pval=i.pval,
      gwas.beta=beta, gwas.se=se, gwas.pval=pval,
      gwas.ea=ea, gwas.oa=oa, pqtl.ea=i.ea, pqtl.oa=i.oa,
      gwas.maf=eaf, pqtl.maf=i.MAF)]
  cat(sprintf("Matched: %d\n", nrow(matched)))
  if (nrow(matched) < 10) next

  # Harmonize
  p_ea <- matched$pqtl.ea; p_oa <- matched$pqtl.oa
  g_ea <- matched$gwas.ea; g_oa <- matched$gwas.oa
  g_beta <- matched$gwas.beta
  direct <- (p_ea==g_ea & p_oa==g_oa)
  flip   <- (p_ea==g_oa & p_oa==g_ea)
  sw_direct <- (p_ea==comp(g_ea) & p_oa==comp(g_oa))
  sw_flip   <- (p_ea==comp(g_oa) & p_oa==comp(g_ea))
  g_beta[flip|sw_flip] <- -g_beta[flip|sw_flip]
  keep <- direct|flip|sw_direct|sw_flip
  matched <- matched[keep]
  matched$gwas.beta <- g_beta[keep]
  cat(sprintf("Harmonized: %d\n", nrow(matched)))
  if (nrow(matched) < 10) next

  # Filter NA MAF (required for sdY estimation in quant coloc)
  maf_ok <- !is.na(matched$pqtl.maf)
  cat(sprintf("With MAF: %d / %d\n", sum(maf_ok), nrow(matched)))
  matched <- matched[maf_ok]
  if (nrow(matched) < 10) next

  # colocalization
  d1 <- list(beta=matched$pqtl.beta, varbeta=matched$pqtl.se^2,
             snp=paste0("chr",g$chr,":",matched$bp.pqtl),
             position=matched$bp.pqtl, type="quant",
             MAF=matched$pqtl.maf,
             N=max(pqtl$n,na.rm=TRUE))
  d2 <- list(beta=matched$gwas.beta, varbeta=matched$gwas.se^2,
             snp=paste0("chr",g$chr,":",matched$bp.gwas),
             position=matched$bp.gwas, type="cc",
             N=185616, s=78473/185616)  # GWAS_N = 78473 cases + 107143 controls

  res <- tryCatch(coloc.abf(d1, d2), error=function(e){ cat("ERROR:", e$message, "\n"); NULL })
  if (is.null(res)) next

  s <- res$summary
  cat(sprintf("PPH4=%.4f PPH3=%.4f (nsnps=%d)\n",
              s["PP.H4.abf"], s["PP.H3.abf"], nrow(matched)))

  coloc_results[[gname]] <- data.table(
    gene=gname, nsnps=nrow(matched),
    PPH0=s["PP.H0.abf"], PPH1=s["PP.H1.abf"], PPH2=s["PP.H2.abf"],
    PPH3=s["PP.H3.abf"], PPH4=s["PP.H4.abf"])
}

cat("\n\n=== deCODE pQTL coloc Summary ===\n")
if (length(coloc_results) > 0) {
  coloc_combined <- rbindlist(coloc_results)
  coloc_combined[, coloc_pass := PPH4 > 0.8]
  print(coloc_combined, digits=4)
  fwrite(coloc_combined, "results/phase2_decode/phase2_decode_pqtl_coloc.csv")
} else {
  cat("No coloc results.\n")
}

cat("\nDone.\n")
