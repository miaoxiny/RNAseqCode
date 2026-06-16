# =============================================================
#  07_DEG_extraction.R — Step 7 V2
#  Strict policy: RAW LFC from DESeq2 Wald test (no apeglm)
#  Uses the existing dds_fitted.rds from Step 6 (no upstream changes)
#  
#  Pipeline:
#    dds_fitted.rds (from Step 6)
#      ↓
#    results() with contrast — RAW LFC (no shrinkage)
#      ↓
#    Annotate (ENSEMBL → SYMBOL / GENENAME / ENTREZID)
#      ↓
#    Classify DEG (STRICT vs RELAXED)
#      ↓
#    Save tables (overwrites previous shrunk versions)
# =============================================================

library(DESeq2)
library(AnnotationDbi)
library(org.Rn.eg.db)
library(dplyr)
library(tibble)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")

# ─── Thresholds ───
padj_cutoff <- 0.05
lfc_cutoff  <- 0.585    # 1.5x fold change

# ─── Load fitted dds from Step 6 ───
dds <- readRDS(file.path(results_dir, "dds_fitted.rds"))

cat("=== Loaded dds_fitted ===\n")
cat("Genes:  ", nrow(dds), "\n")
cat("Samples:", ncol(dds), "\n")
cat("Design: ", deparse(design(dds)), "\n")
cat("Model fitted (dispersions present):",
    !is.null(dispersions(dds)), "\n\n")

# =============================================================
#  Step 7.1: Extract RAW results (no lfcShrink)
# =============================================================

cat("=== Extracting raw DESeq2 results (no shrinkage) ===\n\n")

# Primary: Netrin vs DIV2 control
res_netrin_raw <- results(
  dds,
  contrast = c("condition", "DIV2_Netrin", "DIV2"),
  alpha    = padj_cutoff
)

# Secondary: DIV2 vs DIV1
res_div_raw <- results(
  dds,
  contrast = c("condition", "DIV2", "DIV1"),
  alpha    = padj_cutoff
)

cat("--- Primary: DIV2 + Netrin-1 vs DIV2 ---\n")
summary(res_netrin_raw)

cat("\n--- Secondary: DIV2 vs DIV1 ---\n")
summary(res_div_raw)

# =============================================================
#  Step 7.2: Annotate (ENSEMBL -> SYMBOL / GENENAME / ENTREZID)
# =============================================================

cat("\n=== Annotating with org.Rn.eg.db ===\n")

annotate_results <- function(res, contrast_name) {
  df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_id")
  
  df$symbol    <- mapIds(org.Rn.eg.db,
                         keys      = df$ensembl_id,
                         column    = "SYMBOL",
                         keytype   = "ENSEMBL",
                         multiVals = "first")
  df$gene_name <- mapIds(org.Rn.eg.db,
                         keys      = df$ensembl_id,
                         column    = "GENENAME",
                         keytype   = "ENSEMBL",
                         multiVals = "first")
  df$entrez_id <- mapIds(org.Rn.eg.db,
                         keys      = df$ensembl_id,
                         column    = "ENTREZID",
                         keytype   = "ENSEMBL",
                         multiVals = "first")
  
  df_clean <- df %>%
    select(ensembl_id, symbol, gene_name, entrez_id,
           baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    arrange(padj)
  
  cat(sprintf("  %s: symbol assigned %d/%d (%.1f%%)\n",
              contrast_name,
              sum(!is.na(df_clean$symbol)),
              nrow(df_clean),
              100 * sum(!is.na(df_clean$symbol)) / nrow(df_clean)))
  
  df_clean
}

res_netrin_df <- annotate_results(res_netrin_raw, "Netrin vs DIV2")
res_div_df    <- annotate_results(res_div_raw,    "DIV2 vs DIV1  ")

# =============================================================
#  Step 7.3: Classify DEG (RAW LFC thresholds)
# =============================================================

classify_deg <- function(df) {
  df$direction <- "ns"
  df$direction[!is.na(df$padj) & df$padj < padj_cutoff &
                 df$log2FoldChange >  lfc_cutoff]  <- "up"
  df$direction[!is.na(df$padj) & df$padj < padj_cutoff &
                 df$log2FoldChange < -lfc_cutoff]  <- "down"
  df$direction <- factor(df$direction, levels = c("up", "down", "ns"))
  df
}

res_netrin_df <- classify_deg(res_netrin_df)
res_div_df    <- classify_deg(res_div_df)

# =============================================================
#  Step 7.4: Count DEGs (two tiers)
# =============================================================

cat("\n=== DEG counts ===\n")
cat(sprintf("Thresholds:\n"))
cat(sprintf("  STRICT:  padj < %.2f  &  |log2FC| > %.3f\n",
            padj_cutoff, lfc_cutoff))
cat(sprintf("  RELAXED: padj < %.2f only\n\n", padj_cutoff))

cat("STRICT — Primary contrast (Netrin vs DIV2):\n")
print(table(res_netrin_df$direction))

cat("\nSTRICT — Secondary contrast (DIV2 vs DIV1):\n")
print(table(res_div_df$direction))

relaxed_netrin <- sum(!is.na(res_netrin_df$padj) & 
                        res_netrin_df$padj < padj_cutoff)
relaxed_div    <- sum(!is.na(res_div_df$padj) & 
                        res_div_df$padj < padj_cutoff)

cat(sprintf("\nRELAXED (padj<0.05) summary:\n"))
cat(sprintf("  Netrin vs DIV2:  %d genes\n", relaxed_netrin))
cat(sprintf("  DIV2 vs DIV1:    %d genes\n", relaxed_div))

# =============================================================
#  Step 7.5: Save annotated results (OVERWRITES previous)
# =============================================================

saveRDS(res_netrin_df, file.path(results_dir, "res_netrin_vs_div2.rds"))
saveRDS(res_div_df,    file.path(results_dir, "res_div2_vs_div1.rds"))

cat("\n=== Saved (raw LFC, no shrinkage) ===\n")
cat("  res_netrin_vs_div2.rds   ← overwrote\n")
cat("  res_div2_vs_div1.rds     ← overwrote\n")

# =============================================================
#  Step 7.6: Sanity check — Top 10 DEG by padj (Netrin vs DIV2)
# =============================================================

cat("\n=== Top 10 DEG by padj (Netrin vs DIV2) ===\n")
top10 <- res_netrin_df %>%
  filter(direction != "ns") %>%
  arrange(padj) %>%
  head(10) %>%
  select(symbol, gene_name, baseMean, log2FoldChange, padj, direction) %>%
  mutate(
    baseMean       = round(baseMean, 1),
    log2FoldChange = round(log2FoldChange, 2),
    padj           = signif(padj, 3)
  )
print(top10)