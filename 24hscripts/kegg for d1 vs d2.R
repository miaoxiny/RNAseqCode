library(clusterProfiler)
library(dplyr)

project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
ora_dir     <- file.path(results_dir, "04_ORA_div2_vs_div1")

res_div_df <- readRDS(file.path(results_dir, "res_div2_vs_div1.rds"))

universe_entrez <- res_div_df %>%
  filter(!is.na(entrez_id), !is.na(padj)) %>%
  pull(entrez_id) %>% unique() %>% as.character()

genes_up_strict <- res_div_df %>%
  filter(!is.na(padj), padj < 0.05, log2FoldChange >  0.585, !is.na(entrez_id)) %>%
  pull(entrez_id) %>% unique() %>% as.character()

genes_down_strict <- res_div_df %>%
  filter(!is.na(padj), padj < 0.05, log2FoldChange < -0.585, !is.na(entrez_id)) %>%
  pull(entrez_id) %>% unique() %>% as.character()

cat("KEGG ORA input — Up:", length(genes_up_strict),
    " Down:", length(genes_down_strict),
    " Universe:", length(universe_entrez), "\n\n")

run_KEGG_ORA <- function(gene_set, universe, name) {
  cat(sprintf("=== KEGG ORA · %s ===\n", name))
  if (length(gene_set) < 5) { cat("  Too few genes\n\n"); return(NULL) }
  ek <- enrichKEGG(
    gene          = gene_set,
    universe      = universe,
    organism      = "rno",
    keyType       = "ncbi-geneid",   # NCBI Entrez（规范写法）
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,            # 对齐 GO
    qvalueCutoff  = 0.20,            # 对齐 GO
    minGSSize     = 10,             # 对齐 GO
    maxGSSize     = 500             # 对齐 GO
  )
  if (is.null(ek) || nrow(ek@result) == 0) { cat("  No results\n\n"); return(NULL) }
  cat(sprintf("  significant (p.adjust<0.05): %d / total in @result: %d\n",
              sum(ek@result$p.adjust < 0.05, na.rm = TRUE), nrow(ek@result)))
  cat("  实际 universe:", length(ek@universe), " | 实际 gene 输入:", length(ek@gene), "\n\n")
  ek
}

kegg_up   <- run_KEGG_ORA(genes_up_strict,   universe_entrez, "Up STRICT")
kegg_down <- run_KEGG_ORA(genes_down_strict, universe_entrez, "Down STRICT")

cat("=== Up top 10 ===\n")
if(!is.null(kegg_up))   print(head(kegg_up@result[,c("ID","Description","p.adjust","Count")], 10))
cat("\n=== Down top 10 ===\n")
if(!is.null(kegg_down)) print(head(kegg_down@result[,c("ID","Description","p.adjust","Count")], 10))

# ── 存盘（KEGG 不分 raw/simplified，直接存对象 + 显著子集 csv）──
if(!is.null(kegg_up)) {
  saveRDS(kegg_up, file.path(ora_dir, "ORA_KEGG_div_up.rds"))
  write.csv(kegg_up@result, file.path(ora_dir, "ORA_KEGG_div_up.csv"), row.names = FALSE)
}
if(!is.null(kegg_down)) {
  saveRDS(kegg_down, file.path(ora_dir, "ORA_KEGG_div_down.rds"))
  write.csv(kegg_down@result, file.path(ora_dir, "ORA_KEGG_div_down.csv"), row.names = FALSE)
}
cat("\n=== Saved ===\n")
cat("  ORA_KEGG_div_up.rds + .csv\n  ORA_KEGG_div_down.rds + .csv\n")