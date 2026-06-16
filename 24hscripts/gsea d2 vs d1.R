# =============================================================
#  GSEA · DIV2 vs DIV1 (Contrast 2) — 算+自检，暂不存盘
#  流程照搬 Contrast 1 (GSEA_Clean.R)；数据源换 res_div2_vs_div1
# =============================================================
library(dplyr); library(clusterProfiler); library(org.Rn.eg.db); library(AnnotationDbi)

project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
gsea_dir    <- file.path(results_dir, "03_GSEA")

res_div_df <- readRDS(file.path(results_dir, "res_div2_vs_div1.rds"))

# ── Step 1 — Ranked gene list (Wald stat, dedup by |stat|) ──
ranked_df <- res_div_df %>%
  filter(!is.na(entrez_id), !is.na(stat)) %>%
  group_by(entrez_id) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()
geneList <- sort(setNames(ranked_df$stat, ranked_df$entrez_id), decreasing = TRUE)

cat("=== Ranked gene list ===\n")
cat(sprintf("  Total ranked: %d genes\n", length(geneList)))
cat(sprintf("  Stat range:   %.2f to %.2f\n", min(geneList), max(geneList)))
cat(sprintf("  Median:       %.4f\n\n", median(geneList)))

# ── Step 2 — gseGO (BP) ──
cat("=== Running gseGO (BP) ... ===\n")
set.seed(123)
gsea_go <- gseGO(
  geneList = geneList, OrgDb = org.Rn.eg.db, ont = "BP", keyType = "ENTREZID",
  minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH",
  eps = 0, nPermSimple = 10000, seed = TRUE, verbose = FALSE
)
gsea_go <- setReadable(gsea_go, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
gsea_go@result$direction <- ifelse(gsea_go@result$NES > 0, "activated", "suppressed")
cat(sprintf("  Before post-hoc: %d terms (%d sig)\n",
            nrow(gsea_go@result), sum(gsea_go@result$p.adjust < 0.05, na.rm = TRUE)))
gsea_go@result <- gsea_go@result %>% filter(setSize >= 15)
cat(sprintf("  After  post-hoc: %d terms (%d sig)\n\n",
            nrow(gsea_go@result), sum(gsea_go@result$p.adjust < 0.05, na.rm = TRUE)))

# ── Step 3 — gseKEGG ──
cat("=== Running gseKEGG ... ===\n")
set.seed(123)
gsea_kegg <- gseKEGG(
  geneList = geneList, organism = "rno", keyType = "kegg",
  minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH",
  eps = 0, nPermSimple = 10000, seed = TRUE, verbose = FALSE
)
gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
gsea_kegg@result$direction <- ifelse(gsea_kegg@result$NES > 0, "activated", "suppressed")
cat(sprintf("  Before post-hoc: %d pathways (%d sig)\n",
            nrow(gsea_kegg@result), sum(gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)))
gsea_kegg@result <- gsea_kegg@result %>% filter(setSize >= 15)
cat(sprintf("  After  post-hoc: %d pathways (%d sig)\n\n",
            nrow(gsea_kegg@result), sum(gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)))

# ── Step 4 — Save （暂不存盘，确认自检合理后再解开）──
# saveRDS(gsea_go,   file.path(gsea_dir, "GSEA_GO_BP_div2_vs_div1_final.rds"))
# write.csv(gsea_go@result,   file.path(gsea_dir, "GSEA_GO_BP_div2_vs_div1_final.csv"), row.names = FALSE)
# saveRDS(gsea_kegg, file.path(gsea_dir, "GSEA_KEGG_div2_vs_div1_final.rds"))
# write.csv(gsea_kegg@result, file.path(gsea_dir, "GSEA_KEGG_div2_vs_div1_final.csv"), row.names = FALSE)

# ── Step 5 — 自检（Contrast 2，无预设预期值，现场看）──
cat("══════ Contrast 2 GSEA 自检 ══════\n")
cat("geneList 长度:", length(geneList), "(应≈12110)\n")
cat("GO BP 显著:", sum(gsea_go@result$p.adjust < 0.05, na.rm=TRUE), "\n")
cat("  激活:", sum(gsea_go@result$p.adjust<0.05 & gsea_go@result$NES>0, na.rm=TRUE),
    "| 抑制:", sum(gsea_go@result$p.adjust<0.05 & gsea_go@result$NES<0, na.rm=TRUE), "\n")
cat("KEGG 显著:", sum(gsea_kegg@result$p.adjust < 0.05, na.rm=TRUE), "\n")
cat("  激活:", sum(gsea_kegg@result$p.adjust<0.05 & gsea_kegg@result$NES>0, na.rm=TRUE),
    "| 抑制:", sum(gsea_kegg@result$p.adjust<0.05 & gsea_kegg@result$NES<0, na.rm=TRUE), "\n")

# 看 top 激活/抑制，确认方向符合发育叙事
cat("\n— GO BP top 5 激活 (NES最高) —\n")
print(gsea_go@result %>% filter(p.adjust<0.05) %>% arrange(desc(NES)) %>% head(5) %>% select(Description, NES, p.adjust))
cat("\n— GO BP top 5 抑制 (NES最低) —\n")
print(gsea_go@result %>% filter(p.adjust<0.05) %>% arrange(NES) %>% head(5) %>% select(Description, NES, p.adjust))
cat("\n— KEGG top 5 激活 —\n")
print(gsea_kegg@result %>% filter(p.adjust<0.05) %>% arrange(desc(NES)) %>% head(5) %>% select(Description, NES, p.adjust))
cat("\n— KEGG top 5 抑制 —\n")
print(gsea_kegg@result %>% filter(p.adjust<0.05) %>% arrange(NES) %>% head(5) %>% select(Description, NES, p.adjust))