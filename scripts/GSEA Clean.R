# =============================================================
#  GSEA_final.R — Netrin-1 vs DIV2 · FINAL
#  Wald stat ranking · nPermSimple=10000 · setSize>=15 post-hoc
#  输出: results/03_GSEA/GSEA_*_netrin_final.rds/csv
# =============================================================
library(dplyr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(AnnotationDbi)

project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
gsea_dir    <- file.path(results_dir, "03_GSEA")
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

res_netrin_df <- readRDS(file.path(results_dir, "res_netrin_vs_div2.rds"))

# ── Step 1 — Ranked gene list (Wald stat, dedup by |stat|) ──
ranked_df <- res_netrin_df %>%
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
  geneList      = geneList,
  OrgDb         = org.Rn.eg.db,
  ont           = "BP",
  keyType       = "ENTREZID",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 0,
  nPermSimple   = 10000,
  seed          = TRUE,
  verbose       = FALSE
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
  geneList      = geneList,
  organism      = "rno",
  keyType       = "kegg",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 0,
  nPermSimple   = 10000,
  seed          = TRUE,
  verbose       = FALSE
)
gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
gsea_kegg@result$direction <- ifelse(gsea_kegg@result$NES > 0, "activated", "suppressed")
cat(sprintf("  Before post-hoc: %d pathways (%d sig)\n",
            nrow(gsea_kegg@result), sum(gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)))
gsea_kegg@result <- gsea_kegg@result %>% filter(setSize >= 15)
cat(sprintf("  After  post-hoc: %d pathways (%d sig)\n\n",
            nrow(gsea_kegg@result), sum(gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)))

# ── Step 4 — Save as FINAL ──
saveRDS(gsea_go,   file.path(gsea_dir, "GSEA_GO_BP_netrin_final.rds"))
write.csv(gsea_go@result,   file.path(gsea_dir, "GSEA_GO_BP_netrin_final.csv"), row.names = FALSE)
saveRDS(gsea_kegg, file.path(gsea_dir, "GSEA_KEGG_netrin_final.rds"))
write.csv(gsea_kegg@result, file.path(gsea_dir, "GSEA_KEGG_netrin_final.csv"), row.names = FALSE)
cat("已保存 FINAL: GSEA_GO_BP_netrin_final.rds / GSEA_KEGG_netrin_final.rds\n\n")

# ── Step 5 — 自检 ──
cat("══════ FINAL 自检 ══════\n")
cat("GO BP 显著:", sum(gsea_go@result$p.adjust < 0.05, na.rm=TRUE), "(应 481)\n")
cat("GO BP NA padj:", sum(is.na(gsea_go@result$p.adjust)), "\n")
cat("  激活:", sum(gsea_go@result$p.adjust<0.05 & gsea_go@result$NES>0, na.rm=TRUE),
    "| 抑制:", sum(gsea_go@result$p.adjust<0.05 & gsea_go@result$NES<0, na.rm=TRUE), "\n")
cat("KEGG 显著:", sum(gsea_kegg@result$p.adjust < 0.05, na.rm=TRUE), "(应 42)\n")
cat("KEGG NA padj:", sum(is.na(gsea_kegg@result$p.adjust)), "\n")
ag <- gsea_kegg@result[gsea_kegg@result$ID=="rno04360", ]
cat("axon guidance (KEGG): NES", round(ag$NES[1],3), "padj", signif(ag$p.adjust[1],3), "\n")