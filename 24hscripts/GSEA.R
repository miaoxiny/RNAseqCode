# =============================================================
#  11_GSEA_netrin.R — FINAL CLEAN VERSION
#  GSEA pipeline for Netrin-1 vs DIV2 contrast
#  
#  Input:  12,110 genes (all with valid Wald stat + Entrez ID)
#  Method: Wald statistic ranking, slice_max for duplicates
#  Filter: minGSSize=15 (algorithm) + manual setSize>=15 (post-hoc)
# =============================================================

library(dplyr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
gsea_dir    <- file.path(results_dir, "03_GSEA")
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

# ─── Load DEG results ───
res_netrin_df <- readRDS(file.path(results_dir, "res_netrin_vs_div2.rds"))

# =============================================================
#  Step 1 — Build ranked gene list
# =============================================================
# Ranking: Wald statistic (DESeq2 stat)
# Duplicates: keep gene with largest |stat| per Entrez

ranked_df <- res_netrin_df %>%
  filter(!is.na(entrez_id), !is.na(stat)) %>%
  group_by(entrez_id) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()

geneList <- ranked_df$stat
names(geneList) <- ranked_df$entrez_id
geneList <- sort(geneList, decreasing = TRUE)

cat("=== Ranked gene list ===\n")
cat(sprintf("  Total ranked: %d genes\n", length(geneList)))
cat(sprintf("  Stat range:   %.2f to %.2f\n",
            min(geneList), max(geneList)))
cat(sprintf("  Median:       %.4f\n\n", median(geneList)))

# =============================================================
#  Step 2 — GSEA on GO Biological Process
# =============================================================

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

gsea_go <- setReadable(gsea_go,
                       OrgDb   = org.Rn.eg.db,
                       keyType = "ENTREZID")
gsea_go@result$direction <- ifelse(gsea_go@result$NES > 0,
                                   "activated",
                                   "suppressed")

cat(sprintf("  Before post-hoc filter: %d total terms (%d sig)\n",
            nrow(gsea_go@result),
            sum(gsea_go@result$p.adjust < 0.05, na.rm = TRUE)))

# Post-hoc filter: setSize >= 15
gsea_go@result <- gsea_go@result %>%
  filter(setSize >= 15)

cat(sprintf("  After  post-hoc filter: %d total terms (%d sig)\n\n",
            nrow(gsea_go@result),
            sum(gsea_go@result$p.adjust < 0.05, na.rm = TRUE)))

# =============================================================
#  Step 3 — GSEA on KEGG
# =============================================================

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

gsea_kegg <- setReadable(gsea_kegg,
                         OrgDb   = org.Rn.eg.db,
                         keyType = "ENTREZID")
gsea_kegg@result$direction <- ifelse(gsea_kegg@result$NES > 0,
                                     "activated",
                                     "suppressed")

cat(sprintf("  Before post-hoc filter: %d total pathways (%d sig)\n",
            nrow(gsea_kegg@result),
            sum(gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)))

gsea_kegg@result <- gsea_kegg@result %>%
  filter(setSize >= 15)

cat(sprintf("  After  post-hoc filter: %d total pathways (%d sig)\n\n",
            nrow(gsea_kegg@result),
            sum(gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)))

# =============================================================
#  Step 4 — Save (overwrites previous)
# =============================================================

saveRDS(gsea_go, file.path(gsea_dir, "GSEA_GO_BP_netrin.rds"))
write.csv(gsea_go@result,
          file.path(gsea_dir, "GSEA_GO_BP_netrin.csv"),
          row.names = FALSE)

saveRDS(gsea_kegg, file.path(gsea_dir, "GSEA_KEGG_netrin.rds"))
write.csv(gsea_kegg@result,
          file.path(gsea_dir, "GSEA_KEGG_netrin.csv"),
          row.names = FALSE)

# =============================================================
#  Step 5 — Display top results
# =============================================================

cat("=== TOP 15 ACTIVATED (GO BP) ===\n")
gsea_go@result %>%
  filter(NES > 0, p.adjust < 0.05) %>%
  arrange(p.adjust) %>% head(15) %>%
  select(ID, Description, setSize, NES, p.adjust) %>%
  mutate(NES = round(NES, 2),
         p.adjust = signif(p.adjust, 3)) %>%
  print()

cat("\n=== TOP 15 SUPPRESSED (GO BP) ===\n")
gsea_go@result %>%
  filter(NES < 0, p.adjust < 0.05) %>%
  arrange(p.adjust) %>% head(15) %>%
  select(ID, Description, setSize, NES, p.adjust) %>%
  mutate(NES = round(NES, 2),
         p.adjust = signif(p.adjust, 3)) %>%
  print()

cat("\n=== TOP 10 ACTIVATED (KEGG) ===\n")
gsea_kegg@result %>%
  filter(NES > 0, p.adjust < 0.05) %>%
  arrange(p.adjust) %>% head(10) %>%
  select(ID, Description, setSize, NES, p.adjust) %>%
  mutate(NES = round(NES, 2),
         p.adjust = signif(p.adjust, 3)) %>%
  print()

cat("\n=== TOP 10 SUPPRESSED (KEGG) ===\n")
gsea_kegg@result %>%
  filter(NES < 0, p.adjust < 0.05) %>%
  arrange(p.adjust) %>% head(10) %>%
  select(ID, Description, setSize, NES, p.adjust) %>%
  mutate(NES = round(NES, 2),
         p.adjust = signif(p.adjust, 3)) %>%
  print()

# ─── Verify key thesis pathway ───
cat("\n=== Axon guidance verification ===\n")
gsea_kegg@result %>%
  filter(ID == "rno04360") %>%
  select(ID, Description, setSize, NES, p.adjust) %>%
  mutate(NES = round(NES, 3),
         p.adjust = signif(p.adjust, 3)) %>%
  print()

cat("\n=== Saved to results/03_GSEA/ ===\n")