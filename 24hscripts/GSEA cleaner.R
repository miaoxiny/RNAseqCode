# =============================================================
#  GSEA 全量表重跑 — pvalueCutoff=1, 拿完整结果表
#  全新文件名 *_FULLTABLE, 绝不覆盖任何旧 *_final 文件
#  其余参数与原 GSEA_Clean.R 完全一致(可比)
# =============================================================
library(dplyr); library(clusterProfiler); library(org.Rn.eg.db); library(AnnotationDbi)

project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
gsea_dir    <- file.path(results_dir, "03_GSEA")

run_gsea_fulltable <- function(res_path, tag) {
  cat("\n========== ", tag, " ==========\n")
  res_df <- readRDS(res_path)
  
  # Step 1 — Wald stat 排序去重(同原脚本)
  ranked_df <- res_df %>%
    filter(!is.na(entrez_id), !is.na(stat)) %>%
    group_by(entrez_id) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  geneList <- sort(setNames(ranked_df$stat, ranked_df$entrez_id), decreasing = TRUE)
  cat("geneList:", length(geneList), "\n")
  
  # Step 2 — gseGO BP, pvalueCutoff=1 (⭐ 唯一关键改动)
  set.seed(123)
  go <- gseGO(geneList=geneList, OrgDb=org.Rn.eg.db, ont="BP", keyType="ENTREZID",
              minGSSize=15, maxGSSize=500, pvalueCutoff=1, pAdjustMethod="BH",
              eps=0, nPermSimple=10000, seed=TRUE, verbose=FALSE)
  go <- setReadable(go, OrgDb=org.Rn.eg.db, keyType="ENTREZID")
  go@result$direction <- ifelse(go@result$NES > 0, "activated", "suppressed")
  go@result <- go@result %>% filter(setSize >= 15)   # post-hoc 同原脚本
  
  # Step 3 — gseKEGG, pvalueCutoff=1
  set.seed(123)
  kg <- gseKEGG(geneList=geneList, organism="rno", keyType="kegg",
                minGSSize=15, maxGSSize=500, pvalueCutoff=1, pAdjustMethod="BH",
                eps=0, nPermSimple=10000, seed=TRUE, verbose=FALSE)
  kg <- setReadable(kg, OrgDb=org.Rn.eg.db, keyType="ENTREZID")
  kg@result$direction <- ifelse(kg@result$NES > 0, "activated", "suppressed")
  kg@result <- kg@result %>% filter(setSize >= 15)
  
  # 自检：全量表行数 + 各层级计数
  cat("--- GO BP 全量表 ---\n")
  cat("  总行数(全量):", nrow(go@result), "\n")
  cat("  p.adjust<0.05 (应=旧的显著数):", sum(go@result$p.adjust<0.05, na.rm=TRUE), "\n")
  cat("  p.adjust<0.10:", sum(go@result$p.adjust<0.10, na.rm=TRUE), "\n")
  cat("  nominal pvalue<0.05:", sum(go@result$pvalue<0.05, na.rm=TRUE), "\n")
  cat("--- KEGG 全量表 ---\n")
  cat("  总行数(全量):", nrow(kg@result), "\n")
  cat("  p.adjust<0.05 (应=旧的显著数):", sum(kg@result$p.adjust<0.05, na.rm=TRUE), "\n")
  cat("  p.adjust<0.10:", sum(kg@result$p.adjust<0.10, na.rm=TRUE), "\n")
  cat("  nominal pvalue<0.05:", sum(kg@result$pvalue<0.05, na.rm=TRUE), "\n")
  
  # Step 4 — 存全新文件名 *_FULLTABLE (绝不覆盖旧 *_final)
  saveRDS(go, file.path(gsea_dir, sprintf("GSEA_GO_BP_%s_FULLTABLE.rds", tag)))
  write.csv(go@result, file.path(gsea_dir, sprintf("GSEA_GO_BP_%s_FULLTABLE.csv", tag)), row.names=FALSE)
  saveRDS(kg, file.path(gsea_dir, sprintf("GSEA_KEGG_%s_FULLTABLE.rds", tag)))
  write.csv(kg@result, file.path(gsea_dir, sprintf("GSEA_KEGG_%s_FULLTABLE.csv", tag)), row.names=FALSE)
  cat("已保存:", sprintf("GSEA_GO_BP_%s_FULLTABLE + GSEA_KEGG_%s_FULLTABLE (.rds/.csv)\n", tag, tag))
}

# ⚠️ 安全检查：先确认不会覆盖旧文件
old_files <- c("GSEA_GO_BP_netrin_final.rds","GSEA_KEGG_netrin_final.rds",
               "GSEA_GO_BP_div2_vs_div1_final.rds","GSEA_KEGG_div2_vs_div1_final.rds")
new_files <- c("GSEA_GO_BP_netrin_FULLTABLE.rds","GSEA_KEGG_netrin_FULLTABLE.rds",
               "GSEA_GO_BP_div2_vs_div1_FULLTABLE.rds","GSEA_KEGG_div2_vs_div1_FULLTABLE.rds")
cat("新文件名是否与已存在文件冲突:\n")
for (f in new_files) cat(" ", f, ":", if(file.exists(file.path(gsea_dir,f))) "⚠️已存在!" else "OK(不冲突)", "\n")

# 跑 C1 和 C2
run_gsea_fulltable(file.path(results_dir, "res_netrin_vs_div2.rds"),   "netrin")
run_gsea_fulltable(file.path(results_dir, "res_div2_vs_div1.rds"),     "div2_vs_div1")