## ============================================================================
## acute Netrin-1 RNA-seq · GSEA 重跑(统一规范版 · 最终)
## 与 08/09 同口径:Wald stat 排序、entrez 去重(|stat|最大)、
##   minGSSize=15 maxGSSize=500 BH eps=0 nPermSimple=10000 seed=123、setSize>=15
## 改进:pvalueCutoff=1 完整导出 → 外部分层;新增 GO CC;KEGG 用 ncbi-geneid。
## 方向标签统一 activated/suppressed。新文件名 _v2full,防覆盖,先 print 后存。
## ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(clusterProfiler); library(org.Rn.eg.db); library(AnnotationDbi)
})
out_base <- "D:/Dropbox/Dropbox/RNAseq 2025/Time analysis"
obj_dir  <- file.path(out_base, "results/objects")
gsea_dir <- file.path(out_base, "results/04_GSEA")
dir.create(gsea_dir, recursive=TRUE, showWarnings=FALSE)

## ---- 注释 + 结果对象(同 08) ----
anno <- read.csv(file.path(obj_dir, "gene_annotation.csv"), stringsAsFactors=FALSE)
rownames(anno) <- anno$ensembl
res_all <- readRDS(file.path(obj_dir, "res_all.rds"))

## ---- 构造 preranked geneList(逐行同 08) ----
make_ranks <- function(res){
  df <- as.data.frame(res); df$ensembl <- rownames(df)
  df$entrez_id <- anno[df$ensembl, "entrez"]
  ranked <- df %>% dplyr::filter(!is.na(entrez_id), !is.na(stat)) %>%
    dplyr::group_by(entrez_id) %>%
    dplyr::slice_max(order_by=abs(stat), n=1, with_ties=FALSE) %>% dplyr::ungroup()
  sort(setNames(ranked$stat, ranked$entrez_id), decreasing=TRUE)
}

## ---- GO(BP/CC):pvalueCutoff=1 完整导出 ----
run_go <- function(geneList, ont){
  set.seed(123)
  g <- gseGO(geneList=geneList, OrgDb=org.Rn.eg.db, ont=ont, keyType="ENTREZID",
             minGSSize=15, maxGSSize=500, pvalueCutoff=1, pAdjustMethod="BH",
             eps=0, nPermSimple=10000, seed=TRUE, verbose=FALSE)
  g <- setReadable(g, OrgDb=org.Rn.eg.db, keyType="ENTREZID")
  g@result$direction <- ifelse(is.na(g@result$NES), NA,
                               ifelse(g@result$NES > 0, "activated", "suppressed"))
  g@result <- g@result %>% dplyr::filter(setSize >= 15)
  g
}

## ---- KEGG:ncbi-geneid + 不做 setReadable ----
run_kegg <- function(geneList){
  set.seed(123)
  g <- gseKEGG(geneList=geneList, organism="rno", keyType="ncbi-geneid",
               minGSSize=15, maxGSSize=500, pvalueCutoff=1, pAdjustMethod="BH",
               eps=0, nPermSimple=10000, seed=TRUE, verbose=FALSE)
  g@result$direction <- ifelse(is.na(g@result$NES), NA,
                               ifelse(g@result$NES > 0, "activated", "suppressed"))
  g@result <- g@result %>% dplyr::filter(setSize >= 15)
  g
}

## ---- 分层报告(na.rm 防护) ----
report <- function(g, tag){
  r <- g@result
  cat(sprintf("[%-14s] tested %4d | strict(FDR<0.05) %3d | suggestive(FDR<0.10) %3d | nominal(p<0.05) %3d\n",
              tag, nrow(r),
              sum(r$p.adjust < 0.05, na.rm = TRUE),
              sum(r$p.adjust < 0.10, na.rm = TRUE),
              sum(r$pvalue   < 0.05, na.rm = TRUE)))
}

## ---- strict 子集的 simplify(仅 GO,去冗余,供展示) ----
simp_strict <- function(g){
  gs <- g; gs@result <- g@result[which(g@result$p.adjust < 0.05), ]
  if(nrow(gs@result) == 0) return(gs)
  simplify(gs, cutoff=0.7, by="p.adjust", select_fun=min)
}

## ===================== 三对比 × 三本体 =====================
contrasts <- list("5v0"=res_all$res_5v0, "15v0"=res_all$res_15v0, "15v5"=res_all$res_15v5)
GSEA <- list()
for(cn in names(contrasts)){
  gl <- make_ranks(contrasts[[cn]])
  cat(sprintf("\n=== %s | ranked genes: %d ===\n", cn, length(gl)))
  GSEA[[cn]] <- list(BP=run_go(gl,"BP"), CC=run_go(gl,"CC"), KEGG=run_kegg(gl))
  report(GSEA[[cn]]$BP,   paste(cn,"GO BP"))
  report(GSEA[[cn]]$CC,   paste(cn,"GO CC"))
  report(GSEA[[cn]]$KEGG, paste(cn,"KEGG"))
}

## ---- 新旧对照(strict 应接近旧值) ----
cat("\n--- 新旧对照(strict, FDR<0.05) ---\n")
cat("5v0  BP   =", sum(GSEA[["5v0"]]$BP@result$p.adjust<0.05,  na.rm=TRUE), " (旧 313)\n")
cat("5v0  KEGG =", sum(GSEA[["5v0"]]$KEGG@result$p.adjust<0.05,na.rm=TRUE), " (旧 14, 注:keyType改ncbi可能变)\n")
cat("5v0  CC   =", sum(GSEA[["5v0"]]$CC@result$p.adjust<0.05,  na.rm=TRUE), " (新增)\n")
cat("15v5 BP   =", sum(GSEA[["15v5"]]$BP@result$p.adjust<0.05, na.rm=TRUE), " (旧 295)\n")
cat("15v0 BP   =", sum(GSEA[["15v0"]]$BP@result$p.adjust<0.05, na.rm=TRUE), " (旧 0)\n")

## ---- 预览:5v0 CC 上/下调端(看是否佐证 BP) ----
cc5 <- GSEA[["5v0"]]$CC@result
cat("\n--- 5v0 GO CC 上调端 top12 ---\n")
print(head(cc5[order(-cc5$NES), c("Description","NES","setSize","p.adjust")], 12))
cat("\n--- 5v0 GO CC 下调端 top10 ---\n")
print(head(cc5[order( cc5$NES), c("Description","NES","setSize","p.adjust")], 10))

## ============================================================================
## 存盘:全部新文件名 _v2full,防覆盖,不动任何原 08/09 输出。预览满意后解注。
## ============================================================================
# save_one <- function(g, fname){
#   f <- file.path(gsea_dir, fname); stopifnot(!file.exists(f))
#   write.csv(g@result, f, row.names=FALSE)
# }
# for(cn in names(GSEA)){
#   save_one(GSEA[[cn]]$BP,   sprintf("GSEA_GO_BP_%s_v2full.csv", cn))
#   save_one(GSEA[[cn]]$CC,   sprintf("GSEA_GO_CC_%s_v2full.csv", cn))
#   save_one(GSEA[[cn]]$KEGG, sprintf("GSEA_KEGG_%s_v2full.csv",  cn))
#   save_one(simp_strict(GSEA[[cn]]$BP), sprintf("GSEA_GO_BP_%s_strict_simplified_v2.csv", cn))
#   save_one(simp_strict(GSEA[[cn]]$CC), sprintf("GSEA_GO_CC_%s_strict_simplified_v2.csv", cn))
# }
# saveRDS(GSEA, file.path(gsea_dir, "GSEA_all_v2full.rds"))
# stopifnot(!file.exists(file.path(gsea_dir,"GSEA_all_v2full.rds")))  # (上一行已写,如需严格防覆盖把saveRDS也包stopifnot)
# cat("[saved] GSEA *_v2full.* (原文件未动)\n")