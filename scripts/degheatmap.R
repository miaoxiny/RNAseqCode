library(DESeq2); library(ComplexHeatmap); library(circlize); library(limma); library(dplyr); library(grid)

project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
deg_dir     <- file.path(results_dir, "02_DEG")

# ─── 输入 ───
dds_fitted    <- readRDS(file.path(results_dir, "dds_fitted.rds"))
res_netrin_df <- readRDS(file.path(results_dir, "res_netrin_vs_div2.rds"))

# ─── 66 DEG, 按 padj 排序(复现原图初始行顺序)───
deg_genes <- res_netrin_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)
cat(sprintf("DEG: %d (Up %d / Down %d)\n", nrow(deg_genes),
            sum(deg_genes$log2FoldChange>0), sum(deg_genes$log2FoldChange<0)))

# ─── VST + batch correction (visualization only) ───
vsd     <- vst(dds_fitted, blind = FALSE)
vst_mat <- assay(vsd)
coldata <- as.data.frame(colData(dds_fitted))
mod     <- model.matrix(~ condition, data = coldata)
vst_corrected <- limma::removeBatchEffect(vst_mat, batch = coldata$batch, design = mod)

# ─── A 矩阵 (校正后, z-score) ───
mat_A <- vst_corrected[deg_genes$ensembl_id, ]
rownames(mat_A) <- ifelse(!is.na(deg_genes$symbol), deg_genes$symbol, deg_genes$ensembl_id)
z_A <- t(scale(t(mat_A)))

# ─── 行分块 (up/down) ───
row_split <- factor(ifelse(deg_genes$log2FoldChange > 0, "Upregulated", "Downregulated"),
                    levels = c("Upregulated","Downregulated"))

# ─── 列注释: Replicate 上 / Condition 下 ───
coldata$Condition <- factor(recode(as.character(coldata$condition),
                                   "DIV1"="DIV1","DIV2"="DIV2","DIV2_Netrin"="DIV2 + Netrin-1"),
                            levels = c("DIV1","DIV2","DIV2 + Netrin-1"))
coldata$Replicate <- factor(paste("Replicate", coldata$batch),
                            levels = c("Replicate 1","Replicate 2","Replicate 3"))

col_ha <- HeatmapAnnotation(
  Replicate = coldata$Replicate,
  Condition = coldata$Condition,
  col = list(
    Replicate = c("Replicate 1"="#E8E8E8","Replicate 2"="#888888","Replicate 3"="#2B2B2B"),
    Condition = c("DIV1"="#DCD1BC","DIV2"="#2E5F58","DIV2 + Netrin-1"="#D88040")
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  gap = unit(1, "mm"),
  annotation_legend_param = list(
    Replicate = list(title = "Replicate"),
    Condition = list(title = "Condition")
  )
)

z_col_fun <- colorRamp2(c(-2,0,2), c("#3B7AB8","white","#C73E3E"))

# ─── Heatmap A ───
ht_A <- Heatmap(z_A,
                name = "Z-score",
                col  = z_col_fun,
                top_annotation = col_ha,
                
                # 行: up/down 分块 + 块内 ward.D2 聚类, 不显示行树
                cluster_rows             = TRUE,
                cluster_row_slices       = FALSE,
                clustering_distance_rows = "euclidean",
                clustering_method_rows   = "ward.D2",
                show_row_dend            = FALSE,
                row_split                = row_split,
                row_title                = c("Upregulated","Downregulated"),
                row_title_gp             = gpar(fontsize = 11, fontface = "bold"),
                row_title_rot            = 90,
                row_gap                  = unit(2, "mm"),
                
                # 列: 全局 ward.D2 聚类, 显示列树
                cluster_columns             = TRUE,
                clustering_distance_columns = "euclidean",
                clustering_method_columns   = "ward.D2",
                column_dend_height          = unit(1.5, "cm"),
                
                show_row_names    = TRUE,
                row_names_gp      = gpar(fontsize = 7),
                row_names_side    = "right",
                show_column_names = FALSE,
                
                heatmap_legend_param = list(title = "Z-score", at = c(-2,-1,0,1,2),
                                            legend_height = unit(3, "cm")),
                border = TRUE
)

# ─── 预览 (legend: Replicate→Condition→Z-score) ───
draw(ht_A,
     column_title = "Batch-corrected expression of 66 DEGs (visualization only)",
     column_title_gp = gpar(fontsize = 13, fontface = "bold"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)