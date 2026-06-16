library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)

project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
deg_dir     <- file.path(results_dir, "02_DEG")

res_div_df <- readRDS(file.path(results_dir, "res_div2_vs_div1.rds"))
vst_c      <- readRDS(file.path(results_dir, "vst_mat_corrected.rds"))

strict <- res_div_df %>% filter(!is.na(padj), padj<0.05, abs(log2FoldChange)>0.585) %>%
  mutate(label = ifelse(is.na(symbol) | symbol=="", ensembl_id, symbol))

col_fun <- colorRamp2(c(-2, 0, 2), c("#3B7AB8", "white", "#C73E3E"))

build_col_ha <- function(samp) {
  cond <- recode(sub("^[0-9]+_", "", samp), "DIV2+Netrin-1" = "DIV2 + Netrin-1")
  repl <- paste0("Replicate ", sub("_.*", "", samp))
  HeatmapAnnotation(
    Replicate = repl,
    Condition = cond,
    col = list(
      Condition = c("DIV1"="#DCD1BC","DIV2"="#2E5F58","DIV2 + Netrin-1"="#D88040"),
      Replicate = c("Replicate 1"="#E8E8E8","Replicate 2"="#888888","Replicate 3"="#2B2B2B")),
    annotation_name_side = "left", gap = unit(1,"mm"),
    annotation_legend_param = list(Condition=list(title="Condition"),
                                   Replicate=list(title="Replicate")))
}

cluster_order <- function(ids, mat_z) {
  m <- mat_z[ids,,drop=FALSE]; if(nrow(m)<3) return(ids)
  ids[hclust(dist(m), method="complete")$order]
}

make_ht <- function(ids_up, ids_down, title, show_rows, fontsize=7) {
  all_ids <- c(ids_up, ids_down)
  lab_map <- setNames(strict$label, strict$ensembl_id)
  mat_z <- t(scale(t(vst_c[all_ids, ])))
  up_ord   <- cluster_order(ids_up,   mat_z)
  down_ord <- cluster_order(ids_down, mat_z)
  gene_order <- c(up_ord, down_ord)
  mat_z <- mat_z[gene_order, ]
  if (show_rows) rownames(mat_z) <- lab_map[gene_order]
  row_split <- factor(c(rep("Upregulated", length(up_ord)),
                        rep("Downregulated", length(down_ord))),
                      levels=c("Upregulated","Downregulated"))
  cond_rank <- c("DIV1"=1, "DIV2"=2, "DIV2+Netrin-1"=3)
  col_w <- as.numeric(cond_rank[ sub("^[0-9]+_", "", colnames(mat_z)) ])
  
  Heatmap(mat_z, name="Z-score", col=col_fun,
          top_annotation = build_col_ha(colnames(mat_z)),
          row_split = row_split, row_title_rot = 90, row_gap = unit(2.5,"mm"),
          cluster_rows = FALSE,
          cluster_columns = TRUE,
          column_dend_reorder = col_w,
          show_column_names = FALSE,
          show_row_names = show_rows, row_names_side = "right",
          row_names_gp = gpar(fontsize = fontsize),
          column_title = title, column_title_gp = gpar(fontsize=13, fontface="bold"),
          heatmap_legend_param = list(title="Z-score", at=c(-2,-1,0,1,2),
                                      legend_height=unit(3,"cm")))
}

# ═══ 图1: top60 ═══
top_up   <- strict %>% filter(log2FoldChange>0) %>% arrange(padj) %>% head(30)
top_down <- strict %>% filter(log2FoldChange<0) %>% arrange(padj) %>% head(30)
ht_top60 <- make_ht(top_up$ensembl_id, top_down$ensembl_id,
                    "Top 60 DEGs (DIV2 vs DIV1, padj top30 each direction)\nbatch-corrected, visualization only",
                    show_rows = TRUE, fontsize = 7)

png(file.path(deg_dir, "07_heatmap_div2_vs_div1_top60.png"), width=10, height=11, units="in", res=300)
draw(ht_top60, heatmap_legend_side="right", annotation_legend_side="right", merge_legend=TRUE); dev.off()
pdf(file.path(deg_dir, "07_heatmap_div2_vs_div1_top60.pdf"), width=10, height=11)
draw(ht_top60, heatmap_legend_side="right", annotation_legend_side="right", merge_legend=TRUE); dev.off()
cat("已保存: 07_heatmap_div2_vs_div1_top60.png/pdf\n")

# ═══ 图2: all1495 ═══
up_all   <- strict %>% filter(log2FoldChange>0) %>% pull(ensembl_id)
down_all <- strict %>% filter(log2FoldChange<0) %>% pull(ensembl_id)
ht_all <- make_ht(up_all, down_all,
                  "All 1495 strict DEGs (DIV2 vs DIV1)\nbatch-corrected, visualization only",
                  show_rows = FALSE)

png(file.path(deg_dir, "08_heatmap_div2_vs_div1_all1495.png"), width=10, height=11, units="in", res=300)
draw(ht_all, heatmap_legend_side="right", annotation_legend_side="right", merge_legend=TRUE); dev.off()
pdf(file.path(deg_dir, "08_heatmap_div2_vs_div1_all1495.pdf"), width=10, height=11)
draw(ht_all, heatmap_legend_side="right", annotation_legend_side="right", merge_legend=TRUE); dev.off()
cat("已保存: 08_heatmap_div2_vs_div1_all1495.png/pdf\n")