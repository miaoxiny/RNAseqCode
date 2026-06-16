# =============================================================
#  03_QC.R
#  Step 3: Exploratory QC - Filtering, DESeq2 object,
#          sample correlation
# =============================================================

library(DESeq2)
library(ComplexHeatmap)
library(circlize)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
qc_dir      <- file.path(results_dir, "01_QC")

# ─── Load imported data ───
txi     <- readRDS(file.path(results_dir, "txi.rds"))
coldata <- readRDS(file.path(results_dir, "coldata.rds"))

cat("=== Loaded data ===\n")
cat("Genes:  ", nrow(txi$counts), "\n")
cat("Samples:", ncol(txi$counts), "\n\n")

# =============================================================
#  3.1 Build DESeq2 object
# =============================================================

cat("=== Building DESeq2 object ===\n")
dds <- DESeqDataSetFromTximport(
  txi     = txi,
  colData = coldata,
  design  = ~ batch + condition    # account for replicate batch
)

cat("Initial dds dimensions:", dim(dds), "\n\n")

# =============================================================
#  3.2 Filter low-expressed genes
# =============================================================
# Rule: keep genes with >= 10 counts in at least 3 samples
# (3 = the smallest group size in our design)

keep_idx <- rowSums(counts(dds) >= 10) >= 3
cat("=== Gene filtering ===\n")
cat("Before filtering:", nrow(dds), "\n")
dds <- dds[keep_idx, ]
cat("After filtering: ", nrow(dds), "\n")
cat("Removed:         ", sum(!keep_idx),
    sprintf("(%.1f%%)\n\n", 100 * sum(!keep_idx) / length(keep_idx)))

# =============================================================
#  3.3 VST transformation for visualization
# =============================================================
# blind = FALSE so transformation uses the design

cat("=== Running VST (variance-stabilizing transformation) ===\n")
vsd <- vst(dds, blind = FALSE)
cat("VST done.\n\n")

# =============================================================
#  3.4 Sample-to-sample correlation heatmap
# =============================================================

# Project palette
condition_colors <- c(
  "DIV1"        = "#9CA5A6",
  "DIV2"        = "#3D6B8C",
  "DIV2_Netrin" = "#D97B43"
)
batch_colors <- c("1" = "#E8E8E8", "2" = "#B8B8B8", "3" = "#888888")

# Compute sample-to-sample correlation on VST data
cor_matrix <- cor(assay(vsd), method = "pearson")
cat("=== Pearson correlation summary ===\n")
cat("Min:   ", round(min(cor_matrix[upper.tri(cor_matrix)]), 4), "\n")
cat("Mean:  ", round(mean(cor_matrix[upper.tri(cor_matrix)]), 4), "\n")
cat("Max:   ", round(max(cor_matrix[upper.tri(cor_matrix)]), 4), "\n\n")

# Annotation strip (top)
ha <- HeatmapAnnotation(
  Condition = coldata$condition,
  Batch     = coldata$batch,
  col       = list(Condition = condition_colors,
                   Batch     = batch_colors),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 8),
  simple_anno_size     = unit(3, "mm"),
  border               = TRUE
)

# Color scale: white -> dark navy for correlation values
col_fun <- colorRamp2(
  c(min(cor_matrix), (min(cor_matrix) + 1) / 2, 1),
  c("#F7F7F7", "#92B4C8", "#1F3A5F")
)

# Plot heatmap
hm <- Heatmap(
  cor_matrix,
  name              = "Pearson r",
  col               = col_fun,
  top_annotation    = ha,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  row_names_gp      = gpar(fontsize = 9),
  column_names_gp   = gpar(fontsize = 9),
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  rect_gp           = gpar(col = "white", lwd = 0.4),
  border            = TRUE,
  column_title      = "Sample-to-sample correlation (VST, Pearson)",
  column_title_gp   = gpar(fontsize = 11, fontface = "plain"),
  heatmap_legend_param = list(
    title_gp     = gpar(fontsize = 9),
    labels_gp    = gpar(fontsize = 8),
    legend_height = unit(3, "cm")
  )
)

# Save to PDF
pdf(file.path(qc_dir, "01_sample_correlation_heatmap.pdf"),
    width = 7, height = 6)
draw(hm)
dev.off()

cat("Saved: 01_sample_correlation_heatmap.pdf\n\n")

# =============================================================
#  3.5 Save filtered DESeq2 object for next step
# =============================================================

saveRDS(dds, file.path(results_dir, "dds_filtered.rds"))
saveRDS(vsd, file.path(results_dir, "vsd.rds"))
cat("Saved: dds_filtered.rds, vsd.rds\n")