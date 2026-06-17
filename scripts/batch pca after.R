# step3_batch_correction.R

# 如未安装：
# install.packages(c("readr", "dplyr", "tibble", "ggplot2", "pheatmap"))
# install.packages("BiocManager")
# BiocManager::install("limma")

library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(limma)

root <- "D:/Dropbox/Dropbox/RNAseq 2025/Codex acute"

step2_dir <- file.path(root, "RNAseq_stepwise_tximport", "step2_qc")
out_dir <- file.path(root, "RNAseq_stepwise_tximport", "step3_batch_correction")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Read Step 2 outputs
# -----------------------------

sample_metadata <- read_csv(
  file.path(root, "RNAseq_stepwise_tximport", "step1_count_matrix", "sample_metadata.csv"),
  show_col_types = FALSE
)

sample_metadata <- as.data.frame(sample_metadata)

sample_metadata$condition <- factor(
  sample_metadata$condition,
  levels = c("0min", "5min", "15min")
)

sample_metadata$replicate <- factor(sample_metadata$replicate)

rownames(sample_metadata) <- sample_metadata$sample_id

vst_df <- read_csv(
  file.path(step2_dir, "step2_vst_matrix_filtered_for_qc.csv"),
  show_col_types = FALSE
)

vst_mat <- vst_df %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

vst_mat <- vst_mat[, sample_metadata$sample_id]

# -----------------------------
# 2. Batch correction
# -----------------------------
# Batch = biological replicate / culture prep
# Preserve condition/time-point biological effect

design_condition <- model.matrix(~ condition, data = sample_metadata)

vst_batch_corrected <- removeBatchEffect(
  vst_mat,
  batch = sample_metadata$replicate,
  design = design_condition
)

write_csv(
  data.frame(
    gene_id = rownames(vst_batch_corrected),
    vst_batch_corrected,
    check.names = FALSE
  ),
  file.path(out_dir, "step3_vst_batch_corrected_replicate_for_visualization.csv")
)

# -----------------------------
# 3. Correlation before and after correction
# -----------------------------

cor_before <- cor(vst_mat, method = "pearson")
cor_after <- cor(vst_batch_corrected, method = "pearson")

write_csv(
  data.frame(
    sample_id = rownames(cor_before),
    cor_before,
    check.names = FALSE
  ),
  file.path(out_dir, "step3_sample_correlation_before_batch_correction.csv")
)

write_csv(
  data.frame(
    sample_id = rownames(cor_after),
    cor_after,
    check.names = FALSE
  ),
  file.path(out_dir, "step3_sample_correlation_after_batch_correction.csv")
)

annotation_col <- data.frame(
  condition = sample_metadata$condition,
  replicate = sample_metadata$replicate
)

rownames(annotation_col) <- sample_metadata$sample_id

cor_min_before <- min(cor_before[upper.tri(cor_before)])
cor_min_after <- min(cor_after[upper.tri(cor_after)])

cor_floor <- max(
  0.8,
  floor(min(cor_min_before, cor_min_after) * 100) / 100
)

pdf(file.path(out_dir, "step3_correlation_before_batch_correction.pdf"), width = 6.5, height = 6)
pheatmap(
  cor_before,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  color = colorRampPalette(c("white", "#4C78A8"))(100),
  breaks = seq(cor_floor, 1, length.out = 101),
  main = "Before batch correction, VST"
)
dev.off()

png(file.path(out_dir, "step3_correlation_before_batch_correction.png"), width = 1800, height = 1800, res = 300)
pheatmap(
  cor_before,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  color = colorRampPalette(c("white", "#4C78A8"))(100),
  breaks = seq(cor_floor, 1, length.out = 101),
  main = "Before batch correction, VST"
)
dev.off()

pdf(file.path(out_dir, "step3_correlation_after_batch_correction.pdf"), width = 6.5, height = 6)
pheatmap(
  cor_after,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  color = colorRampPalette(c("white", "#54A24B"))(100),
  breaks = seq(cor_floor, 1, length.out = 101),
  main = "After replicate correction, VST"
)
dev.off()

png(file.path(out_dir, "step3_correlation_after_batch_correction.png"), width = 1800, height = 1800, res = 300)
pheatmap(
  cor_after,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  color = colorRampPalette(c("white", "#54A24B"))(100),
  breaks = seq(cor_floor, 1, length.out = 101),
  main = "After replicate correction, VST"
)
dev.off()

# -----------------------------
# 4. PCA function
# -----------------------------

# -----------------------------
# 4. PCA using the same top 500 genes before and after correction
# -----------------------------

top_n <- 500

# Select top 500 variable genes from BEFORE-correction VST matrix
gene_var_before <- apply(vst_mat, 1, var)

top_genes_fixed <- names(sort(gene_var_before, decreasing = TRUE))[1:min(top_n, length(gene_var_before))]

write_csv(
  data.frame(
    gene_id = top_genes_fixed,
    variance_before = gene_var_before[top_genes_fixed],
    rank = seq_along(top_genes_fixed)
  ),
  file.path(out_dir, "step3_fixed_top500_variable_genes_selected_before_correction.csv")
)

run_pca_fixed_genes <- function(mat, sample_metadata, genes) {
  pca_input <- mat[genes, , drop = FALSE]
  
  pca <- prcomp(t(pca_input), scale. = FALSE)
  
  percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)
  
  pca_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    condition = sample_metadata[rownames(pca$x), "condition", drop = TRUE],
    time_min = sample_metadata[rownames(pca$x), "time_min", drop = TRUE],
    replicate = sample_metadata[rownames(pca$x), "replicate", drop = TRUE],
    stringsAsFactors = FALSE
  )
  
  list(
    pca_df = pca_df,
    percent_var = percent_var
  )
}

pca_before <- run_pca_fixed_genes(
  mat = vst_mat,
  sample_metadata = sample_metadata,
  genes = top_genes_fixed
)

pca_after <- run_pca_fixed_genes(
  mat = vst_batch_corrected,
  sample_metadata = sample_metadata,
  genes = top_genes_fixed
)


# -----------------------------
# 5. PCA plots before / after
# -----------------------------

condition_colors <- c(
  "0min" = "#4C78A8",
  "5min" = "#F58518",
  "15min" = "#54A24B"
)

p_before <- ggplot(
  pca_before$pca_df,
  aes(x = PC1, y = PC2, color = condition, shape = replicate)
) +
  geom_point(size = 4) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "PCA before batch correction",
    subtitle = "Top 500 variable genes, VST",
    x = paste0("PC1: ", pca_before$percent_var[1], "% variance"),
    y = paste0("PC2: ", pca_before$percent_var[2], "% variance")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank()
  )

p_after <- ggplot(
  pca_after$pca_df,
  aes(x = PC1, y = PC2, color = condition, shape = replicate)
) +
  geom_point(size = 4) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "PCA after replicate correction",
    subtitle = "Top 500 variable genes, batch-corrected VST",
    x = paste0("PC1: ", pca_after$percent_var[1], "% variance"),
    y = paste0("PC2: ", pca_after$percent_var[2], "% variance")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank()
  )

ggsave(
  file.path(out_dir, "step3_pca_before_batch_correction_top500.png"),
  p_before,
  width = 5.5,
  height = 4.5,
  dpi = 300
)

ggsave(
  file.path(out_dir, "step3_pca_before_batch_correction_top500.pdf"),
  p_before,
  width = 5.5,
  height = 4.5
)

ggsave(
  file.path(out_dir, "step3_pca_after_batch_correction_top500.png"),
  p_after,
  width = 5.5,
  height = 4.5,
  dpi = 300
)

ggsave(
  file.path(out_dir, "step3_pca_after_batch_correction_top500.pdf"),
  p_after,
  width = 5.5,
  height = 4.5
)

# -----------------------------
# 6. Combined PCA comparison plot
# -----------------------------

pca_before_plot <- pca_before$pca_df
pca_before_plot$status <- "Before correction"

pca_after_plot <- pca_after$pca_df
pca_after_plot$status <- "After replicate correction"

pca_compare <- bind_rows(pca_before_plot, pca_after_plot)

p_compare <- ggplot(
  pca_compare,
  aes(x = PC1, y = PC2, color = condition, shape = replicate)
) +
  geom_point(size = 3.5) +
  facet_wrap(~ status, scales = "free") +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "PCA before and after replicate correction",
    subtitle = "Top 500 variable genes",
    x = "PC1",
    y = "PC2"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank()
  )

ggsave(
  file.path(out_dir, "step3_pca_before_after_batch_correction_top500.png"),
  p_compare,
  width = 9,
  height = 4.5,
  dpi = 300
)

ggsave(
  file.path(out_dir, "step3_pca_before_after_batch_correction_top500.pdf"),
  p_compare,
  width = 9,
  height = 4.5
)

# -----------------------------
# 7. Summary table
# -----------------------------

batch_correction_summary <- data.frame(
  metric = c(
    "genes_in_vst_matrix",
    "top_variable_genes_for_pca",
    "batch_variable",
    "biological_variable_preserved",
    "min_sample_correlation_before",
    "min_sample_correlation_after",
    "pc1_percent_before",
    "pc2_percent_before",
    "pc1_percent_after",
    "pc2_percent_after"
  ),
  value = c(
    nrow(vst_mat),
    500,
    "replicate",
    "condition",
    cor_min_before,
    cor_min_after,
    pca_before$percent_var[1],
    pca_before$percent_var[2],
    pca_after$percent_var[1],
    pca_after$percent_var[2]
  )
)

write_csv(
  batch_correction_summary,
  file.path(out_dir, "step3_batch_correction_summary.csv")
)

cat("\nStep 3 batch correction complete.\n")
cat("Output folder:\n")
cat(out_dir, "\n\n")

cat("Batch correction summary:\n")
print(batch_correction_summary)

cat("\nImportant note:\n")
cat("Batch-corrected VST matrix is for visualization/clustering only.\n")
cat("Do not use batch-corrected expression values for DESeq2 differential expression.\n")
cat("For DESeq2, use design = ~ replicate + condition.\n")