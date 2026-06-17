# =============================================================
#  04_PCA.R · FINAL VERSION
#  PCA before & after batch correction
#  - Shape = biological replicate
#  - Color = condition (Palette F)
#  - Legend: Replicate ABOVE Condition
# =============================================================

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(matrixStats)   # for rowVars (loaded via DESeq2 usually)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
qc_dir      <- file.path(results_dir, "01_QC")

# ─── Load saved data ───
vsd               <- readRDS(file.path(results_dir, "vsd.rds"))
vst_mat_corrected <- readRDS(file.path(results_dir, "vst_mat_corrected.rds"))
coldata           <- as.data.frame(colData(vsd))

# ─── Palette F ───
condition_display <- c(
  "DIV1"        = "DIV1",
  "DIV2"        = "DIV2",
  "DIV2_Netrin" = "DIV2 + Netrin-1"
)
condition_pretty <- factor(
  condition_display[as.character(coldata$condition)],
  levels = c("DIV1", "DIV2", "DIV2 + Netrin-1")
)

condition_colors <- c(
  "DIV1"            = "#DCD1BC",   # cream
  "DIV2"            = "#2E5F58",   # dark teal
  "DIV2 + Netrin-1" = "#D88040"    # warm orange
)
replicate_shapes <- c("1" = 16, "2" = 17, "3" = 15)   # ● ▲ ■

# =============================================================
#  PCA before batch correction
# =============================================================
n_top  <- 500
rv_u   <- rowVars(assay(vsd))
sel_u  <- order(rv_u, decreasing = TRUE)[seq_len(min(n_top, length(rv_u)))]
pca_u  <- prcomp(t(assay(vsd)[sel_u, ]), scale. = FALSE)
pct_u  <- round(100 * pca_u$sdev^2 / sum(pca_u$sdev^2), 1)

df_u <- data.frame(
  sample    = rownames(pca_u$x),
  PC1       = pca_u$x[, 1],
  PC2       = pca_u$x[, 2],
  Condition = condition_pretty,
  Replicate = factor(coldata$batch, levels = c("1", "2", "3"))
)

# =============================================================
#  PCA after batch correction
# =============================================================
rv_c   <- rowVars(vst_mat_corrected)
sel_c  <- order(rv_c, decreasing = TRUE)[seq_len(min(n_top, length(rv_c)))]
pca_c  <- prcomp(t(vst_mat_corrected[sel_c, ]), scale. = FALSE)
pct_c  <- round(100 * pca_c$sdev^2 / sum(pca_c$sdev^2), 1)

df_c <- data.frame(
  sample    = rownames(pca_c$x),
  PC1       = pca_c$x[, 1],
  PC2       = pca_c$x[, 2],
  Condition = condition_pretty,
  Replicate = factor(coldata$batch, levels = c("1", "2", "3"))
)

# =============================================================
#  Plot helper — Replicate above Condition in legend
# =============================================================
plot_pca <- function(df, pct_var, title) {
  ggplot(df, aes(x = PC1, y = PC2,
                 color = Condition, shape = Replicate,
                 label = sample)) +
    geom_point(size = 5, stroke = 0.6, alpha = 0.95) +
    geom_text_repel(size = 3.2,
                    max.overlaps  = Inf,
                    box.padding   = 0.6,
                    color         = "grey30",
                    segment.color = "grey60",
                    segment.size  = 0.3) +
    scale_color_manual(values = condition_colors) +
    scale_shape_manual(values = replicate_shapes) +
    guides(
      shape = guide_legend(order = 1, title = "Replicate"),
      color = guide_legend(order = 2, title = "Condition")
    ) +
    labs(title    = title,
         subtitle = "Top 500 variable genes  ·  Shape = biological replicate",
         x = sprintf("PC1: %.1f%% variance", pct_var[1]),
         y = sprintf("PC2: %.1f%% variance", pct_var[2])) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 13, face = "bold"),
          plot.subtitle    = element_text(size = 10, color = "grey40"),
          legend.position  = "right",
          legend.title     = element_text(size = 11, face = "bold"),
          legend.text      = element_text(size = 10))
}

# ─── Build individual plots ───
p_before <- plot_pca(df_u, pct_u, "PCA before batch correction")
p_after  <- plot_pca(df_c, pct_c, "PCA after batch correction")

# ─── Display individually ───
print(p_before)
print(p_after)