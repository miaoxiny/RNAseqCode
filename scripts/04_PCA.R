# =============================================================
#  04_PCA.R — Step 4: PCA before batch correction
#  Inspect-only mode (display only, no auto-save)
# =============================================================

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")

# ─── Load data ───
vsd     <- readRDS(file.path(results_dir, "vsd.rds"))
coldata <- as.data.frame(colData(vsd))

# ─── Project palette (final, locked) ───
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
  "DIV1"             = "#DCD1BC",   # cream
  "DIV2"             = "#2E5F58",   # dark teal
  "DIV2 + Netrin-1"  = "#D88040"    # warm orange
)
replicate_shapes <- c("1" = 16, "2" = 17, "3" = 15)   # circle / triangle / square
replicate_colors <- c("1" = "#9A9A9A", "2" = "#5A5A5A", "3" = "#1A1A1A")

# ─── Compute PCA on top 500 most variable genes ───
n_top      <- 5