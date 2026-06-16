# =============================================================
#  05_batch_correction.R — Step 5: Batch correction
#  (for visualization only — does NOT affect DEG analysis)
#  Inspect-only mode
# =============================================================

library(DESeq2)
library(limma)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")

# ─── Load data ───
vsd     <- readRDS(file.path(results_dir, "vsd.rds"))
coldata <- as.data.frame(colData(vsd))

# ─── Project palette ───
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
replicate_shapes <- c("1" = 16, "2" = 17, "3" = 15)
replicate_colors <- c("1" = "#9A9A9A", "2" = "#5A5A5A", "3" = "#1A1A1A")

# =============================================================
#  Step 5.1: Apply batch correction (for visualization only)
# =============================================================

# Extract VST matrix and metadata
vst_mat <- assay(vsd)
batch   <- coldata$batch

# Build design matrix that PRESERVES condition information
# (this tells removeBatchEffect: remove batch variation, keep condition)
design_to_keep <- model.matrix(~ condition_pretty)

# Remove batch effect for visualization
vst_mat_corrected <- limma::removeBatchEffect(
  vst_mat,
  batch  = batch,
  design = design_to_keep
)

cat("=== Batch correction summary ===\n")
cat("Original VST matrix dimensions:  ", dim(vst_mat), "\n")
cat("Corrected VST matrix dimensions: ", dim(vst_mat_corrected), "\n")
cat("Mean correlation between before/after (should be < 1):",
    round(mean(diag(cor(vst_mat, vst_mat_corrected))), 4), "\n\n")

# Build a "vsd-like" object with corrected data for downstream functions
vsd_corrected <- vsd
assay(vsd_corrected) <- vst_mat_corrected

# =============================================================
#  Step 5.2: PCA on corrected data
# ================
