# =============================================================
#  06_DESeq2_fit.R — Step 6: DESeq2 model fitting
#  Inspect-only mode (display dispersion plot, no auto-save)
# =============================================================

library(DESeq2)
library(ggplot2)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")

# ─── Load filtered DESeq2 object from Step 3 ───
dds <- readRDS(file.path(results_dir, "dds_filtered.rds"))

cat("=== Loaded dds ===\n")
cat("Genes:  ", nrow(dds), "\n")
cat("Samples:", ncol(dds), "\n")
cat("Design: ", deparse(design(dds)), "\n\n")

# ─── Verify design includes batch + condition ───
stopifnot(identical(as.character(design(dds)),
                    c("~", "batch", "+", "condition")))

# =============================================================
#  Step 6.1: Run DESeq2
# =============================================================
# This does: size factor estimation, dispersion estimation, Wald test

cat("=== Running DESeq2... ===\n")
dds <- DESeq(dds)
cat("DESeq2 fit complete.\n\n")

# =============================================================
#  Step 6.2: Inspect size factors
# =============================================================

cat("=== Size factors ===\n")
sf <- sizeFactors(dds)
print(round(sf, 3))
cat("\nSize factor summary:\n")
print(summary(sf))
cat("\n")

# =============================================================
#  Step 6.3: Inspect dispersion estimates
# =============================================================

cat("=== Dispersion estimates ===\n")
disp_summary <- list(
  n_genes        = nrow(dds),
  median_disp    = median(dispersions(dds), na.rm = TRUE),
  mean_disp      = mean(dispersions(dds),   na.rm = TRUE),
  range_disp     = range(dispersions(dds),  na.rm = TRUE)
)
cat(sprintf("  Genes:           %d\n", disp_summary$n_genes))
cat(sprintf("  Median dispersion: %.4f\n", disp_summary$median_disp))
cat(sprintf("  Mean dispersion:   %.4f\n", disp_summary$mean_disp))
cat(sprintf("  Range:             %.4f to %.4f\n",
            disp_summary$range_disp[1], disp_summary$range_disp[2]))
cat("\n")

# =============================================================
#  Step 6.4: Display dispersion plot
# =============================================================
# This is the most important QC plot for the DESeq2 model

plotDispEsts(dds,
             main = "Dispersion estimates",
             genecol = "grey60",
             fitcol  = "#D88040",   # warm orange (project palette)
             finalcol = "#2E5F58")  # dark teal (project palette)

# =============================================================
#  Step 6.5: Inspect available contrasts
# =============================================================

cat("=== Available results / coefficients ===\n")
print(resultsNames(dds))
cat("\n")

# =============================================================
#  Step 6.6: Save fitted dds object
# =============================================================

saveRDS(dds, file.path(results_dir, "dds_fitted.rds"))
cat("Saved: dds_fitted.rds (ready for Step 7 results extraction)\n")