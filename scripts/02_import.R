# =============================================================
#  02_import.R
#  Step 2: Data Import - Build sample sheet and import .sf files
#  Project: Netrin-1 RNA-seq (Group A)
# =============================================================

library(tximport)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
data_dir    <- file.path(project_dir, "data", "sf_files")
results_dir <- file.path(project_dir, "results")

# ─── List .sf files ───
sf_files <- list.files(data_dir, pattern = "\\.sf$", full.names = TRUE)
cat("Found", length(sf_files), "sf files\n")

# ─── Build sample sheet ───
coldata <- data.frame(
  sample = c("1D1", "1D2", "1D2N",
             "2D1", "2D2", "2D2N",
             "3D1", "3D2", "3D2N"),
  condition = factor(
    c("DIV1",  "DIV2",  "DIV2_Netrin",
      "DIV1",  "DIV2",  "DIV2_Netrin",
      "DIV1",  "DIV2",  "DIV2_Netrin"),
    levels = c("DIV2", "DIV2_Netrin", "DIV1")
  ),
  batch = factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
  stringsAsFactors = FALSE
)
rownames(coldata) <- coldata$sample

cat("\n=== Sample sheet ===\n")
print(coldata)

# ─── Match files to samples ───
expected_pattern <- c("1D1_",  "1D2_",  "1D2N_",
                      "2D1_",  "2D2_",  "2D2N_",
                      "3D1_",  "3D2_",  "3D2N_")

ordered_files <- sapply(expected_pattern, function(p) {
  matched <- sf_files[grepl(paste0("^", p), basename(sf_files))]
  if (length(matched) != 1) stop("Could not match: ", p)
  matched
})
names(ordered_files) <- coldata$sample

cat("\n=== File-to-sample mapping ===\n")
for (i in seq_along(ordered_files)) {
  cat(sprintf("  %-5s -> %s\n", names(ordered_files)[i],
              basename(ordered_files[i])))
}

# ─── Import with tximport ───
cat("\n=== Importing with tximport... ===\n")
txi <- tximport(
  files       = ordered_files,
  type        = "salmon",
  txIn        = FALSE,
  txOut       = FALSE,
  dropInfReps = TRUE
)

# ─── Inspect ───
cat("\n=== Import successful ===\n")
cat("Number of genes:  ", nrow(txi$counts), "\n")
cat("Number of samples:", ncol(txi$counts), "\n\n")

cat("=== Count matrix preview ===\n")
print(round(head(txi$counts, 5), 2))

cat("\n=== Per-sample library size (millions) ===\n")
lib_sizes <- colSums(txi$counts) / 1e6
print(round(lib_sizes, 2))

# ─── Save intermediate objects ───
saveRDS(txi,     file.path(results_dir, "txi.rds"))
saveRDS(coldata, file.path(results_dir, "coldata.rds"))
cat("\nSaved txi.rds and coldata.rds to results/\n")