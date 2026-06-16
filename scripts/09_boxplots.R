# =============================================================
#  09_boxplots.R — Step 9: Boxplots of 8 key DEGs
#    2x4 grid (top: 4 down, bottom: 4 up)
#    DESeq2 normalized counts, per-sample dots
# =============================================================

library(DESeq2)
library(dplyr)
library(ggplot2)
library(patchwork)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
deg_dir     <- file.path(results_dir, "02_DEG")

# ─── Load data ───
dds <- readRDS(file.path(results_dir, "dds_fitted.rds"))
res_netrin_df <- readRDS(file.path(results_dir, "res_netrin_vs_div2.rds"))

# ─── Get normalized counts ───
norm_counts <- counts(dds, normalized = TRUE)
coldata     <- as.data.frame(colData(dds))

# ─── Condition labels (project palette) ───
condition_display <- c(
  "DIV1"        = "DIV1",
  "DIV2"        = "DIV2",
  "DIV2_Netrin" = "DIV2 + Netrin-1"
)
coldata$Condition <- factor(
  condition_display[as.character(coldata$condition)],
  levels = c("DIV1", "DIV2", "DIV2 + Netrin-1")
)

condition_colors <- c(
  "DIV1"             = "#DCD1BC",
  "DIV2"             = "#2E5F58",
  "DIV2 + Netrin-1"  = "#D88040"
)

# ─── Selected genes ───
genes_down <- c("Olig1", "Ccnd1", "Nbl1", "Itgb5")
genes_up   <- c("C3", "Acod1", "Adgre1", "Apoe")
all_genes  <- c(genes_down, genes_up)

# Get full names for subtitle
gene_info <- res_netrin_df %>%
  filter(symbol %in% all_genes) %>%
  select(symbol, gene_name, ensembl_id, log2FoldChange, padj)

cat("=== Selected genes ===\n")
print(gene_info %>%
        mutate(log2FoldChange = round(log2FoldChange, 2),
               padj           = signif(padj, 3)))
cat("\n")

# ─── Plot helper ───
make_boxplot <- function(gene_symbol, gene_info, norm_counts, coldata) {
  # Look up ensembl ID and gene name
  info  <- gene_info %>% filter(symbol == gene_symbol)
  ens   <- info$ensembl_id[1]
  gname <- info$gene_name[1]
  
  # Build plot data
  plot_df <- data.frame(
    Condition  = coldata$Condition,
    expression = norm_counts[ens, ]
  )
  
  # Trim gene_name if too long (keep subtitle readable)
  if (nchar(gname) > 38) {
    gname <- paste0(substr(gname, 1, 35), "...")
  }
  
  ggplot(plot_df, aes(x = Condition, y = expression, fill = Condition)) +
    geom_boxplot(width = 0.55, alpha = 0.85,
                 outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(width = 0.15, size = 2.0, alpha = 0.95,
                stroke = 0.3, color = "grey20") +
    scale_fill_manual(values = condition_colors) +
    labs(title    = gene_symbol,
         subtitle = gname,
         x = NULL,
         y = "Normalized counts") +
    theme_bw(base_size = 11) +
    theme(legend.position  = "none",
          plot.title       = element_text(size = 13, face = "italic"),
          plot.subtitle    = element_text(size = 9, color = "grey40"),
          panel.grid.minor = element_blank(),
          axis.text.x      = element_text(size = 9))
}

# ─── Build all 8 plots ───
plot_list <- lapply(all_genes, make_boxplot,
                    gene_info  = gene_info,
                    norm_counts = norm_counts,
                    coldata     = coldata)
names(plot_list) <- all_genes

# ─── Combine into 2x4 grid ───
# Top row: 4 down-regulated
# Bottom row: 4 up-regulated
combined <- (plot_list$Olig1 | plot_list$Ccnd1 | plot_list$Nbl1 | plot_list$Itgb5) /
  (plot_list$C3    | plot_list$Acod1 | plot_list$Adgre1 | plot_list$Apoe) +
  plot_annotation(
    title    = "Key DEGs: expression across conditions",
    subtitle = "Top row: Netrin-1 down-regulated  ·  Bottom row: Netrin-1 up-regulated  ·  n = 3 per condition",
    theme    = theme(plot.title    = element_text(size = 15, face = "bold"),
                     plot.subtitle = element_text(size = 11, color = "grey40"))
  )

print(combined)