# =============================================================
#  08_DEG_visualize.R — Step 8 V2 (FINAL)
#  CSV export + Volcano plot — all RAW LFC
#  Latest adjustments:
#    - Legend: "Upregulated" / "Downregulated" (no hyphen)
#    - Point sizes: ns=1.4, relaxed=2.4, strict=3.2
#    - X-axis breaks: -6, -4, -2, 0, 2, 4, 6
#    - Legend dots sized to match plot points
# =============================================================

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
deg_dir     <- file.path(results_dir, "02_DEG")
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

# ─── Thresholds ───
padj_cutoff <- 0.05
lfc_cutoff  <- 0.585

# ─── Load results (raw LFC from Step 7 V2) ───
res_netrin_df <- readRDS(file.path(results_dir, "res_netrin_vs_div2.rds"))
res_div_df    <- readRDS(file.path(results_dir, "res_div2_vs_div1.rds"))

cat("=== Loaded annotated DEG tables (raw LFC) ===\n")
cat("Netrin vs DIV2:", nrow(res_netrin_df), "rows\n")
cat("DIV2 vs DIV1:  ", nrow(res_div_df),    "rows\n\n")

# =============================================================
#  Step 8.1: Define tier subsets
# =============================================================

relaxed_netrin <- res_netrin_df %>%
  filter(!is.na(padj), padj < padj_cutoff) %>%
  arrange(padj)
strict_netrin <- relaxed_netrin %>%
  filter(abs(log2FoldChange) > lfc_cutoff)

relaxed_div <- res_div_df %>%
  filter(!is.na(padj), padj < padj_cutoff) %>%
  arrange(padj)
strict_div <- relaxed_div %>%
  filter(abs(log2FoldChange) > lfc_cutoff)

cat("=== DEG tier summary ===\n")
cat(sprintf("Netrin vs DIV2 — STRICT:  %d genes\n", nrow(strict_netrin)))
cat(sprintf("Netrin vs DIV2 — RELAXED: %d genes\n", nrow(relaxed_netrin)))
cat(sprintf("DIV2  vs DIV1  — STRICT:  %d genes\n", nrow(strict_div)))
cat(sprintf("DIV2  vs DIV1  — RELAXED: %d genes\n\n", nrow(relaxed_div)))

# =============================================================
#  Step 8.2: Export CSV tables (6 files)
# =============================================================

write.csv(res_netrin_df,
          file.path(deg_dir, "DEG_full_DIV2-Netrin_vs_DIV2.csv"),
          row.names = FALSE)
write.csv(relaxed_netrin,
          file.path(deg_dir, "DEG_significant_DIV2-Netrin_vs_DIV2.csv"),
          row.names = FALSE)
write.csv(strict_netrin,
          file.path(deg_dir, "DEG_highlighted_DIV2-Netrin_vs_DIV2.csv"),
          row.names = FALSE)

write.csv(res_div_df,
          file.path(deg_dir, "DEG_full_DIV2_vs_DIV1.csv"),
          row.names = FALSE)
write.csv(relaxed_div,
          file.path(deg_dir, "DEG_significant_DIV2_vs_DIV1.csv"),
          row.names = FALSE)
write.csv(strict_div,
          file.path(deg_dir, "DEG_highlighted_DIV2_vs_DIV1.csv"),
          row.names = FALSE)

cat("=== CSV files saved (6 files) ===\n")
cat("Primary contrast (Netrin vs DIV2):\n")
cat("  DEG_full_DIV2-Netrin_vs_DIV2.csv          (", nrow(res_netrin_df),    "rows)\n")
cat("  DEG_significant_DIV2-Netrin_vs_DIV2.csv   (", nrow(relaxed_netrin),   "rows)\n")
cat("  DEG_highlighted_DIV2-Netrin_vs_DIV2.csv   (", nrow(strict_netrin),    "rows)\n")
cat("Secondary contrast (DIV2 vs DIV1):\n")
cat("  DEG_full_DIV2_vs_DIV1.csv                 (", nrow(res_div_df),       "rows)\n")
cat("  DEG_significant_DIV2_vs_DIV1.csv          (", nrow(relaxed_div),      "rows)\n")
cat("  DEG_highlighted_DIV2_vs_DIV1.csv          (", nrow(strict_div),       "rows)\n\n")

# =============================================================
#  Step 8.3: Volcano plot — primary contrast
# =============================================================

# Build plot data — all raw LFC
volcano_df <- res_netrin_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    neg_log10_padj = pmin(neg_log10_padj, 25),
    direction = case_when(
      padj < padj_cutoff & log2FoldChange >  0 ~ "up",
      padj < padj_cutoff & log2FoldChange <  0 ~ "down",
      TRUE                                      ~ "ns"
    ),
    is_strict = padj < padj_cutoff &
      abs(log2FoldChange) > lfc_cutoff
  )

cat("=== Volcano data composition ===\n")
print(table(volcano_df$direction, volcano_df$is_strict))
cat("\n")

# Colors (red/blue convention)
direction_colors <- c(
  "up"   = "#C73E3E",
  "down" = "#3B7AB8",
  "ns"   = "#CFCFCC"
)

# Symmetric x-axis
x_limit <- max(abs(volcano_df$log2FoldChange), na.rm = TRUE) * 1.05

# Build plot
p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = neg_log10_padj)) +
  # Layer 1: ns points (gray cloud, larger than before)
  geom_point(data = filter(volcano_df, direction == "ns"),
             aes(color = direction),
             size = 1.4, alpha = 0.5) +
  # Layer 2: relaxed significant (medium points, no label)
  geom_point(data = filter(volcano_df,
                           direction != "ns" & !is_strict),
             aes(color = direction),
             size = 2.4, alpha = 0.85) +
  # Layer 3: strict (large points, will be labeled)
  geom_point(data = filter(volcano_df, is_strict),
             aes(color = direction),
             size = 3.2, alpha = 0.95) +
  # Threshold lines
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  # Labels — only on strict genes with symbol
  geom_text_repel(
    data = filter(volcano_df, is_strict, !is.na(symbol)),
    aes(label = symbol),
    size = 3.3, max.overlaps = Inf,
    box.padding = 0.45, point.padding = 0.25,
    segment.color = "grey55", segment.size = 0.3,
    color = "grey20", min.segment.length = 0.2,
    force = 2
  ) +
  scale_color_manual(
    values = direction_colors,
    labels = c("up"   = "Upregulated",
               "down" = "Downregulated",
               "ns"   = "Not significant"),
    breaks = c("up", "down", "ns"),
    name   = NULL
  ) +
  # Symmetric x-axis with custom breaks
  scale_x_continuous(limits = c(-x_limit, x_limit),
                     breaks = seq(-6, 6, by = 2)) +
  labs(
    title    = "Transcriptional response to Netrin-1",
    subtitle = sprintf("DIV2 + Netrin-1 (24 h) vs DIV2 control · %d highlighted DEGs (%d up + %d down)",
                       nrow(strict_netrin),
                       sum(strict_netrin$log2FoldChange > 0),
                       sum(strict_netrin$log2FoldChange < 0)),
    x        = expression(log[2]~"fold change"),
    y        = expression(-log[10]~"adjusted"~italic(p))
  ) +
  # Override legend dot sizes to match plot points
  guides(color = guide_legend(
    override.aes = list(size = c(3.2, 3.2, 1.6), alpha = 0.9)
  )) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 13, face = "bold"),
        plot.subtitle    = element_text(size = 10, color = "grey40"),
        legend.position  = "right",
        legend.text      = element_text(size = 10))

# ─── DISPLAY ───
print(p_volcano)