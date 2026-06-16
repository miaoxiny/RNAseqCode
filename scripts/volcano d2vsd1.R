# =============================================================
#  14_volcano_div2_vs_div1.R — PREVIEW MODE
#  点大小对齐 Contrast 1: 显著=3.0, 不显著=1.3
#  标签: 每方向先滤无 symbol，再取 padj top 10（保证各 10 个）
# =============================================================

library(ggplot2)
library(ggrepel)
library(dplyr)

# ─── Paths ───
project_dir <- "D:/Dropbox/Dropbox/RNAseq 2025/Claude analysis"
results_dir <- file.path(project_dir, "results")
deg_dir     <- file.path(results_dir, "02_DEG")

# ─── Load DEG results ───
res_div_df <- readRDS(file.path(results_dir, "res_div2_vs_div1.rds"))

# =============================================================
#  Step 14.1: Prepare data (3-class)
# =============================================================

volcano_df <- res_div_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    neglog10_padj = -log10(padj),
    category = case_when(
      padj < 0.05 & abs(log2FoldChange) > 0.585 & log2FoldChange > 0  ~ "Upregulated",
      padj < 0.05 & abs(log2FoldChange) > 0.585 & log2FoldChange < 0  ~ "Downregulated",
      TRUE                                                              ~ "Not significant"
    )
  )

cat("=== Category counts ===\n")
print(volcano_df %>% count(category))

# =============================================================
#  Step 14.2: Select labels — 先滤无 symbol，再取 top 10
# =============================================================

top_up <- volcano_df %>%
  filter(category == "Upregulated", !is.na(symbol), symbol != "") %>%
  arrange(padj) %>%
  head(10) %>%
  pull(symbol)

top_down <- volcano_df %>%
  filter(category == "Downregulated", !is.na(symbol), symbol != "") %>%
  arrange(padj) %>%
  head(10) %>%
  pull(symbol)

label_genes <- unique(c(top_up, top_down))
label_genes <- label_genes[!is.na(label_genes) & label_genes != ""]

cat(sprintf("\n=== Labels: %d genes ===\n", length(label_genes)))
cat("UP:   ", paste(top_up,   collapse = ", "), "\n")
cat("DOWN: ", paste(top_down, collapse = ", "), "\n")

volcano_df <- volcano_df %>%
  mutate(label_this = symbol %in% label_genes & !is.na(symbol))

# =============================================================
#  Step 14.3: Color, size, layer order
# =============================================================

volcano_colors <- c(
  "Upregulated"     = "#C73E3E",
  "Downregulated"   = "#3B7AB8",
  "Not significant" = "#CFCFCC"
)

volcano_df$point_size <- case_when(
  volcano_df$category != "Not significant" ~ 3.0,
  TRUE                                      ~ 1.3
)

# Layer order: ns at bottom, significant on top
volcano_df <- volcano_df %>%
  mutate(plot_order = ifelse(category == "Not significant", 1, 2)) %>%
  arrange(plot_order)

# Subtitle
n_strict <- sum(volcano_df$category != "Not significant")
n_up     <- sum(volcano_df$category == "Upregulated")
n_down   <- sum(volcano_df$category == "Downregulated")

subtitle_text <- sprintf(
  "DIV2 vs DIV1 · %d highlighted DEGs (%d up + %d down)",
  n_strict, n_up, n_down
)

# =============================================================
#  Step 14.4: Build plot
# =============================================================

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange,
                        y = neglog10_padj,
                        color = category,
                        size = point_size)) +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.585, 0.585),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  
  geom_point(alpha = 0.85, stroke = 0) +
  
  scale_color_manual(
    values = volcano_colors,
    breaks = c("Upregulated", "Downregulated", "Not significant"),
    name = NULL
  ) +
  
  scale_size_identity() +
  
  ggrepel::geom_text_repel(
    data = subset(volcano_df, label_this),
    aes(label = symbol),
    size = 3.3, color = "grey15", max.overlaps = Inf,
    box.padding = 0.5, point.padding = 0.3,
    segment.color = "grey60", segment.size = 0.3,
    min.segment.length = 0, force = 5, force_pull = 0.5,
    seed = 42, bg.color = "white", bg.r = 0.15
  ) +
  
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 2)) +
  scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, by = 50)) +
  
  labs(
    title    = "Transcriptional changes: DIV1 \u2192 DIV2",
    subtitle = subtitle_text,
    x = expression(log[2]~"fold change"),
    y = expression(-log[10]~"adjusted "*italic(p))
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    legend.position  = "right",
    legend.text      = element_text(size = 13),
    legend.key.size  = unit(0.8, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.key       = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))

# =============================================================
#  Step 14.5: Display (PREVIEW)
# =============================================================

print(p_volcano)

# =============================================================
#  Step 14.6: Save — COMMENTED OUT FOR PREVIEW
#  解开注释前请先确认预览满意。文件名 05_volcano_div2_vs_div1，不覆盖 Contrast 1。
# =============================================================

# pdf(file.path(deg_dir, "05_volcano_div2_vs_div1.pdf"), width = 11, height = 9)
# print(p_volcano); dev.off()
#
# png(file.path(deg_dir, "05_volcano_div2_vs_div1.png"), width = 11, height = 9, units = "in", res = 300)
# print(p_volcano); dev.off()
#
# cat("\n=== Saved ===\n")
# cat("  results/02_DEG/05_volcano_div2_vs_div1.pdf\n")
# cat("  results/02_DEG/05_volcano_div2_vs_div1.png\n")