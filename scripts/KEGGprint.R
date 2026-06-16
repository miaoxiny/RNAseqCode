library(ggplot2)
library(dplyr)
library(stringr)

gsea_dir <- file.path(results_dir, "03_GSEA")
gsea_kegg <- readRDS(file.path(gsea_dir, "GSEA_KEGG_netrin_final.rds"))

top_n <- 10
kegg_bar <- gsea_kegg@result %>%
  filter(!is.na(p.adjust), p.adjust < 0.05) %>%
  mutate(direction = ifelse(NES > 0, "Enriched in DIV2 + Netrin-1", "Enriched in DIV2")) %>%
  group_by(direction) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  mutate(Description_wrapped = str_wrap(Description, width = 45)) %>%
  arrange(NES) %>%
  mutate(Description_wrapped = factor(Description_wrapped, levels = Description_wrapped))

n_sig <- sum(gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)
cat("=== KEGG plot data ===\n")
cat(sprintf("  显著总数: %d\n", n_sig))
cat(sprintf("  图中激活: %d | 抑制: %d\n",
            sum(kegg_bar$direction == "Enriched in DIV2 + Netrin-1"),
            sum(kegg_bar$direction == "Enriched in DIV2")))

RED <- "#C73E3E"; BLUE <- "#3B7AB8"

p_kegg_bar <- ggplot(kegg_bar, aes(x = NES, y = Description_wrapped, fill = direction)) +
  geom_col(width = 0.72, alpha = 0.92) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = c("Enriched in DIV2" = BLUE,
                               "Enriched in DIV2 + Netrin-1" = RED), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title    = "KEGG Pathway GSEA",
       subtitle = sprintf("DIV2 + Netrin-1 vs DIV2 control · top %d pathways per direction · %d significant KEGG pathways (padj < 0.05)", top_n, n_sig),
       x = "Normalized enrichment score (NES)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(size = 14, face = "bold", color = "black"),
        plot.subtitle = element_text(size = 9.5, color = "grey25"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text   = element_text(color = "black"),
        axis.title  = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        axis.text.y = element_text(size = 9, color = "black"))

print(p_kegg_bar)