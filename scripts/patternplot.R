library(DESeq2); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
df <- readRDS(file.path(results_dir, "pattern_classification_tau115.rds"))

# 拆 P5→降/升,再重映射到新编号
df <- df %>% mutate(
  pat_old = case_when(
    pattern=="P5" & step2 <  0 ~ "P5old_down",
    pattern=="P5" & step2 >= 0 ~ "P5old_up",
    TRUE ~ pattern),
  pattern = recode(pat_old,
                   "P1"="P1", "P4"="P2", "P5old_up"="P3",
                   "P2"="P4", "P3"="P5", "P5old_down"="P6"))

cat("=== 新编号计数 ===\n")
print(table(factor(df$pattern, levels=paste0("P",1:6))))
cat("合计:", nrow(df), " | P3(Netrin-up):",
    paste(na.omit(df$symbol[df$pattern=="P3"]), collapse=", "), "\n\n")

cond_levels <- c("DIV1","DIV2","DIV2 + Netrin-1")
long <- df %>%
  mutate(fc_DIV1=1, fc_DIV2=DIV2/DIV1, fc_Netrin=Netrin/DIV1) %>%
  select(symbol, ensembl, pattern, fc_DIV1, fc_DIV2, fc_Netrin) %>%
  pivot_longer(starts_with("fc_"), names_to="cond", values_to="fc") %>%
  mutate(cond=factor(recode(cond, fc_DIV1="DIV1", fc_DIV2="DIV2", fc_Netrin="DIV2 + Netrin-1"),
                     levels=cond_levels))

# 配色: P1灰(空), P2/P3红(Netrin上调), P4/P5/P6蓝(Netrin下调)
RED <- "#C73E3E"; BLUE <- "#3B7AB8"; GREY <- "#9AA0A6"
pat_meta <- list(
  P1 = list(title="P1 · Monotonic up",   col=RED),
  P2 = list(title="P2 · Dip at DIV2",    col=RED),
  P3 = list(title="P3 · Netrin-up",      col=RED),
  P4 = list(title="P4 · Monotonic down", col=BLUE),
  P5 = list(title="P5 · Peak at DIV2",   col=BLUE),
  P6 = list(title="P6 · Netrin-down",    col=BLUE)
)
panel_order <- paste0("P",1:6)
Y_LIM <- c(0.11, 14); Y_BREAKS <- c(0.125,0.25,0.5,1,2,4,8)

make_panel <- function(pcode) {
  meta <- pat_meta[[pcode]]
  sub  <- long %>% filter(pattern == pcode)
  n_genes <- length(unique(df$ensembl[df$pattern==pcode]))
  base <- ggplot(mapping=aes(x=cond, y=fc)) +
    geom_hline(yintercept=1, linetype="dashed", color="grey60", linewidth=0.4) +
    scale_x_discrete(limits=cond_levels, drop=FALSE) +
    scale_y_continuous(trans="log2", breaks=Y_BREAKS) +
    coord_cartesian(ylim=Y_LIM) +
    labs(title=meta$title, subtitle=sprintf("n = %d", n_genes),
         x=NULL, y="Fold change vs DIV1") +
    theme_bw(base_size=11) +
    theme(plot.title=element_text(size=12, face="bold", color=meta$col),  # 标题也用该色(P1灰)
          plot.subtitle=element_text(size=9, color="grey40"),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=8),
          plot.margin=margin(5,10,5,5))
  if (n_genes == 0)
    return(base + annotate("text",x=2,y=1.7,label="n = 0",size=5,color="grey50") +
             annotate("text",x=2,y=0.9,label="(no gene)",size=3.5,color="grey60"))
  mean_traj <- sub %>% group_by(cond) %>% summarise(fc=mean(fc), .groups="drop")
  base +
    geom_line(data=sub, aes(group=ensembl), alpha=0.35, linewidth=0.5, color=meta$col) +
    geom_line(data=mean_traj, aes(group=1), linewidth=1.1, color=meta$col) +
    geom_point(data=mean_traj, aes(group=1), size=2.4, color=meta$col)
}

panels <- lapply(panel_order, make_panel)
combined <- wrap_plots(panels, nrow=2) +
  plot_annotation(
    title    = "Expression trajectory patterns across DIV1, DIV2 and DIV2 + Netrin-1",
    subtitle = "66 Netrin-responsive DEGs · top row up-regulated by Netrin-1 (P1-P3) · bottom row down-regulated (P4-P6) · 1.15-fold threshold · P1=0",
    theme = theme(plot.title=element_text(size=14, face="bold"),
                  plot.subtitle=element_text(size=10, color="grey40")))
print(combined)