## ============================================================================
## GSEA 可视化 v2 — 全部从 _v2full 出(不用 simplified)；风格同 10_gsea_plots.R
## ============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(ggrepel) })
out_base <- "D:/Dropbox/Dropbox/RNAseq 2025/Time analysis"
gsea_dir <- file.path(out_base, "results/04_GSEA")

gsea_bar <- function(csv, title, sub, topn=15, outfile, w=10, h=8) {
  sr <- read.csv(file.path(gsea_dir, csv), stringsAsFactors=FALSE)
  s  <- sr %>% dplyr::filter(p.adjust<0.05)
  up <- s %>% dplyr::filter(NES>0) %>% dplyr::arrange(desc(NES)) %>% head(topn)
  dn <- s %>% dplyr::filter(NES<0) %>% dplyr::arrange(NES)       %>% head(topn)
  d  <- dplyr::bind_rows(up, dn)
  d$Description <- factor(d$Description, levels=d$Description[order(d$NES)])
  p <- ggplot(d, aes(NES, Description, fill=NES>0)) +
    geom_col(width=0.7) + geom_vline(xintercept=0, color="grey40", linewidth=0.4) +
    scale_fill_manual(values=c(`TRUE`="#C73E3E",`FALSE`="#3B7AB8"),
                      labels=c(`TRUE`="Activated (NES>0)",`FALSE`="Suppressed (NES<0)"), name=NULL) +
    labs(title=title, subtitle=sub, x="NES", y=NULL) +
    theme_bw(base_size=11) +
    theme(panel.grid.minor=element_blank(), axis.text.y=element_text(size=11),
          plot.title=element_text(size=13, face="bold"),
          plot.subtitle=element_text(size=9, color="grey40"), legend.position="bottom")
  fp <- file.path(gsea_dir, paste0(outfile,".pdf"))
  fg <- file.path(gsea_dir, paste0(outfile,".png"))
  stopifnot(!file.exists(fp), !file.exists(fg))
  ggsave(fp, p, width=w, height=h); ggsave(fg, p, width=w, height=h, dpi=300)
  cat("[fig]", outfile, "\n")
}

## ===== GO BP（全谱）=====
gsea_bar("GSEA_GO_BP_5v0_v2full.csv",
         "GSEA GO BP * 5min Netrin-1 vs Vehicle",
         "Top 15 activated + top 15 suppressed (by NES); FDR<0.05",
         15, "08_GSEA_GO_BP_5v0_barplot_v2")
gsea_bar("GSEA_GO_BP_15v5_v2full.csv",
         "GSEA GO BP * 15min vs 5min Netrin-1 (time evolution)",
         "Top 15 each (by NES); FDR<0.05  Mirror of 5v0",
         15, "10_GSEA_GO_BP_15v5_barplot_v2")

## ===== GO CC（全谱，新增）=====
gsea_bar("GSEA_GO_CC_5v0_v2full.csv",
         "GSEA GO CC * 5min Netrin-1 vs Vehicle",
         "Top 15 each (by NES); FDR<0.05  Cellular components",
         15, "14_GSEA_GO_CC_5v0_barplot_v2")
gsea_bar("GSEA_GO_CC_15v5_v2full.csv",
         "GSEA GO CC * 15min vs 5min Netrin-1 (time evolution)",
         "Top 15 each (by NES); FDR<0.05", 15, "15_GSEA_GO_CC_15v5_barplot_v2")

## ===== KEGG（全谱）=====
gsea_bar("GSEA_KEGG_5v0_v2full.csv", "GSEA KEGG * 5min Netrin-1 vs Vehicle",
         "All significant pathways (by NES); FDR<0.05", 50, "09_GSEA_KEGG_5v0_barplot_v2", h=6)
gsea_bar("GSEA_KEGG_15v5_v2full.csv", "GSEA KEGG * 15min vs 5min Netrin-1 (time evolution)",
         "Top 15 each of significant (by NES); FDR<0.05", 15, "11_GSEA_KEGG_15v5_barplot_v2", h=7)
gsea_bar("GSEA_KEGG_15v0_v2full.csv", "GSEA KEGG * 15min Netrin-1 vs Vehicle",
         "All significant; FDR<0.05  (GO BP: 0 sig = returned to baseline)", 50,
         "12_GSEA_KEGG_15v0_barplot_v2", h=4.5)

## ---- 镜像散点：完整谱 _v2full ----
a <- read.csv(file.path(gsea_dir,"GSEA_GO_BP_5v0_v2full.csv"))
b <- read.csv(file.path(gsea_dir,"GSEA_GO_BP_15v5_v2full.csv"))
m <- inner_join(a %>% dplyr::select(ID, Description, NES_5v0=NES, padj_5v0=p.adjust),
                b %>% dplyr::select(ID, NES_15v5=NES, padj_15v5=p.adjust), by="ID") %>%
  dplyr::filter(!is.na(NES_5v0), !is.na(NES_15v5))
m$pattern <- dplyr::case_when(
  m$NES_5v0>0 & m$NES_15v5<0 ~ "transient up (5min up, 15min down)",
  m$NES_5v0<0 & m$NES_15v5>0 ~ "transient down (5min down, 15min up)",
  TRUE ~ "same direction")
r <- cor(m$NES_5v0, m$NES_15v5)
cat(sprintf("镜像: 共有通路 %d, r = %.3f\n", nrow(m), r))
sig <- m %>% dplyr::filter(padj_5v0<0.05, padj_15v5<0.05) %>%
  dplyr::mutate(score=abs(NES_5v0)+abs(NES_15v5))
key <- dplyr::bind_rows(
  sig %>% dplyr::filter(NES_5v0>0) %>% dplyr::arrange(desc(score)) %>% head(5),
  sig %>% dplyr::filter(NES_5v0<0) %>% dplyr::arrange(desc(score)) %>% head(5))
pal <- c("transient up (5min up, 15min down)"="#C73E3E",
         "transient down (5min down, 15min up)"="#3B7AB8","same direction"="grey70")
ps <- ggplot(m, aes(NES_5v0, NES_15v5)) +
  geom_hline(yintercept=0, color="grey80", linewidth=0.3) +
  geom_vline(xintercept=0, color="grey80", linewidth=0.3) +
  geom_abline(slope=-1, intercept=0, linetype="dashed", color="grey60") +
  geom_point(aes(color=pattern), size=2, alpha=0.6) +
  ggrepel::geom_text_repel(data=key, aes(label=Description, color=pattern), size=2.8,
                           max.overlaps=20, min.segment.length=0, show.legend=FALSE, box.padding=0.5, seed=1) +
  scale_color_manual(values=pal, name=NULL) +
  labs(title="GSEA pathway mirror: 5v0 vs 15v5 (GO BP)",
       subtitle=sprintf("Each point = one shared pathway (n=%d). NES r = %.2f. Dashed = perfect mirror (y=-x)", nrow(m), r),
       x="NES, 5min vs Vehicle (5v0)", y="NES, 15min vs 5min (15v5)") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(size=13, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"), legend.position="bottom")
fp <- file.path(gsea_dir,"13_GSEA_mirror_scatter_5v0_vs_15v5_v2.pdf")
fg <- file.path(gsea_dir,"13_GSEA_mirror_scatter_5v0_vs_15v5_v2.png")
stopifnot(!file.exists(fp), !file.exists(fg))
ggsave(fp, ps, width=9, height=8); ggsave(fg, ps, width=9, height=8, dpi=300)
cat(sprintf("[fig] 13_mirror_v2  r=%.3f\n", r))