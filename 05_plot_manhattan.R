# create Manhattan plots from summary stats (Figure 5A), magma.genes.out (Figure 5B), and FINDOR (Figure S4)

# libraries
library(qqman)
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(ggpubr)

# read summary stats
gwas.snp <- read.table(file = "results/04_GWAS/manhattan.var", header = TRUE)
gwas.snp$SNP <- as.character(gwas.snp$ID)
gwas.snp$P <- as.numeric(as.character(gwas.snp$P))
gwas.snp$BP <- as.numeric(as.character(gwas.snp$BP))
gwas.snp$CHR <- as.numeric(as.character(gwas.snp$CHR))

# read gene-based output from magma
gwas.gene <- read.table(file = "results/05_FUMA_VA/magma.genes.out", header = TRUE)
gwas.gene <- gwas.gene %>%
  select(SYMBOL, P, START, CHR) %>%
  rename(SNP = SYMBOL) %>%
  rename(BP = START)
gwas.gene$SNP <- as.character(gwas.gene$SNP)
gwas.gene$P <- as.numeric(as.character(gwas.gene$P))
gwas.gene$BP <- as.numeric(as.character(gwas.gene$BP))
gwas.gene$CHR <- as.numeric(as.character(gwas.gene$CHR))

# read sumstats from FINDOR
gwas.findor <- read.table(file = "results/04_GWAS/vision.data.reweighted", header = TRUE)
gwas.findor <- gwas.findor %>%
  select(SNP, P_weighted, BP, CHR) %>%
  rename(P = P_weighted)
gwas.findor$SNP <- as.character(gwas.findor$SNP)
gwas.findor$P <- as.numeric(as.character(gwas.findor$P))
gwas.findor$BP <- as.numeric(as.character(gwas.findor$BP))
gwas.findor$CHR <- as.numeric(as.character(gwas.findor$CHR))

mypalette <- c("#FEE8C8", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000") 
sig.snp = 5e-8 # significant threshold line
sugg.snp = 1e-5 # suggestive threshold line
sig.gene = 2.7e-6 # significant threshold line (GWAS: SNPs were mapped to 18360 protein-coding genes)
sugg.gene = 1e-4 # suggestive threshold line

# Define Function ====
gg.manhattan <- function(df, sig, sugg, threshold, hlight, col, ylims){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by = c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate(BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate(is_highlight = ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate(is_annotate = ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 2) +
    scale_color_manual(values = rep(col, 22 )) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(sig)) +
    geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    
    # Custom the theme:
    theme_bw(base_size = 16) +
    theme(
      text = element_text(size = 16),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.margin = unit(c(1,0,2,2), "lines")
    )
}

# plot figures
Fig5a <- gg.manhattan(gwas.snp, sig.snp, sugg.snp, NA, NA, ylims = c(0,8), col = mypalette)
Fig5b <- gg.manhattan(gwas.gene, sig.gene, sugg.gene, sugg.gene, NA, ylims = c(0,7), col = mypalette)

FigS4 <- gg.manhattan(gwas.findor, sig.snp, sugg.snp, NA, NA, ylims = c(0,8), col = mypalette)

# tiff in publication quality
tiff("outputs/plots/Vision_Fig5_manhattan_plots.tiff", units = "in", width = 15, height = 12, res = 300)
ggarrange(Fig5a, Fig5b,
          labels = c("A", "B"),
          ncol = 1, nrow = 2,
          font.label = list(size = 18))
dev.off()

tiff("outputs/plots/Vision_FigS3_qq_snp.tiff", units = "in", width = 7, height = 7, res = 300)
qq(gwas.snp$P)
dev.off()

tiff("outputs/plots/Vision_FigS4_manhattan_plot_findor.tiff", units = "in", width = 15, height = 6, res = 300)
FigS4
dev.off()

tiff("outputs/plots/Vision_FigS5_qq_gene.tiff", units = "in", width = 7, height = 7, res = 300)
qq(gwas.gene$P)
dev.off()
