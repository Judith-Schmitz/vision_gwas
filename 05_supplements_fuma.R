# create supplementary tables and structure FUMA output

# libraries
library(readr)
library(dplyr)
library(tidyverse)
library(xlsx)
library(summarytools)


### Table S6: Top SNPs (association p < .00001) for VA GWAS.
tophits.var <- read.table(file = "results/04_GWAS/tophits.var", header = TRUE)
annov.var <- read.table(file = "results/05_FUMA_VA/annov.txt", header = TRUE, sep = "\t")
annov.var <- within(annov.var, uniqID <- data.frame(do.call('rbind', strsplit(as.character(uniqID), ':', fixed = TRUE))))
annov.var$ID <- paste0(annov.var$uniqID$X1, ":", annov.var$uniqID$X2)
table.s6 <- left_join(tophits.var, annov.var, by = "ID") %>%
  select(SNP, ID, ALLELE0, ALLELE1, A1FREQ, CHISQ_BOLT_LMM_INF, BETA, SE, P_BOLT_LMM_INF, symbol, annot, dist)
table.s6 <- table.s6[order(table.s6$P_BOLT_LMM_INF),] 
#write.xlsx(table.s6, file = "results/05_FUMA_VA/Table_S6.xlsx")


###  Table S7: Results of replication analysis for VA. 
# Input file for FUMA
all.snps <- read.xlsx(file = "data/replication.xlsx", sheetIndex = 1) %>%
  filter(P < 5e-08) %>%
  filter(Study == "Fan_2016" | Study == "Hysi_2020" | Study == "Tedja_2018" | Study == "Verhoefen_2013" | Study == "Kiefer_2013")
all.snps <- all.snps[!duplicated(all.snps$SNP),]
top.snps <- all.snps %>%
  select(SNP, CHR, BP)
#write.table(top.snps, file = "outputs/replication.snps.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
# top.snps (1033 SNPs) has been used as input file for pre-defined lead SNPs for FUMA
loci <- read.table(file = "results/05_FUMA_VA_replication/snps.txt", header = TRUE)
loci <- loci %>%
  filter(!is.na(gwasP)) # overlap with ALSPAC: 743 SNPs
table.s7 <- loci %>%
  filter(gwasP < 0.05/length(unique(loci$GenomicLocus))) %>% 
  rename(SNP = rsID)
table.s7 <- table.s7[order(table.s7$gwasP),]
table.s7 <- left_join(table.s7, all.snps, by = 'SNP') %>%
  select(GenomicLocus, SNP, chr, pos, non_effect_allele, effect_allele, MAF, beta, se, gwasP, nearestGene, dist, func, Beta, P, Phenotype, Study)
rm(all.snps, top.snps, loci)
#write.xlsx(table.s7, file = "results/05_FUMA_VA_replication/Table_S7.xlsx")


### Table S8: Associations of independent significant SNPs (rs11656126) and correlated SNPs in previous GWAS obtained from GWAS catalog (FUMA). 
table.s8 <- read.table(file = "results/05_FUMA_VA/gwascatalog.txt", header = TRUE, sep = "\t", fill = TRUE)
table.s8 <- table.s8 %>%
  select(IndSigSNP, GenomicLocus, chr, bp, snp, FirstAuth, Trait, MappedGene, P, OrBeta)
#write.xlsx(table.s8, file = "results/05_FUMA_VA/Table_S8.xlsx")


### Table S9: eQTLs of top hits on chromosome 17 (FUMA), sorted in increasing order of VA GWAS association p value. 
eqtl.var <- read.table(file = "results/05_FUMA_VA/eqtl.txt", header = TRUE)
eqtl.var <- eqtl.var %>%
  mutate(ID = paste0(chr, ":", pos))
tophits.var <- tophits.var %>% select(ID, SNP, P_BOLT_LMM_INF)
table.s9 <- left_join(eqtl.var, tophits.var, by = "ID") %>%
  filter(!is.na(SNP))
table.s9 <- table.s9[order(table.s9$P_BOLT_LMM_INF),] 
table.s9 <- table.s9 %>%
  select(SNP, ID, db, tissue, symbol, RiskIncAllele, signed_stats, FDR)
#write.xlsx(table.s9, file = "results/05_FUMA_VA/Table_S9.xlsx")

### Table S10: FINDOR: Top SNPs (FINDOR re-weighted association p < .00001) for VA GWAS.
gwas.findor <- read.table(file = "results/04_GWAS/vision.data.reweighted", header = TRUE) %>%
  mutate(ID = paste0(CHR, ":", BP))
table.s10 <- left_join(gwas.findor, annov.var, by = "ID") %>%
  select(SNP, ID, Allele1, Allele2, A1FREQ, Z, P_weighted) %>%
  filter(P_weighted < .00001)
table.s10 <- table.s10[order(table.s10$P_weighted),] 
rm(annov.var, gwas.findor)
#write.xlsx(table.s10, file = "results/04_GWAS/Table_S10.xlsx")


### Table S11: Full results of PGS analysis for VA including sex, age, two PCs, and SES as covariates. See 08_pgs_plot.r




# No tables, to report in paper:

# Gene-based results
gene.var <- read.table(file = "results/05_FUMA_VA/magma.genes.out", header = TRUE)
gene.var <- gene.var %>%
  filter(P < 0.05/18360)
gene.var <- gene.var[order(gene.var$P),] 

# Gene set
gsa.var <- read.table(file = "results/05_FUMA_VA/magma.gsa.out", header = TRUE)
gsa.var <- gsa.var %>%
  filter(str_detect(VARIABLE, "GO_bp:go_")) %>%
  filter(P < 0.05/7343)
gsa.var <- gsa.var[order(gsa.var$P),] 







