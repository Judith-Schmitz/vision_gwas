# PGS vision

library(tidyverse)
library(summarytools)
library(ggplot2)
library(ggpubr)
library(xlsx)
library(corrplot)
library(ggcorrplot)
library(RColorBrewer)

color.new <- brewer.pal(n = 8, name = 'OrRd')

thresholds <- c(0.00100005, 0.0500001, 0.1, 0.2, 0.3, 0.4, 0.5, 1)

pgs.ADHD.thresholds <- read.table(file = "results/07_PGS/pgs.covar.ses/ADHD.var.prsice", header = TRUE) %>%
  filter(Threshold %in% thresholds)
pgs.ADHD.best <- read.table(file = "results/07_PGS/pgs.covar.ses/ADHD.var.prsice", header = TRUE) %>%
  filter(R2 == max(R2))
pgs.ADHD <- rbind(pgs.ADHD.best, pgs.ADHD.thresholds)
pgs.ADHD <- pgs.ADHD[, 3:8]
rm(pgs.ADHD.best, pgs.ADHD.thresholds)
#write.xlsx(pgs.ADHD, file = "results/07_PGS/Table_S11.xlsx", sheetName = "ADHD", append = TRUE)

pgs.ASD.thresholds <- read.table(file = "results/07_PGS/pgs.covar.ses/ASD.var.prsice", header = TRUE) %>%
  filter(Threshold %in% thresholds)
pgs.ASD.best <- read.table(file = "results/07_PGS/pgs.covar.ses/ASD.var.prsice", header = TRUE) %>%
  filter(R2 == max(R2))
pgs.ASD <- rbind(pgs.ASD.best, pgs.ASD.thresholds)
pgs.ASD <- pgs.ASD[, 3:8]
rm(pgs.ASD.best, pgs.ASD.thresholds)
#write.xlsx(pgs.ASD, file = "results/07_PGS/Table_S11.xlsx", sheetName = "ASD", append = TRUE)

pgs.BIP.thresholds <- read.table(file = "results/07_PGS/pgs.covar.ses/BIP.var.prsice", header = TRUE) %>%
  filter(Threshold %in% thresholds)
pgs.BIP.best <- read.table(file = "results/07_PGS/pgs.covar.ses/BIP.var.prsice", header = TRUE) %>%
  filter(R2 == max(R2))
pgs.BIP <- rbind(pgs.BIP.best, pgs.BIP.thresholds)
pgs.BIP <- pgs.BIP[, 3:8]
rm(pgs.BIP.best, pgs.BIP.thresholds)
#write.xlsx(pgs.BIP, file = "results/07_PGS/Table_S11.xlsx", sheetName = "BIP", append = TRUE)

pgs.DYS.thresholds <- read.table(file = "results/07_PGS/pgs.covar.ses/DYS.var.prsice", header = TRUE) %>%
  filter(Threshold %in% thresholds)
pgs.DYS.best <- read.table(file = "results/07_PGS/pgs.covar.ses/DYS.var.prsice", header = TRUE) %>%
  filter(R2 == max(R2))
pgs.DYS <- rbind(pgs.DYS.best, pgs.DYS.thresholds)
pgs.DYS <- pgs.DYS[, 3:8]
rm(pgs.DYS.best, pgs.DYS.thresholds)
#write.xlsx(pgs.DYS, file = "results/07_PGS/Table_S11.xlsx", sheetName = "RD", append = TRUE)

pgs.SCZ.thresholds <- read.table(file = "results/07_PGS/pgs.covar.ses/SCZ.var.prsice", header = TRUE) %>%
  filter(Threshold %in% thresholds)
pgs.SCZ.best <- read.table(file = "results/07_PGS/pgs.covar.ses/SCZ.var.prsice", header = TRUE) %>%
  filter(R2 == max(R2))
pgs.SCZ <- rbind(pgs.SCZ.best, pgs.SCZ.thresholds)
pgs.SCZ <- pgs.SCZ[, 3:8]
rm(pgs.SCZ.best, pgs.SCZ.thresholds)
#write.xlsx(pgs.SCZ, file = "results/07_PGS/Table_S11.xlsx", sheetName = "SCZ", append = TRUE)

pgs.EA.thresholds <- read.table(file = "results/07_PGS/pgs.covar.ses/EA.var.prsice", header = TRUE) %>%
  filter(Threshold %in% thresholds)
pgs.EA.best <- read.table(file = "results/07_PGS/pgs.covar.ses/EA.var.prsice", header = TRUE) %>%
  filter(R2 == max(R2))
pgs.EA <- rbind(pgs.EA.best, pgs.EA.thresholds)
pgs.EA <- pgs.EA[, 3:8]
rm(pgs.EA.best, pgs.EA.thresholds)
#write.xlsx(pgs.EA, file = "results/07_PGS/Table_S11.xlsx", sheetName = "EA", append = TRUE)

pgs.INT.thresholds <- read.table(file = "results/07_PGS/pgs.covar.ses/INT.var.prsice", header = TRUE) %>%
  filter(Threshold %in% thresholds)
pgs.INT.best <- read.table(file = "results/07_PGS/pgs.covar.ses/INT.var.prsice", header = TRUE) %>%
  filter(R2 == max(R2))
pgs.INT <- rbind(pgs.INT.best, pgs.INT.thresholds)
pgs.INT <- pgs.INT[, 3:8]
rm(pgs.INT.best, pgs.INT.thresholds)
#write.xlsx(pgs.INT, file = "results/07_PGS/Table_S11.xlsx", sheetName = "INT", append = TRUE)

pgs.MYO.thresholds <- read.table(file = "results/07_PGS/pgs.covar.ses/MYO.var.prsice", header = TRUE) %>%
  filter(Threshold %in% thresholds)
pgs.MYO.best <- read.table(file = "results/07_PGS/pgs.covar.ses/MYO.var.prsice", header = TRUE) %>%
  filter(R2 == max(R2))
pgs.MYO <- rbind(pgs.MYO.best, pgs.MYO.thresholds)
pgs.MYO <- pgs.MYO[, 3:8]
rm(pgs.MYO.best, pgs.MYO.thresholds)
#write.xlsx(pgs.MYO, file = "results/07_PGS/Table_S11.xlsx", sheetName = "MYO", append = TRUE)

pgs.REF.thresholds <- read.table(file = "results/07_PGS/pgs.covar.ses/REF.var.prsice", header = TRUE) %>%
  filter(Threshold %in% thresholds)
pgs.REF.best <- read.table(file = "results/07_PGS/pgs.covar.ses/REF.var.prsice", header = TRUE) %>%
  filter(R2 == max(R2))
pgs.REF <- rbind(pgs.REF.best, pgs.REF.thresholds)
pgs.REF <- pgs.REF[, 3:8]
rm(pgs.REF.best, pgs.REF.thresholds)
#write.xlsx(pgs.REF, file = "results/07_PGS/Table_S11.xlsx", sheetName = "REF", append = TRUE)

pgs.r2.data <- cbind(pgs.ADHD$R2, pgs.ASD$R2, pgs.BIP$R2,  pgs.DYS$R2, pgs.SCZ$R2, pgs.EA$R2, pgs.INT$R2, pgs.MYO$R2, pgs.REF$R2)
rownames(pgs.r2.data) <- c("best ", "0.001 ", "0.05 ", "0.1 ", "0.2 ", "0.3 ", "0.4 ", "0.5 ", "1 ")
colnames(pgs.r2.data) <- c(" ADHD", " ASD", " BIP", " RD", " SCZ", " EA", " INT", " MYO", " REF")
pgs.r2.data <- round(pgs.r2.data*100, 2)

pgs.p.data <- cbind(pgs.ADHD$P, pgs.ASD$P, pgs.BIP$P, pgs.DYS$P, pgs.SCZ$P, pgs.EA$P, pgs.INT$P, pgs.MYO$P, pgs.REF$P)
rownames(pgs.p.data) <- c("best ", "0.001 ", "0.05 ", "0.1 ", "0.2 ", "0.3 ", "0.4 ", "0.5 ", "1 ")
colnames(pgs.p.data) <- c(" ADHD", " ASD", " BIP", " RD", " SCZ", " EA", " INT", " MYO", " REF")

tiff("outputs/plots/Vision_FigS6_PGS.tiff", units = "cm", width = 11, height = 12, res = 600)
corrplot(pgs.r2.data, 
         tl.col = "black",
         method = "color",
         addCoef.col = "white",
         addgrid.col = "dimgrey",
         is.corr = FALSE,
         cl.pos = "n",
         cl.lim = c(0,1),
         col = color.new,
         tl.srt = 45,
         p.mat = pgs.p.data,
         insig = "pch",
         sig.level = 0.05/(nrow(pgs.r2.data)*ncol(pgs.r2.data)),
         number.cex = 0.8, 
         tl.cex = 0.8,
         cl.cex = 0.8,
         mar = c(0, 1, 0, 3))
colorlegend(colbar = color.new, 
            labels = c(0,0.05,0.1,0.15,2,0.25,0.3), 
            xlim = c(0.5, 9.5), 
            ylim = c(-1, 0), 
            align = "l", 
            cex = 0.8,
            vertical = FALSE)
mtext(text = expression(bold(paste(italic("p "), "value threshold"))), side = 2, line = 3, at = 5, las = 0, cex = 0.8)
mtext(text = expression(bold("training GWAS")), side = 3, line = 2, at = 5.5, las = 1, cex = 0.8)
dev.off()

T2ADHD <- read.table(file = "results/07_PGS/pgs.covar.ses/ADHD.var.summary", header = TRUE)
T2ASD <- read.table(file = "results/07_PGS/pgs.covar.ses/ASD.var.summary", header = TRUE)
T2BIP <- read.table(file = "results/07_PGS/pgs.covar.ses/BIP.var.summary", header = TRUE)
T2DYS <- read.table(file = "results/07_PGS/pgs.covar.ses/DYS.var.summary", header = TRUE)
T2SCZ <- read.table(file = "results/07_PGS/pgs.covar.ses/SCZ.var.summary", header = TRUE)
T2EA <- read.table(file = "results/07_PGS/pgs.covar.ses/EA.var.summary", header = TRUE)
T2INT <- read.table(file = "results/07_PGS/pgs.covar.ses/INT.var.summary", header = TRUE)
T2MYO <- read.table(file = "results/07_PGS/pgs.covar.ses/MYO.var.summary", header = TRUE)
T2REF <- read.table(file = "results/07_PGS/pgs.covar.ses/REF.var.summary", header = TRUE)

T2 <- rbind(T2ADHD, T2ASD, T2BIP, T2DYS, T2SCZ, T2EA, T2INT, T2MYO, T2REF) %>%
  select(Threshold, Num_SNP, Coefficient, Standard.Error, PRS.R2, Full.R2, P)
#write.xlsx(T2, file = "results/07_PGS/Table_2.xlsx")

