
source(file = 'scripts/01_behavioural_analysis_functions.r')

##########################
#
# PART 1: READ ALSPAC

# Read SPSS file
raw.alspac <- read.alspac()
alspac.gcse <- read.alspac.newfile()
raw.alspac <- left_join(raw.alspac, alspac.gcse, by = 'ID_1')

# Apply exclusion criteria (sensory disabilities, full phenotype data)
alspac.clean <- clean.alspac(mydata = raw.alspac)

# age-correct reading and wold scores
alspac.clean$read.7[!is.na(alspac.clean$read.7)] <- scale(resid(glm(read.7 ~ age.weeks.7, data = alspac.clean)))[,]
alspac.clean$nonword.rep[!is.na(alspac.clean$nonword.rep)] <- scale(resid(glm(nonword.rep ~ age.weeks.8, data = alspac.clean)))[,]
alspac.clean$wold.comp[!is.na(alspac.clean$wold.comp)] <- scale(resid(glm(wold.comp ~ age.weeks.8, data = alspac.clean)))[,]
alspac.clean$sum.ccc[!is.na(alspac.clean$sum.ccc)] <- scale(resid(glm(sum.ccc ~ age.weeks.9, data = alspac.clean)))[,]

# Normalize data
pheno.data <- normalize.alspac(dataset = alspac.clean)

# Read in data.sample file to identify individuals to be included in GWAS
gwas.id <- read.table('data/data.sample', header = TRUE) %>%
  filter(ID_1 != 0) %>%
  select("ID_1", "ID_2", "sex") %>%
  mutate(ID_1 = as.character(ID_1))

# Merge GWAS IDs with pheno.data to obtain gwas.data
gwas.data <- pheno.data %>% select(-sex)
gwas.data <- left_join(gwas.id, gwas.data, by = 'ID_1') %>%
  filter(!is.na(var))

# Unrelated individuals
unrelated <- read.table("data/eye_rel.grm.id", header = FALSE)
unrel.data <- gwas.data %>%
  filter(ID_1 %in% unrelated$V1) 

# Clean environment
rm(raw.alspac, alspac.gcse, alspac.clean, gwas.id)


##########################
#
# PART 2: DESCRIPTIVE STATISTICS

# PHENO.DATA (N = 6807)
# GWAS.DATA (N = 5571)
# UNREL.DATA (N = 5153)

data <- pheno.data # pheno.data or gwas.data or unrel.data

## sex and age descriptives
freq(data$sex)

round(mean(data$age.weeks.11)/52,2)
round(sd(data$age.weeks.11)/52,2)

## VA descriptives
describe(data$var)

## sex and age effects on VA
t.test(var ~ sex, data = data)

group_by(subset(data, !is.na(sex)), sex) %>%
  summarise(
    count = n(),
    mean = round(mean(var),2),
    sd = round(sd(var),2)
  )

summary(lm(var ~ age.weeks.11, data = data))


##########################
#
# PART 3: BEHAVIOURAL ANALYSIS

## VA and cognitive traits
parcor.cog <- data %>% 
  select(sex, age.weeks.11, var,
         read.7, sum.ccc, wold.comp, nonword.rep, wisc.total.iq, wisc.verbal.iq, wisc.performance.iq, GCSE) %>%
  mutate(sex = as.numeric(sex))

# partial correlation adjusting for sex and age
parcor.cog.r <- matrix(partial.r(parcor.cog, c(3:ncol(parcor.cog)), c(1:2)), ncol = ncol(parcor.cog) - 2)
rownames(parcor.cog.r) <- c("VAR ", "reading skill ",  "communication skills ", "listening comprehension ", "short term memory ", "total IQ ", "verbal IQ ", "performance IQ ", "GCSE ")
colnames(parcor.cog.r) <- c(" VAR", " reading skill",  " communication skills", " listening comprehension", " short term memory", " total IQ", " verbal IQ", " performance IQ", " GCSE")

parcor.cog.p <- corr.p(parcor.cog.r, adjust = "none", n = freq(!is.na(parcor.cog$GCSE))[2, 1] - 2)$p


## VA and hearing
parcor.hear <- data %>%
  filter(hear.imp == 1 ) %>%
  filter(hear.loss == 1 ) %>%
  filter(hf.hear.loss == 1 ) %>%
  filter(!is.na(best.ear)) %>%
  mutate(var = as.numeric(var)) %>%
  mutate(best.ear = as.numeric(best.ear)) %>%
  mutate(sex = as.numeric(sex))

print(nrow(parcor.hear))

freq(parcor.hear$sex)

# partial correlation adjusting for sex and age
parcor.hear.r <- matrix(partial.r(parcor.hear, c(10, 12), c(4, 5)), ncol = 2)
parcor.hear.r
rownames(parcor.hear.r) <- c("VAR ", "HT ")
colnames(parcor.hear.r) <- c(" VAR", " HT")

corr.p(parcor.hear.r, adjust = "none", n = ncol(parcor.hear))$p


## VA and cognitive traits, correcting for SES, ...

parcor.ses <- data %>% 
  select(sex, age.weeks.11, mother.education, read.7, wisc.verbal.iq, wisc.performance.iq, var, GCSE) %>%
  mutate(mother.education = as.numeric(mother.education)) %>%
  mutate(sex = as.numeric(sex))

# partial correlation adjusting for sex, age, SES, reading and IQ
parcor.ses.r <- matrix(partial.r(parcor.ses, c(7:ncol(parcor.ses)), c(1:6)), ncol = ncol(parcor.ses) - 6)
parcor.ses.r
rownames(parcor.ses.r) <- c("VAR ", "GCSE ")
colnames(parcor.ses.r) <- c(" VAR", " GCSE")

corr.p(parcor.ses.r, adjust = "none", n = freq(!is.na(parcor.ses$GCSE))[2, 1] - 6)$p

# Clean environment
rm(parcor.hear, parcor.hear.r, parcor.ses, parcor.ses.r)


## VA and SES
data$mother.education <- factor(data$mother.education, 
                                levels = c('1', '2', '3', '4', '5'),
                                labels = c('CSE', 'Vocational', 'O level', 'A level', 'Degree'))

print(nrow(subset(data, !is.na(mother.education))))

# Unadjusted ANOVA
ses.aov <- aov(var ~ mother.education, data = subset(data, !is.na(mother.education)))
summary(ses.aov)

# ANOVA adjusted for sex and age
ses.aov.sex.age <- aov(var ~ mother.education + sex + age.weeks.11, data = subset(data, !is.na(mother.education)))
summary(ses.aov.sex.age)

# Posthoc tests unadjusted ANOVA
TukeyHSD(ses.aov)

group_by(subset(data, !is.na(mother.education)), mother.education) %>%
  summarise(
    count = n(),
    mean = sprintf("%0.2f", mean(var)),
    sd = round(sd(var, na.rm = TRUE),2)
  )

# Clean environment
rm(ses.aov, ses.aov.sex.age)


## VA and neurodevelopmental subgroups

id.link <- read.xlsx(file = "data/scerri_subgroups.xlsx", sheetName = "Subgroups")
id.link[id.link == "NV"] <- NA
colnames(id.link) <- c("ID_1", "Status")

# Control 0
# SLI 1
# RD	2
# SLI + RD	3
# ADHD	4
# SLI + ADHD	5
# RD + ADHD	6
# RD + SLI + ADHD	7

data.neurodev <- left_join(data, id.link, by = 'ID_1') %>%
  filter(!is.na(Status) | asd == 2) %>%
  filter(!is.na(sex)) %>%
  mutate(group = ifelse(test = Status == 3 | Status == 6 | Status == 7, yes = 2, no = Status)) %>% # comorbid with rd --> rd group
  mutate(group = ifelse(test = asd == 2 & Status == 0, yes = 3, no = group)) %>%
  mutate(group = ifelse(test = is.na(asd) & Status == 0, yes = Status, no = group)) %>% 
  mutate(group = ifelse(test = asd == 2 & is.na(Status), yes = 3, no = group)) %>%
  filter(group != 5) # excluding the 4 cases with sli + adhd

data.cases <- data.neurodev %>%
  filter(group != 0)

data.controls <- data.neurodev %>%
  filter(group == 0)

# Control	0
# SLI	1
# RD (including comorbid cases)	2
# ASD 3
# ADHD 4

# M:F ratio ####
# Gender ratio between DLD cases and RD cases shold be very similar
gender_ratio_all_cases <- round(table(data.cases %>% select(sex))[[2]] / table(data.cases %>% select(sex))[[3]], 2)
print(paste0('M:F in cases: ',gender_ratio_all_cases))

gender_controls <- round(table(data.controls %>% select(sex))[[2]] / table(data.controls %>% select(sex))[[3]], 2)
print(paste0('M:F in controls: ',gender_controls))

set.seed(10000)

# Sex-match all controls ####
if (gender_controls < gender_ratio_all_cases){
  print('M:F ratio higher in cases, removing random females from controls')
  while (gender_controls < gender_ratio_all_cases) {
    data.controls <- data.controls[-sample(which(data.controls$sex == 2 & data.controls$group == 0), 1), ]
    print(paste0('Total sample: ', nrow(data.controls)))
    gender_controls <- round(table(data.controls %>% filter(group==0) %>% select(sex))[[2]] / table(data.controls %>% filter(group==0) %>% select(sex))[[3]], 2)
    print(paste0('M:F in controls: ',gender_controls))
  }
} else {
  print('M:F ratio higher in controls, removing random males from controls')
  while (gender_controls > gender_ratio_all_cases) {
    data.controls <- data.controls[-sample(which(data.controls$sex == 1 & data.controls$group == 0), 1), ]
    print(paste0('Total sample: ', nrow(data.controls)))
    gender_controls <- round(table(data.controls %>% filter(group==0) %>% select(sex))[[2]] / table(data.controls %>% filter(group==0) %>% select(sex))[[3]], 2)
    print(paste0('M:F in controls: ',gender_controls))
  }
}

data.neurodev <- rbind(data.cases, data.controls)

data.neurodev$group <- factor(data.neurodev$group,
                              levels = c('0', '1', '2', '3', '4'),
                              labels = c('control', 'LI', 'RD', 'ASD', 'ADHD'))

# Unadjusted ANOVA
neurodev.aov <- aov(var ~ group, data = data.neurodev)
summary(neurodev.aov)

# ANOVA adjusted for sex and age
neurodev.aov.sex.age <- aov(var ~ group + sex + age.weeks.11, data = data.neurodev)
summary(neurodev.aov.sex.age)

# Posthoc tests unadjusted ANOVA
TukeyHSD(neurodev.aov)

group_by(data.neurodev, group) %>%
  summarise(
    count = n(),
    mean = sprintf("%0.2f", mean(var)),
    sd = round(sd(var, na.rm = TRUE),2)
  )


## Filtering outliers
mean.neurodev <- mean(data.neurodev$var)
sd.neurodev <- sd(data.neurodev$var)
lower.neurodev <- mean.neurodev - 3*sd.neurodev
upper.neurodev <- mean.neurodev + 3*sd.neurodev

filtered.neurodev <- data.neurodev %>%
  filter(var <= upper.neurodev) %>%
  filter(var >= lower.neurodev)

# Who is excluded?
group_by(data.neurodev, group) %>%
  filter(var > upper.neurodev | var < lower.neurodev) %>%
  summarise(
    count = n()
  )

# Unadjusted ANOVA
neurodev.aov.filter <- aov(var ~ group, data = filtered.neurodev)
summary(neurodev.aov.filter)

# ANOVA adjusted for sex and age
neurodev.aov.filter.sex.age <- aov(var ~ group + sex + age.weeks.11, data = filtered.neurodev)
summary(neurodev.aov.filter.sex.age)

# Posthoc tests unadjusted ANOVA
TukeyHSD(neurodev.aov.filter)

group_by(filtered.neurodev, group) %>% 
  summarise(
    count = n(),
    mean = sprintf("%0.2f", mean(var)),
    sd = round(sd(var, na.rm = TRUE),2)
  )

# Clean environment
rm(data.cases, data.controls, mean.neurodev, sd.neurodev, lower.neurodev, upper.neurodev, id.link)
rm(neurodev.aov, neurodev.aov.sex.age, neurodev.aov.filter, neurodev.aov.filter.sex.age)
rm(gender_controls, gender_ratio_all_cases)


##########################
#
# PART 4: FIGURES

## FIGURE 2 (FOR FIGURE S1: data <- gwas.data)

# plot densities
data.sex <- data %>% filter(!is.na(sex))
data.sex$sex <- factor(data.sex$sex,
                       levels = c('1', '2'),
                       labels = c('male', 'female'))

# var
Fig2a <- ggplot(aes(x = var), data = data) +
  geom_density(fill = "#7F0000", alpha = 0.8) +
  scale_y_continuous(limits = c(0,0.12), breaks = seq(0, 0.12, by = 0.02)) + 
  scale_x_continuous(limits = c(80,120), breaks = seq(80, 120, by = 10)) + 
  labs(y = "density", x = "VA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"))

Fig2b <- ggplot(aes(x = var, group = sex, fill = sex), data = data.sex) +
  geom_density(alpha = 0.8) +
  scale_y_continuous(limits = c(0,0.12), breaks = seq(0, 0.12, by = 0.02)) + 
  scale_x_continuous(limits = c(80,120), breaks = seq(80, 120, by = 10)) + 
  labs(y = "density", x = "VA") +
  scale_fill_manual(values = c("#EF6548", "#FEE8C8")) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12),
        legend.position = c(.3, .5),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1))

tiff("outputs/plots/Vision_Fig2_pheno_dist.tiff", units = "in", width = 9, height = 3, res = 300)
ggarrange(Fig2a, Fig2b,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()


## FIGURE 3

fig3a.r <- parcor.cog.r[2:nrow(parcor.cog.r),1]
fig3a.p <- parcor.cog.p[2:nrow(parcor.cog.p),1]

fig3b.data <- as.matrix(read.xlsx(file = 'results/03_SNPh2/GCTA_genetic_correlation.xlsx', sheetName = 'figure', row.names = T, header = T))
fig3b.r <- as.matrix(fig2b.data[1,])
fig3b.p <- as.matrix(fig2b.data[2,])

fig3.r <- cbind(fig2b.r, fig2a.r)
fig3.p <- cbind(fig2b.p, fig2a.p)

colnames(fig3.r) <- c("B ", "A ")
rownames(fig3.r) <- c(" reading skill",  " communication skill", " listening comprehension", " short term memory", " total IQ", " verbal IQ", " performance IQ", " GCSE")
colnames(fig3.p) <- c("B ", "A ")
rownames(fig3.p) <- c(" reading skill",  " communication skill", " listening comprehension", " short term memory", " total IQ", " verbal IQ", " performance IQ", " GCSE")

tiff("outputs/plots/Vision_Fig3_correlation_cognition.tiff", units = "cm", width = 16, height = 10, res = 600)
ggcorrplot(fig3.r,
           lab = TRUE,
           p.mat = fig3.p,
           sig.level = 0.05/8,
           insig = "blank",
           colors = c("white", "white", "white"),
           outline.color = "black",
           show.legend = FALSE)
dev.off()


## FIGURE S2

options(scipen = 1)

FigS2a <- ggplot(aes(x = var, y = read.7), data = subset(parcor.cog, !is.na(read.7))) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(y = "reading skill", x = "VA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  stat_cor(method = "pearson", label.x = 90, label.y = 3)

FigS2b <- ggplot(aes(x = var, y = sum.ccc), data = subset(parcor.cog, !is.na(sum.ccc))) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(y = "communication skills", x = "VA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  stat_cor(method = "pearson", label.x = 90, label.y = 2)

FigS2c <- ggplot(aes(x = var, y = wold.comp), data = subset(parcor.cog, !is.na(wold.comp))) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(y = "listening comprehension", x = "VA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  stat_cor(method = "pearson", label.x = 90, label.y = 4.2)

FigS2d <- ggplot(aes(x = var, y = nonword.rep), data = subset(parcor.cog, !is.na(nonword.rep))) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(y = "short term memory", x = "VA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  stat_cor(method = "pearson", label.x = 90, label.y = 2.2)

FigS2e <- ggplot(aes(x = var, y = wisc.total.iq), data = subset(parcor.cog, !is.na(wisc.verbal.iq))) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(y = "total IQ", x = "VA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  stat_cor(method = "pearson", label.x = 90, label.y = 165)

FigS2f <- ggplot(aes(x = var, y = wisc.verbal.iq), data = subset(parcor.cog, !is.na(wisc.verbal.iq))) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(y = "verbal IQ", x = "VA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  stat_cor(method = "pearson", label.x = 90, label.y = 165)

FigS2g <- ggplot(aes(x = var, y = wisc.performance.iq), data = subset(parcor.cog, !is.na(wisc.performance.iq))) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(y = "performance IQ", x = "VA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  stat_cor(method = "pearson", label.x = 90, label.y = 160)

FigS2h <- ggplot(aes(x = var, y = GCSE), data = subset(parcor.cog, !is.na(GCSE))) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(y = "GCSE score", x = "VA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  stat_cor(method = "pearson", label.x = 90, label.y = 600)

tiff("outputs/plots/Vision_FigS2_cognitive_measures.tiff", units="in", width=8, height=12, res=300)
ggarrange(FigS2a, FigS2b, FigS2c, FigS2d, FigS2e, FigS2f, FigS2g, FigS2h,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 2, nrow = 4)
dev.off()


## FIGURE 4

options(scipen = 999)

color.new <- brewer.pal(n = 5, name = 'OrRd')

Fig4a <- ggplot(subset(data, !is.na(mother.education)), 
                aes(x = mother.education, y = var)) +
  geom_violin(aes(fill = mother.education), trim = FALSE) +  
  scale_y_continuous(limits = c(80,140), breaks = round(seq(80, 140, by = 20),2)) + 
  scale_fill_manual(values = color.new) + 
  stat_summary(
    fun.data = "mean_sdl",  fun.args = list(mult = 1), 
    geom = "pointrange", color = "black") + 
  labs(x = "SES (mother's highest education)", y = "VA") + 
  theme_bw() + 
  theme(legend.position = "none") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  geom_signif(comparisons = list(c("CSE", "O level")), vjust = 0.5, map_signif_level = TRUE, y_position = 138, annotations = "*") + 
  geom_signif(comparisons = list(c("CSE", "A level")), vjust = 0.5, map_signif_level = TRUE, y_position = 136, annotations = "***") + 
  geom_signif(comparisons = list(c("CSE", "Degree")), vjust = 0.5, map_signif_level = TRUE, y_position = 134, annotations = "***") + 
  geom_signif(comparisons = list(c("Vocational", "A level")), vjust = 0.5, map_signif_level = TRUE, y_position = 132, annotations = "**") + 
  geom_signif(comparisons = list(c("Vocational", "Degree")), vjust = 0.5, map_signif_level = TRUE, y_position = 130, annotations = "***") + 
  geom_signif(comparisons = list(c("O level", "A level")), vjust = 0.5, map_signif_level = TRUE, y_position = 128, annotations = "*") + 
  geom_signif(comparisons = list(c("O level", "Degree")), vjust = 0.5, map_signif_level = TRUE, y_position = 126, annotations = "***") + 
  geom_signif(comparisons = list(c("A level", "Degree")), vjust = 0.5, map_signif_level = TRUE, y_position = 124, annotations = "*")

Fig4b <- ggplot(data.neurodev,
                aes(x = group, y = var)) +
  geom_violin(aes(fill = group), trim = FALSE) +  
  scale_fill_manual(values = color.new) + 
  scale_y_continuous(limits = c(80,140), breaks = round(seq(80, 140, by = 20),2)) + 
  stat_summary(
    fun.data = "mean_sdl",  fun.args = list(mult = 1), 
    geom = "pointrange", color = "black") +
  labs(x = "neurodevelopmental condition", y = "VA") + 
  theme_bw() + 
  theme(legend.position = "none") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  geom_signif(comparisons = list(c("control", "RD")), vjust = 0.5, map_signif_level = TRUE, y_position = 124, annotations = "***")

tiff("outputs/plots/Vision_Fig4_SES_neurodev.tiff", units = "in", width = 10, height = 5, res = 300)
ggarrange(Fig4a, Fig4b,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()



