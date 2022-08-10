
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

# Clean environment
rm(raw.alspac, alspac.gcse, alspac.clean)



##########################
#
# PART 2: CREATE PHENOTYPE FILE FOR PGS

# Read in data.sample file to identify individuals to be included in PGS analysis
gwas.id <- read.table('data/data.sample', header = TRUE) %>%
  filter(ID_1 != 0) %>%
  select("ID_1", "ID_2", "sex") %>%
  mutate(ID_1 = as.character(ID_1))

# Merge GWAS IDs with pheno.data to obtain pgs.data
pgs.data <- pheno.data %>% select(-sex)
pgs.data <- left_join(gwas.id, pgs.data, by = 'ID_1') %>%
  filter(!is.na(var)) 

pgs.data$var <- RankNorm(pgs.data$var)

sample.eye <- left_join(gwas.id[,1:3], pgs.data[,c(1,13)], by = 'ID_1') 

colnames(sample.eye) <- c("FID", "IID", "sex", "var")

#write.table(sample.eye, 'outputs/pgs/sample.eye', quote = F, row.names = F, col.names = T)



##########################
#
# PART 3: CREATE COVARIATE FILES FOR PGS

# Read in PCs
pcs <- read.table('data/eye.child.eigenvec') %>%
  select(V1, V3, V4) %>%
  rename(ID_1 = V1,
         pc1 = V3,
         pc2 = V4)

pgs.cov <- pgs.data %>%
  select(-c(ID_2,sex))

pgs.cov$GCSE[pgs.cov$GCSE == "-10"] <- NA
pgs.cov$GCSE[pgs.cov$GCSE == "-1"] <- NA

ea.pgs <- read.table('data/pgs/EDU.var.best', h = T)
ea.pgs <- ea.pgs %>%
  select(IID, PRS)
colnames(ea.pgs) <- c("ID_1", "EA_PRS")

int.pgs <- read.table('data/pgs/INT.var.best', h = T)
int.pgs <- int.pgs %>%
  select(IID, PRS)
colnames(int.pgs) <- c("ID_1", "INT_PRS")

adhd.pgs <- read.table('data/pgs/ADHD.var.best', h = T)
adhd.pgs <- adhd.pgs %>%
  select(IID, PRS)
colnames(adhd.pgs) <- c("ID_1", "ADHD_PRS")

dys.pgs <- read.table('data/pgs/DYS.var.best', h = T)
dys.pgs <- dys.pgs %>%
  select(IID, PRS)
colnames(dys.pgs) <- c("ID_1", "DYS_PRS")

cov.eye <- left_join(gwas.id, pgs.cov, by = 'ID_1') %>%
  select(ID_1, ID_2, sex, age.weeks.11, read.7, sum.ccc, wold.comp, nonword.rep, wisc.total.iq, wisc.verbal.iq, wisc.performance.iq, GCSE, mother.education)
cov.eye <- left_join(cov.eye, pcs, by = 'ID_1') %>%
  select(ID_1, ID_2, sex, age.weeks.11, pc1, pc2, read.7, sum.ccc, wold.comp, nonword.rep, wisc.total.iq, wisc.verbal.iq, wisc.performance.iq, GCSE, mother.education)
cov.eye <- left_join(cov.eye, ea.pgs, by = 'ID_1')
cov.eye <- left_join(cov.eye, int.pgs, by = 'ID_1') 
cov.eye <- left_join(cov.eye, adhd.pgs, by = 'ID_1') 
cov.eye <- left_join(cov.eye, dys.pgs, by = 'ID_1') 

colnames(cov.eye) <- c("FID", "IID", "sex", "age", "pc1", "pc2", "reading.ability", "communication.skills", "listening.comprehension", "short.term.memory", "total.iq", "verbal.iq", "performance.iq", "GCSE", "SES", "EA_PRS", "INT_PRS", "ADHD_PRS", "DYS_PRS")

#write.table(cov.eye, 'outputs/pgs/cov.eye', quote = F, row.names = F, col.names=T)

