
source(file = 'scripts/02_prepare_bolt_gcta_functions.r')

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

# transform phenotypes
pheno.data <- normalize.alspac(dataset = alspac.clean)

# Clean environment
rm(raw.alspac, alspac.gcse, alspac.clean)


##########################
#
# PART 2: CREATE PHENO FILE FOR BOLT.LMM GWAS

# Read in data.sample file to identify individuals to be included in GWAS
gwas.id <- read.table('data/data.sample', header = TRUE) %>%
  filter(ID_1 != 0) %>%
  select("ID_1", "ID_2", "sex") %>%
  mutate(ID_1 = as.character(ID_1))

# Merge GWAS IDs with pheno.data
gwas.data <- left_join(gwas.id, pheno.data, by = 'ID_1') %>%
  filter(!is.na(var))

# create keep file
keep.eye <- gwas.data %>%
  select(ID_1, ID_2)

#write.table(keep.eye, 'outputs/keep files/keep.eye', sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Normalise VA for GWAS
gwas.data$var <- RankNorm(gwas.data$var)

# Read in PCs
pcs <- read.table('data/eye.child.eigenvec') %>%
  select(V1, V3, V4) %>%
  rename(ID_1 = V1,
         pc1 = V3,
         pc2 = V4)

bolt.pheno <- left_join(gwas.data, pcs, by = 'ID_1') %>%
  select(-hear.imp, -hear.loss, -hf.hear.loss) %>%
  select(ID_1, ID_2, sex, pc1, pc2, age.weeks.11, var, mother.education, read.7, sum.ccc, wold.comp, nonword.rep, wisc.total.iq, wisc.verbal.iq, wisc.performance.iq, GCSE) %>%
  rename(FID = ID_1,
         IID = ID_2,
         SEX = sex,
         read7 = read.7,
         SES = mother.education,
         age = age.weeks.11)

#write.table(bolt.pheno, 'outputs/bolt.files/child.eye.pheno', sep = "\t", quote = F, row.names = F)

bolt.pheno.ear <- left_join(gwas.data, pcs, by = 'ID_1') %>%
  filter(!is.na(best.ear)) %>%
  filter(hear.imp == 1) %>%
  filter(hear.loss == 1) %>%
  filter(hf.hear.loss == 1) %>%
  select(ID_1, ID_2, sex, pc1, pc2, age.weeks.11, var, best.ear) %>%
  rename(FID = ID_1,
         IID = ID_2,
         SEX = sex,
         age = age.weeks.11)

bolt.pheno.ear$best.ear <- RankNorm(bolt.pheno.ear$best.ear)

#write.table(bolt.pheno.ear, 'outputs/bolt.files/child.eye.ear.pheno', sep = "\t", quote = F, row.names = F)

# Clean environment
rm(gwas.id, pcs)



##########################
#
# PART 3: CREATE PHENO FILE FOR GCTA-GREML

unrelated <- read.table("data/eye_rel.grm.id", header = FALSE)

gcta.unrelated <- bolt.pheno %>%
  filter(FID %in% unrelated$V1)

gcta.ear.unrelated <- bolt.pheno.ear %>%
  filter(FID %in% unrelated$V1) 
  
gcta.var.pheno <- gcta.unrelated %>%
  select(FID, IID, var)
write.table(gcta.var.pheno, 'outputs/gcta/gcta.eye.pheno', sep = "\t", quote = F, row.names = F, col.names = F)

gcta.var.cognition.pheno <- gcta.unrelated %>%
  select(FID, IID, var, read7, sum.ccc, wold.comp, nonword.rep, wisc.total.iq, wisc.verbal.iq, wisc.performance.iq, GCSE)
write.table(gcta.var.cognition.pheno, 'outputs/gcta/gcta.eye.cognition.pheno', sep = "\t", quote = F, row.names = F, col.names = F)

gcta.var.covar <- gcta.unrelated %>%
  select(FID, IID, SEX)
write.table(gcta.var.covar, 'outputs/gcta/gcta.eye.cognition.covar', sep = "\t", quote = F, row.names = F, col.names = F)

gcta.var.qcovar <- gcta.unrelated %>%
  select(FID, IID, pc1, pc2, age)
write.table(gcta.var.qcovar, 'outputs/gcta/gcta.eye.cognition.qcovar', sep = "\t", quote = F, row.names = F, col.names = F)

gcta.ear.pheno <- gcta.ear.unrelated %>%
  select(FID, IID, var, best.ear)
write.table(gcta.ear.pheno, 'outputs/gcta/gcta.eye.ear.pheno', sep = "\t", quote = F, row.names = F, col.names = F)

gcta.ear.covar <- gcta.ear.unrelated %>%
  select(FID, IID, SEX)
write.table(gcta.ear.covar, 'outputs/gcta/gcta.eye.ear.covar', sep = "\t", quote = F, row.names = F, col.names = F)

gcta.ear.qcovar <- gcta.ear.unrelated %>%
  select(FID, IID, pc1, pc2, age)
write.table(gcta.ear.qcovar, 'outputs/gcta/gcta.eye.ear.qcovar', sep = "\t", quote = F, row.names = F, col.names = F)


