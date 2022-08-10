#!/bin/bash

pheno=var
phenofile=eye

cd results.gwas/

# Make file for FINDOR

# 2: CHR
# 3: BP
# 5: ALLELE1
# 6: ALLELE2
# 7: A1FREQ
# 11: BETA
# 12: SE
# 14: P_BOLT_LMM_INF
# 15: SNP

# N: 5571
# Z: BETA/SE

echo "Preparing file for FINDOR"
cut -d$'\t' -f 2,3,5,6,7,11,12,14,15 sumstats.${pheno} > findor.sumstats.${pheno}1.temp

# Calculate Z
awk '$0=$0"\t"(NR==1?"Z":$6/$7)' findor.sumstats.${pheno}1.temp > findor.sumstats.${pheno}2.temp

cut -d$'\t' -f 1,2,3,4,5,8,9,10 findor.sumstats.${pheno}2.temp > findor.sumstats.${pheno}

# Add N column
sed -i "s/$/\t5771/" findor.sumstats.${pheno}

# Header for sumstats
sed -i s/ALLELE1/Allele1/ findor.sumstats.${pheno}
sed -i s/ALLELE0/Allele2/ findor.sumstats.${pheno}
sed -i s/P_BOLT_LMM_INF/P/ findor.sumstats.${pheno}
#sed -i s/5571/N/ findor.sumstats.${pheno}

rm  *.temp

cd ../