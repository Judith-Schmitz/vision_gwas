#!/bin/bash

# This script includes sex and the first 2 pcs as covariates

pheno=var
phenofile=eye
cohort=child

cd results.gwas/

for i in {1..22}

do

cat > $pheno.chr${i}.sh <<EOF1
#!/bin/bash

#SBATCH --job-name="gwas"
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

~/apps/conda/envs/geno_utils/bin/bolt \
--bfile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_data/new_alspac_external_500G_4/genetics/ARRAY_2019-05-13/all1/data/data \
--bgenFile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_data/new_alspac_external_500G_4/genetics/HRC_2019-05-13/all1/data/bgen_8bit/data_8bit_$i.bgen \
--sampleFile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_data/new_alspac_external_500G_4/genetics/HRC_2019-05-13/all1/data/data.sample \
--bgenMinMAF=0.05 \
--bgenMinINFO=0.8 \
--phenoFile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_resources/bolt.lmm/phenos/${cohort}.${phenofile}.pheno \
--phenoCol=$pheno \
--covarFile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_resources/bolt.lmm/phenos/${cohort}.${phenofile}.pheno \
--covarCol=SEX \
--qCovarCol=pc1 \
--qCovarCol=pc2 \
--qCovarCol=age \
--lmm \
--LDscoresFile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_resources/bolt.lmm/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_resources/bolt.lmm/tables/genetic_map_hg19_withX.txt.gz \
--numThreads=8 \
--statsFile=$pheno.chr${i}.stats.gz \
--statsFileBgenSnps=$pheno.chr${i}.bgen.stats.gz \
--verboseStats

EOF1

chmod +x $pheno.chr${i}.sh

sbatch $pheno.chr${i}.sh

done

cat > ${pheno}.combine.results.sh <<EOF2
#!/bin/bash

#SBATCH --job-name="combine.results"
#SBATCH --mem=16G

echo "Opening files"

for i in {1..22}
do
  gunzip ${pheno}.chr\${i}.bgen.stats.gz
	awk 'NR > 1' ${pheno}.chr\${i}.bgen.stats > ${pheno}.chr\${i}.temp
  gzip ${pheno}.chr\${i}.bgen.stats
done

echo "Concatenating files"
cat ${pheno}.chr1.temp ${pheno}.chr2.temp ${pheno}.chr3.temp ${pheno}.chr4.temp ${pheno}.chr5.temp ${pheno}.chr6.temp ${pheno}.chr7.temp ${pheno}.chr8.temp ${pheno}.chr9.temp ${pheno}.chr10.temp ${pheno}.chr11.temp ${pheno}.chr12.temp  ${pheno}.chr13.temp ${pheno}.chr14.temp ${pheno}.chr15.temp ${pheno}.chr16.temp ${pheno}.chr17.temp ${pheno}.chr18.temp ${pheno}.chr19.temp ${pheno}.chr20.temp ${pheno}.chr21.temp ${pheno}.chr22.temp > sumstats.${pheno}

rm  *.temp

# Make file for Manhattan plot
# 1: SNP
# 2: CHROM
# 3: BP
# 14: P_BOLT_LMM_INF

echo "Preparing file for Manhattan plot"
cut -d$'\t' -f 1,2,3,14 sumstats.${pheno} > manhattan.${pheno}
sed -i $'1 i\\\ID\tCHR\tBP\tP'  manhattan.${pheno}

# Header for sumstats
sed -i $'1 i\\\ID\tCHR\tBP\tGENPOS\tALLELE1\tALLELE0\tA1FREQ\tINFO\tCHISQ_LINREG\tP_LINREG\tBETA\tSE\tCHISQ_BOLT_LMM_INF\tP_BOLT_LMM_INF' sumstats.${pheno}

# Update SNP IDs
echo "Recode SNP IDs"
Rscript ../chrpos2rsid.bolt.R sumstats.${pheno}

echo "Finished"

EOF2

chmod +x ${pheno}.combine.results.sh

cd ../


