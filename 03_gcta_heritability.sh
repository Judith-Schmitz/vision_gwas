#!/bin/bash

pheno=var
phenofile=eye

cd results.her/

cat > $phenofile.heritability.sh <<EOF1
#!/bin/bash

#SBATCH --job-name="heritability"
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

mkdir data/
mkdir her/
mkdir bivar/

# Making a GRM from all SNPs
~/apps/conda/envs/geno_utils/bin/gcta64 \
--bfile /mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_data/new_alspac_external_500G_4/genetics/ARRAY_2019-05-13/all1/data/data \
--keep ../data/keep.$phenofile \
--make-grm \
--thread-num 8 \
--out data/$phenofile

# Apply GRM cutoff
~/apps/conda/envs/geno_utils/bin/gcta64 \
--grm data/$phenofile \
--grm-cutoff 0.05 \
--make-grm \
--thread-num 8 \
--out data/${phenofile}_rel

# Running a REML analysis
~/apps/conda/envs/geno_utils/bin/gcta64 \
--reml \
--grm data/${phenofile}_rel \
--pheno ../data/gcta.$phenofile.pheno \
--covar ../data/gcta.$phenofile.cognition.covar \
--qcovar ../data/gcta.$phenofile.cognition.qcovar \
--thread-num 8 \
--out her/${pheno}

# Running bivariate REML
~/apps/conda/envs/geno_utils/bin/gcta64 \
--reml-bivar 1 2 \
--reml-bivar-lrt-rg 0 \
--grm data/${phenofile}_rel \
--pheno ../data/gcta.${phenofile}.cognition.pheno \
--covar ../data/gcta.$phenofile.cognition.covar \
--qcovar ../data/gcta.$phenofile.cognition.qcovar \
--thread-num 8 \
--out bivar/${pheno}_reading_ability

# Running bivariate REML
~/apps/conda/envs/geno_utils/bin/gcta64 \
--reml-bivar 1 3 \
--reml-bivar-lrt-rg 0 \
--grm data/${phenofile}_rel \
--pheno ../data/gcta.${phenofile}.cognition.pheno \
--covar ../data/gcta.$phenofile.cognition.covar \
--qcovar ../data/gcta.$phenofile.cognition.qcovar \
--thread-num 8 \
--out bivar/${pheno}_communication_skills

# Running bivariate REML
~/apps/conda/envs/geno_utils/bin/gcta64 \
--reml-bivar 1 4 \
--reml-bivar-lrt-rg 0 \
--grm data/${phenofile}_rel \
--pheno ../data/gcta.${phenofile}.cognition.pheno \
--covar ../data/gcta.$phenofile.cognition.covar \
--qcovar ../data/gcta.$phenofile.cognition.qcovar \
--thread-num 8 \
--out bivar/${pheno}_listening_comprehension

# Running bivariate REML
~/apps/conda/envs/geno_utils/bin/gcta64 \
--reml-bivar 1 5 \
--reml-bivar-lrt-rg 0 \
--grm data/${phenofile}_rel \
--pheno ../data/gcta.${phenofile}.cognition.pheno \
--covar ../data/gcta.$phenofile.cognition.covar \
--qcovar ../data/gcta.$phenofile.cognition.qcovar \
--thread-num 8 \
--out bivar/${pheno}_short_term_memory

# Running bivariate REML
~/apps/conda/envs/geno_utils/bin/gcta64 \
--reml-bivar 1 6 \
--reml-bivar-lrt-rg 0 \
--grm data/${phenofile}_rel \
--pheno ../data/gcta.${phenofile}.cognition.pheno \
--covar ../data/gcta.$phenofile.cognition.covar \
--qcovar ../data/gcta.$phenofile.cognition.qcovar \
--thread-num 8 \
--out bivar/${pheno}_total_IQ

# Running bivariate REML
~/apps/conda/envs/geno_utils/bin/gcta64 \
--reml-bivar 1 7 \
--reml-bivar-lrt-rg 0 \
--grm data/${phenofile}_rel \
--pheno ../data/gcta.${phenofile}.cognition.pheno \
--covar ../data/gcta.$phenofile.cognition.covar \
--qcovar ../data/gcta.$phenofile.cognition.qcovar \
--thread-num 8 \
--out bivar/${pheno}_verbal_IQ

# Running bivariate REML
~/apps/conda/envs/geno_utils/bin/gcta64 \
--reml-bivar 1 8 \
--reml-bivar-lrt-rg 0 \
--grm data/${phenofile}_rel \
--pheno ../data/gcta.${phenofile}.cognition.pheno \
--covar ../data/gcta.$phenofile.cognition.covar \
--qcovar ../data/gcta.$phenofile.cognition.qcovar \
--thread-num 8 \
--out bivar/${pheno}_performance_IQ

# Running bivariate REML
~/apps/conda/envs/geno_utils/bin/gcta64 \
--reml-bivar 1 9 \
--reml-bivar-lrt-rg 0 \
--grm data/${phenofile}_rel \
--pheno ../data/gcta.${phenofile}.cognition.pheno \
--covar ../data/gcta.$phenofile.cognition.covar \
--qcovar ../data/gcta.$phenofile.cognition.qcovar \
--thread-num 8 \
--out bivar/${pheno}_GCSE

# Running bivariate REML
~/apps/conda/envs/geno_utils/bin/gcta64 \
--reml-bivar 1 2 \
--reml-bivar-lrt-rg 0 \
--grm data/${phenofile}_rel \
--pheno ../data/gcta.${phenofile}.ear.pheno \
--covar ../data/gcta.$phenofile.ear.covar \
--qcovar ../data/gcta.$phenofile.ear.qcovar \
--thread-num 8 \
--out bivar/${pheno}_ear

EOF1

chmod +x $phenofile.heritability.sh

sbatch $phenofile.heritability.sh

cd ../


