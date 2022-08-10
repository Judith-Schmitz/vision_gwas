#! /bin/bash

pheno=var
phenofile=eye
binary=F
stat=BETA

# OR = ADHD ASD
# BETA = BIP EA INT SCZ REF

for base in REF

do

cd results.pgs/

cat > pgs.${base}.${pheno}.sh <<EOF1
#!/bin/bash

#SBATCH --job-name="PGS"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

Rscript ~/apps/PRSice/PRSice.R --dir . \
    --prsice ~/apps/PRSice/PRSice_linux \
    --base /mnt/shared/scratch/jschmitz/private/vision_gwas/sumstats.new/${base}.chrpos.gz \
    --extract /mnt/shared/scratch/jschmitz/private/vision_gwas/sumstats.new/${base}.valid \
    --stat ${stat} \
    --binary-target ${binary} \
    --type bgen \
    --target-list ~/projects/shared/shared_resources/data_PRS/target.bgen.files.gruffalo \
    --geno 0.1 \
    --info 0.9 \
    --maf 0.05 \
    --pheno ~/projects/shared/shared_resources/data_PRS/phenos/sample.${phenofile} \
    --pheno-col ${pheno} \
    --cov ~/projects/shared/shared_resources/data_PRS/phenos/cov.${phenofile} \
    --cov-col sex,age,pc1,pc2,SES \
    --cov-factor sex,SES \
    --out ${base}.${pheno} \
    --allow-inter \
    --thread 8     
    
EOF1

chmod +x pgs.${base}.${pheno}.sh

sbatch pgs.${base}.${pheno}.sh

cd ../

done
