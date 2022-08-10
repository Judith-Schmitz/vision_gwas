#! /bin/bash

pheno=var
phenofile=eye
binary=F

cd results.pgs.23andMe/

for base in DYS MYO

do

cat > pgs.${base}.${pheno}.sh <<EOF1
#!/bin/bash

#SBATCH --job-name="PGS"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

Rscript ~/apps/PRSice/PRSice.R --dir . \
    --prsice ~/apps/PRSice/PRSice_linux \
    --base ~/projects/shared/shared_resources/data_PRS/base/${base}.chrpos.gz \
    --snp ID \
    --beta \
    --stat effect \
    --binary-target ${binary} \
    --type bgen \
    --target-list ~/projects/shared/shared_resources/data_PRS/target.bgen.files.gruffalo \
    --extract ~/projects/shared/shared_resources/data_PRS/valid/${base}.valid \
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

done

cd ../