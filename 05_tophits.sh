#!/bin/bash

pheno=var

cd results.gwas/

awk 'NR == 1' sumstats.${pheno} > header.${pheno}.temp

awk '{ if ($14 < 0.00001) { print } }' sumstats.${pheno} > tophits.${pheno}.temp

cat header.${pheno}.temp tophits.${pheno}.temp > tophits.${pheno}

rm *.temp

cd ../