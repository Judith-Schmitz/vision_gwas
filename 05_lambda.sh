#!/bin/bash

pheno=var

~/apps/conda/envs/geno_utils/bin/python 05_compute_lambda.py \
-i results.gwas/sumstats.${pheno} \
-w \
-f CHISQ_BOLT_LMM_INF \
--chi2 \
--snp-field SNP
