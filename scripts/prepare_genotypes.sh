#!/bin/bash

# download CeNDR VCF
wget https://storage.googleapis.com/caendr-site-public-bucket/dataset_release/c_elegans/20210121/variation/WI.20210121.hard-filter.vcf.gz

# generate plink files
for chrom in I II III IV V X
do
plink --allow-extra-chr \
  --biallelic-only \
  --make-bed \
  --out $chrom \
  --output-missing-genotype 9 \
  --recode 01 \
  --set-missing-var-ids @:# \
  --snps-only \
  --chr $chrom \
  --vcf ../WI.20210121.hard-filter.isotype_with_cM.vcf.gz 
done