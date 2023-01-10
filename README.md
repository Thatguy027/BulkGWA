# BulkGWA

### Included files

* data/
  * `strain_isotype_lookup.tsv`
    * file that has the CeNDR strain-to-isotype mapping to convert strain to isotype names if needed
  * `strains_used.txt`
    * file that contains strain names used in the experiment
  * `WI.20210121.hard-filter_allele_counts.tsv.gz`
    * allele counts for variants present in the CeNDR population, used for genotype pruning
  * `aser/ASER_test.table.gz` - test GATK ASER output
  * `aser/ASER_test2.table.gz` - test GATK ASER output
* scripts/
  * `prepare_plink_genotypes.sh`
    * download CeNDR hard-filtered VCF and convert to plink file format
  * `process_genotypes_n_counts.R`
    * required inputs are PLINK formatted genotypes and an aser directory that contains the output of GATK ASER
    * output is the input for `prune_genotypes.R`
  * `prune_genotypes.R`
    * input is the output of `process_genotypes_n_counts.R`
    * filters variant sites to only be those with high genotyping rates in the CeNDR population
    * output is the input for `strain_frequency_inference.R`

### How to use

* generate input files:
 * run `prepare_plink_genotypes.sh`
 * generate GATK ASER output (not included in repo yet)
 * run `process_genotypes_n_counts.R`
 * run `prune_genotypes.R`
 * run `strain_frequency_inference.R`
