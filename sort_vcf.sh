#!/bin/bash
#SBATCH --mem=4g
#SBATCH --time=20:20:0
##SBATCH --array=1-23

outvcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/LIJMC_all_chrs_no_ambiguous_sorted.vcf

bcftools sort /vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/LIJMC_all_chrs_no_ambiguous.vcf.gz > $outvcf
bgzip $outvcf
tabix ${outvcf}.gz

