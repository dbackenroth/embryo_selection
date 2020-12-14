#!/bin/bash
#SBATCH --mem=6g
#SBATCH --time=2:20:0
##SBATCH --array=1-23


outvcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/LIJMC_all_chrs.vcf

bcftools concat /vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/*vcf.gz --allow-overlaps -o $outvcf
bgzip $outvcf
tabix ${outvcf}.gz


