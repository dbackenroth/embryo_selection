#!/bin/bash
#SBATCH --mem=6g
#SBATCH --time=2:20:0
##SBATCH --array=1-23

vcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/LIJMC_all_chrs.vcf.gz

outvcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/LIJMC_all_chrs_no_ambiguous.vcf

gunzip -c $vcf | awk '{if( !(($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="C" && $5=="G") || ($4=="G" && $5=="C"))) print }' > $outvcf

bgzip $outvcf
tabix ${outvcf}.gz


