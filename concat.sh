#!/bin/bash
#SBATCH --mem=6g
#SBATCH --time=2:20:0
##SBATCH --array=1-23

invcf=$1
outvcf=$2

#invcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/crohns/phasing/phased*vcf.gz
#outvcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/crohns/phased.vcf


#invcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/*vcf.gz
#outvcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/LIJMC_all_chrs.vcf

bcftools concat $invcf --allow-overlaps -o $outvcf
bgzip $outvcf
tabix ${outvcf}.gz
