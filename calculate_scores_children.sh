#!/bin/bash
#SBATCH --mem=6g
#SBATCH --time=2:20:0
##SBATCH --array=1-23

vcf=$1

dir=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/
#vcf=${dir}/Peds/3000_5.vcf
selected_snps=${dir}/daner.valid.snp
gwas_info=${dir}/daner_flipped.tsv
out=${vcf}.scores
#freq=${dir}/LIJMC_score_snps.recode_freq.afreq

plink2 --vcf $vcf --score $gwas_info list-variants 3 5 9 header --extract $selected_snps --out $out
#--read-freq $freq

