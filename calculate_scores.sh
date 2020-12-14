#!/bin/bash
#SBATCH --mem=6g
#SBATCH --time=2:20:0
##SBATCH --array=1-23

dir=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/
vcf=${dir}/LIJMC/LIJMC_all_chrs_no_ambiguous.vcf.gz
selected_snps=${dir}/daner.valid.snp
gwas_info=${dir}/daner_flipped.tsv
out=${dir}/LIJMC_scores

plink2 --vcf $vcf --score $gwas_info list-variants 3 5 9 header --extract $selected_snps --out $out

