#!/bin/bash
#SBATCH --mem=6g
#SBATCH --time=2:20:0
##SBATCH --array=1-23

invcfgz=$1
outvcf=$2

#invcfgz=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/crohns/phased.vcf.gz
#outvcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/crohns/phased_no_ambiguous.vcf

#invcfgz=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/LIJMC_all_chrs.vcf.gz
#outvcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/LIJMC_all_chrs_no_ambiguous.vcf

gunzip -c $invcfgz | awk '{if( !(($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="C" && $5=="G") || ($4=="G" && $5=="C") || ($4=="D" || $4 =="I" || $4 =="." || $5=="D" || $5=="I" || $5=="."))) print }' > $outvcf

bgzip $outvcf
tabix ${outvcf}.gz


