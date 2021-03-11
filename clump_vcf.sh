#!/bin/bash
#SBATCH --mem=15g
#SBATCH --time=2:20:0
##SBATCH --array=1-23

#sbatch count_bases.sh

vcf=$1
daner=$2
outfile=$3
#clumpfield=$4

#vcf="phased_no_ambiguous.vcf.gz"
#daner="CD_trans_ethnic_association_summ_stats_b37.txt.gz"
#outfile="CD_no_ambiguous"

#parentdir=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/
#vcf=${parentdir}/LIJMC/LIJMC_all_chrs_no_ambiguous.vcf.gz
#daner=${parentdir}/daner_PGC_SCZ_w3_90_0518d__lencz.gz
#outfile=${parentdir}/daner_no_ambiguous

plink1.9 --vcf $vcf --maf 0.01 --clump $daner --clump-kb 250 --clump-r2 0.1 --clump-p1 0.05 --out $outfile
