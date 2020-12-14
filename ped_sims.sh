#!/bin/bash
#SBATCH --mem=10g
#SBATCH --time=24:20:0
##SBATCH --array=1-2

mapfile=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/Maps
pedsim=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/ped-sim/ped-sim
indef=/cs/icore/db2175/embryos/children.def
invcf=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/LIJMC_all_chrs_no_ambiguous_sorted.vcf
outprefix=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/Peds/sample

$pedsim --bp --fam --pois -d $indef -m $mapfile -o $outprefix

# ./ped-sim -d <in.def> -m <map file> -i <in.vcf/in.vcf.gz> -o <out_prefix> --intf <filename>
