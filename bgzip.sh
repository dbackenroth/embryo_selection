#!/bin/bash
#SBATCH --mem=6g
#SBATCH --time=2:20:0
##SBATCH --array=1-23

parentdir=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/
dir=${parentdir}/LIJMC

for i in {1..22}
do
  bgzip ${dir}/LIJMC37_phased_${i}.vcf
  tabix ${dir}/LIJMC37_phased_${i}.vcf.gz
done

