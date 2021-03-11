#!/bin/bash

#SBATCH --mem=5g
#SBATCH --time=7:0:0
#SBATCH --array=1-22
#SBATCH --cpus-per-task=2

gmapfile=/cs/icore/db2175/bin/Eagle_v2.4/tables/genetic_map_hg19_withX.txt.gz
dir=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/crohns/
phasingdir=${dir}/phasing/
vcf=${dir}/ToddSamples.vcf.gz
chr=$SLURM_ARRAY_TASK_ID
#chr=11

unphased=${phasingdir}/unphased_${chr}.vcf
phased=${phasingdir}/phased_${chr}

eagle=/cs/icore/db2175/bin/Eagle_v2.4/eagle

# eagle accepts bed/bim/fam input but these bed/bim/fams have data on chromosome 0 which throws off eagle, so instead just convert to vcf and then select the right vcf
bcftools view $vcf --regions $chr -O v -o ${unphased} 
bgzip $unphased
$eagle --numThreads 2 --chrom $chr --geneticMapFile $gmapfile --vcf ${unphased}.gz --outPrefix ${phased}
tabix ${phased}.vcf.gz
