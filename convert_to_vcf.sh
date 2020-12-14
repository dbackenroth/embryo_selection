#!/bin/bash
#SBATCH --mem=4g
#SBATCH --time=20:20:0
#SBATCH --array=1-23


chr=$SLURM_ARRAY_TASK_ID

parentdir=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/
dir=${parentdir}/LIJMC
samplefile=${dir}/LIJMC37_.${chr}
vcf=${dir}/LIJMC37_phased_${chr}.vcf

shapeit=/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit

$shapeit -convert --input-haps $samplefile --output-vcf $vcf

