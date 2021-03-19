#!/bin/bash
#SBATCH --mem=6g
#SBATCH --time=0:40:0
#SBATCH --array=3-1500

##1-300

disease=$1

R -q -e "DISEASE=$disease; source('recombine.R'); DoRun(1500, $SLURM_ARRAY_TASK_ID)"

