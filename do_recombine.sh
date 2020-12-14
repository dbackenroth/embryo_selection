#!/bin/bash
#SBATCH --mem=6g
#SBATCH --time=0:40:0
#SBATCH --array=1-300
0
##1-300

R -q -e "source('recombine.R'); DoRun(300, $SLURM_ARRAY_TASK_ID)"

