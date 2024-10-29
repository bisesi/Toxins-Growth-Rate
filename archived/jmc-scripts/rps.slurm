#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=10g
#SBATCH --mail-type=NONE

module load gurobi/9.0.2
module load python
unset PYTHONPATH
unset PYTHONHOME
unset PYTHONSTARTUP


cd comets_RPS

python growth_rate_vs_toxins_commandline.py $SLURM_ARRAY_TASK_ID
