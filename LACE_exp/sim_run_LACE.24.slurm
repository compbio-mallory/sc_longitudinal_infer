#!/bin/bash
#SBATCH --job-name=LACE24
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp3/rep5/output24.out
#SBATCH --error=sim_output/bp3/rep5/error24.out
module load R
Rscript run_LACE.R ./sim_input/bp3/rep5/bp3_rep5.RData sim_output/bp3/rep5/bp3_rep5_inf.rdata sim_output/bp3/rep5/bp3_rep5_lgTree.xml
