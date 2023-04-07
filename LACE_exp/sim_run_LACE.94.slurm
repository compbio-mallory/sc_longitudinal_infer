#!/bin/bash
#SBATCH --job-name=LACE94
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u3/rep5/output94.out
#SBATCH --error=sim_output/u3/rep5/error94.out
module load R
Rscript run_LACE.R ./sim_input/u3/rep5/u3_rep5.RData sim_output/u3/rep5/u3_rep5_inf.rdata sim_output/u3/rep5/u3_rep5_lgTree.xml
