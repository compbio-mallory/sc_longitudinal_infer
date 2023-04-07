#!/bin/bash
#SBATCH --job-name=LACE22
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp3/rep3/output22.out
#SBATCH --error=sim_output/bp3/rep3/error22.out
module load R
Rscript run_LACE.R ./sim_input/bp3/rep3/bp3_rep3.RData sim_output/bp3/rep3/bp3_rep3_inf.rdata sim_output/bp3/rep3/bp3_rep3_lgTree.xml
