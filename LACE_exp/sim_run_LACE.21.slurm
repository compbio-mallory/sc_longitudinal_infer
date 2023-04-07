#!/bin/bash
#SBATCH --job-name=LACE21
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp3/rep2/output21.out
#SBATCH --error=sim_output/bp3/rep2/error21.out
module load R
Rscript run_LACE.R ./sim_input/bp3/rep2/bp3_rep2.RData sim_output/bp3/rep2/bp3_rep2_inf.rdata sim_output/bp3/rep2/bp3_rep2_lgTree.xml
