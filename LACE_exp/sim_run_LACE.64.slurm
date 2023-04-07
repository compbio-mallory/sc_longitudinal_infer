#!/bin/bash
#SBATCH --job-name=LACE64
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/t4/rep5/output64.out
#SBATCH --error=sim_output/t4/rep5/error64.out
module load R
Rscript run_LACE.R ./sim_input/t4/rep5/t4_rep5.RData sim_output/t4/rep5/t4_rep5_inf.rdata sim_output/t4/rep5/t4_rep5_lgTree.xml
