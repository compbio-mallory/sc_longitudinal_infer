#!/bin/bash
#SBATCH --job-name=LACE61
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/t4/rep2/output61.out
#SBATCH --error=sim_output/t4/rep2/error61.out
module load R
Rscript run_LACE.R ./sim_input/t4/rep2/t4_rep2.RData sim_output/t4/rep2/t4_rep2_inf.rdata sim_output/t4/rep2/t4_rep2_lgTree.xml
