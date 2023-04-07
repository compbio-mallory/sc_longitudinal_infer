#!/bin/bash
#SBATCH --job-name=LACE37
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mc12/rep3/output37.out
#SBATCH --error=sim_output/mc12/rep3/error37.out
module load R
Rscript run_LACE.R ./sim_input/mc12/rep3/mc12_rep3.RData sim_output/mc12/rep3/mc12_rep3_inf.rdata sim_output/mc12/rep3/mc12_rep3_lgTree.xml
