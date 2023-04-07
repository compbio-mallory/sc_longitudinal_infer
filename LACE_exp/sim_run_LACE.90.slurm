#!/bin/bash
#SBATCH --job-name=LACE90
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u3/rep1/output90.out
#SBATCH --error=sim_output/u3/rep1/error90.out
module load R
Rscript run_LACE.R ./sim_input/u3/rep1/u3_rep1.RData sim_output/u3/rep1/u3_rep1_inf.rdata sim_output/u3/rep1/u3_rep1_lgTree.xml
