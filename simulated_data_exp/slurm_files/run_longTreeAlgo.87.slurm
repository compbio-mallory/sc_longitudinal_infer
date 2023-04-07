#!/bin/bash
#SBATCH --job-name=lgTree87
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u1/rep3/lgOutput87.out
#SBATCH --error=sim_output/u1/rep3/lgError87.out
time python gen_bnpc_getGmatrix.py -cc sim_output/u1/rep3/m1/assignment.txt -gp sim_output/u1/rep3/m1/genotypes_posterior_mean.tsv -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -doublet false -op sim_consensus_genotype/u1_rep3_m1.G
time python inputForTree.py -cc sim_output/u1/rep3/m1/assignment.txt -cg sim_consensus_genotype/u1_rep3_m1.G.tsv -sim false -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -op sim_input_ForTree/u1_rep3_m1.Gprime.csv -ccop sim_input_ForTree/u1_rep3_m1.npy
time python getEstimatedFPFN.py -G sim_consensus_genotype/u1_rep3_m1.G.tsv -D sim_input/u1/rep3/input_u1_rep3.D.csv -op sim_output/u1/rep3/m1/u1_rep3_m1_FPFN.npy
time python longitudinalTree.py -input sim_input_ForTree/u1_rep3_m1.Gprime.csv -cg sim_consensus_genotype/u1_rep3_m1.G.tsv -cc sim_input_ForTree/u1_rep3_m1.npy -D sim_input/u1/rep3/input_u1_rep3.D.csv -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -FPFN sim_output/u1/rep3/m1/u1_rep3_m1_FPFN.npy -sim false -op_tree all_trees/u1_rep3_m1_tree.csv -cloneNode all_trees/u1_rep3_m1_cloneNode.npy -cloneCells all_trees/u1_rep3_m1_cloneCells.npy
time python gen_bnpc_getGmatrix.py -cc sim_output/u1/rep3/m2/assignment.txt -gp sim_output/u1/rep3/m2/genotypes_posterior_mean.tsv -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -doublet false -op sim_consensus_genotype/u1_rep3_m2.G
time python inputForTree.py -cc sim_output/u1/rep3/m2/assignment.txt -cg sim_consensus_genotype/u1_rep3_m2.G.tsv -sim false -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -op sim_input_ForTree/u1_rep3_m2.Gprime.csv -ccop sim_input_ForTree/u1_rep3_m2.npy
time python getEstimatedFPFN.py -G sim_consensus_genotype/u1_rep3_m2.G.tsv -D sim_input/u1/rep3/input_u1_rep3.D.csv -op sim_output/u1/rep3/m2/u1_rep3_m2_FPFN.npy
time python longitudinalTree.py -input sim_input_ForTree/u1_rep3_m2.Gprime.csv -cg sim_consensus_genotype/u1_rep3_m2.G.tsv -cc sim_input_ForTree/u1_rep3_m2.npy -D sim_input/u1/rep3/input_u1_rep3.D.csv -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -FPFN sim_output/u1/rep3/m2/u1_rep3_m2_FPFN.npy -sim false -op_tree all_trees/u1_rep3_m2_tree.csv -cloneNode all_trees/u1_rep3_m2_cloneNode.npy -cloneCells all_trees/u1_rep3_m2_cloneCells.npy
time python gen_bnpc_getGmatrix.py -cc sim_output/u1/rep3/m3/assignment.txt -gp sim_output/u1/rep3/m3/genotypes_posterior_mean.tsv -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -doublet false -op sim_consensus_genotype/u1_rep3_m3.G
time python inputForTree.py -cc sim_output/u1/rep3/m3/assignment.txt -cg sim_consensus_genotype/u1_rep3_m3.G.tsv -sim false -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -op sim_input_ForTree/u1_rep3_m3.Gprime.csv -ccop sim_input_ForTree/u1_rep3_m3.npy
time python getEstimatedFPFN.py -G sim_consensus_genotype/u1_rep3_m3.G.tsv -D sim_input/u1/rep3/input_u1_rep3.D.csv -op sim_output/u1/rep3/m3/u1_rep3_m3_FPFN.npy
time python longitudinalTree.py -input sim_input_ForTree/u1_rep3_m3.Gprime.csv -cg sim_consensus_genotype/u1_rep3_m3.G.tsv -cc sim_input_ForTree/u1_rep3_m3.npy -D sim_input/u1/rep3/input_u1_rep3.D.csv -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -FPFN sim_output/u1/rep3/m3/u1_rep3_m3_FPFN.npy -sim false -op_tree all_trees/u1_rep3_m3_tree.csv -cloneNode all_trees/u1_rep3_m3_cloneNode.npy -cloneCells all_trees/u1_rep3_m3_cloneCells.npy
time python gen_bnpc_getGmatrix.py -cc sim_output/u1/rep3/m4/assignment.txt -gp sim_output/u1/rep3/m4/genotypes_posterior_mean.tsv -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -doublet false -op sim_consensus_genotype/u1_rep3_m4.G
time python inputForTree.py -cc sim_output/u1/rep3/m4/assignment.txt -cg sim_consensus_genotype/u1_rep3_m4.G.tsv -sim false -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -op sim_input_ForTree/u1_rep3_m4.Gprime.csv -ccop sim_input_ForTree/u1_rep3_m4.npy
time python getEstimatedFPFN.py -G sim_consensus_genotype/u1_rep3_m4.G.tsv -D sim_input/u1/rep3/input_u1_rep3.D.csv -op sim_output/u1/rep3/m4/u1_rep3_m4_FPFN.npy
time python longitudinalTree.py -input sim_input_ForTree/u1_rep3_m4.Gprime.csv -cg sim_consensus_genotype/u1_rep3_m4.G.tsv -cc sim_input_ForTree/u1_rep3_m4.npy -D sim_input/u1/rep3/input_u1_rep3.D.csv -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -FPFN sim_output/u1/rep3/m4/u1_rep3_m4_FPFN.npy -sim false -op_tree all_trees/u1_rep3_m4_tree.csv -cloneNode all_trees/u1_rep3_m4_cloneNode.npy -cloneCells all_trees/u1_rep3_m4_cloneCells.npy
time python gen_bnpc_getGmatrix.py -cc sim_output/u1/rep3/m5/assignment.txt -gp sim_output/u1/rep3/m5/genotypes_posterior_mean.tsv -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -doublet false -op sim_consensus_genotype/u1_rep3_m5.G
time python inputForTree.py -cc sim_output/u1/rep3/m5/assignment.txt -cg sim_consensus_genotype/u1_rep3_m5.G.tsv -sim false -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -op sim_input_ForTree/u1_rep3_m5.Gprime.csv -ccop sim_input_ForTree/u1_rep3_m5.npy
time python getEstimatedFPFN.py -G sim_consensus_genotype/u1_rep3_m5.G.tsv -D sim_input/u1/rep3/input_u1_rep3.D.csv -op sim_output/u1/rep3/m5/u1_rep3_m5_FPFN.npy
time python longitudinalTree.py -input sim_input_ForTree/u1_rep3_m5.Gprime.csv -cg sim_consensus_genotype/u1_rep3_m5.G.tsv -cc sim_input_ForTree/u1_rep3_m5.npy -D sim_input/u1/rep3/input_u1_rep3.D.csv -celltp sim_input/u1/rep3/input_u1_rep3.SNVcell.csv -FPFN sim_output/u1/rep3/m5/u1_rep3_m5_FPFN.npy -sim false -op_tree all_trees/u1_rep3_m5_tree.csv -cloneNode all_trees/u1_rep3_m5_cloneNode.npy -cloneCells all_trees/u1_rep3_m5_cloneCells.npy
time python selectBestCluster.py -tree all_trees/u1_rep3_m1_tree.csv all_trees/u1_rep3_m2_tree.csv all_trees/u1_rep3_m3_tree.csv all_trees/u1_rep3_m4_tree.csv all_trees/u1_rep3_m5_tree.csv -bnpcOp sim_output/u1/rep3
