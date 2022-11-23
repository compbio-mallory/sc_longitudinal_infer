
python getGT_ConsensusGenotype.py -cg ../simulation/default/rep5/default_rep5.G.csv -cell ../simulation/default/rep5/default_rep5.SNVcell.csv -op ground_truth/default/rep5/rep5

# Evaluation script to compare the consensus genotypes at diff timepoints.
python evaluateMetrics.py -cg sim_consensus_genotype/default_rep1_t2.csv -gtG ground_truth/default/rep1/rep1_t2.csv -h1 true
python ../evaluation.py -i "bnpc:"./sim_output/default/rep1/assignment.txt -G ./sim_input/default/rep1/default_rep1.G.csv -v 

# Calculate precision recall
python calculate_PrecisionRecall.py -mut ../simulation/default/rep1/default_rep1.mut.csv -gtTree ../simulation/default/rep1/tree_default_rep1.csv -lgTree inferred_tree/rep1_tree.csv 
#-cg sim_input_ForTree/default_rep1.Gprime.csv -clone clone_node.npy 

# Evaluate the unobserved subclones
python evaluate_unobservedSubclones.py -mut ../simulation/default/rep1/default_rep1.mut.csv -gtTree ../simulation/default/rep1/tree_default_rep1.csv -lgTree inferred_tree/rep1_tree.csv 
