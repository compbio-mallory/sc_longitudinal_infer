# scLongTree
Computational tool to infer longitudinal tree for scDNAseq data

## Table of Contents
- [Commands to run BnpC, scLongTree, getting SNV placements for all tools, evaluation and simulation.](#commands_3methods)
    * [Clustering](#clustering)
    * [scLongTree](#scLongTree)
    * [SNV placements](#snvAccuracy)
    * [Evaluation](#evaluation)
    * [Simulation](#simulation)

# <a name="usage_of_scDNAseq_clustering"></a>Usage of scDNAseq Clustering.
## <a name="software_requirements"></a>Software Requirements ##

1. Python 2.7.15 or up.

# <a name="commands_3methods"></a>Commands to run BnpC, scLongTree, getting SNV placements for all tools, evaluation and simulation. #

## <a name="clustering"></a>Clustering ##

Since we used BnpC as a clustering algorithm in our pipeline for the evaluation we will briefly list the steps we followed:

1. We process the noisy input matrix provided in ``` sample/input.D.csv ``` with the following scripts provided in ``` Clustering ```

	``` python processSimInput.py -input $j -output outputFile.csv ```

2. We run Bnpc on the processed input data and get the consensus genotype using the following script:

	``` python gen_bnpc_getGmatrix.py -cc assignment.txt -gp BnpCs_genotypes_posterior_mean.tsv -celltp sample/cellsInSubclones.csv -op consensus_genotype.tsv ```

3. You can also run your own clustering algorithm and get a cluster assignment file similar to assignment.txt and a consensus genotype similar to one provided in ``` Clustering```

## <a name="scLongTree"></a>scLongTree ## 

To run the algorithm we have all the necessary scripts in ``` algorithm ``` and sample inputs in ``` sample ```. The below steps will be explained with scripts provided in these folders.

1. Once we have the clustering results we will stratify it across timepoints with the following script:

	``` python algorithm/inputForTree.py -cc Clustering/assignment.txt -cg Clustering/consensus_genotype.tsv -sim false -celltp sample/cellsInSubclones.csv -op consensus_genotype_withTimepoints.csv -ccop cluster_cell.npy ```

2. Next we run our algorithm to infer the longitudinal tree:

	``` python algorithm/longitudinalTree.py -input consensus_genotype_withTimepoints.csv -cg Clustering/consensus_genotype.tsv -cc cluster_cell.npy -D sample/input.D.csv -celltp sample/cellsInSubclones.csv -sim false -op_tree inferredTree.csv -cloneNode cloneNode.npy -cloneCells cloneCells.npy ```

3. Now, given that we run the algorithm multiple times we select the best tree with following script. You can pass n no. of trees as an argument. I passed three trees below as arguments:

	``` python selectBestClusterAndTree.py -tree inferredTree1.csv inferredTree2.csv inferredTree3.csv ```

## <a name="snvAccuracy"></a>Calculating SNV accuracy ##

1. Given the results we evaluate based on placement of SNVs on the tree. Since the tree inferred is different for each tool we used the following scripts for each respective trees to get the pair of SNVS for ancestral, on same branch and incomparable/parallel:

	scLongTree: ``` python algorithm/getSNVrelations.py -tree inferredTree.csv -output outputFileLoc ```

	LACE:	``` python LACE_exp/getSNVrelations.py -B lace_B.json -output outputFileLoc ```

	SCITE: ``` python SCITE_exp/getSNVpairs.py -tree sciteTree.gv -output outputFileLoc ```

	SiCloneFit: ``` python SiCloneFit_exp/getSNVpairs.py -tree best_MAP_tree.txt -genotype best_MAP_predicted_genotype.txt -output outputFileLoc ```
	
## <a name="evaluation"></a>Evaluation ##

1. Given the pair of SNVs we get the SNV pair accuracy by running the following script:

	``` python evaluation/evaluateSNVpairs.py -gt groundTruth_locationOf_SavedSNVpairs -scLongTree locationOf_SavedSNVpairs -siclonefit locationOf_SavedSNVpairs -scite locationOf_SavedSNVpairs -lace locationOf_SavedSNVpairs ```

2. We calculate the precision and recall for LACE and scLongTree:

	``` python calculate_PrecisionRecall.py -mut sample/clone_mut.csv -gtTree groundTruthTree -lgTree inferredTree -lace true/false ```

## <a name="simulation"></a>Simulation ##

### Steps to run the simulator using different variables. ###

1. To generate the tree with desired parameters use the following script:
   
   ``` python gen_tree.py -T 3 -t_v 15 7 0 -v1 0.2 -v2 0.4 -u 0.2 -B 0.2 -theta 0.2 -o output_tree.csv ```
   
   Parameters:
   
   * ```-T``` pass the number of timepoints.
   * ```-t_v``` years representing the timepoints.
   * ```-v1``` condition variable to decide distribution of mutations.
   * ```-v2``` condition variable to decide distribution of mutations.
   * ```-u``` condition variable deciding unobserved subclones.
   * ```-B``` beta split variable.
   * ```-theta``` condition variable to decide if a subclone should have new mutations.
   * ```-o``` output file.

3. Once we have the tree we can assign cells, mutations and include false positives, false negatives, missing data using the following script:

   	``` python sim_par.py -a 0.01 -b 0.2 -m 0.2 -t 3 -mc 1.7 -f output_tree.csv -P outputFilePrefix ```
   
   Parameters:
   * ```-a``` false positive rate.
   * ```-b``` false negative rate.
   * ```-m``` missing rate.
   * ```-t``` timepoints.
   * ```-mc``` mutation rate.
   * ```-f``` the tree generated earlier using number of leaves and beta splitting variable.
   * ```-P``` the prefix and directory for output files.

