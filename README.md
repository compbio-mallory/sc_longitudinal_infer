# scLongTree
Computational tool to infer longitudinal tree for scDNAseq data

## Table of Contents
- [Commands to run BnpC, scLongTree, getting SNV placements for all tools, evaluation and simulation.](#commands_3methods)
    * [BnpC](#bnpc)
    * [scLongTree](#scLongTree)
    * [SNV placements](#snvAccuracy)
    * [Evaluation](#evaluation)
    * [Simulation](#simulation)

# <a name="usage_of_scDNAseq_clustering"></a>Usage of scDNAseq Clustering.
## <a name="software_requirements"></a>Software Requirements ##

1. Python 2.7.15 or up.

# <a name="commands_3methods"></a>Commands to run BnpC, scLongTree, getting SNV placements for all tools, evaluation and simulation. #

## <a name="bnpc"></a>BnpC ##

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

1. To generate the tree with beta-splitting variable and desired number of leaves use the following script:
        * ``` python gen_tree.py -F $i -B 0.2 -o output_tree.csv ```
        * Parameters:
                * ```-F``` pass the number of leaves or desired number of clones.
                * ```-B``` pass the beta split variable value.
                * ```-o``` output tree name to save.

2. Once we have the tree we can assign cells, mutations and include false positives, false negatives, missing data, doublets using the following script:
        * ``` python sim_par.py -a 0.01 -b 0.2 -m 0.2 -c 500 -n 200 -e 0.1 -f output_tree.csv -P outputFilePrefix ```
        * Parameters:
                * ```-a``` false positive rate.
                * ```-b``` false negative rate.
                * ```-m``` missing rate.
                * ```-c``` number of cells.
                * ```-n``` number of mutations.
                * ```-e``` doublet rate.
                * ```-f``` the tree generated earlier using number of leaves and beta splitting variable.
                * ```-P``` the prefix and directory for output files.

