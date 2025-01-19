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

We used BnpC as a clustering algorithm and used their results from multiple runs as an input to our algorithm.

## <a name="scLongTree"></a>scLongTree ## 

To run the algorithm we have included all scripts in ``` algorithm ``` and sample inputs are in ``` sample ```. The below steps will be explained with scripts provided in these folders.

1. We can run scLongTree using the below command:

	``` python algorithm/selectBestTree.py -m 5 -t t1 t2 t3 -loc sample/bnpc_runs -cells sample/cell_timepoints.csv -D sample/input.D.csv -k 0 -op output ```
	
	Parameters:

	* ```-m``` number of BnpC runs.
	* ```-t``` different timepoints.
	* ```-loc``` file path having BnpC runs.
	* ```-cells``` file describing cells belonging to each timepoint.
	* ```-D``` input genotype matrix with cells as rows and mutations as columns.
	* ```-k``` value of k for k-Dollo model.
	* ```-op``` output file path to save the results.

## <a name="snvAccuracy"></a>Calculating SNV accuracy ##

1. Given the results we evaluate based on placement of SNVs on the tree. Since the tree inferred is different for each tool we used the following scripts for each respective trees to get the pair of SNVS for ancestral, on same branch and parallel:

	scLongTree: ``` python algorithm/scLongTree_getSNVrelations.py -tree output/tree.csv -output outputFilePath ```

	LACE:	``` python LACE_exp/getSNVrelations.py -B lace_B.csv -output outputFilePath ```

	SCITE: ``` python SCITE_exp/getSNVrelations.py -tree sciteTree.gv -output outputFilePath ```

	SiCloneFit: ``` python SiCloneFit_exp/getSNVrelations.py -tree best_MAP_tree.txt -genotype best_MAP_predicted_genotype.txt -output outputFileLoc ```
	
## <a name="evaluation"></a>Evaluation ##

1. Given the pair of SNVs we get the SNV pair accuracy by running the following script:

	``` python evaluation/evaluateSNVpairs.py -gt sample/gtSNVpairs -scLongTree sample/scLongTreeSNVpairs -siclonefit sample/siclonefitSNVpairs -scite sample/sciteSNVpairs -lace sample/laceSNVpairs ```

	Parameters:

	* ```-gt``` File path for SNV pairs from groud truth Tree.
	* ```-scLongTree``` File path for SNV pairs from scLongTree's inferred tree.
	* ```-siclonefit``` File path for SNV pairs from SiCloneFit's inferred tree.
	* ```-scite``` File path for SNV pairs from SCITE's inferred tree.
	* ```-lace``` File path for SNV pairs from LACE's inferred tree.

2. We calculate the precision and recall for LACE and scLongTree:

	``` python evaluation/calculate_PrecisionRecall.py -mut sample/clone_mut.csv -gtTree sample/groundTruthTree.csv -lgTree sample/inferredTree.csv ```
	
	Parameters:

	* ```-mut``` ground truth mutations of each clones.
	* ```-gtTree``` ground truth Tree.
	* ```-lgTree``` longitudinal tree inferred by the algorithms.

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
   * ```-theta``` condition variable deciding new mutations for subclones.
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
   * ```-P``` the prefix and file path for output files.

