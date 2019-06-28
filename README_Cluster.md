README for Dual Threshold Optimization (DTO)

**I. Prepping the data:**
Every dataset must be organized in a single directory with the following files:
1. TFNames.csv - a CSV file with a list of the transcription factors (or conditions) for which data is available
2. GeneNames.csv - a CSV file containing a matrix of the genes that correspond to each score in the data matrix
   1. Each row in this matrix corresponds to the appropriate index in the list of TFNames
   2. Even if the data is organized with the same order of genes in each condition, the matrix must be filled with one row per TF/condition
3. Data.csv - a CSV file containing the values that are going to be used to rank-order the genes in GeneNames.csv
   1. Each row in this matrix corresponds to the appropriate index in the list of TFNames
   2. The data can be of any form (Pvalue, log-fold change, etc.)
4. Using NetProphet data:
	The output network from NetProphet can be used easily with DTO. regulators.txt corresponds to TFNames.csv and genes.txt corresponds to GeneNames.csv after it has been copied to fill the full matrix of data. The adjcency matrix can either be converted directly to Data.csv, but often times it is desirable to restrict the number of edges to the top N edges. This can be done using computeTopNEdges() in loadData.py. Simply modify the command in the main method to adjust the input file, output directory, and the desired number of edges. This will produce the Data.csv file that can be used directly in DTO, maintaining the order of TFs and genes in the matrix. 

**II. Overview of the analysis**
1. Scripts
   1. loadData.py - contains functions both to parse data into the appropriate format and functions used by other scripts to load saved data for analysis
   2. thresholdSearch.py - the script used to prepare the data and sbatch scripts for the dual threshold analysis (this is the script you will run)
   3. runDualThreshold.py - the script that perfoms that actual dual thresholding (run by the sbatch scripts created by thresholdSearch.py)
   4. statistics.py - contains the funtions used by runDualThreshold.py to compute various statistical measures on the subsets
   5. compileResults.py - compiles the results from the many csv files produced by either standard and randomized dual threshold analyses into a single file
   6. computeAcceptableTFs.py - parses the results from a corresponding pair of compiled standard and randomized results and produces the following outputs:
    1. Results for the acceptable TFs 
    2. Target genes for the acceptable TFs
    3. Edges for the acceptable TFs
    4. Cutoffs for each TF from the randomizations used to determine the set of acceptable TFs
   7. createBinaryEdgeFile.py - takes in the set of acceptable TFs from two different analyses and produces a combined set of edges for those two analyses with binary flags indicated which analyses each edge participates in
2. The following are the parameters necessary to run a DTO analysis
   * `-d/--de_dir`: the path to the directory containing differential expression data
   * `-b/--bin_dir`: the path to the directory containing binding data
   * `-j/--DE_decreasing`: a boolean indicating if the DE data should be ranked in increasing or decreasing order
   * `-k/--Bin_decreasing`: a boolean indicating if the binding data should be ranked in increasing or decreasing order
   * `-r/--random`: a boolean indicating whether this analysis should perform the standard optimization or randomized optimizations
   * `-g/--geneNames_file`: the path to the conversion file between common and systematic names of genes for yeast
   * `-f/--sbatch_loc (default = ".")`: the path to the output directory (or the directory from which the sbatch script will run)
   * `-w/--rank_width (default = "1.01")`: a parameter used to define the resolution of the search (it is recommended to keep this at the default value)
   * `-o/--opt_crit (default = "pval")`: a parameter used to define the paramater used for optimization (it is recommended to keep this at the default value)
   * `-u/--genes_universe (default = "")`: the path to a file that defines the universe of genes that should be included in the analysis
   * `-a/--organism (default = "yeast")`: a string indicating the organism that is being studied
3. Running a DTO requires two different analyses to be run, one of which will be randomized. After each one has been run, they must both be compiled (in a sense, organizing the results into a single file). Finally, the number of acceptable TFs must be computed and a final set of results produced.

**III. Running the analysis** 
(The following steps show how to walk through each of these steps for an analysis. In this example, we will run an analysis that intersects the Harbison ChIP data with a ZEV-based NetProphet network.)
1. Parse Harbison ChIP data into the proper format for DTO:
	The parseData() function in loadData.py will convert a directory that contains a distinct file for each condition into the TF list and gene/data matrices necessary for DTO. To use this function, modify the skeleton in the main method, and run the script on the command line. The function takes in the following parameters:
   1. locIn - The path to the directory of input data (each file must have the TF/condition in its name, and must be formatted as a list of genes with scores for each)
   2. locOut - The path to the directory to output the data to
   3. colToKeep - The column of the input files that contains the scores (it is assumed that the gene names are in column 0)
   4. charsToDrop - The number of characters to drop from the end of the file name to describe each tf (Ex. if the file name was TFName.peaks.txt, you would enter 10 for this parameter)
This function can handle many but not all formats of raw data. You may need to define your own function for this to ensure it matches the form described in Part I.
2. Parse the NetProphet network into the proper form and restrict the network to the top N edges:
	The computeTopNEdges() function in loadData.py will select the top N edges from the network and replace the rest of the edges with a 0. To use this function, modify the skeleton in the main method, and run the script on the command line. The function takes in the following parameters:
   1. locIn - The path to the directory of input data (this is most likely the .adjmtr output file from the NetProphet output)
   2. locOut - The path to the directory to output the data to
   3. numEdges - The number of edges to restrict the network to
3. Create the necessary directories/files on the cluster:
   1. Create a directory for the standard analysis (Ex. Harbison_NP_ZEV)
   2. Create a directory for the randomized analysis (Ex. Harbison_NP_ZEV_Rand)
   3. Prepare a text file with a list of all genes that are to be included in the "universe" of genes for the analysis
   4. If running an analysis with yeast, prepare a gene conversion file with common names in column 0 and systematic names in column 1
4. Run the analyses:
   1. Run the standard analysis - modify the file paths of the following command to agree with your file structure and run it on the command line:
	`python thresholdSearch.py --de_dir home/Data/NetProphet_Z_15_45_90_150K/ --bin_dir home/Data/HarbisonYPD/ --sbatch_loc home/Analyses/Harbison_NetProphet_ZEV_15_45_90/ --genes_universe home/ExtraFiles/Harbison_NetProphet_Z_Universe.txt --geneNames_file home/ExtraFiles/YeastCommonAndSystematicGeneNames.csv --DE_decreasing True --Bin_decreasing False --random False`
   2. Run the randomized analysis - modify the file paths of the following command to agree with your file structure and run it on the command line:
	`python thresholdSearch.py --de_dir home/Data/NetProphet_Z_15_45_90_150K/ --bin_dir home/Data/HarbisonYPD/ --sbatch_loc home/Analyses/Harbison_NetProphet_ZEV_15_45_90_Rand/ --genes_universe home/ExtraFiles/Harbison_NetProphet_Z_Universe.txt --geneNames_file home/ExtraFiles/YeastCommonAndSystematicGeneNames.csv --DE_decreasing True --Bin_decreasing False --random True`
   3. Compile the results - modify the file paths of the following command to agree with your file structure and run it on the command line:
	`python compileResults.py --csv_dir home/Analyses/Harbison_NetProphet_ZEV_15_45_90/`
   4. Compile the results - modify the file paths of the following command to agree with your file structure and run it on the command line:
	`python compileResults.py --csv_dir home/Analyses/Yeast/Harbison_NetProphet_ZEV_15_45_90_Rand/ --rand_type "split"`
   5. Identify the acceptable TFs - modify the file paths of the following command to agree with your file structure and run it on the command line:
	`python computeAcceptableTFs.py --data_file home/Analyses/Harbison_NetProphet_ZEV_15_45_90/Harbison_NetProphet_ZEV_15_45_90.csv --rand_file home/Analyses/Harbison_NetProphet_ZEV_15_45_90_Rand/Harbison_NetProphet_ZEV_15_45_90_Rand.csv --output_dir home/Analyses/Harbison_NetProphet_ZEV_15_45_90/`
