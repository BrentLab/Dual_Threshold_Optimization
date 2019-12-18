# Dual Threshold Optimzation (DTO)
Dual Threshold Optimization (DTO) is a method that sets the thresholds for TF binding and TF-perturbation response by considering both data sets together. DTO chooses, for each TF, the pair of (binding, response) thresholds that minimizes the probability that the overlap between the bound and responsive sets results from random gene selection

### Package requirement
Install the following Python libraries if not installed, using your favorite package management tool such `pip` or `conda`:
```
pip install scipy numpy pandas
```

### [Optional] Job queuing system requirement
Job queuing system [SLURM](http://slurm.schedmd.com/documentation.html) is preferred to enable high-throughput computing. SLURM version tested: `slurm-wlm 17.11.7`.  


### Input data

### Code usage
####1. Run DTO on authentic datasets
```
python thresholdSearch.py -d <response_csv> --DE_decreasing <True/False> -b <binding_csv> --Bin_decreasing <True/False> --sbatch_loc <output_dir>/authentic_model/ --genes_universe <gene_universe_file> [--geneNames_file <gene_name_file>] [--run_local]
```
Required arguments:
- `-d <response_csv>` CSV file of transcriptional response levels of TF-perturbation. The response levels are based on differential expression or DE analysis on TF-perturbed expression profile vs. pre-perturbation profile. Those levels can be represented as absolute of log fold-changes (LFCs) or p-values.
- `-b <binding_csv>` CSV file of TF binding strengths. The binding strengths can be represented as occupancy levels of the peaks (e.g. heights of ChIP peaks) or statistical significances (e.g. p-values).
- `--DE_decreasing <True/False>` Use 'True' if the response levels should be raneked in descending order. Use 'False' if ranked in ascending order. For example: If using LFCs, set 'True' as higher absolute LFC represents stronger response. If using p-values, set 'False' as lower p-value represents stronger response.
- `--Bin_decreasing <True/False>` Use 'True' if the binding strengths should be ranked in descending order. Use 'False' if ranked in ascending order. For example: If using occupancy levels, set 'True' as higher occupancy level represents stronger TF binding. If using p-values, set 'False' as lower p-value represents stronger TF binding.
- `--sbatch_loc <output_dir>/authentic_model/` Directory for output data.
- `--genes_universe <gene_universe_file>` Gene universe list file, which contains all possible genes to be considered as the universe for calculating overlap statistics.

Optional arguments: 
- `--geneNames_file <gene_name_file>` CSV file of gene name conversion, if the systematic gene names in the data files are preferred to be converted into common gene names.
- `--run_local` If no SLURM is available, set this flag to run DTO in serial fashion on your local machine.

####2. Run DTO on randomized datasets
```
python thresholdSearch.py --random True -d <response_csv> --DE_decreasing <True/False> -b <binding_csv> --Bin_decreasing <True/False> --sbatch_loc <output_dir>/random_models/ --genes_universe <gene_universe_file> [--geneNames_file <gene_name_file>] [--run_local]
```
Required arguments:
- Same as above, except:
- `--random True` Use 'True' for running DTO on randomized input data. The DTO results on randomized data is used to create an empirical null distribution.
- `--sbatch_loc <output_dir>/random_models/` Make a new output directory for null model output.

####3. Summarize DTO results
```
python summarizeFinalResults.py -i <output_dir> [--run local]
```
Required arguments:
- `-i <output_dir>` Directory of DTO results including both authentic model and randomized models.

Optional arguments:
- `--run_local` If no SLURM is available, set this flag to run DTO in serial fashion on your local machine.

### Output data
