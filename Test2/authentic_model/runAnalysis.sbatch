#!/bin/bash
#SBATCH -D /scratch/mblab/acharyas/Human_TF_Project/DTO/Dual_Threshold_Optimization/Code
#SBATCH --mem=2G
#SBATCH -J runDTO_%A_%a
#SBATCH -o ../Test2/authentic_model//log/runDTO_%A_%a.out
#SBATCH -e ../Test2/authentic_model//log/runDTO_%A_%a.err
#SBATCH --array=0-1%50
ID=${SLURM_ARRAY_TASK_ID}
python runDualThreshold.py --de_file ../../../Overlap_Remap_ENCODE_Perturbation/DTO_Input/DTO_run_25_percent_threshold_Perturbation_And_Binding_padj/test_perturbation.txt --bin_file ../../../Overlap_Remap_ENCODE_Perturbation/DTO_Input/DTO_run_25_percent_threshold_Perturbation_And_Binding_padj/test_binding.txt --DE_decreasing False --Bin_decreasing True --TF_num ${ID} --rank_width 1.01 --opt_crit pval --genes_universe ../../../Overlap_Remap_ENCODE_Perturbation/gene_universe_remap_encodePerturbation.txt --geneNames_file "" --DE_pval_lower_bound 1 --Bin_pval_lower_bound 0.1 --output_dir ../Test2/authentic_model/
