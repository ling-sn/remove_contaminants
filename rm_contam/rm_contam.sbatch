#!/usr/bin/env bash
#SBATCH --job-name=remove_contaminants
#SBATCH --mail-user=lingsn@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=rm_contam_%u_%A_%a.out
#SBATCH --array=0-11%4
#SBATCH --account=cweidman0
#SBATCH --time=4:00:00
#SBATCH --mem=20000m
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=1
################################################################################
# Edit the strings under 'declare -a tasks=(' to match your experiments.
#
# The #SBATCH --array variable above creates an array [0,1,2,3]. Change it so that the length
# is how many jobs you need (same as number of strings under $tasks).
#
# This script is submitted that many times, but only one line from $tasks is
# evaluated each time.
#
# For more info on #SBATCH variables, see https://arc.umich.edu/greatlakes/slurm-user-guide/
# and https://slurm.schedmd.com/sbatch.html
#
# This requires a conda environment with samtools and bowtie2 (RNA-SEQ-2), run create_env.sbatch 
# before running this file 
# 
# To call this script:
# sbatch remove_contaminants.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda activate ~/miniconda3/envs/RNA-SEQ-2

declare -a tasks=(
"python3 rm_contam.py -u --input 7KO-Cyto-BS_processed_fastqs --output 7KO-Cyto-BS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input 7KO-Cyto-NBS_processed_fastqs --output 7KO-Cyto-NBS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input 7KO-Nuc-BS_processed_fastqs --output 7KO-Nuc-BS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input 7KO-Nuc-NBS_processed_fastqs --output 7KO-Nuc-NBS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input 7LKO-Cyto-BS_processed_fastqs --output 7LKO-Cyto-BS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input 7LKO-Cyto-NBS_processed_fastqs --output 7LKO-Cyto-NBS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input 7LKO-Nuc-BS_processed_fastqs --output 7LKO-Nuc-BS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input 7LKO-Nuc-NBS_processed_fastqs --output 7LKO-Nuc-NBS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input WT-Cyto-BS_processed_fastqs --output WT-Cyto-BS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input WT-Cyto-NBS_processed_fastqs --output WT-Cyto-NBS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input WT-Nuc-BS_processed_fastqs --output WT-Nuc-BS_filtered_processed_fastqs"
"python3 rm_contam.py -u --input WT-Nuc-NBS_processed_fastqs --output WT-Nuc-NBS_filtered_processed_fastqs"
)
eval ${tasks[$SLURM_ARRAY_TASK_ID]}