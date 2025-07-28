#!/usr/bin/env bash
#SBATCH --job-name=create_env
#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/dev/null
#SBATCH --account=cweidman0
#SBATCH --time=10:00
#SBATCH --mem=16000m
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
# This requires a conda environment for genome alignment, edit to your named version in the activate command
# 
# To call this script:
# sbatch remove_contaminants.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda create -n RNA-SEQ-2 python=3.10 r-base=4.3.1 hisat2=2.2.1 samtools=1.17 Bioconductor-deseq2=1.40.2 subread=2.0.6 matplotlib=3.7.2 numpy=1.25.1 pandas=2.0.3 scipy=1.11.1 seaborn=0.12.2 scikit-learn=1.3.0 rpy2=3.5.11 bowtie2=2.5.1
