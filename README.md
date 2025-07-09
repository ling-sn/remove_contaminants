# **Removing contaminant human RNAs from processed fastp files**
## Necessary files
<img src="https://github.com/user-attachments/assets/36c3a2b7-2e37-4377-b131-e5fba40acf9e" width="600"/>

* bowtie2 index
  * **Option 1 (Manual):** First, copy `contaminants.fa` to the working directory containing your data folders. Then, copy over `build_index.py` and `build_index.sbatch` and run the SBATCH file.
  * **Option 2 (Pre-Built):** Skip bowtie2 index creation by navigating into the `contaminants_index` folder and moving all `*.bt2` files into the working directory.
* `create_env.sbatch`
* `rm_contam.py` and `rm_contam.sbatch`
## Instructions
1. Run `create_env.sbatch` to create the RNA-SEQ-2 conda environment
3. Activate conda environment via `conda activate RNA-SEQ-2`
4. Create bowtie2 index (choose manual vs. pre-built method)
5. Edit `rm_contam.sbatch` to match your experiments
6. Run `rm_contam.sbatch` to separate out contaminant RNAs from each fastq file
## Tools used in contaminant removal script
* **bowtie2** is used to align single-end reads (merged/unpaired) and paired reads (unmerged) according to the contaminants index
* **samtools** is used to convert the .sam outputs from bowtie2 into compressed .bam and .bai files
## When do I use this pipeline?
This is applied after running the fastp script on your raw data (fastq files). 
## Additional information about contaminants FASTA
`contaminants.fa` contains human rRNA, tRNA, snoRNA, and snRNA sequences sourced from the anRNAlab genomic database and various online sources.
 * rRNA online source: [fallerlab](https://github.com/fallerlab/ARF/blob/main/rRNAs/4V6X_human_rRNAs.fa)
 * tRNA online source: [GtRNAdb](https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-seq.html)
 * snoRNA online source: [snoRNABase](https://www-snorna.biotoul.fr/browse.php)
