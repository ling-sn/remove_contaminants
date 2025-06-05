# **Extracting contaminant human RNAs from processed fastp files**
## Instructions
1. Run `create_env.sbatch` to create the RNA-SEQ-2 conda environment
3. Activate conda environment via `conda activate RNA-SEQ-2`
4. Run `build_index.sbatch` one time to build bowtie2 index for contaminants
5. Run `rm_contam.sbatch` to separate out contaminant RNAs from each fastq file
## Tools used in contaminant removal script
* **bowtie2** is used to align single-end reads (merged/unpaired) and paired reads (unmerged) according to the contaminants index
* **samtools** is used to convert the .sam outputs from bowtie2 into compressed .bam and .bai files
## When do I use this pipeline?
This is applied after running the fastp script on your raw data (fastq files). 
## Additional information about contaminants FASTA
* Contains human rRNA, tRNA, snoRNA, and snRNA sequences sourced from anRNAlab genomic database and various online databases
  * rRNA online database: [GitHub](https://github.com/fallerlab/ARF/blob/main/rRNAs/4V6X_human_rRNAs.fa)
  * tRNA online database: [GtRNAdb](https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-seq.html)
  * snoRNA online database: [snoRNABase](https://www-snorna.biotoul.fr/browse.php)
