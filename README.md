# **Removing contaminant human RNAs from processed fastp files**
## Necessary files
* `create_env.sbatch`
* bowtie2 index
  * OPTION 1: Copy over `build_index.py` and `build_index.sbatch` and run the SBATCH file.
  * OPTION 2: Skip the bowtie2 index creation step by navigating into the `contaminants_index` folder and moving all `*.bt2` files into the parent directory.
* `contaminants.fa`
* `rm_contam.py` and `rm_contam.sbatch`
## Instructions
1. Run `create_env.sbatch` to create the RNA-SEQ-2 conda environment
3. Activate conda environment via `conda activate RNA-SEQ-2`
4. Run `build_index.sbatch` one time to build contaminant index with bowtie2
5. Run `rm_contam.sbatch` to separate out contaminant RNAs from each fastq file
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
## Using the test data
This allows you to test the script on example data.
1. Navigate to `test` folder
2. Run `mv 7KO-Cyto-BS_processed_fastqs ../` and navigate back to parent folder
3. Follow Steps 1-3 from "Instructions" section (see above)
   * If running code in local Linux environment (_i.e._, WSL), copy/paste code from .sbatch files instead of directly running them
4. Run `python3 rm_contam.py -u --input 7KO-Cyto-BS_processed_fastqs --output 7KO-Cyto-BS_filtered_processed_fastqs`
