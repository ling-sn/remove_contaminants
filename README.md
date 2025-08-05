## <ins>**Removing contaminant human RNAs from processed fastqs**</ins>
### Necessary files
<img src="https://github.com/user-attachments/assets/54af626d-27a8-49e1-9822-62c617c1f937" width="400"/>

* bowtie2 index
  * **Option 1 (Manual):** First, copy `contaminants.fa` to the working directory containing your data folders. Then, copy over `build_index.py` and `build_index.sbatch` and run the SBATCH file.
  * **Option 2 (Pre-Built; Recommended):** Skip bowtie2 index creation by navigating into the `contaminants_index` folder and moving all `*.bt2` files into the working directory.
* `create_env.sbatch`
* `rm_contam.py` and `rm_contam.sbatch`
### Instructions
1. Run `create_env.sbatch` to create the RNA-SEQ-2 conda environment
3. Activate conda environment via `conda activate RNA-SEQ-2`
4. Create bowtie2 index (choose manual vs. pre-built method)
5. Edit `rm_contam.sbatch` to match your experiments
   * Change the following:
     * `#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu`
     * `#SBATCH --array=0-11%4`
     * `#SBATCH --time=4:00:00`
     * Strings under `declare -a tasks=(`
6. Run `rm_contam.sbatch` to separate out contaminant RNAs from each fastq file
### Tools used in contaminant removal script
* **bowtie2** is used to align single-end reads (merged/unpaired) and paired reads (unmerged) according to the contaminants index
* **samtools** is used to convert the .sam outputs from bowtie2 into compressed .bam and .bai files
### When do I use this script?
* This script is used after running `run_fastp` on your raw data.
### What does this script do?
* The resulting merged, unmerged, and unpaired fastqs from `run_fastp` are aligned to a contaminants FASTA using bowtie2.
  * The output folder from this process is split into the subfolders `mapped_contam` and `rm_contam`, which respectively group (i) the fastqs mapped to the contaminants and (ii) the fastqs with removed contaminants. The latter will be used for downstream processing.
### Understanding the rm_contam SBATCH
```
python3 rm_contam.py -u --folder_name 7KO-Cyto-BS_processed_fastqs -B
```
* **--folder_name:** Name of folder containing merged, paired, and unpaired fastqs. DO NOT INPUT A PATH.
* **-B:** Creates a .bam file for mapped contaminants, which can be used for verification/sanity checks in IGV. The default option is `--no-samfile`, which is the same as not including this option (_i.e._, if you only specify `--folder_name` and not `-B`).
### Additional information about contaminants FASTA
`contaminants.fa` contains human rRNA, tRNA, snoRNA, and snRNA sequences sourced from the anRNAlab genomic database and various online sources.
 * rRNA online source: [fallerlab](https://github.com/fallerlab/ARF/blob/main/rRNAs/4V6X_human_rRNAs.fa)
 * tRNA online source: [GtRNAdb](https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-seq.html)
 * snoRNA online source: [snoRNABase](https://www-snorna.biotoul.fr/browse.php)
### Citations
* Zhang et al. BID-seq for transcriptome-wide quantitative sequencing of mRNA pseudouridine at base resolution. _Nature Protocols_ 19, 517–538 (2024). https://doi.org/10.1038/s41596-023-00917-5
