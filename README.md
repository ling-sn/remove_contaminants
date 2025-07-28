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
6. Run `rm_contam.sbatch` to separate out contaminant RNAs from each fastq file
### Tools used in contaminant removal script
* **bowtie2** is used to align single-end reads (merged/unpaired) and paired reads (unmerged) according to the contaminants index
* **samtools** is used to convert the .sam outputs from bowtie2 into compressed .bam and .bai files
### When do I use this pipeline?
This is applied after running the fastp script on your raw data (fastq files). 
### Understanding the rm_contam SBATCH
```
python3 rm_contam.py -u --input 7KO-Cyto-BS_processed_fastqs --output 7KO-Cyto-BS_filtered_processed_fastqs
```
* **--input:** Name of folder containing merged, paired, and unpaired fastqs. DO NOT INPUT A PATH.
* **--output:** Name of desired output folder.
### Additional information about contaminants FASTA
`contaminants.fa` contains human rRNA, tRNA, snoRNA, and snRNA sequences sourced from the anRNAlab genomic database and various online sources.
 * rRNA online source: [fallerlab](https://github.com/fallerlab/ARF/blob/main/rRNAs/4V6X_human_rRNAs.fa)
 * tRNA online source: [GtRNAdb](https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-seq.html)
 * snoRNA online source: [snoRNABase](https://www-snorna.biotoul.fr/browse.php)
### Citations
* Zhang et al. BID-seq for transcriptome-wide quantitative sequencing of mRNA pseudouridine at base resolution. _Nature Protocols_ 19, 517–538 (2024). https://doi.org/10.1038/s41596-023-00917-5
