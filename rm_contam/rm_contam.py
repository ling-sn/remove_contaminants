import subprocess
from pathlib import Path
import traceback
import argparse

class Bowtie2Aligner:
    def __init__(self, current_path):
        self.bowtie2_index = current_path/"contaminants_index"
        self.r1_filename = None
        self.r2_filename = None
        
    def single_reads(self, bamfile, file, mapped_folder, unmapped_folder, samtools_folder):
        """
        Align single-end reads (merged/unpaired)
        """
        filename = file.name.split(".")[0]
        
        try:
            sam_output = samtools_folder/f"{filename}_mapped.sam"
            contam_output = mapped_folder/f"{filename}_mapped.fastq.gz"
            rmcontam_output = unmapped_folder/f"{filename}_unmapped.fastq.gz"

            cmd = ["bowtie2", 
                    "-x", str(self.bowtie2_index), ## provides index that we created earlier
                    "-U", str(file), ## specifies single-read
                    "--un-gz", str(rmcontam_output), ## merged/unpaired reads that failed to align (mrna)
                    "--al-gz", str(contam_output)] ## merged/unpaired reads that aligned ≥1 times
            if bamfile:
                cmd.extend(["-S", str(sam_output)]) ## .sam file of contam reads
                
            result = subprocess.run(cmd, 
                                    check = True, 
                                    capture_output = True, 
                                    text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to align merged or unpaired fastq file {file.name}: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise
        return result

    def detect_reps(self, subfolder):
        """
        Defines unmerged R1 and R2 filename variables
        """
        for r1_file in subfolder.glob("*unmerged_R1*"):
            self.r1_filename = r1_file
            self.r2_filename = r1_file.with_name(r1_file.name.replace("_R1_", "_R2_"))

    def paired_reads(self, bamfile, file, mapped_folder, unmapped_folder, samtools_folder):
        """
        Align paired-end reads (unmerged)
        """
        filename = file.name.split(".")[0]

        try:
            sam_output = samtools_folder/f"{filename}_mapped.sam"
            contam_output = mapped_folder/f"{filename}_mapped.fastq.gz"
            rmcontam_output = unmapped_folder/f"{filename}_unmapped.fastq.gz"
            self.detect_reps(file.parent)

            cmd = ["bowtie2",
                    "-x", str(self.bowtie2_index),
                    "-1", str(self.r1_filename),
                    "-2", str(self.r2_filename),
                    "--un-conc-gz", str(rmcontam_output), ## unmerged reads that failed to align (mrna)
                    "--al-conc-gz", str(contam_output)] ## unmerged reads that aligned ≥1 times
            if bamfile:
                cmd.extend(["-S", str(sam_output)]) ## .sam file of contam reads

            result = subprocess.run(cmd, 
                                    check = True, 
                                    capture_output = True, 
                                    text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to align unmerged fastq file {file.name}: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise
        return result
    
    def convert_sam(self, samtools_folder, file):
        """
        Converts .sam output from bowtie2 into .bam
        """
        filename = file.name.split(".")[0]
        sam_input = samtools_folder/f"{filename}_mapped.sam"
        sam_filename = Path(sam_input).stem
        bam_output = samtools_folder/f"{sam_filename}_out.bam"

        try:
            subprocess.run(["samtools", "sort", "-O", "BAM", ## convert .sam to .bam and sort
                            "-o", str(bam_output), ## output file name
                            str(sam_input)], ## input file name
                            check = True,
                            capture_output = True,
                            text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to convert {sam_input.name} to .bam: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise

    def merge_bam(self, samtools_folder, file):
        """
        Merges all .bam files, then 
        sorts and indexes into .bai
        """
        base_name = file.stem.split("_")[0] # create general filestem
        merged_bam = samtools_folder/f"{base_name}_mapped.bam" # merged output
        bam_list = [*samtools_folder.glob("*.bam")] # detect .bam files
        rm_list = [*samtools_folder.glob("*out.bam"), *samtools_folder.glob("*.sam")]

        try:
            subprocess.run(["samtools", "merge", ## merge all .bam files into one
                            str(merged_bam), *map(str, bam_list)], ## asterisk = expand list & iterate through
                            check = True, 
                            capture_output = True,
                            text = True)
            subprocess.run(["samtools", "index", str(merged_bam)], ## create .bai from .bam
                            check = True,
                            capture_output = True,
                            text = True)
            subprocess.run(["rm", *map(str, rm_list)], ## remove original .sam and .bam files
                            check = True, ## ensures that this block only runs if previous 2 were successful
                            capture_output = True,
                            text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to create {merged_bam.name} and convert to .bai: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise

def rmcontam_pipeline(folder_name, bamfile):
    ## define input directory
    current_path = Path.cwd()
    input_dir = current_path/folder_name

    ## initialize class
    aligner = Bowtie2Aligner(current_path)

    for subfolder in input_dir.iterdir(): ## amount of subfolders = number of replicates
        if subfolder.is_dir():
            mapped_folder = current_path/"filtered_processed_fastqs"/"mapped_contam"/folder_name/f"{subfolder.name}"
            mapped_folder.mkdir(exist_ok=True, parents=True)

            unmapped_folder = current_path/"filtered_processed_fastqs"/"removed_contam"/folder_name/f"{subfolder.name}"
            unmapped_folder.mkdir(exist_ok=True, parents=True)

            if bamfile:
                samtools_folder = mapped_folder/"samtools"
                samtools_folder.mkdir(exist_ok=True, parents=True)

            for file in subfolder.glob("*.fastq.gz"): ## iterate through indiv. files in subfolder
                try:
                    ## run bowtie2 alignment functions
                    if "_merged" in file.name or "_unpaired" in file.name:
                        aligner.single_reads(bamfile, file, mapped_folder, unmapped_folder, samtools_folder)
                    elif "_unmerged" in file.name:
                        aligner.paired_reads(bamfile, file, mapped_folder, unmapped_folder, samtools_folder)

                    ## if user requests it, then a BAM file will be outputted
                    if bamfile:
                        ## run samtools function
                        aligner.convert_sam(samtools_folder, file)
                except Exception as e:
                    print(f"Failed to align {file.name} with bowtie2: {e}")
                    traceback.print_exc()
                    continue
            
            if bamfile:
                ## merge bam files, convert to bai, & remove files
                aligner.merge_bam(samtools_folder, file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Run contaminant removal pipeline.")
    parser.add_argument("--folder_name", required = True, help = "Input processed_fastqs folder name")
    parser.add_argument("-B", "--bamfile", action = argparse.BooleanOptionalAction, default = False, help = "Creates BAM file for mapped contaminants; disabled by default")

    args = parser.parse_args()

    print("Starting contaminant removal pipeline...")
    rmcontam_pipeline(args.folder_name, args.bamfile)
    print("Pipeline finished.")