import os
import subprocess
from pathlib import Path
import traceback
import argparse
import re

class Bowtie2Aligner:
    def __init__(self):
        self.current_path = Path.cwd()
        self.contaminants_dir = self.current_path/"contaminants.fa"
        self.bowtie2_index = self.current_path/"contaminants_index"
        self.r1_filename = None
        self.r2_filename = None
        
    def single_reads(self, file, processed_folder, samtools_folder):
        """
        Align single-end reads (merged/unpaired)
        """
        try:
            sam_output = samtools_folder/f"{file.stem}_mapped.sam" ## file.stem = filename w/o both extensions
            rmcontam_output = processed_folder/f"{file.stem}_unmapped.fastq.gz"
            contam_output = processed_folder/f"{file.stem}_mapped.fastq.gz"
            cmd = ["bowtie2", 
                    "-x", str(self.bowtie2_index), ## provides index that we created earlier
                    "-U", str(file), ## specifies single-read
                    "-S", str(sam_output), ## .sam file of contam reads
                    "--un-gz", str(rmcontam_output), ## merged/unpaired reads that failed to align (mrna)
                    "--al-gz", str(contam_output)] ## merged/unpaired reads that aligned ≥1 times
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

    def paired_reads(self, file, processed_folder, samtools_folder):
        """
        Align paired-end reads (unmerged)
        """
        try:
            sam_output = samtools_folder/f"{file.stem}_mapped.sam" ## file.stem = filename w/o both extensions
            rmcontam_output = processed_folder/f"{file.stem}_unmapped.fastq.gz"
            contam_output = processed_folder/f"{file.stem}_mapped.fastq.gz"
            self.detect_reps(file.parent)
            cmd = ["bowtie2",
                    "-x", str(self.bowtie2_index),
                    "-1", str(self.r1_filename),
                    "-2", str(self.r2_filename),
                    "-S", str(sam_output), ## .sam file of contam reads
                    "--un-conc-gz", str(rmcontam_output), ## unmerged reads that failed to align (mrna)
                    "--al-conc-gz", str(contam_output)] ## unmerged reads that aligned ≥1 times
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
    
class SamtoolsConversion:
    def __init__(self, samtools_folder, subfolder):
        self.base_name = subfolder.stem.split("_")[0] # create general filestem
        self.merged_bam = samtools_folder/f"{self.base_name}_mapped.bam" # merged output
        self.open_file = ["samtools", "view", "-H", str(self.merged_bam)]

    def convert_sam(self, samtools_folder, file):
        """
        Converts .sam output from bowtie2 into .bam
        """
        sam_input = samtools_folder/f"{file.stem}_mapped.sam"
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

    def check_regex(self):
        """
        Checks if headers in .bam have correct regex
        """
        header = subprocess.run(self.open_file,
                                check = True,
                                capture_output = True,
                                text = True).stdout
        invalid_regex = re.compile(r'@SQ.*SN:[^A-Za-z0-9!#$%&*+./:;=?@^_|~-]')
        return bool(invalid_regex.search(header))
    
    def sanitize_header(self):
        """
        Sanitize .bam headers if incorrect regex detected
        """
        if not self.check_regex():
            print("No invalid characters detected in .bam headers.")
            return
        
        else:
            try:
                p1 = subprocess.Popen(self.open_file, 
                                      stdout=subprocess.PIPE)
                p2 = subprocess.Popen(["sed", "-E", 
                                      r's/SN:([^A-Za-z0-9!#$%&*+.\/:;=?@^_|~-])//g'], 
                                      stdin=p1.stdout, 
                                      stdout=subprocess.PIPE)
                p3 = subprocess.Popen(["samtools", "reheader", "-", 
                                       str(self.merged_bam)], 
                                       stdin=p2.stdout)
                p3.communicate()
            except subprocess.CalledProcessError as e: ## error handling
                print(f"Failed to modify regex in {self.merged_bam.name}: {e}")
                print("STDERR:", e.stderr)
                print("STDOUT:", e.stdout)
                traceback.print_exc()
                raise

    def merge_bam(self, samtools_folder):
        """
        Merges all .bam files, then 
        sorts and indexes into .bai
        """
        bam_list = [*samtools_folder.glob("*.bam")] # detect .bam files
        rm_list = [*samtools_folder.glob("*out.bam"), *samtools_folder.glob("*.sam")]

        try:
            subprocess.run(["samtools", "merge", ## merge all .bam files into one
                            str(self.merged_bam), *map(str, bam_list)], ## asterisk = expand list & iterate through
                            check = True, 
                            capture_output = True,
                            text = True)
            self.sanitize_header(samtools_folder)
            subprocess.run(["samtools", "index", str(self.merged_bam)], ## create .bai from .bam
                            check = True,
                            capture_output = True,
                            text = True)
            subprocess.run(["rm", *map(str, rm_list)], ## remove original .sam and .bam files
                            check = True, ## ensures this only runs if previous blocks were successful
                            capture_output = True,
                            text = True)

        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to create {self.merged_bam.name} and convert to .bai: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise

def rmcontam_pipeline(folder_path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path) ## if output_path doesn't already exist, create new empty output folder

    input_dir = Path(folder_path) ## used Path to improve readability of code below
    output_dir = Path(output_path)

    ## initialize class
    aligner = Bowtie2Aligner()

    for subfolder in input_dir.iterdir(): ## amount of subfolders = number of replicates
        if subfolder.is_dir():
            processed_folder = output_dir/f"{subfolder.name}_bowtie2" ## processed folders for bowtie2 outputs
            processed_folder.mkdir(exist_ok=True) ## if directory already exists, suppress OSError
            samtools_folder = processed_folder/"samtools"
            samtools_folder.mkdir(exist_ok=True)
            samtools = SamtoolsConversion(samtools_folder, subfolder)

            for file in subfolder.glob("*.fastq.gz"): ## iterate through indiv. files in subfolder
                try:
                    ## run bowtie2 alignment functions
                    if "_merged" in file.name or "_unpaired" in file.name:
                        aligner.single_reads(file, processed_folder, samtools_folder)
                    elif "_unmerged" in file.name:
                        aligner.paired_reads(file, processed_folder, samtools_folder)

                    ## run samtools function
                    samtools.convert_sam(samtools_folder, file)
                    
                except Exception as e:
                    print(f"Failed to align {file.name} with bowtie2 and produce .bam files: {e}")
                    traceback.print_exc()
                    continue
            
            ## merge bam files, convert to bai, & remove files
            samtools.merge_bam(samtools_folder)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Run contaminant removal pipeline.")
    parser.add_argument("--input", required = True, help = "Input processed_fastqs folder name")
    parser.add_argument("--output", required = True, help = "Desired output folder name for processed files")
    parser.add_argument("-u", "--unzip", action = "store_true", help = "Unzips fastq files")

    args = parser.parse_args()

    print("Starting contaminant removal pipeline...")
    rmcontam_pipeline(args.input, args.output)
    print("Pipeline finished.")