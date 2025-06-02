import os
import subprocess
import shutil
from pathlib import Path
import traceback

class Bowtie2Aligner:
    def __init__(self, folder_path):
        self.parent_path = Path(folder_path).parent
        self.contaminants_dir = self.parent_path/"contaminants.fa"
        self.bowtie2_index = self.parent_path/"contaminants_index"
        self.grouped_files = {}

    def build_bowtie2_index(self):
        """
        Must have contaminants.fa in parent directory of folder_path.
            folder_path = Folder name that ends with 'processed_fastqs'
        """
        if not self.bowtie2_index.exists():
            try:
                cmd = ["bowtie2-build", "-f", 
                       str(self.contaminants_dir), 
                       str(self.bowtie2_index)]
                result = subprocess.run(cmd, 
                                        check = True, ## if command returns non-zero exit status, raise error
                                        capture_output = True, 
                                        text = True)
            except subprocess.CalledProcessError as e: ## error handling
                    print(f"Failed to build bowtie2 index: {e}")
                    print("STDERR:", e.stderr)
                    print("STDOUT:", e.stdout)
                    traceback.print_exc()
                    raise
        else:
            pass
        return result
        
    def single_reads(self, file, output_dir):
        """
        Align single-end reads (merged/unpaired)
        """
        sam_output = output_dir/f"{file.stem}_aligned_contam.sam" ## file.stem = og filename w/o extension
        rmcontam_output = output_dir/f"{file.stem}_removed_contam.fastq.gz"
        contam_output = output_dir/f"{file.stem}_aligned_contam.fastq.gz"

        try:
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
        grouped_files = dictionary used in paired_reads()
        base_name = replicates are grouped by sample name
        """
        for r1_file in subfolder.glob("*_R1_*"):
            base_name = r1_file.name.split("_R1_")[0] 
            r1_filename = r1_file
            r2_filename = r1_file.replace("_R1_", "_R2_")
            
            self.grouped_files[base_name] = (r1_filename, r2_filename) ## define dictionary

    def paired_reads(self, file, output_dir):
        """
        Align paired-end reads (unmerged)
        """
        sam_output = output_dir/f"{file.stem}_aligned_contam.sam" ## file.stem = og filename w/o extension
        rmcontam_output = output_dir/f"{file.stem}_removed_contam_%.fastq.gz" ## include % for R1/R2
        contam_output = output_dir/f"{file.stem}_aligned_contam_%.fastq.gz" ## include % for R1/R2

        self.detect_reps(file.parent) ## group files according to dictionary
        for base_name, (r1_filename, r2_filename) in self.grouped_files.items(): ## for [key, value] in dictionary
            try: 
                cmd = ["bowtie2",
                        "-x", str(self.bowtie2_index),
                        "-1", str(r1_filename),
                        "-2", str(r2_filename),
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

def rmcontam_pipeline(folder_path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path) ## if output_path doesn't already exist, create new empty output folder

    input_dir = Path(folder_path) ## imported Path to improve readability of code below
    output_dir = Path(output_path)
    aligner = Bowtie2Aligner(input_dir) ## initialize class
    
    aligner.build_bowtie2_index() ## build bowtie2 index

    for subfolder in input_dir.iterdir(): ## amount of subfolders = number of replicates
        if subfolder.is_dir():
            processed_folder = output_dir/f"{subfolder.name}_bowtie2_out" ## bowtie2 outputs
            processed_folder.mkdir(exist_ok=True) ## if directory already exists, suppress OSError

            for file in subfolder.glob("*.fastq.gz"):
                try:
                    if "_merged" in file.name or "_unpaired" in file.name:
                        aligner.single_reads(file, processed_folder)
                    elif "_unmerged" in file.name:
                        aligner.paired_reads(file, processed_folder)
                    
                except Exception as e:
                    print(f"Failed to align {file.name} with bowtie2: {e}")
                    traceback.print_exc()
                   
