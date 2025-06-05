import subprocess
from pathlib import Path
import traceback

def build_bowtie2_index(folder_path):
    """
    For this script, you only need to list one (1) folder_path
    because all folders use the same bowtie2 index
    
    folder_path = Folder name that ends with 'processed_fastqs'
    """
    parent_path = Path(folder_path).parent
    contaminants_dir = parent_path/"contaminants.fa"
    bowtie2_index = parent_path/"contaminants_index"
    bt2_files = list(parent_path.glob("*.bt2")) ## produces list of files
    
    if not bt2_files: ## checks if list is empty; if so, proceed
        try:
            cmd = ["bowtie2-build",
                    str(contaminants_dir),
                    str(bowtie2_index)]
            result = subprocess.run(cmd, 
                                    check = True, ## if command returns non-zero exit status, raise error
                                    capture_output = True, 
                                    text = True)
            return result
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to build bowtie2 index: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise