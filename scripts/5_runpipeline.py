#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path
import logging
import datetime

# --- Configuration Constants ---
DEFAULT_OUTPUT_ROOT = "PairTCR_results"
DEFAULT_INPUT_DIR = "raw"
DEFAULT_PREFIX = "TCR_TSO_18"
DEFAULT_READ_LIMIT = 100000
DEFAULT_THREADS = 90

# Detect MiXCR JAR in scripts directory by default
DEFAULT_MIXCR_JAR = os.path.join("scripts", "mixcr.jar")

def get_scripts_directory():
    """Get the directory containing the scripts, trying multiple possible locations."""
    # Try to get it from the package installation first
    try:
        import pairtcr
        package_dir = os.path.dirname(pairtcr.__file__)
        scripts_dir_package = os.path.join(package_dir, 'scripts')
        if os.path.isdir(scripts_dir_package):
            return scripts_dir_package
    except ImportError:
        pass
    
    # Fallback 1: relative to current script location
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    if os.path.basename(current_script_dir) == 'scripts':
        # We're running from the scripts directory
        return current_script_dir
    
    # Fallback 2: scripts directory relative to current working directory
    scripts_dir_relative = os.path.join(os.getcwd(), 'scripts')
    if os.path.isdir(scripts_dir_relative):
        return scripts_dir_relative
    
    # Fallback 3: assume scripts directory in same directory as this script
    scripts_dir_sibling = os.path.join(os.path.dirname(current_script_dir), 'scripts')
    if os.path.isdir(scripts_dir_sibling):
        return scripts_dir_sibling
    
    # If all else fails, return the relative path and let the caller handle the error
    return "scripts"

class PipelineRunner:
    def __init__(self, input_dir, output_root, prefix, read_limit, threads, mixcr_jar, force_restart=False, use_c_version=False):
        self.input_dir = input_dir
        self.output_root = output_root
        self.prefix = prefix
        self.read_limit = read_limit
        self.threads = threads
        self.mixcr_jar = mixcr_jar
        self.force_restart = force_restart
        self.use_c_version = use_c_version
        
        # Get the correct scripts directory
        self.scripts_dir = get_scripts_directory()
        
        # Define all output directories
        self.step1_output = os.path.join(output_root, "1_preprocess_and_trim_output")
        self.step2_output = os.path.join(output_root, "2_create_umi_pairs_output")
        self.matched_fastq_output = os.path.join(output_root, "matched_fastq_output")
        self.step3_output = os.path.join(output_root, "3_run_mixcr_and_export_output")
        self.step4_output = os.path.join(output_root, "4_pair_and_filter_clones_output")
        self.logs_output = os.path.join(output_root, "logs")
        
        # Define key output files
        self.umi_pairs_file = os.path.join(self.step2_output, "umi_pairs.tsv")
        self.final_output = os.path.join(self.step4_output, "final_paired_clones_filtered.tsv")
        
        # Setup logging
        self.setup_logging()
        
        # Define step completion markers
        self.step_markers = {
            'step1': os.path.join(self.step1_output, f"{self.prefix}_TRA_1.fq.gz"),
            'step2': self.umi_pairs_file,
            'step2.5': os.path.join(self.matched_fastq_output, f"{self.prefix}_matched_TRA_matched_1.fq.gz"),
            'step3': os.path.join(self.step3_output, "TRA_alignments_export_with_headers.tsv"),
            'step4': self.final_output
        }

    def setup_logging(self):
        """Setup logging configuration."""
        # Create logs directory
        os.makedirs(self.logs_output, exist_ok=True)
        
        # Create timestamp for this run
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Setup main pipeline log
        self.pipeline_log = os.path.join(self.logs_output, f"pipeline_{timestamp}.log")
        
        # Create individual step log files
        self.step_logs = {
            'step1': os.path.join(self.logs_output, f"step1_preprocess_{timestamp}.log"),
            'step2': os.path.join(self.logs_output, f"step2_umi_pairs_{timestamp}.log"),
            'step2.5': os.path.join(self.logs_output, f"step2.5_matched_fastq_{timestamp}.log"),
            'step3': os.path.join(self.logs_output, f"step3_mixcr_{timestamp}.log"),
            'step4': os.path.join(self.logs_output, f"step4_pair_filter_{timestamp}.log")
        }
        
        # Configure main logger
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.pipeline_log),
                logging.StreamHandler(sys.stdout)
            ]
        )
        
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Pipeline logging initialized. Main log: {self.pipeline_log}")

    def check_step_completion(self, step_key):
        """Check if a step has been completed successfully."""
        if step_key not in self.step_markers:
            return False
        
        marker_file = self.step_markers[step_key]
        completed = os.path.exists(marker_file) and os.path.getsize(marker_file) > 0
        
        if completed:
            self.logger.info(f"Step {step_key} already completed. Marker file: {marker_file}")
        else:
            self.logger.info(f"Step {step_key} not completed. Marker file missing or empty: {marker_file}")
        
        return completed

    def check_pipeline_completion(self):
        """Check if the entire pipeline has been completed successfully."""
        all_completed = all(self.check_step_completion(step) for step in self.step_markers.keys())
        
        if all_completed:
            self.logger.info("Pipeline appears to have been completed successfully")
        else:
            self.logger.info("Pipeline is incomplete")
        
        return all_completed

    def cleanup_incomplete_run(self):
        """Clean up any incomplete pipeline run."""
        self.logger.info("Cleaning up incomplete pipeline run...")
        
        # List of directories to clean up
        cleanup_dirs = [
            self.step1_output,
            self.step2_output,
            self.matched_fastq_output,
            self.step3_output,
            self.step4_output
        ]
        
        for cleanup_dir in cleanup_dirs:
            if os.path.exists(cleanup_dir):
                try:
                    shutil.rmtree(cleanup_dir)
                    self.logger.info(f"Removed directory: {cleanup_dir}")
                    print(f"Cleaned up directory: {cleanup_dir}")
                except Exception as e:
                    self.logger.warning(f"Failed to remove directory {cleanup_dir}: {e}")
                    print(f"Warning: Failed to remove directory {cleanup_dir}: {e}")
        
        self.logger.info("Cleanup completed")

    def run_command(self, cmd, step_name, shell=False, step_key=None, show_progress=True):
        """Run a command and handle errors."""
        self.logger.info(f"Starting {step_name}")
        self.logger.info(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
        
        print(f"\n{'='*50}")
        print(f"Running {step_name}")
        print(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
        print(f"{'='*50}")
        
        # Determine log file for this step
        log_file = None
        if step_key and step_key in self.step_logs:
            log_file = self.step_logs[step_key]
        
        try:
            if log_file and not show_progress:
                # Redirect output to log file only (no console output)
                with open(log_file, 'w') as log_f:
                    if shell:
                        result = subprocess.run(cmd, shell=True, check=True, 
                                              stdout=log_f, stderr=subprocess.STDOUT)
                    else:
                        result = subprocess.run(cmd, check=True, 
                                              stdout=log_f, stderr=subprocess.STDOUT)
                self.logger.info(f"{step_name} completed successfully. Log saved to: {log_file}")
            else:
                # Show output on console (for progress bars) and also log to file
                if log_file:
                    # Use tee-like behavior: show on console and save to log
                    import subprocess
                    if shell:
                        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                                                 stderr=subprocess.STDOUT, universal_newlines=True)
                    else:
                        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                                                 stderr=subprocess.STDOUT, universal_newlines=True)
                    
                    # Read output line by line and display + log
                    with open(log_file, 'w') as log_f:
                        for line in iter(process.stdout.readline, ''):
                            print(line.rstrip())  # Display on console
                            log_f.write(line)     # Write to log file
                            log_f.flush()
                        
                        process.stdout.close()
                        return_code = process.wait()
                        
                        if return_code != 0:
                            raise subprocess.CalledProcessError(return_code, cmd)
                else:
                    # No log file, just run normally
                    if shell:
                        result = subprocess.run(cmd, shell=True, check=True, 
                                              capture_output=False)
                    else:
                        result = subprocess.run(cmd, check=True, 
                                              capture_output=False)
                
                if log_file:
                    self.logger.info(f"{step_name} completed successfully. Log saved to: {log_file}")
                else:
                    self.logger.info(f"{step_name} completed successfully.")
            
            print(f"{step_name} completed successfully.")
            return True
        except subprocess.CalledProcessError as e:
            error_msg = f"{step_name} failed with return code {e.returncode}"
            self.logger.error(error_msg)
            self.logger.error(f"Command that failed: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
            print(f"Error: {error_msg}")
            print(f"Command that failed: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
            if log_file:
                print(f"Check log file for details: {log_file}")
            sys.exit(1)

    def step1_preprocess_and_trim(self):
        """Step 1: Preprocess and trim FASTQ files."""
        version_info = "C version" if self.use_c_version else "Python version"
        self.logger.info("="*50)
        self.logger.info(f"STEP 1: Starting preprocess and trim ({version_info})")
        self.logger.info("="*50)
        
        # Check if step is already completed
        if self.check_step_completion('step1'):
            print(f"Step 1: Preprocess and Trim ({version_info}) - SKIPPED (already completed)")
            return True
        
        if self.use_c_version:
            # Try to use C version
            c_executable = os.path.join(self.scripts_dir, "1_preprocess_and_trim")
            if not os.path.exists(c_executable):
                print(f"Warning: C executable not found at {c_executable}")
                print("Falling back to Python version for preprocessing...")
                self.logger.warning("C executable not found, falling back to Python version")
                self.use_c_version = False  # Switch to Python version
        
        if self.use_c_version:
            # Use C version (confirmed to exist)
            c_executable = os.path.join(self.scripts_dir, "1_preprocess_and_trim")
            cmd = [
                c_executable,
                self.input_dir,
                "-n", str(self.read_limit),
                "-o", self.prefix,
                "-d", self.step1_output
            ]
            step_name = "Step 1: Preprocess and Trim (C version)"
        else:
            # Use Python version
            script_path = os.path.join(self.scripts_dir, "1_preprocess_and_trim.py")
            cmd = [
                "python3", script_path,
                self.input_dir,
                "-n", str(self.read_limit),
                "-o", self.prefix,
                "-d", self.step1_output
            ]
            step_name = "Step 1: Preprocess and Trim (Python version)"
        
        return self.run_command(cmd, step_name, step_key='step1')

    def step2_create_umi_pairs(self):
        """Step 2: Create UMI pairs."""
        self.logger.info("="*50)
        self.logger.info("STEP 2: Starting UMI pairs creation")
        self.logger.info("="*50)
        
        # Check if step is already completed
        if self.check_step_completion('step2'):
            print("Step 2: Create UMI Pairs - SKIPPED (already completed)")
            return True
        
        script_path = os.path.join(self.scripts_dir, "2_create_umi_pairs.py")
        cmd = [
            "python3", script_path,
            "-i", self.step1_output,
            "-p", self.prefix,
            "-o", self.umi_pairs_file
        ]
        return self.run_command(cmd, "Step 2: Create UMI Pairs", step_key='step2')

    def step2_5_create_matched_fastq(self):
        """Step 2.5: Create matched FASTQ files from UMI pairs."""
        import pandas as pd
        import gzip
        from collections import defaultdict
        from tqdm import tqdm
        
        self.logger.info("="*50)
        self.logger.info("STEP 2.5: Starting matched FASTQ files creation")
        self.logger.info("="*50)
        
        # Check if step is already completed
        if self.check_step_completion('step2.5'):
            print("Step 2.5: Create Matched FASTQ Files - SKIPPED (already completed)")
            return True
        
        print("\n" + "="*50)
        print("Running Step 2.5: Create Matched FASTQ Files")
        print("="*50)
        
        # Setup step-specific logging
        step_log_file = self.step_logs['step2.5']
        step_logger = logging.getLogger('step2.5')
        step_handler = logging.FileHandler(step_log_file)
        step_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        step_logger.addHandler(step_handler)
        step_logger.setLevel(logging.INFO)
        
        # Create output directory
        os.makedirs(self.matched_fastq_output, exist_ok=True)
        step_logger.info(f"Created output directory: {self.matched_fastq_output}")
        
        # Read UMI pairs
        print(f"Reading UMI pairs from: {self.umi_pairs_file}")
        step_logger.info(f"Reading UMI pairs from: {self.umi_pairs_file}")
        umi_pairs_df = pd.read_csv(self.umi_pairs_file, sep='\t')
        
        # Create sets of read IDs for quick lookup
        tra_read_ids = set(umi_pairs_df['TRA_Read_ID_Base'])
        trb_read_ids = set(umi_pairs_df['TRB_Read_ID_Base'])
        
        print(f"Found {len(tra_read_ids)} unique TRA read IDs and {len(trb_read_ids)} unique TRB read IDs")
        step_logger.info(f"Found {len(tra_read_ids)} unique TRA read IDs and {len(trb_read_ids)} unique TRB read IDs")
        
        # Process TRA files
        tra_r1_in = os.path.join(self.step1_output, f"{self.prefix}_TRA_1.fq.gz")
        tra_r2_in = os.path.join(self.step1_output, f"{self.prefix}_TRA_2.fq.gz")
        tra_r1_out = os.path.join(self.matched_fastq_output, f"{self.prefix}_matched_TRA_matched_1.fq.gz")
        tra_r2_out = os.path.join(self.matched_fastq_output, f"{self.prefix}_matched_TRA_matched_2.fq.gz")
        
        self._filter_fastq_by_read_ids(tra_r1_in, tra_r1_out, tra_read_ids, "TRA R1")
        self._filter_fastq_by_read_ids(tra_r2_in, tra_r2_out, tra_read_ids, "TRA R2")
        
        # Process TRB files
        trb_r1_in = os.path.join(self.step1_output, f"{self.prefix}_TRB_1.fq.gz")
        trb_r2_in = os.path.join(self.step1_output, f"{self.prefix}_TRB_2.fq.gz")
        trb_r1_out = os.path.join(self.matched_fastq_output, f"{self.prefix}_matched_TRB_matched_1.fq.gz")
        trb_r2_out = os.path.join(self.matched_fastq_output, f"{self.prefix}_matched_TRB_matched_2.fq.gz")
        
        self._filter_fastq_by_read_ids(trb_r1_in, trb_r1_out, trb_read_ids, "TRB R1")
        self._filter_fastq_by_read_ids(trb_r2_in, trb_r2_out, trb_read_ids, "TRB R2")
        
        print("Step 2.5: Create Matched FASTQ Files completed successfully.")
        step_logger.info("Step 2.5: Create Matched FASTQ Files completed successfully.")
        self.logger.info(f"Step 2.5 completed. Log saved to: {step_log_file}")
        
        # Clean up step logger
        step_logger.removeHandler(step_handler)
        step_handler.close()
        
        return True

    def _filter_fastq_by_read_ids(self, input_file, output_file, target_read_ids, file_desc):
        """Filter FASTQ file to keep only reads with IDs in target_read_ids."""
        import re
        import gzip
        from itertools import islice
        
        def get_base_read_id(header):
            if not header:
                return None
            base_id_part = header.split(maxsplit=1)[0]
            if base_id_part.startswith('@'):
                base_id_part = base_id_part[1:]
            base_id_part = re.sub(r'/[12]$', '', base_id_part)
            return base_id_part
        
        def read_fastq_record(handle):
            lines = list(islice(handle, 4))
            if not lines or len(lines) < 4:
                return None
            return [line.strip() for line in lines]
        
        print(f"Filtering {file_desc}: {os.path.basename(input_file)} -> {os.path.basename(output_file)}")
        
        matched_reads = 0
        total_reads = 0
        
        with gzip.open(input_file, 'rt') as in_f, \
             gzip.open(output_file, 'wt') as out_f:
            
            while True:
                record = read_fastq_record(in_f)
                if record is None:
                    break
                
                total_reads += 1
                header, seq, plus, qual = record
                base_read_id = get_base_read_id(header)
                
                if base_read_id in target_read_ids:
                    out_f.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
                    matched_reads += 1
        
        print(f"  {file_desc}: {matched_reads}/{total_reads} reads matched and written")
        # Log to step logger if available
        step_logger = logging.getLogger('step2.5')
        if step_logger.handlers:
            step_logger.info(f"{file_desc}: {matched_reads}/{total_reads} reads matched and written")

    def step3_run_mixcr(self):
        """Step 3: Run MiXCR analysis."""
        self.logger.info("="*50)
        self.logger.info("STEP 3: Starting MiXCR analysis")
        self.logger.info("="*50)
        
        # Check if step is already completed
        if self.check_step_completion('step3'):
            print("Step 3: Run MiXCR - SKIPPED (already completed)")
            return True
        
        # Determine how to call MiXCR (jar vs executable)
        if str(self.mixcr_jar).endswith('.jar'):
            mixcr_call = f"java -jar \"{self.mixcr_jar}\""
        else:
            mixcr_call = f"\"{self.mixcr_jar}\""
        
        # Build the shell script with the resolved MiXCR command
        mixcr_script = f"""#!/bin/bash

# --- MiXCR Paired Analysis Shell Script (Auto-generated) ---
set -euo pipefail

# How to invoke MiXCR (either jar or executable)
MIXCR_CALL="{mixcr_call}"

# I/O settings
INPUT_DIR="{self.matched_fastq_output}"
PREFIX="{self.prefix}_matched"
OUTPUT_DIR="{self.step3_output}"
THREADS={self.threads}

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Define convenience variables
TRA_R1="$INPUT_DIR/${{PREFIX}}_TRA_matched_1.fq.gz"
TRA_R2="$INPUT_DIR/${{PREFIX}}_TRA_matched_2.fq.gz"
TRB_R1="$INPUT_DIR/${{PREFIX}}_TRB_matched_1.fq.gz"
TRB_R2="$INPUT_DIR/${{PREFIX}}_TRB_matched_2.fq.gz"

# --------------------- TRA chain ---------------------
echo "[MiXCR] TRA analysis (analyze amplicon)"
$MIXCR_CALL analyze amplicon \
  -s hsa \
  --starting-material RNA \
  --5-end no-v-primers \
  --3-end j-primers \
  --adapters no-adapters \
  --report "$OUTPUT_DIR/TRA_analyze.report.log" \
  -t $THREADS \
  --align "-OsaveOriginalReads=true" \
  "$TRA_R1" "$TRA_R2" \
  "$OUTPUT_DIR/TRA.vdjca"

# --------------------- TRB chain ---------------------
echo "[MiXCR] TRB analysis (analyze amplicon)"
$MIXCR_CALL analyze amplicon \
  -s hsa \
  --starting-material RNA \
  --5-end no-v-primers \
  --3-end j-primers \
  --adapters no-adapters \
  --report "$OUTPUT_DIR/TRB_analyze.report.log" \
  -t $THREADS \
  --align "-OsaveOriginalReads=true" \
  "$TRB_R1" "$TRB_R2" \
  "$OUTPUT_DIR/TRB.vdjca"

# --------------------- Assemble clones ---------------------
echo "[MiXCR] Assemble TRA clones"
$MIXCR_CALL assemble \
  -f \
  --report "$OUTPUT_DIR/TRA_assemble.report.log" \
  "$OUTPUT_DIR/TRA.vdjca" \
  "$OUTPUT_DIR/TRA.vdjca.clns"

echo "[MiXCR] Assemble TRB clones"
$MIXCR_CALL assemble \
  -f \
  --report "$OUTPUT_DIR/TRB_assemble.report.log" \
  "$OUTPUT_DIR/TRB.vdjca" \
  "$OUTPUT_DIR/TRB.vdjca.clns"

# --------------------- Export alignments ---------------------
echo "[MiXCR] Export alignments"
$MIXCR_CALL exportAlignments -f -descrsR1 -vGene -jGene -nFeature CDR3 -aaFeature CDR3 \
  "$OUTPUT_DIR/TRA.vdjca" \
  "$OUTPUT_DIR/TRA_alignments_export_with_headers.tsv"

$MIXCR_CALL exportAlignments -f -descrsR1 -vGene -jGene -nFeature CDR3 -aaFeature CDR3 \
  "$OUTPUT_DIR/TRB.vdjca" \
  "$OUTPUT_DIR/TRB_alignments_export_with_headers.tsv"

echo "--- MiXCR Analysis and Export Steps Completed ---"
"""
        
        # Write the script to a temporary file
        script_path = "temp_mixcr_script.sh"
        with open(script_path, 'w') as f:
            f.write(mixcr_script)
        
        # Make it executable
        os.chmod(script_path, 0o755)
        
        try:
            # Run the script
            result = self.run_command(f"bash {script_path}", "Step 3: Run MiXCR", shell=True, step_key='step3')
            return result
        finally:
            # Clean up temporary script
            if os.path.exists(script_path):
                os.remove(script_path)

    def step4_pair_and_filter(self):
        """Step 4: Pair and filter clones."""
        self.logger.info("="*50)
        self.logger.info("STEP 4: Starting pair and filter clones")
        self.logger.info("="*50)
        
        # Check if step is already completed
        if self.check_step_completion('step4'):
            print("Step 4: Pair and Filter Clones - SKIPPED (already completed)")
            return True
        
        tra_export = os.path.join(self.step3_output, "TRA_alignments_export_with_headers.tsv")
        trb_export = os.path.join(self.step3_output, "TRB_alignments_export_with_headers.tsv")
        
        script_path = os.path.join(self.scripts_dir, "4_pair_and_filter_clones.py")
        cmd = [
            "python3", script_path,
            "--umi-pairs", self.umi_pairs_file,
            "--tra-export", tra_export,
            "--trb-export", trb_export,
            "-o", self.final_output
        ]
        return self.run_command(cmd, "Step 4: Pair and Filter Clones", step_key='step4')

    def run_pipeline(self):
        """Run the complete pipeline."""
        print("="*60)
        print("STARTING PAIRTCR PIPELINE")
        print("="*60)
        print(f"Input directory: {self.input_dir}")
        print(f"Output root: {self.output_root}")
        print(f"Prefix: {self.prefix}")
        print(f"Read limit: {self.read_limit}")
        print(f"Threads: {self.threads}")
        print(f"MiXCR JAR: {self.mixcr_jar}")
        print(f"Preprocessor: {'C version (fast)' if self.use_c_version else 'Python version'}")
        print(f"Logs directory: {self.logs_output}")
        print("="*60)
        
        # Log pipeline start
        self.logger.info("="*60)
        self.logger.info("STARTING PAIRTCR PIPELINE")
        self.logger.info("="*60)
        self.logger.info(f"Input directory: {self.input_dir}")
        self.logger.info(f"Output root: {self.output_root}")
        self.logger.info(f"Prefix: {self.prefix}")
        self.logger.info(f"Read limit: {self.read_limit}")
        self.logger.info(f"Threads: {self.threads}")
        self.logger.info(f"MiXCR JAR: {self.mixcr_jar}")
        self.logger.info(f"Preprocessor: {'C version (fast)' if self.use_c_version else 'Python version'}")
        self.logger.info(f"Logs directory: {self.logs_output}")
        self.logger.info("="*60)
        
        # Validate inputs
        if not os.path.exists(self.input_dir):
            error_msg = f"Input directory '{self.input_dir}' does not exist."
            self.logger.error(error_msg)
            print(f"Error: {error_msg}")
            sys.exit(1)
        
        # Resolve MiXCR location: accept jar path or fallback to 'mixcr' executable in PATH
        if not os.path.exists(self.mixcr_jar):
            # Try environment variable first
            env_mixcr = os.environ.get("MIXCR_JAR")
            if env_mixcr and os.path.exists(env_mixcr):
                warn_msg = (f"MiXCR JAR '{self.mixcr_jar}' not found. "
                            f"Falling back to environment variable MIXCR_JAR='{env_mixcr}'.")
                self.logger.warning(warn_msg)
                print(f"Warning: {warn_msg}")
                self.mixcr_jar = env_mixcr
            else:
                alt_mixcr_path = shutil.which("mixcr")
                if alt_mixcr_path:
                    warn_msg = (f"MiXCR JAR '{self.mixcr_jar}' not found. "
                                f"Falling back to 'mixcr' executable found at: {alt_mixcr_path}")
                    self.logger.warning(warn_msg)
                    print(f"Warning: {warn_msg}")
                    self.mixcr_jar = alt_mixcr_path
                else:
                    error_msg = (f"MiXCR JAR file '{self.mixcr_jar}' does not exist and "
                                 f"'mixcr' executable could not be found in PATH.")
                    self.logger.error(error_msg)
                    print(f"Error: {error_msg}")
                    sys.exit(1)
        
        # Create output root directory
        os.makedirs(self.output_root, exist_ok=True)
        self.logger.info(f"Created/verified output root directory: {self.output_root}")
        
        # Check if pipeline is already completed
        if self.check_pipeline_completion() and not self.force_restart:
            print("\n" + "="*60)
            print("PIPELINE ALREADY COMPLETED!")
            print("="*60)
            print(f"Final output file: {self.final_output}")
            print("Use --force to restart from beginning")
            self.logger.info("Pipeline already completed, exiting")
            return
        elif self.force_restart:
            # Force restart: clean up everything and start fresh
            print("\nForce restart requested. Cleaning up all previous results...")
            self.logger.info("Force restart requested")
            self.cleanup_incomplete_run()
        else:
            # If not completely finished, clean up partial results and restart
            print("\nDetected incomplete pipeline run. Cleaning up and restarting...")
            self.cleanup_incomplete_run()
        
        try:
            # Run all pipeline steps with proper error checking
            if not self.step1_preprocess_and_trim():
                error_msg = "Step 1 (Preprocess and Trim) failed"
                self.logger.error(error_msg)
                print(f"Error: {error_msg}")
                sys.exit(1)
                
            if not self.step2_create_umi_pairs():
                error_msg = "Step 2 (Create UMI Pairs) failed"
                self.logger.error(error_msg)
                print(f"Error: {error_msg}")
                sys.exit(1)
                
            if not self.step2_5_create_matched_fastq():
                error_msg = "Step 2.5 (Create Matched FASTQ) failed"
                self.logger.error(error_msg)
                print(f"Error: {error_msg}")
                sys.exit(1)
                
            if not self.step3_run_mixcr():
                error_msg = "Step 3 (Run MiXCR) failed"
                self.logger.error(error_msg)
                print(f"Error: {error_msg}")
                sys.exit(1)
                
            if not self.step4_pair_and_filter():
                error_msg = "Step 4 (Pair and Filter Clones) failed"
                self.logger.error(error_msg)
                print(f"Error: {error_msg}")
                sys.exit(1)
            
            print("\n" + "="*60)
            print("PIPELINE COMPLETED SUCCESSFULLY!")
            print("="*60)
            print(f"Final output file: {self.final_output}")
            print(f"All intermediate files are saved in: {self.output_root}")
            
            # Log successful completion
            self.logger.info("="*60)
            self.logger.info("PIPELINE COMPLETED SUCCESSFULLY!")
            self.logger.info("="*60)
            self.logger.info(f"Final output file: {self.final_output}")
            self.logger.info(f"All intermediate files are saved in: {self.output_root}")
            
            # Print summary of outputs
            print("\nOutput Summary:")
            print(f"  Step 1 outputs: {self.step1_output}")
            print(f"  Step 2 outputs: {self.step2_output}")
            print(f"  Matched FASTQ: {self.matched_fastq_output}")
            print(f"  Step 3 outputs: {self.step3_output}")
            print(f"  Step 4 outputs: {self.step4_output}")
            print(f"  Logs directory: {self.logs_output}")
            print(f"  Final result: {self.final_output}")
            
            # Log output summary
            self.logger.info("Output Summary:")
            self.logger.info(f"  Step 1 outputs: {self.step1_output}")
            self.logger.info(f"  Step 2 outputs: {self.step2_output}")
            self.logger.info(f"  Matched FASTQ: {self.matched_fastq_output}")
            self.logger.info(f"  Step 3 outputs: {self.step3_output}")
            self.logger.info(f"  Step 4 outputs: {self.step4_output}")
            self.logger.info(f"  Logs directory: {self.logs_output}")
            self.logger.info(f"  Final result: {self.final_output}")
            
            # Print log files summary
            print("\nLog Files:")
            print(f"  Main pipeline log: {self.pipeline_log}")
            for step, log_file in self.step_logs.items():
                print(f"  {step.upper()} log: {log_file}")
            
            self.logger.info("Log Files:")
            self.logger.info(f"  Main pipeline log: {self.pipeline_log}")
            for step, log_file in self.step_logs.items():
                self.logger.info(f"  {step.upper()} log: {log_file}")
            
        except Exception as e:
            error_msg = f"Pipeline failed with error: {e}"
            self.logger.error(error_msg)
            print(f"\n{error_msg}")
            print(f"Check the main log file for details: {self.pipeline_log}")
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="PairTCR Pipeline: Complete workflow for paired TCR analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("input_dir", nargs='?', default=DEFAULT_INPUT_DIR,
                        help="Path to directory containing input FASTQ files")
    parser.add_argument("-o", "--output-root", default=DEFAULT_OUTPUT_ROOT,
                        help="Root directory for all pipeline outputs")
    parser.add_argument("-p", "--prefix", default=DEFAULT_PREFIX,
                        help="Prefix for output files")
    parser.add_argument("-n", "--read-limit", type=int, default=DEFAULT_READ_LIMIT,
                        help="Maximum number of read pairs to process")
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS,
                        help="Number of threads for MiXCR")
    parser.add_argument("--mixcr-jar", default=DEFAULT_MIXCR_JAR,
                        help="Path to MiXCR JAR file")
    parser.add_argument("--force", action="store_true",
                        help="Force restart pipeline from beginning, even if already completed")
    parser.add_argument("--use-c", action="store_true",
                        help="Use C version of preprocessor for faster processing (requires compilation)")
    
    args = parser.parse_args()
    
    # Create pipeline runner and execute
    pipeline = PipelineRunner(
        input_dir=args.input_dir,
        output_root=args.output_root,
        prefix=args.prefix,
        read_limit=args.read_limit,
        threads=args.threads,
        mixcr_jar=args.mixcr_jar,
        force_restart=args.force,
        use_c_version=args.use_c
    )
    
    pipeline.run_pipeline()


if __name__ == "__main__":
    main() 