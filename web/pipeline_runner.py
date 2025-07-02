#!/usr/bin/env python3
"""
Web-enhanced PairTCR Pipeline Runner
Extended version with real-time logging for web interface
"""

import os
import sys
import subprocess
from pathlib import Path

# Import the original pipeline runner from scripts
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
scripts_dir = os.path.join(parent_dir, 'scripts')
sys.path.insert(0, parent_dir)

# Import from the scripts directory
sys.path.insert(0, scripts_dir)
from scripts.pipeline_runner_base import PipelineRunner

class WebPipelineRunner(PipelineRunner):
    """Enhanced pipeline runner for web interface with real-time logging"""
    
    def __init__(self, input_dir, output_root, prefix, read_limit, threads, mixcr_jar, 
                 force_restart=False, use_c_version=False, job_id=None, job_manager=None):
        # Initialize parent class
        super().__init__(input_dir, output_root, prefix, read_limit, threads, mixcr_jar, 
                        force_restart, use_c_version)
        
        # Web-specific attributes
        self.job_id = job_id
        self.job_manager = job_manager
        
        # Override logging setup for web
        self.setup_web_logging()
    
    def setup_web_logging(self):
        """Setup web-specific logging that sends messages to job manager"""
        # Keep the original file-based logging
        super().setup_logging()
        
        # Add custom log handler for web interface
        if self.job_manager and self.job_id:
            self.web_log = lambda msg: self.job_manager.add_log_message(self.job_id, msg)
        else:
            self.web_log = lambda msg: print(f"[Web Log] {msg}")
    
    def run_command(self, cmd, step_name, shell=False, step_key=None, show_progress=True):
        """Enhanced command runner with web logging"""
        self.web_log(f"Starting {step_name}")
        self.web_log(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
        
        # Update job progress if possible
        if self.job_manager and self.job_id:
            progress_map = {
                'step1': 20,
                'step2': 40,
                'step2.5': 50,
                'step3': 70,
                'step4': 90
            }
            if step_key in progress_map:
                self.job_manager.update_job(self.job_id, 
                                          current_step=step_name,
                                          progress=progress_map[step_key])
        
        try:
            # For web interface, we want to capture output and send it live
            if shell:
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                                         stderr=subprocess.STDOUT, universal_newlines=True)
            else:
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                                         stderr=subprocess.STDOUT, universal_newlines=True)
            
            # Read output line by line and send to web log
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                line = line.rstrip()
                if line:  # Only log non-empty lines
                    self.web_log(f"[{step_name}] {line}")
            
            process.stdout.close()
            return_code = process.wait()
            
            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, cmd)
            
            self.web_log(f"{step_name} completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            error_msg = f"{step_name} failed with return code {e.returncode}"
            self.web_log(f"ERROR: {error_msg}")
            if self.job_manager and self.job_id:
                self.job_manager.update_job(self.job_id, 
                                          status='failed',
                                          error_message=error_msg)
            raise e
    
    def step1_preprocess_and_trim(self):
        """Step 1 with web logging"""
        version_info = "C version" if self.use_c_version else "Python version"
        self.web_log(f"="*50)
        self.web_log(f"STEP 1: Starting preprocess and trim ({version_info})")
        self.web_log(f"="*50)
        
        return super().step1_preprocess_and_trim()
    
    def step2_create_umi_pairs(self):
        """Step 2 with web logging"""
        self.web_log("="*50)
        self.web_log("STEP 2: Starting UMI pairs creation")
        self.web_log("="*50)
        
        return super().step2_create_umi_pairs()
    
    def step2_5_create_matched_fastq(self):
        """Step 2.5 with enhanced web logging"""
        self.web_log("="*50)
        self.web_log("STEP 2.5: Starting matched FASTQ files creation")
        self.web_log("="*50)
        
        # Check if step is already completed
        if self.check_step_completion('step2.5'):
            self.web_log("Step 2.5: Create Matched FASTQ Files - SKIPPED (already completed)")
            return True
        
        try:
            import pandas as pd
            import gzip
            from collections import defaultdict
            
            # Create output directory
            os.makedirs(self.matched_fastq_output, exist_ok=True)
            self.web_log(f"Created output directory: {self.matched_fastq_output}")
            
            # Read UMI pairs
            self.web_log(f"Reading UMI pairs from: {self.umi_pairs_file}")
            umi_pairs_df = pd.read_csv(self.umi_pairs_file, sep='\t')
            
            # Create sets of read IDs for quick lookup
            tra_read_ids = set(umi_pairs_df['TRA_Read_ID_Base'])
            trb_read_ids = set(umi_pairs_df['TRB_Read_ID_Base'])
            
            self.web_log(f"Found {len(tra_read_ids)} unique TRA read IDs and {len(trb_read_ids)} unique TRB read IDs")
            
            # Process TRA files
            tra_r1_in = os.path.join(self.step1_output, f"{self.prefix}_TRA_1.fq.gz")
            tra_r2_in = os.path.join(self.step1_output, f"{self.prefix}_TRA_2.fq.gz")
            tra_r1_out = os.path.join(self.matched_fastq_output, f"{self.prefix}_matched_TRA_matched_1.fq.gz")
            tra_r2_out = os.path.join(self.matched_fastq_output, f"{self.prefix}_matched_TRA_matched_2.fq.gz")
            
            self._filter_fastq_by_read_ids_web(tra_r1_in, tra_r1_out, tra_read_ids, "TRA R1")
            self._filter_fastq_by_read_ids_web(tra_r2_in, tra_r2_out, tra_read_ids, "TRA R2")
            
            # Process TRB files
            trb_r1_in = os.path.join(self.step1_output, f"{self.prefix}_TRB_1.fq.gz")
            trb_r2_in = os.path.join(self.step1_output, f"{self.prefix}_TRB_2.fq.gz")
            trb_r1_out = os.path.join(self.matched_fastq_output, f"{self.prefix}_matched_TRB_matched_1.fq.gz")
            trb_r2_out = os.path.join(self.matched_fastq_output, f"{self.prefix}_matched_TRB_matched_2.fq.gz")
            
            self._filter_fastq_by_read_ids_web(trb_r1_in, trb_r1_out, trb_read_ids, "TRB R1")
            self._filter_fastq_by_read_ids_web(trb_r2_in, trb_r2_out, trb_read_ids, "TRB R2")
            
            self.web_log("Step 2.5: Create Matched FASTQ Files completed successfully")
            return True
            
        except Exception as e:
            self.web_log(f"ERROR in Step 2.5: {str(e)}")
            raise e
    
    def _filter_fastq_by_read_ids_web(self, input_file, output_file, target_read_ids, file_desc):
        """Filter FASTQ file with web logging"""
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
        
        self.web_log(f"Filtering {file_desc}: {os.path.basename(input_file)} -> {os.path.basename(output_file)}")
        
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
                
                # Log progress every 10000 reads
                if total_reads % 10000 == 0:
                    self.web_log(f"  {file_desc}: Processed {total_reads} reads, {matched_reads} matched")
        
        self.web_log(f"  {file_desc}: {matched_reads}/{total_reads} reads matched and written")
    
    def step3_run_mixcr(self):
        """Step 3 with web logging"""
        self.web_log("="*50)
        self.web_log("STEP 3: Starting MiXCR analysis")
        self.web_log("="*50)
        
        return super().step3_run_mixcr()
    
    def step4_pair_and_filter(self):
        """Step 4 with web logging"""
        self.web_log("="*50)
        self.web_log("STEP 4: Starting pair and filter clones")
        self.web_log("="*50)
        
        return super().step4_pair_and_filter()
    
    def run_pipeline(self):
        """Run the complete pipeline with enhanced web logging"""
        self.web_log("="*60)
        self.web_log("STARTING PAIRTCR PIPELINE")
        self.web_log("="*60)
        self.web_log(f"Input directory: {self.input_dir}")
        self.web_log(f"Output root: {self.output_root}")
        self.web_log(f"Prefix: {self.prefix}")
        self.web_log(f"Read limit: {self.read_limit}")
        self.web_log(f"Threads: {self.threads}")
        self.web_log(f"MiXCR JAR: {self.mixcr_jar}")
        self.web_log(f"Preprocessor: {'C version (fast)' if self.use_c_version else 'Python version'}")
        self.web_log("="*60)
        
        # Run the pipeline using parent method but with enhanced error handling
        try:
            return super().run_pipeline()
        except Exception as e:
            self.web_log(f"PIPELINE FAILED: {str(e)}")
            if self.job_manager and self.job_id:
                self.job_manager.update_job(self.job_id, 
                                          status='failed',
                                          error_message=str(e))
            raise e 