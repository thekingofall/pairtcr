#!/usr/bin/env python3
"""
PairTCR Command Line Interface

This module provides command-line entry points that call the existing scripts
in the scripts directory.
"""

import os
import sys
import subprocess
import argparse
from . import get_version, get_package_info

def get_scripts_dir():
    """Get the path to the scripts directory"""
    package_dir = os.path.dirname(__file__)
    # The scripts directory is typically installed inside the package directory
    scripts_dir_internal = os.path.join(package_dir, 'scripts')

    # Fallback #1: scripts directory as a *sibling* of the package directory (older layout)
    scripts_dir_sibling = os.path.join(os.path.dirname(package_dir), 'scripts')

    # Fallback #2: repository-style layout (two levels up from this file)
    scripts_dir_repo = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'scripts')

    # Pick the first existing directory
    for candidate in (scripts_dir_internal, scripts_dir_sibling, scripts_dir_repo):
        if os.path.isdir(candidate):
            return candidate

    # If nothing exists, return the internal path (default) – the caller will handle the error
    return scripts_dir_internal

def run_script(script_name, args=None):
    """Run a script from the scripts directory"""
    scripts_dir = get_scripts_dir()
    script_path = os.path.join(scripts_dir, script_name)
    
    if not os.path.exists(script_path):
        print(f"Error: Script not found: {script_path}", file=sys.stderr)
        sys.exit(1)
    
    # Prepare command
    if script_name.endswith('.py'):
        cmd = [sys.executable, script_path]
    elif script_name.endswith('.sh'):
        cmd = ['bash', script_path]
    else:
        cmd = [script_path]
    
    if args:
        cmd.extend(args)
    
    # Run the script
    try:
        result = subprocess.run(cmd, check=False)
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        print("\nInterrupted by user", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error running script: {e}", file=sys.stderr)
        sys.exit(1)

def run_pipeline():
    """Entry point for pairtcr-pipeline command"""
    # Get the correct scripts directory path
    scripts_dir = get_scripts_dir()
    mixcr_jar_path = os.path.join(scripts_dir, 'mixcr.jar')
    
    # Parse arguments to check if --mixcr-jar is already specified
    args = sys.argv[1:] if len(sys.argv) > 1 else []
    
    # Check if --mixcr-jar is already in the arguments
    mixcr_jar_specified = False
    for i, arg in enumerate(args):
        if arg == '--mixcr-jar' or arg.startswith('--mixcr-jar='):
            mixcr_jar_specified = True
            break
    
    # If --mixcr-jar is not specified and the package mixcr.jar exists, add it
    if not mixcr_jar_specified and os.path.exists(mixcr_jar_path):
        args.extend(['--mixcr-jar', mixcr_jar_path])
    
    run_script('5_runpipeline.py', args)

def run_preprocess():
    """Entry point for pairtcr-preprocess command"""
    args = sys.argv[1:] if len(sys.argv) > 1 else []
    run_script('1_preprocess_and_trim.py', args)

def run_umi_pairs():
    """Entry point for pairtcr-umi-pairs command"""
    args = sys.argv[1:] if len(sys.argv) > 1 else []
    run_script('2_create_umi_pairs.py', args)

def run_pair_filter():
    """Entry point for pairtcr-pair-filter command"""
    args = sys.argv[1:] if len(sys.argv) > 1 else []
    run_script('4_pair_and_filter_clones.py', args)

def main():
    """Main entry point for pairtcr command"""
    parser = argparse.ArgumentParser(
        description="PairTCR - Paired T-cell Receptor Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available commands:
  pairtcr-pipeline      Run the complete pipeline
  pairtcr-preprocess    Preprocess and trim FASTQ files  
  pairtcr-umi-pairs     Create UMI pairs
  pairtcr-pair-filter   Pair and filter clones

Examples:
  pairtcr --version
  pairtcr --help
  pairtcr-pipeline input_dir/ -n 50000 --threads 8
  pairtcr-preprocess input_dir/ -n 100000
"""
    )
    
    parser.add_argument('--version', action='version', 
                        version=f'PairTCR {get_version()}')
    parser.add_argument('--info', action='store_true',
                        help='Show package information')
    
    args = parser.parse_args()
    
    if args.info:
        info = get_package_info()
        print(f"Package: {info['name']}")
        print(f"Version: {info['version']}")
        print(f"Description: {info['description']}")
        print(f"Author: {info['author']}")
        print(f"License: {info['license']}")
        return
    
    # If no arguments, show help
    parser.print_help()

if __name__ == "__main__":
    main() 