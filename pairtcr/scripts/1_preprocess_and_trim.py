#!/usr/bin/env python3
import gzip
import re
import argparse
import sys
import os
import glob
from itertools import islice
from tqdm import tqdm # Import tqdm

# --- Configuration Constants ---
# Default output directory placed under the common results folder
DEFAULT_OUTPUT_DIR = os.path.join("PairTCR_results", "1_preprocess_and_trim_output")
UMI1_LEN = 7 # Adjust if your UMI lengths are different
UMI2_LEN = 7 # Adjust if your UMI lengths are different

# Define the components for the required structures
PRE_UMI1_TRA = "GACTCTGATGACGACGCACA"
LINKER_FWD_TRA = "GTACACGCTGGATCCGACTTGTAGA"
FLANK_TRA_SEQ = "TACTCTGCTGATACCGATGC"

PRE_UMI1_TRB = "GCATCGGTATCAGCAGAGTA"
LINKER_REV_TRB = "TCTACAAGTCGGATCCAGCGTGTAC"
FLANK_TRB_SEQ = "TGTGCGTCGTCATCAGAGTC"

# Compile Regular Expressions for the *entire* structure search
TRA_STRUCTURE_PATTERN = re.compile(
    f"({re.escape(PRE_UMI1_TRA)})"
    f"(.{{{UMI1_LEN}}})"
    f"({re.escape(LINKER_FWD_TRA)})"
    f"(.{{{UMI2_LEN}}})"
    f"({re.escape(FLANK_TRA_SEQ)}[AT])"
)

TRB_STRUCTURE_PATTERN = re.compile(
    f"({re.escape(PRE_UMI1_TRB)})"
    f"(.{{{UMI1_LEN}}})"
    f"({re.escape(LINKER_REV_TRB)})"
    f"(.{{{UMI2_LEN}}})"
    f"({re.escape(FLANK_TRB_SEQ)}[AT])"
)
# ---

def find_fastq_pair(directory):
    """Finds the first matching R1/R2 FASTQ file pair in the directory."""
    r1_pattern = os.path.join(directory, '*_1.fq.gz')
    r1_files = glob.glob(r1_pattern)

    if not r1_files:
        print(f"Error: No R1 file (*_1.fq.gz) found in directory: {directory}", file=sys.stderr)
        return None, None, None

    r1_file = sorted(r1_files)[0]
    base_name = os.path.basename(r1_file).replace('_1.fq.gz', '')
    r2_file = os.path.join(directory, f"{base_name}_2.fq.gz")

    if not os.path.exists(r2_file):
        print(f"Error: Corresponding R2 file ({os.path.basename(r2_file)}) not found for R1 file: {os.path.basename(r1_file)}", file=sys.stderr)
        return None, None, None

    print(f"Found FASTQ pair: {os.path.basename(r1_file)}, {os.path.basename(r2_file)}")
    return r1_file, r2_file, base_name

def read_fastq_record(handle):
    """Reads one FASTQ record (4 lines) from a file handle."""
    lines = list(islice(handle, 4))
    if not lines or len(lines) < 4:
        return None
    return [line.strip() for line in lines]

def reverse_complement(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    try:
        return "".join(complement[base] for base in reversed(seq))
    except KeyError as e:
        print(f"Warning: Unknown base '{e}' found in sequence. Treating as 'N'. Seq: {seq}", file=sys.stderr)
        return "".join(complement.get(base, 'N') for base in reversed(seq))


def process_reads(input_dir, out_prefix, output_dir, read_limit):
    """Processes paired FASTQ files, identifies TRA/TRB structure, extracts UMIs,
       adds UMI to header, *keeps only sequence downstream of the structure* in the
       identified read, handles reverse complements, limits reads, shows progress,
       and writes output to specified directory."""

    r1_file, r2_file, base_name = find_fastq_pair(input_dir)
    if not r1_file or not r2_file:
        sys.exit(1)

    if not out_prefix:
        out_prefix = base_name
        print(f"Using '{base_name}' as output prefix.")

    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as e:
        print(f"Error creating output directory {output_dir}: {e}", file=sys.stderr)
        sys.exit(1)

    tra_r1_out_file = os.path.join(output_dir, f"{out_prefix}_TRA_1.fq.gz")
    tra_r2_out_file = os.path.join(output_dir, f"{out_prefix}_TRA_2.fq.gz")
    trb_r1_out_file = os.path.join(output_dir, f"{out_prefix}_TRB_1.fq.gz")
    trb_r2_out_file = os.path.join(output_dir, f"{out_prefix}_TRB_2.fq.gz")

    processed_pairs = 0
    tra_pairs = 0
    trb_pairs = 0

    try:
        with gzip.open(r1_file, 'rt') as r1_in, \
             gzip.open(r2_file, 'rt') as r2_in, \
             gzip.open(tra_r1_out_file, 'wt') as tra_r1_out, \
             gzip.open(tra_r2_out_file, 'wt') as tra_r2_out, \
             gzip.open(trb_r1_out_file, 'wt') as trb_r1_out, \
             gzip.open(trb_r2_out_file, 'wt') as trb_r2_out, \
             tqdm(total=read_limit, desc="Processing Reads", unit="pair", ascii=True) as pbar:

            while processed_pairs < read_limit:
                r1_record = read_fastq_record(r1_in)
                r2_record = read_fastq_record(r2_in)

                if r1_record is None or r2_record is None:
                    if r1_record is not None or r2_record is not None:
                        print("\nWarning: Input files ended prematurely or have different numbers of reads.", file=sys.stderr)
                    pbar.total = processed_pairs
                    pbar.refresh()
                    break

                processed_pairs += 1
                pbar.update(1)

                r1_header_orig, r1_seq, r1_plus, r1_qual = r1_record
                r2_header_orig, r2_seq, r2_plus, r2_qual = r2_record

                read_type = None
                umi_str = None
                r1_header_mod = r1_header_orig
                r2_header_mod = r2_header_orig
                # Initialize processed sequences with the original ones
                processed_r1_seq = r1_seq
                processed_r1_qual = r1_qual
                processed_r2_seq = r2_seq
                processed_r2_qual = r2_qual
                found_in_rc = False

                # --- Check R1 forward and reverse complement for TRA structure ---
                tra_match = TRA_STRUCTURE_PATTERN.search(r1_seq)
                r1_seq_rc = None # Calculate only if needed
                if not tra_match:
                    r1_seq_rc = reverse_complement(r1_seq)
                    tra_match = TRA_STRUCTURE_PATTERN.search(r1_seq_rc)
                    if tra_match:
                        found_in_rc = True

                if tra_match:
                    read_type = "TRA"
                    umi1 = tra_match.group(2)
                    umi2 = tra_match.group(4)
                    umi_str = f"{umi1}_{umi2}"
                    rc_tag = ":RC" if found_in_rc else ""
                    r1_header_mod = f"{r1_header_orig} UMI:TRA:{umi_str}{rc_tag}"
                    r2_header_mod = f"{r2_header_orig} UMI:TRA:{umi_str}{rc_tag}" # Modify R2 header as well

                    # Trim R1: Keep only sequence *after* the match end
                    match_end = tra_match.end()

                    if not found_in_rc:
                        # Keep sequence after the match in original R1
                        processed_r1_seq = r1_seq[match_end:]
                        processed_r1_qual = r1_qual[match_end:]
                    else:
                        # Match was in RC. The sequence *after* the structure in the RC
                        # corresponds to the sequence *before* the structure's start
                        # position in the original R1 read.
                        L = len(r1_seq)
                        orig_structure_start_pos = L - match_end
                        processed_r1_seq = r1_seq[:orig_structure_start_pos]
                        processed_r1_qual = r1_qual[:orig_structure_start_pos]

                    # R2 sequence remains unchanged for TRA pairs
                    processed_r2_seq = r2_seq
                    processed_r2_qual = r2_qual


                # --- Check R2 forward and reverse complement for TRB structure (only if not TRA) ---
                if not read_type:
                    found_in_rc = False # Reset for TRB check
                    trb_match = TRB_STRUCTURE_PATTERN.search(r2_seq)
                    r2_seq_rc = None # Calculate only if needed
                    if not trb_match:
                        r2_seq_rc = reverse_complement(r2_seq)
                        trb_match = TRB_STRUCTURE_PATTERN.search(r2_seq_rc)
                        if trb_match:
                            found_in_rc = True

                    if trb_match:
                        read_type = "TRB"
                        umi1 = trb_match.group(2)
                        umi2 = trb_match.group(4)
                        umi_str = f"{umi1}_{umi2}"
                        rc_tag = ":RC" if found_in_rc else ""
                        r1_header_mod = f"{r1_header_orig} UMI:TRB:{umi_str}{rc_tag}" # Modify R1 header
                        r2_header_mod = f"{r2_header_orig} UMI:TRB:{umi_str}{rc_tag}"

                        # Trim R2: Keep only sequence *after* the match end
                        match_end = trb_match.end()

                        if not found_in_rc:
                             # Keep sequence after the match in original R2
                            processed_r2_seq = r2_seq[match_end:]
                            processed_r2_qual = r2_qual[match_end:]
                        else:
                             # Match was in RC. Keep sequence before structure start in original R2
                            L = len(r2_seq)
                            orig_structure_start_pos = L - match_end
                            processed_r2_seq = r2_seq[:orig_structure_start_pos]
                            processed_r2_qual = r2_qual[:orig_structure_start_pos]

                        # R1 sequence remains unchanged for TRB pairs
                        processed_r1_seq = r1_seq
                        processed_r1_qual = r1_qual


                # --- Write output (potentially trimmed sequence in one read, modified headers) ---
                # Write only if *both* reads have sequence content after potential trimming
                if read_type == "TRA":
                    if len(processed_r1_seq) > 0 and len(processed_r2_seq) > 0:
                         tra_r1_out.write(f"{r1_header_mod}\n{processed_r1_seq}\n{r1_plus}\n{processed_r1_qual}\n")
                         tra_r2_out.write(f"{r2_header_mod}\n{processed_r2_seq}\n{r2_plus}\n{processed_r2_qual}\n") # R2 is original seq
                         tra_pairs += 1
                elif read_type == "TRB":
                    if len(processed_r1_seq) > 0 and len(processed_r2_seq) > 0:
                        trb_r1_out.write(f"{r1_header_mod}\n{processed_r1_seq}\n{r1_plus}\n{processed_r1_qual}\n") # R1 is original seq
                        trb_r2_out.write(f"{r2_header_mod}\n{processed_r2_seq}\n{r2_plus}\n{processed_r2_qual}\n")
                        trb_pairs += 1

    except FileNotFoundError as e:
        print(f"\nError: File not found - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\nAn error occurred during processing: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    print("\n--- Processing Summary ---")
    print(f"Processed {processed_pairs} read pairs (limit was {read_limit}).")
    print(f"TRA pairs identified (UMI added, R1 trimmed to downstream): {tra_pairs}")
    print(f"TRB pairs identified (UMI added, R2 trimmed to downstream): {trb_pairs}")
    print(f"Output files written to directory: {output_dir}")
    print(f"  TRA R1: {os.path.basename(tra_r1_out_file)}")
    print(f"  TRA R2: {os.path.basename(tra_r2_out_file)}")
    print(f"  TRB R1: {os.path.basename(trb_r1_out_file)}")
    print(f"  TRB R2: {os.path.basename(trb_r2_out_file)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process paired-end FASTQ files to identify TRA/TRB structures, "
                    "extract UMIs, add UMI info to headers, trim identified read "
                    "to keep only sequence downstream of the structure, and write filtered/trimmed pairs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("input_dir", help="Path to the directory containing input FASTQ files (e.g., 'raw/'). "
                                          "Script will look for the first '*_1.fq.gz' and corresponding '*_2.fq.gz' pair.")
    parser.add_argument("-n", "--limit", type=int, default=100000,
                        help="Maximum number of read pairs to process.")
    parser.add_argument("-o", "--output_prefix", default=None,
                        help="Prefix for the output files. If not provided, it's derived from the input filenames.")
    parser.add_argument("-d", "--outdir", default=None,
                        help=f"Directory to write the output FASTQ files. Defaults to '{DEFAULT_OUTPUT_DIR}' if not specified.")

    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    if args.outdir is None:
        output_directory = DEFAULT_OUTPUT_DIR
        print(f"Output directory not specified, using default: '{output_directory}'")
    else:
        output_directory = args.outdir
        print(f"Using specified output directory: '{output_directory}'")

    print("Starting processing...")
    print(f"Input Directory: {args.input_dir}")
    print(f"Read Limit: {args.limit}")
    print(f"Output Prefix: {args.output_prefix if args.output_prefix else '(derived from input)'}")
    print(f"Output Directory: {output_directory}")
    print(f"Absolute Output Path: {os.path.abspath(output_directory)}")
    print(f"TRA Structure Pattern Prefix: {PRE_UMI1_TRA}...UMI...{LINKER_FWD_TRA}...UMI...{FLANK_TRA_SEQ}[AT]")
    print(f"TRB Structure Pattern Prefix: {PRE_UMI1_TRB}...UMI...{LINKER_REV_TRB}...UMI...{FLANK_TRB_SEQ}[AT]")
    print(f"UMI lengths: {UMI1_LEN} + {UMI2_LEN}")
    print("Trimming strategy: Keep only sequence downstream of identified structure in R1 (TRA) or R2 (TRB).")
    print("-" * 20)

    process_reads(args.input_dir, args.output_prefix, output_directory, args.limit)

    print("\nProcessing finished.")
