#!/usr/bin/env python3
import gzip
import re
import argparse
import sys
import os
# No longer need glob as prefix is defaulted or provided
from collections import defaultdict
from tqdm import tqdm

# --- Configuration ---
# Default input directory (output directory of script 1)
DEFAULT_INPUT_DIR = os.path.join("PairTCR_results", "1_preprocess_and_trim_output")
# Default prefix (assuming this common prefix from the pipeline)
DEFAULT_PREFIX = "TCR_TSO_18"
# Default output TSV placed inside the common results folder
DEFAULT_OUTPUT_FILE = os.path.join("PairTCR_results", "2_create_umi_pairs_output", "umi_pairs.tsv")

# --- Helper Functions --- (Keep these as they are)
def read_fastq_record(handle):
    lines = [handle.readline().strip() for _ in range(4)]
    if not lines[0]: return None
    if len(lines) < 4 or not lines[3]:
         if any(lines): pass
         return None
    return lines

def extract_umi_from_header(header):
    match = re.search(r'UMI:(?:TRA|TRB):([A-Za-z0-9_]+)(?::RC)?', header)
    if match: return match.group(1)
    return None

def get_base_read_id(header):
    if not header: return None
    base_id_part = header.split(maxsplit=1)[0]
    base_id_part = re.sub(r'/[12]$', '', base_id_part)
    return base_id_part[1:] if base_id_part.startswith('@') else base_id_part

def reverse_complement(seq):
    if seq is None: return None
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '_': '_'}
    try:
        return "".join(complement.get(base, 'N') for base in reversed(seq))
    except TypeError: return None

def reverse_complement_umi(umi_seq):
    if umi_seq is None or '_' not in umi_seq: return None
    parts = umi_seq.split('_')
    if len(parts) != 2: return None
    umi1, umi2 = parts
    rc_umi1 = reverse_complement(umi1)
    rc_umi2 = reverse_complement(umi2)
    if rc_umi1 is None or rc_umi2 is None: return None
    return f"{rc_umi2}_{rc_umi1}"

# --- Main Processing Function --- (Keep as is)
def find_umi_pairs(input_dir, prefix, output_file):
    # This function now assumes prefix is correctly determined *before* calling it
    tra_r1_file = os.path.join(input_dir, f"{prefix}_TRA_1.fq.gz")
    trb_r1_file = os.path.join(input_dir, f"{prefix}_TRB_1.fq.gz")

    # Input dir validation is still useful
    if not os.path.isdir(input_dir):
        print(f"Error: Input directory not found or is not a directory: {input_dir}", file=sys.stderr)
        print(f"Please ensure the directory exists or specify a different one using -i.", file=sys.stderr)
        sys.exit(1)

    # Validate input file existence based on the provided/defaulted prefix
    missing_files = []
    for f in [tra_r1_file, trb_r1_file]:
        if not os.path.exists(f):
            missing_files.append(os.path.basename(f))
    if missing_files:
        print(f"Error: Input file(s) not found in '{input_dir}' for prefix '{prefix}': {', '.join(missing_files)}", file=sys.stderr)
        print(f"Please ensure files exist or specify the correct directory/prefix using -i / -p.", file=sys.stderr)
        sys.exit(1)

    # --- Steps 1, 2, 3: Read TRA, RC, Read TRB, Find Pairs ---
    # (Code for these steps remains unchanged - shortened for brevity)
    # Step 1: Read TRA Data
    tra_data = defaultdict(set); processed_tra_count = 0; extracted_tra_umi_count = 0
    print(f"Reading TRA UMIs from {os.path.basename(tra_r1_file)}...")
    try:
        with gzip.open(tra_r1_file, 'rt') as tra_in, tqdm(desc="Processing TRA R1", unit="record", ascii=True, mininterval=1.0, leave=False) as pbar:
             while True:
                record = read_fastq_record(tra_in);
                if record is None: break
                pbar.update(1); processed_tra_count += 1
                umi = extract_umi_from_header(record[0]); read_id = get_base_read_id(record[0])
                if umi and read_id: tra_data[umi].add(read_id); extracted_tra_umi_count += 1
    except Exception as e: print(f"\nError reading/processing {tra_r1_file}: {e}", file=sys.stderr); sys.exit(1)
    if not tra_data: print(f"Warning: No TRA UMIs extracted from {processed_tra_count} records. Cannot find pairs.", file=sys.stderr); sys.exit(0) # Exit gracefully
    print(f"Processed {processed_tra_count} TRA records, extracted {extracted_tra_umi_count} UMIs ({len(tra_data)} unique).")

    # Step 2: Calculate RC TRA UMIs
    tra_rc_map = {}; tra_rc_to_ids = defaultdict(set); rc_failed_count = 0
    print("Calculating reverse complements of TRA UMIs...")
    for tra_umi, read_ids in tqdm(tra_data.items(), desc="RC TRA UMIs", unit="UMI", ascii=True, leave=False):
        rc_tra_umi = reverse_complement_umi(tra_umi)
        if rc_tra_umi is not None: tra_rc_map[rc_tra_umi] = tra_umi; tra_rc_to_ids[rc_tra_umi].update(read_ids)
        else: rc_failed_count += 1
    if rc_failed_count > 0: print(f"Warning: Failed to reverse complement {rc_failed_count} unique TRA UMIs.", file=sys.stderr)
    if not tra_rc_map: print("Warning: No valid reverse complement TRA UMIs generated. Cannot find pairs.", file=sys.stderr); sys.exit(0) # Exit gracefully
    print(f"Generated reverse complements for {len(tra_rc_map)} unique TRA UMI patterns.")

    # Step 3: Read TRB Data and Find Pairs
    print(f"Reading TRB UMIs from {os.path.basename(trb_r1_file)} and finding pairs...")
    pairs_found = 0; processed_trb_count = 0; extracted_trb_umi_count = 0; matched_trb_umi_count = 0
    try:
        with gzip.open(trb_r1_file, 'rt') as trb_in, open(output_file, 'w') as out_f:
            out_f.write("TRA_UMI\tTRB_UMI\tTRA_Read_ID_Base\tTRB_Read_ID_Base\n")
            with tqdm(desc="Processing TRB R1", unit="record", ascii=True, mininterval=1.0, leave=False) as pbar:
                while True:
                    record = read_fastq_record(trb_in)
                    if record is None: break
                    pbar.update(1); processed_trb_count += 1
                    trb_umi = extract_umi_from_header(record[0]); trb_read_id = get_base_read_id(record[0])
                    if trb_umi and trb_read_id:
                        extracted_trb_umi_count += 1
                        if trb_umi in tra_rc_to_ids:
                            matched_trb_umi_count +=1; original_tra_umi = tra_rc_map[trb_umi]
                            associated_tra_read_ids = tra_rc_to_ids[trb_umi]
                            for tra_id in associated_tra_read_ids:
                                out_f.write(f"{original_tra_umi}\t{trb_umi}\t{tra_id}\t{trb_read_id}\n")
                                pairs_found += 1
    except IOError as e: print(f"Error writing to output file {output_file}: {e}", file=sys.stderr); sys.exit(1)
    except Exception as e: print(f"\nError reading/processing {trb_r1_file}: {e}", file=sys.stderr); sys.exit(1)

    print(f"\n--- Pairing Summary ---")
    print(f"Processed {processed_trb_count} TRB records, extracted {extracted_trb_umi_count} UMIs.")
    print(f"Found {matched_trb_umi_count} TRB records with UMIs matching a reverse-complemented TRA UMI.")
    print(f"Wrote {pairs_found} total TRA-TRB UMI pairing lines (combinations of read IDs).")
    print(f"Output written to: {os.path.abspath(output_file)}")


# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Finds paired TRA/TRB UMIs from processed FASTQ files where RC(TRA_UMI) == TRB_UMI. "
                    "Defaults input directory and prefix based on a common workflow.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Shows defaults in help
        )
    # Input directory HAS a default
    parser.add_argument("-i", "--input-dir",
                        default=DEFAULT_INPUT_DIR,
                        help="Path to the directory containing the processed FASTQ files.")
    # Prefix HAS A DEFAULT now
    parser.add_argument("-p", "--prefix",
                        default=DEFAULT_PREFIX,
                        help="Prefix for the processed FASTQ files.")
    # Output file has a default
    parser.add_argument("-o", "--output-file", default=DEFAULT_OUTPUT_FILE,
                        help="Path to the output TSV file containing the matched UMI pairs.")

    args = parser.parse_args()

    # No more prefix auto-detection needed, just use args.prefix directly

    # Validate output directory
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            print(f"Error: Output directory '{output_dir}' could not be created: {e}", file=sys.stderr)
            sys.exit(1)

    # --- Run the main function ---
    print("Starting UMI pairing process...")
    # Use the default or user-specified values directly from args
    print(f"Input Directory: {os.path.abspath(args.input_dir)}")
    print(f"File Prefix: {args.prefix}") # Will be DEFAULT_PREFIX if not overridden
    print(f"Output File: {os.path.abspath(args.output_file)}")
    print("-" * 20)

    # Pass the arguments directly to the function
    find_umi_pairs(args.input_dir, args.prefix, args.output_file)

    print("\nPairing finished.")
