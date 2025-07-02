# PairTCR - Paired T-cell Receptor Analysis Pipeline

A high-performance toolkit for processing paired TCR sequencing data with UMI-based error correction and comprehensive quality control.

## Features

- **High-performance preprocessing**: Fast identification of TRA/TRB structures with UMI extraction
- **UMI-based pairing**: Sophisticated UMI matching algorithm for paired TCR analysis  
- **MiXCR integration**: Seamless integration with MiXCR for sequence alignment
- **Quality filtering**: Advanced filtering to remove cross-alignments and ensure data quality
- **Complete pipeline**: End-to-end automation with resume/checkpoint functionality
- **C optimization**: Optional C-based preprocessing for maximum speed

## Installation

### From PyPI (coming soon)
```bash
pip install pairtcr
```

### From source
```bash
git clone https://github.com/yourusername/pairtcr.git
cd pairtcr
pip install -e .
```

### With C extensions (for maximum performance)
```bash
pip install pairtcr[c_extensions]
```

## Quick Start

### Complete Pipeline
```bash
# Run the complete pipeline
pairtcr-pipeline input_directory/ -n 50000 --threads 8

# With custom output directory
pairtcr-pipeline input_directory/ -o my_results/ -n 100000
```

### Individual Steps
```bash
# Step 1: Preprocess and trim FASTQ files
pairtcr-preprocess input_directory/ -n 100000

# Step 2: Create UMI pairs  
pairtcr-umi-pairs -i preprocessed_output/ -p sample_prefix

# Step 3: Run MiXCR (manual step)
# See documentation for MiXCR commands

# Step 4: Pair and filter clones
pairtcr-pair-filter --umi-pairs umi_pairs.tsv --tra-export tra_export.tsv --trb-export trb_export.tsv
```

## Commands

- `pairtcr`: Main command with version and help information
- `pairtcr-pipeline`: Run the complete analysis pipeline
- `pairtcr-preprocess`: Preprocess and trim FASTQ files to identify TRA/TRB structures
- `pairtcr-umi-pairs`: Create UMI pairs from preprocessed data
- `pairtcr-pair-filter`: Pair and filter TCR clones using MiXCR output

## Requirements

- Python 3.7+
- pandas >= 1.0.0
- tqdm >= 4.50.0  
- numpy >= 1.18.0
- MiXCR (for sequence alignment step)

## Input Data Format

Input should be paired-end FASTQ files in gzipped format:
- `*_1.fq.gz` (R1 files)
- `*_2.fq.gz` (R2 files)

## Output

The pipeline generates comprehensive results including:
- Preprocessed and trimmed FASTQ files
- UMI pair mappings
- MiXCR alignment results  
- Final paired and filtered TCR clones

## Performance

For optimal performance:
- Use SSD storage for input/output
- Allocate sufficient RAM (8GB+ recommended)
- Use multiple CPU cores with `--threads` option
- Consider C extensions for preprocessing-heavy workloads

## License

MIT License - see LICENSE file for details.

## Citation

If you use PairTCR in your research, please cite:
```
[Citation information to be added]
```

## Support

- Issues: https://github.com/yourusername/pairtcr/issues
- Documentation: https://pairtcr.readthedocs.io/
- Email: your-email@example.com
