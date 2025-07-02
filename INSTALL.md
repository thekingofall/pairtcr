# PairTCR Installation and Packaging Guide

## Quick Installation

### 1. Install from source (recommended for development)

```bash
# Clone or navigate to the PairTCR directory
cd PairTCR

# Install in development mode (editable installation)
pip install -e .

# Verify installation
pairtcr --version
pairtcr --help
```

### 2. Install from a built package

```bash
# Build the package
python setup.py sdist bdist_wheel

# Install from built package
pip install dist/pairtcr-0.1.0-py3-none-any.whl
```

## Available Commands After Installation

- `pairtcr` - Main command with version and help
- `pairtcr-pipeline` - Run complete pipeline
- `pairtcr-preprocess` - Preprocess FASTQ files  
- `pairtcr-umi-pairs` - Create UMI pairs
- `pairtcr-pair-filter` - Pair and filter clones

## Example Usage

```bash
# Run complete pipeline
pairtcr-pipeline input_directory/ -n 50000 --threads 8

# Run individual steps
pairtcr-preprocess input_directory/ -n 100000
pairtcr-umi-pairs -i preprocessed_output/ -p sample_prefix
```

## Building and Distributing

### Build source and wheel distributions

```bash
# Install build dependencies
pip install build

# Build packages
python -m build

# This creates:
# dist/pairtcr-0.1.0.tar.gz (source distribution)
# dist/pairtcr-0.1.0-py3-none-any.whl (wheel distribution)
```

### Upload to PyPI (when ready)

```bash
# Install twine for uploading
pip install twine

# Upload to test PyPI first
twine upload --repository testpypi dist/*

# Upload to real PyPI
twine upload dist/*
```

## Dependencies

The package automatically installs these Python dependencies:
- pandas >= 1.0.0
- tqdm >= 4.50.0
- numpy >= 1.18.0

## Optional C Extensions

For maximum performance, you can install the C extensions:

```bash
# Install build tools
pip install pybind11

# Compile C extensions manually
cd scripts
make

# Use C version in pipeline
pairtcr-pipeline input_dir/ --use-c
```

## Development Setup

```bash
# Clone repository
git clone <repository-url>
cd PairTCR

# Install in development mode with extra dependencies
pip install -e .[dev]

# Run tests (if available)
pytest

# Code formatting
black pairtcr/
flake8 pairtcr/
```

## Package Structure

```
PairTCR/
├── setup.py              # Package configuration
├── README.md             # Package documentation
├── LICENSE               # MIT license
├── MANIFEST.in           # Include additional files
├── requirements.txt      # Dependencies
├── pairtcr/              # Python package
│   ├── __init__.py       # Package initialization
│   └── cli.py           # Command-line interface
└── scripts/              # Original scripts (included in package)
    ├── 1_preprocess_and_trim.py
    ├── 2_create_umi_pairs.py
    ├── 3_run_mixcr_and_export.sh
    ├── 4_pair_and_filter_clones.py
    ├── 5_runpipeline.py
    ├── 1_preprocess_and_trim.c
    └── Makefile
```

## Troubleshooting

### Command not found after installation

If commands like `pairtcr` are not found after installation:

```bash
# Check if ~/.local/bin is in PATH
echo $PATH

# Add to PATH if needed (add to ~/.bashrc for permanent)
export PATH=$PATH:~/.local/bin

# Or use full path
~/.local/bin/pairtcr --version
```

### Permission errors during installation

```bash
# Use user installation
pip install --user -e .

# Or use virtual environment
python -m venv pairtcr-env
source pairtcr-env/bin/activate
pip install -e .
```

### Python version compatibility

The package requires Python 3.6 or higher. Check your version:

```bash
python --version
```

If you have multiple Python versions, use the specific version:

```bash
python3.7 -m pip install -e .
``` 