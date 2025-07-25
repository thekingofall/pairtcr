#!/usr/bin/env python3
"""
PairTCR - Paired T-cell Receptor Analysis Pipeline
A high-performance toolkit for processing paired TCR sequencing data
"""

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import os
import shutil

# Read the contents of README file
def read_file(fname):
    """Read file contents"""
    try:
        with open(os.path.join(os.path.dirname(__file__), fname), encoding='utf-8') as f:
            return f.read()
    except FileNotFoundError:
        return ""

# Get version
def get_version():
    """Get version from __init__.py"""
    version_file = os.path.join(os.path.dirname(__file__), 'pairtcr', '__init__.py')
    try:
        with open(version_file) as f:
            for line in f:
                if line.startswith('__version__'):
                    return line.split('"')[1]
    except FileNotFoundError:
        pass
    return "0.1.0"

class CustomBuildPy(build_py):
    """Custom build command to copy scripts into package directory"""
    def run(self):
        # Copy scripts directory into pairtcr package if it doesn't exist
        source_scripts = os.path.join(self.build_lib, '..', '..', 'scripts')
        target_scripts = os.path.join(self.build_lib, 'pairtcr', 'scripts')
        
        # Create the target directory if it doesn't exist
        os.makedirs(os.path.dirname(target_scripts), exist_ok=True)
        
        # Copy scripts if source exists and target doesn't exist
        if os.path.exists(source_scripts) and not os.path.exists(target_scripts):
            print(f"Copying scripts from {source_scripts} to {target_scripts}")
            shutil.copytree(source_scripts, target_scripts)
        
        # Run the standard build
        build_py.run(self)

# Requirements
install_requires = [
    "pandas>=1.0.0",
    "tqdm>=4.50.0",
    "numpy>=1.18.0",
]

# Extra requirements for development
extras_require = {
    'dev': [
        'pytest>=6.0',
        'pytest-cov',
        'flake8',
        'black',
        'mypy',
    ],
    'c_extensions': [
        'pybind11>=2.6.0',
    ]
}

# Console scripts entry points that call existing scripts
entry_points = {
    'console_scripts': [
        'pairtcr-pipeline=pairtcr.cli:run_pipeline',
        'pairtcr-preprocess=pairtcr.cli:run_preprocess',
        'pairtcr-umi-pairs=pairtcr.cli:run_umi_pairs',
        'pairtcr-pair-filter=pairtcr.cli:run_pair_filter',
        'pairtcr=pairtcr.cli:main',
    ],
}

setup(
    name="pairtcr",
    version=get_version(),
    author="PairTCR Development Team",
    author_email="your-email@example.com",
    description="Paired T-cell Receptor Analysis Pipeline",
    long_description=read_file("README.md"),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/pairtcr",
    project_urls={
        "Bug Tracker": "https://github.com/yourusername/pairtcr/issues",
        "Documentation": "https://pairtcr.readthedocs.io/",
        "Source Code": "https://github.com/yourusername/pairtcr",
    },
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
    ],
    python_requires=">=3.6",
    install_requires=install_requires,
    extras_require=extras_require,
    entry_points=entry_points,
    include_package_data=True,
    package_data={
        'pairtcr': [
            'scripts/*.py',
            'scripts/*.sh',
            'scripts/*.c',
            'scripts/Makefile',
            'scripts/mixcr',
            'scripts/mixcr.jar',
            'scripts/1_preprocess_and_trim',
        ],
    },
    cmdclass={
        'build_py': CustomBuildPy,
    },
    zip_safe=False,
    keywords=['tcr', 'sequencing', 'bioinformatics', 'immunology', 'paired-end'],
    platforms=['any'],
) 