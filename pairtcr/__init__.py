"""
PairTCR - Paired T-cell Receptor Analysis Pipeline

A high-performance toolkit for processing paired TCR sequencing data with UMI-based 
error correction and comprehensive quality control.

Author: PairTCR Development Team
License: MIT
"""

__version__ = "0.1.0"
__author__ = "PairTCR Development Team"
__email__ = "your-email@example.com"
__license__ = "MIT"

# Package metadata
__all__ = [
    '__version__',
    '__author__',
    '__email__',
    '__license__',
]

# Configuration constants that can be used across modules
DEFAULT_UMI1_LEN = 7
DEFAULT_UMI2_LEN = 7
DEFAULT_READ_LIMIT = 100000
DEFAULT_THREADS = 4

# Package-level constants
PACKAGE_NAME = "pairtcr"
PACKAGE_DESCRIPTION = "Paired T-cell Receptor Analysis Pipeline"

def get_version():
    """Return the package version"""
    return __version__

def get_package_info():
    """Return package information"""
    return {
        'name': PACKAGE_NAME,
        'version': __version__,
        'description': PACKAGE_DESCRIPTION,
        'author': __author__,
        'email': __email__,
        'license': __license__,
    } 