"""
North Atlantic SST Pattern Analysis Package

This package provides tools for analyzing SST patterns in the North Atlantic,
specifically for comparing JJA 2023 SST patterns with MPI_GE model outputs.
"""

__version__ = "0.1.0"
__author__ = "Quan Liu"

from .pattern_correlation import calculate_pattern_correlation

__all__ = ["calculate_pattern_correlation"]
