"""
Dose Response Analyzer (DRA) - A comprehensive Python package for dose-response curve analysis.

This package provides Python equivalent functionality to R's drc package for
dose-response curve analysis with comprehensive model fitting and visualization capabilities.

Classes:
    DoseResponseAnalyzer: Main analysis class for fitting dose-response models.
    DoseResponsePlotter: Comprehensive plotting class for visualization.

Functions:
    example_usage: Demonstrates module usage with synthetic data.
    reconstruct_curve_from_parameters: Reconstruct curves from known parameters.
    reconstruct_curve_from_results: Reconstruct curves from analysis results.

Examples:
    >>> import pandas as pd
    >>> from dra import DoseResponseAnalyzer, DoseResponsePlotter
    >>>
    >>> # Using default column names (Compound, Conc, Rab10)
    >>> analyzer = DoseResponseAnalyzer()
    >>> results = analyzer.fit_best_models(df)
    >>> plotter = DoseResponsePlotter()
    >>> plotter.plot_dose_response_curves(results, analyzer, df)
    >>>
    >>> # Reconstruct curve from known parameters
    >>> from dra import reconstruct_curve_from_parameters
    >>> conc, resp = reconstruct_curve_from_parameters(
    ...     model_name='LL.4', top=1.0, bottom=0.1, ic50=100.0, hillslope=1.5
    ... )
"""

from .dose_response_analyzer import DoseResponseAnalyzer, DoseResponsePlotter, example_usage

__version__ = "0.1.0"
__email__ = "toan.phung@proteo.info"

__all__ = [
    "DoseResponseAnalyzer",
    "DoseResponsePlotter",
    "example_usage",
    "reconstruct_curve_from_parameters",
    "reconstruct_curve_from_results",
]

# Convenience imports for curve reconstruction
reconstruct_curve_from_parameters = DoseResponseAnalyzer.reconstruct_curve_from_parameters
reconstruct_curve_from_results = DoseResponseAnalyzer.reconstruct_curve_from_results
