#!/usr/bin/env python3
"""Surrogate Hessian accelerated parallel line-search: line-search"""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

# from .ErrorSurface import ErrorSurface
from .fitting_function import FittingFunction
from .fitting_result import FittingResult
from .linesearch import LineSearch
from .linesearch_base import LineSearchBase
from .linesearch_grid import LineSearchGrid
from .ls_settings import LsSettings
from .morse_fit import MorseFit
from .morse_result import MorseResult
from .polynomial_fit import PolynomialFit
from .polynomial_result import PolynomialResult
from .spline_fit import SplineFit
from .spline_result import SplineResult
from .target_linesearch import TargetLineSearch
from .target_linesearch_base import TargetLineSearchBase
from .tls_settings import TlsSettings

__all__ = [
    'FittingFunction',
    'FittingResult',
    'LineSearch',
    'LineSearchBase',
    'LineSearchGrid',
    'LsSettings',
    'MorseFit',
    'MorseResult',
    'PolynomialFit',
    'PolynomialResult',
    'SplineFit',
    'SplineResult',
    'TargetLineSearch',
    'TargetLineSearchBase',
    'TlsSettings',
]
