#!/usr/bin/env python3
"""Fitting utilities for line-searches"""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from .fitting_function import FittingFunction
from .fitting_result import FittingResult
from .morse_fit import MorseFit
from .morse_result import MorseResult
from .polynomial_fit import PolynomialFit
from .polynomial_result import PolynomialResult
from .spline_fit import SplineFit
from .spline_result import SplineResult

__all__ = [
    'FittingFunction',
    'FittingResult',
    'MorseFit',
    'MorseResult',
    'PolynomialFit',
    'PolynomialResult',
    'SplineFit',
    'SplineResult',
]
