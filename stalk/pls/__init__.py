#!/usr/bin/env python3
"""Surrogate Hessian accelerated parallel line-search: parallel line-search"""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from stalk.pls.parallel_linesearch import ParallelLineSearch
from stalk.pls.surrogate import Surrogate

__all__ = [
    'ParallelLineSearch',
    'Surrogate',
]
