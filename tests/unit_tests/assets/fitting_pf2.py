#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import argmin, linspace, polyder, polyfit, polyval, roots, where

from stalk.ls.linesearch_grid import LineSearchGrid
from stalk.ls.fitting_result import FittingResult
from stalk.ls.morse_result import MorseResult
from stalk.util.util import morse


def generate_exact_pf2(x_offset, y_offset, N, h=1.0, error=None):
    if N % 2 == 0:
        R = N / 2.0
    else:
        R = (N - 1) / 2.0
    # end if
    offsets = linspace(-R, R, N)
    offsets = offsets
    values = 0.5 * h * (offsets - x_offset)**2 + y_offset

    grid = LineSearchGrid(offsets)
    grid.values = values
    if error is not None:
        grid.errors = len(grid) * [error]
    # end if

    # TODO: x0_err, y0_err not considered if errors > 0
    # TODO: exact fit not considered
    ref = FittingResult(x0=x_offset, y0=y_offset)

    return grid, ref
# end def


def generate_exact_pf3(x_offset, y_offset, N, c=0.01, error=None):
    if N % 2 == 0:
        R = N / 2.0
    else:
        R = (N - 1) / 2.0
    # end if
    offsets = linspace(-R, R, N)
    offsets = offsets
    values = c * (offsets - x_offset)**3 + (offsets - x_offset)**2 + y_offset

    grid = LineSearchGrid(offsets)
    grid.values = values
    if error is not None:
        grid.errors = len(grid) * [error]
    # end if

    ref = FittingResult(x0=x_offset, y0=y_offset)

    return grid, ref
# end def


def minimize_pf(offsets, values, pfn=2):
    # Basically redo PolynomiaFit minimum finder for testing purposes
    pf = polyfit(offsets, values, pfn)
    pfd = polyder(pf)
    r = roots(pfd)
    d = polyval(polyder(pfd), r)
    # filter real minima (maxima with sgn < 0)
    x_mins = r[where((r.imag == 0) & (d > 0))].real
    if len(x_mins) > 0:
        y_mins = polyval(pf, x_mins)
        imin = argmin(abs(x_mins))
    else:
        x_mins = [min(offsets), max(offsets)]
        y_mins = polyval(pf, x_mins)
        imin = argmin(y_mins)  # pick the lowest/highest energy
    # end if
    y0 = y_mins[imin]
    x0 = x_mins[imin]
    return x0, y0, pf
# end def


def generate_exact_morse(x0, x_offset, y0, N, h=1.0):
    offsets = x_offset + linspace(0.5 * x0, 1.5 * x0, N)
    Einf = y0 + 0.5
    a = (2 * (Einf - y0) / h)**0.5
    p = [x0, a, Einf - y0, Einf]
    values = morse(p, offsets)

    grid = LineSearchGrid(offsets)
    grid.values = values
    ref = MorseResult(x0=x0, y0=y0, fit=p)

    return grid, ref
# end def
