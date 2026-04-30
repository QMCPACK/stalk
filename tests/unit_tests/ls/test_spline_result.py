#!/usr/bin/env python

from numpy import linspace
from scipy.interpolate import CubicSpline

from stalk.ls.spline_result import SplineResult
from stalk.util.util import match_to_tol

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test PolynomialResult class
def test_PolynomialResult():

    h = 1.23
    x0 = 0.2
    y0 = -0.4
    grid = linspace(-2, 2, 7)
    values = y0 + 0.5 * h * (grid - x0)**2

    fit = CubicSpline(grid, values)

    res = SplineResult(x0, y0, fit=fit)
    assert match_to_tol(res.get_force(x0), 0.0)
    assert match_to_tol(res.get_hessian(x0), h)
    d = 0.1
    assert match_to_tol(res.get_force(x0 + d), -h * d)
    assert match_to_tol(res.get_values(grid + d), fit(grid + d))

# end def
