#!/usr/bin/env python

from numpy import linspace
from stalk.ls.morse_result import MorseResult
from stalk.util.util import match_to_tol, morse

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test MorseResult class
def test_MorseResult():

    x0 = 1.2
    a = 2.4
    y0 = 0.5
    Einf = 0.1
    fit = [x0, a, y0, Einf]
    res = MorseResult(x0, y0, fit=fit)

    grid = linspace(1.0, 3.0, 5)
    assert match_to_tol(res.get_values(grid), morse(fit, grid))
    assert match_to_tol(res.get_force(x0), 0.0, tol=1e-5)
    H = 2 * y0 / a**2
    assert match_to_tol(res.get_hessian(x0, dx=1e-6), H, tol=1e-4)

# end def
