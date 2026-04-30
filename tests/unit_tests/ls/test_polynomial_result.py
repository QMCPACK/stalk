#!/usr/bin/env python

from numpy import linspace, polyval
from stalk.ls.polynomial_result import PolynomialResult
from stalk.util.util import match_to_tol

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test PolynomialResult class
def test_PolynomialResult():

    # test default init
    a = 1.23
    b = 0.3
    c = 0.0
    fit = [a, b, c]
    # x0 = -b/2a
    x0 = -b / 2 / a
    y0 = -b**2 / 4 / a

    res = PolynomialResult(x0, y0, fit=fit)
    assert match_to_tol(res.get_force(x0), 0.0)
    assert match_to_tol(res.get_hessian(x0), 2 * a)
    d = 0.2
    assert match_to_tol(res.get_force(x0 + d), -2 * a * d)
    grid = linspace(-1.0, 3.0, 5)
    assert match_to_tol(res.get_values(grid), polyval(fit, grid))

# end def
