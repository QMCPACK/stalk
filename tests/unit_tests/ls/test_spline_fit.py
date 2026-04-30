#!/usr/bin/env python

from pytest import raises

from stalk.ls.spline_fit import SplineFit
from stalk.ls.spline_result import SplineResult
from stalk.util.util import match_to_tol

from ..assets.fitting_pf2 import generate_exact_pf2

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test SplineFit class
def test_SplineFit():

    # Test nominal using regular 2-degree polynomial fit
    h = 4.5
    grid, ref = generate_exact_pf2(1.23, 2.34, h=h, N=5, error=0.1)
    fit = SplineFit()
    # Find minimum, noise not requested
    fit_res = fit.find_minimum(grid)
    assert isinstance(fit_res, SplineResult)
    assert fit_res.analyzed
    assert match_to_tol(fit_res.x0, ref.x0)
    assert match_to_tol(fit_res.y0, ref.y0)
    assert fit_res.x0_err == 0.0
    assert fit_res.y0_err == 0.0

    # Test too small grid
    with raises(ValueError):
        grid_small, ref_small = generate_exact_pf2(1.23, 2.34, N=2)
        SplineFit().find_minimum(grid_small)
    # end with

    # Test __eq__
    assert fit == SplineFit()

# end def
