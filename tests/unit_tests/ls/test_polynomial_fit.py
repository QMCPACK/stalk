#!/usr/bin/env python

from numpy import ones
from pytest import raises

from stalk.ls.polynomial_fit import PolynomialFit
from stalk.ls.fitting_function import FittingFunction
from stalk.ls.polynomial_result import PolynomialResult
from stalk.util.util import match_to_tol

from ..assets.fitting_pf2 import generate_exact_pf2, minimize_pf

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test FittingFunction class
def test_PolynomialFit():

    # Degree must be provided
    with raises(TypeError):
        PolynomialFit(None)
    # end with
    # Degree must be > 1
    with raises(ValueError):
        PolynomialFit(1)
    # end with

    # Test nominal using regular 2-degree polynomial fit
    h = 4.5
    grid, ref = generate_exact_pf2(1.23, 2.34, h=h, N=5, error=0.1)
    fit = PolynomialFit(2)
    # Find minimum, noise not requested
    fit_res = fit.find_minimum(grid)
    assert isinstance(fit_res, PolynomialResult)
    assert fit_res.analyzed
    assert match_to_tol(fit_res.x0, ref.x0)
    assert match_to_tol(fit_res.y0, ref.y0)
    assert fit_res.x0_err == 0.0
    assert fit_res.y0_err == 0.0

    # Test too small grid
    with raises(ValueError):
        grid_small, ref_small = generate_exact_pf2(1.23, 2.34, N=4)
        PolynomialFit(4).find_minimum(grid_small)
    # end with

    # Find noisy minimum (elevate by y_offset using Gs)
    y_offset = 2.0
    Gs = y_offset * ones((20, 5))
    fit_noisy = fit.find_noisy_minimum(grid, Gs=Gs)
    assert isinstance(fit_noisy, PolynomialResult)
    assert fit_noisy.analyzed
    assert match_to_tol(fit_noisy.x0, ref.x0)
    assert match_to_tol(fit_noisy.y0, ref.y0)
    # The fluctuation is numerically zero
    assert fit_noisy.x0_err == 0.0
    assert fit_noisy.y0_err == 0.0

    # Test __eq__
    assert fit == PolynomialFit(2)
    assert fit != PolynomialFit(3)
    assert fit != []
    assert fit != FittingFunction(minimize_pf, {'pfn': 3})

# end def
