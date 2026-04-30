#!/usr/bin/env python

from numpy import ones

from pytest import raises
from stalk.ls.fitting_function import FittingFunction
from stalk.ls.fitting_result import FittingResult
from stalk.util.util import match_to_tol

from ..assets.fitting_pf2 import generate_exact_pf2, minimize_pf

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test FittingFunction class
def test_FittingFunction():

    # Function must be callable
    with raises(ValueError):
        FittingFunction(None, None)
    # end with

    # Test nominal using regular 2-degree polynomial fit
    grid, ref = generate_exact_pf2(1.23, 2.34, N=5, error=0.1)
    # Add point that remains disabled and ignored in proper implementation
    grid.add_point(10.0)
    fit = FittingFunction(minimize_pf, {'pfn': 2})
    # Find minimum, noise not requested
    fit_res = fit.find_minimum(grid)
    assert isinstance(fit_res, FittingResult)
    assert fit_res.analyzed
    assert match_to_tol(fit_res.x0, ref.x0)
    assert match_to_tol(fit_res.y0, ref.y0)
    assert fit_res.x0_err == 0.0
    assert fit_res.y0_err == 0.0

    # Find noisy minimum (elevate by y_offset using Gs)
    y_offset = 2.0
    Gs = y_offset * ones((20, 5))
    fit_noisy = fit.find_noisy_minimum(grid, Gs=Gs)
    assert isinstance(fit_noisy, FittingResult)
    assert fit_noisy.analyzed
    assert match_to_tol(fit_noisy.x0, ref.x0)
    assert match_to_tol(fit_noisy.y0, ref.y0)
    # The fluctuation is numerically zero
    assert fit_noisy.x0_err == 0.0
    assert fit_noisy.y0_err == 0.0

    # Test getting distribution (TODO: incomprehensive in statistical sense!)
    N = 10
    x0s, y0s = fit.get_distribution(grid, N=N)
    assert len(x0s) == N
    assert len(y0s) == N
    x0s2 = fit.get_x0_distribution(grid, N=N)
    y0s2 = fit.get_y0_distribution(grid, N=N)
    assert len(x0s2) == N
    assert len(y0s2) == N

    # Test __eq__
    assert fit == FittingFunction(minimize_pf, {'pfn': 2})
    assert fit != []
    assert fit != FittingFunction(minimize_pf, {'pfn': 3})

# end def
