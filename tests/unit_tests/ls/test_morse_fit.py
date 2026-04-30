#!/usr/bin/env python


from pytest import raises

from stalk.ls.linesearch_grid import LineSearchGrid
from stalk.ls.morse_fit import MorseFit
from stalk.ls.morse_result import MorseResult
from stalk.util.util import match_to_tol

from ..assets.fitting_pf2 import generate_exact_morse

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test MorseFit class
def test_MorseFit():

    # Test nominal using regular 2-degree polynomial fit
    x0 = 1.23
    y0 = 2.1
    grid, ref = generate_exact_morse(x0, 0.1, y0, N=9)
    fit = MorseFit()
    # Find minimum, noise not requested
    fit_res = fit.find_minimum(grid)

    assert isinstance(fit_res, MorseResult)
    assert fit_res.analyzed
    assert match_to_tol(fit_res.x0, ref.x0, 1e-4)
    assert match_to_tol(fit_res.y0, ref.y0, 1e-4)
    assert match_to_tol(fit_res.fit, ref.fit, 1e-4)
    assert fit_res.x0_err == 0.0
    assert fit_res.y0_err == 0.0

    # Find minimum of a grid centered at minimum
    grid_cent = LineSearchGrid(grid.offsets - x0)
    grid_cent.values = grid.values
    fit_res = fit.find_minimum(grid_cent)
    assert match_to_tol(fit_res.x0, 0.0, 1e-4)
    assert match_to_tol(fit_res.y0, ref.y0, 1e-4)

    # Test too small grid
    with raises(ValueError):
        grid_small, ref_small = generate_exact_morse(x0, 0.1, y0, N=2)
        MorseFit().find_minimum(grid_small)
    # end with

    # TODO: test with ultrasmall grid?

    # Test __eq__
    assert fit == MorseFit()

# end def
