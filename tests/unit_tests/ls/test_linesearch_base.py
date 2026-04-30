#!/usr/bin/env python

from stalk.ls import PolynomialFit
from stalk.ls.fitting_result import FittingResult
from stalk.ls.linesearch_base import LineSearchBase
from stalk.util.util import match_to_tol

from ..assets.fitting_pf2 import generate_exact_pf2

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test LineSearchBase class
def test_LineSearchBase():

    # test empty/defaults
    ls = LineSearchBase()
    assert ls.settings.sgn == 1
    assert ls.x0 is None
    assert ls.x0_err is None
    assert ls.y0 is None
    assert ls.y0_err is None
    assert ls.settings.fraction == 0.025
    assert ls.settings.fit_func == PolynomialFit(3)

    # Test initialization with grid and values, no noise
    grid, ref = generate_exact_pf2(1.23, 2.34, N=5)
    ls_val = LineSearchBase(offsets=grid.offsets, values=grid.values, fit_kind='pf2')
    assert ls_val.fit_res.analyzed
    assert match_to_tol(ls_val.x0, ref.x0)
    assert match_to_tol(ls_val.y0, ref.y0)
    assert ls_val.x0_err == 0.0
    assert ls_val.y0_err == 0.0

    # Test initialization with grid, values and noise
    fraction = 0.1
    grid, ref = generate_exact_pf2(1.23, 2.34, N=5, error=0.1)
    ls_noisy = LineSearchBase(
        offsets=grid.offsets,
        values=grid.values,
        errors=grid.errors,
        fraction=fraction,
        fit_kind='pf2',
        N=10
    )
    assert ls_noisy.fit_res.analyzed
    assert match_to_tol(ls_noisy.x0, ref.x0)
    assert match_to_tol(ls_noisy.y0, ref.y0)
    assert ls_noisy.x0_err > 0.0
    assert ls_noisy.y0_err > 0.0
    # Test search
    # TODO: add coverage
    res = ls_noisy.search()
    assert isinstance(res, FittingResult)
    assert match_to_tol(res.x0, ref.x0)
    assert match_to_tol(res.y0, ref.y0)
    assert res.x0_err == 0.0
    assert res.y0_err == 0.0
    # Test search with error
    # TODO: add coverage
    res = ls_noisy.search_with_error()
    assert isinstance(res, FittingResult)
    assert match_to_tol(res.x0, ref.x0)
    assert match_to_tol(res.y0, ref.y0)
    assert res.x0_err > 0.0
    assert res.y0_err > 0.0

# end def
