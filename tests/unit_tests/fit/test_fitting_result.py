#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises
from stalk.fit.fitting_result import FittingResult
from stalk.util.util import match_to_tol


# Test FittingResult class
def test_FittingResult(tmp_path):

    # Test empty init
    res = FittingResult()
    assert res.x0 is None
    assert res.y0 is None
    assert res.x0_err == 0.0
    assert res.y0_err == 0.0
    assert res.fit is None
    assert not res.analyzed
    with raises(AssertionError):
        res.save_result(tmp_path / 'test')
    # end with

    # test default init
    x0 = 1.0
    y0 = 2.0
    res = FittingResult(x0, y0)
    assert res.x0 == x0
    assert res.y0 == y0
    assert res.analyzed

    # test full init
    x0_err = 3.0
    y0_err = 4.0
    fit = [6.0, 7.0]
    fraction = 0.5
    res_full = FittingResult(
        x0,
        y0,
        x0_err=x0_err,
        y0_err=y0_err,
        fraction=fraction,
        fit=fit
    )
    assert res_full.analyzed
    assert res_full.x0 == x0
    assert res_full.y0 == y0
    assert res_full.x0_err == x0_err
    assert res_full.y0_err == y0_err
    assert match_to_tol(res_full.fit, fit)
    assert res_full.fraction == fraction

    # Write to a temporary path
    res_full.save_result(tmp_path / 'test')
    # Load to a new instance
    res_loaded = FittingResult()
    res_loaded.load_result(tmp_path / 'test')
    assert res_loaded.analyzed
    assert res_loaded.x0 == x0
    assert res_loaded.y0 == y0
    assert res_loaded.x0_err == x0_err
    assert res_loaded.y0_err == y0_err
    assert match_to_tol(res_loaded.fit, fit)

    # the closed form is generally not implemented
    with raises(NotImplementedError):
        res_full.get_values(0.0)
    # end with
    with raises(NotImplementedError):
        res_full.get_force(0.0)
    # end with
    with raises(NotImplementedError):
        res_full.get_hessian(0.0)
    # end with

# end def
