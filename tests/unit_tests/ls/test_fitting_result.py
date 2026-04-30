#!/usr/bin/env python

from pytest import raises
from stalk.ls.fitting_result import FittingResult

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test FittingResult class
def test_FittingResult():

    # Must provide x0, y0
    with raises(TypeError):
        FittingResult()
    # end with
    with raises(TypeError):
        FittingResult(0.0)
    # end with

    # test default init
    x0 = 1.0
    y0 = 2.0
    res = FittingResult(x0, y0)
    assert res.x0 == x0
    assert res.y0 == y0
    assert res.analyzed
    # test not analyzed
    res.x0 = None
    assert not res.analyzed

    # test full init
    x0_err = 3.0
    y0_err = 4.0
    fit = ''
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
    assert res_full.fit == fit
    assert res_full.fraction == fraction

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
