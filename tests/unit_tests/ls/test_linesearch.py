#!/usr/bin/env python

from pytest import raises
from numpy import linspace

from stalk import LineSearch
from stalk.params.LineSearchPoint import LineSearchPoint
from stalk.params.ParameterStructure import ParameterStructure
from stalk.util import match_to_tol
from ..assets.h2o import get_structure_H2O, get_hessian_H2O

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# test LineSearch class
def test_LineSearch():

    # Empty init
    ls = LineSearch()
    assert ls.d is None
    assert ls.settings.sgn == 1
    assert ls.direction == 0.0
    assert ls.structure is None
    assert ls.hessian is None
    assert ls.W_max is None
    assert ls.R_max == 0.0
    assert len(ls) == 0
    assert len(ls.grid) == 0
    assert len(ls.offsets) == 0
    assert len(ls.values) == 0
    assert len(ls.errors) == 0
    assert ls.shifted_params is None
    with raises(ValueError):
        # M cannot be negative
        ls.figure_out_offsets(M=-1)
    # end with
    with raises(ValueError):
        # Must characterize grid somehow
        ls.figure_out_offsets(M=5)
    # end with
    with raises(ValueError):
        # Cannot use negative R
        ls.figure_out_offsets(M=5, R=-0.1)
    # end with
    with raises(ValueError):
        # Cannot use W before setting Hessian
        ls.figure_out_offsets(M=5, W=0.1)
    # end with
    ls = LineSearch(d=1)
    assert ls.d == 1

    # Without structure, the grid is only abstract points
    ls.set_grid(offsets=[0.1, 0.0, -0.1])
    for point in ls.grid:
        assert isinstance(point, LineSearchPoint)
    # end for
    ls.add_shift(0.2)
    assert len(ls) == 4

    # Test nominal init using actual structure
    structure = get_structure_H2O()
    d = 1
    R = 0.2
    sigma = 3.0
    M = 5
    offsets_ref = linspace(-R, R, M)
    params_ref = structure.params[d] + offsets_ref
    ls_s = LineSearch(structure=structure, M=M, d=d, sigma=sigma, R=R)
    assert ls_s.structure == structure
    assert len(ls_s) == M
    assert ls_s.d == 1
    assert ls_s.sigma == sigma
    assert ls_s.W_max is None
    assert match_to_tol(ls_s.R_max, R)
    assert match_to_tol(ls_s.direction, [0.0, 1.0])
    params = ls_s.shifted_params
    for point, ref in zip(ls_s.grid, offsets_ref):
        assert isinstance(point, ParameterStructure)
        assert point.offset == ref
        assert match_to_tol(point.params[d] - structure.params[d], ref)
    # end for
    for params, ref in zip(ls_s.shifted_params, params_ref):
        assert match_to_tol(params[d], ref)
    # end for
    with raises(ValueError):
        ls_s.d = 2
    # end with

    # Test nominal init using Hessian
    hessian = get_hessian_H2O()
    W = 0.2
    d = 1
    M = 9
    ls_h = LineSearch(hessian=hessian, M=M, d=d, sigma=sigma, W=W)
    assert len(ls_h) == M
    assert ls_h.structure == hessian.structure
    assert ls_h.hessian == hessian
    assert ls_h.d == 1
    assert ls_h.sgn == 1
    assert ls_h.W_max == W
    assert ls_h.valid_W_max == 0.0
    assert ls_h.Lambda == hessian.lambdas[d]
    with raises(ValueError):
        ls_h.sigma = -1.0
    # end with
    with raises(ValueError):
        ls_h.d = 2
    # end with
    with raises(ValueError):
        ls_h.sigma = []
    # end with

    # TODO: test evaluation, fitting etc

# end def
