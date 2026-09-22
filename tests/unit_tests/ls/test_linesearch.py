#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises
from numpy import linspace, array

from stalk import LineSearch
from stalk.params.parameter_set import ParameterSet
from stalk.params.parameter_structure import ParameterStructure
from stalk.util import match_to_tol
from ..assets.h2o import get_structure_H2O, get_hessian_H2O, h2o_pes


# test LineSearch class
def test_LineSearch():

    # Empty init
    ls = LineSearch()
    assert ls.structure is None
    assert ls.direction is None

    # Test default init with structure, no direction
    p = ParameterSet([1.0, 2.0])
    ls = LineSearch(p)
    assert ls.structure == p
    assert ls.direction is None
    with raises(AssertionError):
        ls.grid = [0.1, 0.2, 0.3]
    # end with

    direction = [0.0, 1.0]
    ls = LineSearch(p, direction=direction)
    assert ls.structure == p
    assert match_to_tol(ls.direction, direction)
    assert ls.d is None
    assert ls.Lambda is None
    assert ls.W_max is None
    assert ls.valid_W_max is None
    assert ls.R_max == 0.0
    assert ls.valid_R_max is None
    assert len(ls) == 0
    assert len(ls.grid) == 0
    assert len(ls.offsets) == 0
    assert len(ls.values) == 0
    assert len(ls.errors) == 0
    assert not ls.shifted
    assert ls.shifted_params is None

    # Test reset of the grid
    with raises(ValueError):
        # M cannot be negative
        ls.figure_out_offsets(M=-1)
    # end with

    # Test with Lambda, W
    Lambda = 0.5
    W = 0.1
    M = 5
    ls.Lambda = Lambda
    with raises(ValueError):
        ls.reset_offsets(M=M, W=-W)
    # end with
    ls.reset_offsets(M=M, W=W)
    assert len(ls) == M
    assert match_to_tol(ls.W_max, W)
    R_ref = (2 * W / Lambda)**0.5
    assert match_to_tol(ls.offsets, linspace(-R_ref, R_ref, 5))
    assert ls.shifted
    assert not ls.evaluated

    # Test init with offset, values, errors
    offsets = [-0.1, 0.0, 0.1, 0.3]
    values = [1.0, 2.0, 3.0, 4.0]
    errors = [0.1, 0.2, 0.3, 0.4]
    ls = LineSearch(
        p,
        direction=direction,
        offsets=offsets,
        values=values,
        errors=errors,
        R=0.4
    )
    assert ls.shifted
    assert ls.evaluated
    assert len(ls) == 4
    assert match_to_tol(ls.offsets, offsets)
    assert match_to_tol(ls.values, values)
    assert match_to_tol(ls.errors, errors)

    # Test nominal init using actual structure
    structure = get_structure_H2O()
    d = 1
    R = 0.2
    sigma = 3.0
    M = 5
    offsets_ref = linspace(-R, R, M)
    direction_ref = array([0.0, 1.0])
    params_ref = structure.params[d] + offsets_ref
    ls_s = LineSearch(
        structure=structure,
        direction=direction_ref,
        M=M,
        d=d,
        sigma=sigma,
        R=R
    )
    assert ls_s.structure == structure
    assert match_to_tol(ls_s.direction, direction_ref)
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

    # Test nominal init using Hessian
    hessian = get_hessian_H2O()
    W = 0.2
    d = 1
    M = 9
    ls_h = LineSearch(hessian=hessian, M=M, d=d, sigma=sigma, W=W)
    assert len(ls_h) == M
    assert ls_h.structure == hessian.structure
    assert match_to_tol(ls_h.direction, hessian.directions[d])
    assert ls_h.d == 1
    assert ls_h.W_max == W
    assert ls_h.valid_W_max is None
    assert ls_h.Lambda == hessian.lambdas[d]

    # test evaluation, fitting etc
    assert ls_h.shifted
    assert not ls_h.evaluated
    # Evaluate
    h2o_pes(ls_h, path=None)
    assert ls_h.evaluated
    assert ls_h.valid

# end def
