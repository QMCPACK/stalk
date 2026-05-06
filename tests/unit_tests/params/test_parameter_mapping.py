#!/usr/bin/env python

from stalk.util import match_to_tol
from stalk.params import ParameterMapping

from ..assets.h2 import pos_H2, forward_H2, backward_H2

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test ParameterMapping class
def test_ParameterMapping():

    # Test empty/default initialization
    mp = ParameterMapping()
    assert mp.forward is None
    assert mp.backward is None
    assert mp.dim == 3
    assert mp.incomplete

    # Test nominal initialization (using H2 1-parameter model, pos init)
    fwd_args = {'fwd': 1}
    bck_args = {'bck': 2}
    mp_H2 = ParameterMapping(
        forward_func=forward_H2,
        backward_func=backward_H2,
        forward_args=fwd_args,
        backward_args=bck_args,
    )
    assert mp_H2.forward.func == forward_H2
    assert mp_H2.forward.args == fwd_args
    assert mp_H2.backward.func == backward_H2
    assert mp_H2.backward.args == bck_args

    # Test mapping
    tol = 1e-7
    params_H2 = [1.4]
    pos, axes = mp_H2.map_backward(params_H2)
    assert match_to_tol(pos_H2, pos, tol)
    assert axes is None
    # Supplying axes should make no difference
    params = mp_H2.map_forward(pos_H2, axes=None)
    assert match_to_tol(params_H2, params, tol)

    assert mp_H2.check_params_consistency([1.6])
    assert mp_H2.check_pos_consistency(pos_H2 * 2)
# end def
