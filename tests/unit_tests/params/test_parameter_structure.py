#!/usr/bin/env python

from pytest import raises
from stalk.params.util import distance
from stalk.util import match_to_tol
from stalk.params import ParameterStructure

from ..assets.h2 import backward_H2_alt, forward_H2_alt, pos_H2, forward_H2, backward_H2, elem_H2
from ..assets.gese import params_GeSe, forward_GeSe, backward_GeSe, elem_GeSe, pos_GeSe, axes_GeSe

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test ParameterStructure class
def test_ParameterStructure_open():

    # Test empty/default initialization
    s_empty = ParameterStructure()
    assert s_empty.forward_func is None
    assert s_empty.forward_args == {}
    assert s_empty.backward_func is None
    assert s_empty.backward_args == {}
    assert s_empty.pos is None
    assert s_empty.axes is None
    assert s_empty.elem is None
    assert s_empty.dim == 3
    assert s_empty.params is None
    assert s_empty.params_err is None
    assert s_empty.value is None
    assert s_empty.error == 0.0
    assert s_empty.label == ''
    assert s_empty.units == 'B'
    assert not s_empty.consistent
    assert not s_empty.periodic
    assert s_empty.tol == 1e-7

    # Test nominal initialization (using H2 1-parameter model, pos init)
    value = 1.0
    error = 2.0
    tol = 0.1
    fwd_args = {'fwd': 1}
    bck_args = {'bck': 2}
    label = 'H2'
    units = 'units'
    s_H2 = ParameterStructure(
        forward=forward_H2,
        backward=backward_H2,
        forward_args=fwd_args,
        backward_args=bck_args,
        pos=pos_H2,
        elem=elem_H2,
        value=value,
        error=error,
        label=label,
        tol=tol,
        units=units
    )
    assert s_H2.forward_func == forward_H2
    assert s_H2.forward_args == fwd_args
    assert s_H2.backward_func == backward_H2
    assert s_H2.backward_args == bck_args
    assert match_to_tol(s_H2.pos, pos_H2, tol)
    assert s_H2.axes is None
    for el, el_ref in zip(s_H2.elem, elem_H2):
        assert el == el_ref
    # end for
    assert s_H2.dim == 3
    assert match_to_tol(s_H2.params, [1.4], tol)
    assert match_to_tol(s_H2.params_err, [0.0], tol)
    assert s_H2.value == value
    assert s_H2.error == error
    assert s_H2.label == label
    assert s_H2.units == units
    assert s_H2.consistent
    assert not s_H2.periodic
    assert s_H2.tol == tol

    # test setting position
    with raises(AssertionError):
        # Dimensions must be right
        s_H2.set_position(pos_H2[:, :-1])
    # end with
    # translate=False
    s_H2.set_position(pos_H2 + 2, translate=False)
    assert match_to_tol(s_H2.params, [1.4], tol)
    assert match_to_tol(s_H2.pos, pos_H2 + 2, tol)
    # translate=True
    s_H2.set_position(pos_H2 + 2, translate=True)
    assert match_to_tol(s_H2.params, [1.4], tol)
    assert match_to_tol(s_H2.pos, pos_H2, tol)

    # setting axes causes TypeError
    with raises(TypeError):
        s_H2.copy().set_axes([1.0, 2.0, 3.0])
    # end with

    # test shifting of position and pos difference
    shift_scalar = 0.2
    shift_array = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    shift_array_bad = [0.1, 0.2]
    # If pos is not set
    with raises(AssertionError):
        s_H2_shift = s_H2.copy()
        s_H2_shift.pos = None
        s_H2_shift.shift_pos(shift_scalar)
    # end with
    s_H2_shift = s_H2.copy()
    s_H2_shift.shift_pos(shift_scalar, translate=False)
    assert match_to_tol(-s_H2_shift.pos_difference(pos_H2), 6 * [shift_scalar], tol)
    # maintain parameter consistency
    assert match_to_tol(s_H2_shift.params, [1.4], tol)
    s_H2_shift = s_H2.copy()
    s_H2_shift.shift_pos(shift_array, translate=False)
    assert match_to_tol(-s_H2_shift.pos_difference(pos_H2), shift_array, tol)
    # maintain parameter consistency
    new_param = [distance(s_H2_shift.pos[0], s_H2_shift.pos[1])]
    assert match_to_tol(s_H2_shift.params, new_param, tol)
    with raises(AssertionError):
        s_H2_shift.shift_pos(shift_array_bad)
    # end with

    # test shifting of parameters
    pshift = 2.0
    s_H2_pshift = s_H2.copy()
    s_H2_pshift.shift_params([pshift])
    assert match_to_tol(s_H2_pshift.params, [1.4 + pshift], tol)
    new_pos = s_H2.pos.copy()
    new_pos[0, 2] += pshift / 2
    new_pos[1, 2] -= pshift / 2
    assert match_to_tol(s_H2_pshift.pos, new_pos, tol)

    # test shifting in dpos mode
    s_H2_dpos = s_H2.copy()
    dpos = 1.23
    pshift = 2.0
    s_H2_dpos.pos += 1.23
    s_H2_dpos.shift_params([pshift], dpos_mode=True)
    assert match_to_tol(s_H2_pshift.params, [1.4 + pshift], tol)
    new_pos = s_H2.pos.copy() + dpos
    new_pos[0, 2] += pshift / 2
    new_pos[1, 2] -= pshift / 2
    assert match_to_tol(s_H2_dpos.pos, new_pos, tol)

    # test setting of alternative forward func (params are divided by factor)
    factor = 3
    s_H2.set_forward_func(forward_H2_alt, {'factor': factor})
    assert match_to_tol(s_H2.params, [1.4 / factor], tol)
    assert match_to_tol(s_H2.params_err, [0.0], tol)
    # Not consistent...
    assert not s_H2.consistent
    # ...until matching mapping is provided
    s_H2.set_backward_func(backward_H2_alt, {'factor': factor})
    assert match_to_tol(s_H2.pos, pos_H2, tol)
    assert s_H2.consistent

# end def


def test_ParameterStructure_periodic():

    value = 1.0
    error = 2.0
    tol = 0.1
    fwd_args = {'fwd': 1}
    bck_args = {'bck': 2}
    label = 'GeSe'
    units = 'units'
    s_GeSe = ParameterStructure(
        forward=forward_GeSe,
        backward=backward_GeSe,
        forward_args=fwd_args,
        backward_args=bck_args,
        params=params_GeSe,
        elem=elem_GeSe,
        value=value,
        error=error,
        label=label,
        tol=tol,
        units=units
    )
    assert s_GeSe.forward_func == forward_GeSe
    assert s_GeSe.forward_args == fwd_args
    assert s_GeSe.backward_func == backward_GeSe
    assert s_GeSe.backward_args == bck_args
    assert match_to_tol(s_GeSe.pos, pos_GeSe, tol)
    assert match_to_tol(s_GeSe.axes, axes_GeSe, tol)
    for el, el_ref in zip(s_GeSe.elem, elem_GeSe):
        assert el == el_ref
    # end for
    assert s_GeSe.dim == 3
    assert match_to_tol(s_GeSe.params, params_GeSe, tol)
    assert match_to_tol(s_GeSe.params_err, 5 * [0.0], tol)
    assert s_GeSe.value == value
    assert s_GeSe.error == error
    assert s_GeSe.label == label
    assert s_GeSe.units == units
    assert s_GeSe.consistent
    assert s_GeSe.periodic
    assert s_GeSe.tol == tol

# end def
