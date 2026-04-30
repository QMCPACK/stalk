#!/usr/bin/env python

from pytest import raises
import numpy as np

from stalk.util.util import match_to_tol

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test ParameterSet class
def test_ParameterSet():
    from stalk import ParameterSet
    from stalk.params import Parameter

    # nominal test empty/defaults
    s = ParameterSet([])
    assert s._param_list == []
    assert len(s) == 0
    assert s.params is None
    assert s.params_err is None
    assert s.value is None
    assert s.error == 0.0
    assert s.label is None

    # nominal test, meaningful values
    params = [1, 2]
    params_err = [3, 4]
    units = ['a', 'b']
    labels = ['c', 'd']
    value = 6.
    error = 7.
    label = 'e'
    s = ParameterSet(
        params,
        params_err,
        units=units,
        value=value,
        error=error,
        label=label,
        labels=labels
    )
    assert len(s) == 2
    assert s.value == value
    assert s.error == error
    assert s.label == label
    for param, param_ref in zip(s.params_list, params):
        assert isinstance(param, Parameter)
        assert param.value == param_ref
    # end for
    for param, param_ref in zip(s.params, params):
        assert param == param_ref
    # end for
    for param_err, param_err_ref in zip(s.params_err, params_err):
        assert param_err == param_err_ref
    # end for

    # Test copy original
    s_copy = s.copy()
    assert len(s_copy) == 2
    assert s_copy.value == value
    assert s_copy.error == error
    assert s_copy.label == label
    for param, param_ref in zip(s_copy.params_list, params):
        assert isinstance(param, Parameter)
        assert param.value == param_ref
    # end for
    for param, param_ref in zip(s_copy.params, params):
        assert param == param_ref
    # end for
    for param_err, param_err_ref in zip(s_copy.params_err, params_err):
        assert param_err == param_err_ref
    # end for

    # Test copy with changes
    copy_params = [-1., -2.]
    copy_params_err = [-3., -4.]
    copy_label = 'g'
    s_copy2 = s.copy(params=copy_params, params_err=copy_params_err, label=copy_label)
    assert s_copy2.label == copy_label
    for param, param_ref in zip(s_copy2.params, copy_params):
        assert param == param_ref
    # end for
    for param_err, param_err_ref in zip(s_copy2.params_err, copy_params_err):
        assert param_err == param_err_ref
    # end for

    # Testing shifting of parameters
    ss = s.copy()
    shifts = [0.1, 0.2]
    with raises(ValueError):
        ss.shift_params([1])
    # end with
    ss.shift_params(shifts)
    for param, param_ref, shift in zip(ss.params, params, shifts):
        assert param == param_ref + shift
        assert ss.value is None
    # end for
    # Test difference and distance
    assert match_to_tol((ss - s).params, shifts)
    assert match_to_tol(s.distance2(ss), sum(np.array(shifts)**2))
    assert match_to_tol(s.distance(ss), sum(np.array(shifts)**2)**0.5)

    # Test setting of parameters (w/o error)
    new_params = [10., 11.]
    s.set_params(new_params, None)
    assert match_to_tol(s.params, new_params)
    assert match_to_tol(s.params_err, [0.0, 0.0])
    assert s.value is None
    # Test setting of parameters (w/ error)
    new_params_err = [12., 13.]
    s.set_params(new_params, new_params_err)
    assert match_to_tol(s.params, new_params)
    assert match_to_tol(s.params_err, new_params_err)
    assert s.value is None

    # Test setting of value
    new_value = 100.0
    new_error = 101.0
    s.value = new_value
    assert s.value == new_value and s.error == 0.0
    s.error = new_error
    assert s.value == new_value and s.error == new_error

    # test (placeholder) consistency check
    assert s.check_consistency()

    # TODO: test minimize

# end def
