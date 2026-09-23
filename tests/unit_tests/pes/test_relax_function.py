#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array

from stalk.params.parameter_set import ParameterSet
from stalk.pes.relax_function import RelaxFunction
from stalk.util.util import match_to_tol


def test_RelaxFunction(tmp_path):

    # Minimal function for testing vs callable
    def relax_func(s: ParameterSet, shift=[0.0], value=0.0, error=0.0):
        s.shift_params(shift)
        s.value = value
        s.error = error
        return s
    # end def

    # Test nominal
    value = 5.0
    error = 0.1
    shift = [0.2, 0.2, 0.2]
    rf = RelaxFunction(relax_func, value=value, error=error, shift=shift)

    params = array([1.0, 2.0, 3.0])
    s = ParameterSet(params)
    rf.evaluate(s)
    assert match_to_tol(s.params, params + shift)
    assert match_to_tol(s.value, value)
    assert match_to_tol(s.error, error)

    # Call and write to disk
    s.label = "test1"
    s1 = rf(structure=s, path=tmp_path)
    # Call again and load
    scp = s.copy(label="test1")
    scp.reset_value()
    s2 = rf(structure=scp, path=tmp_path)
    assert match_to_tol(s2.params, s1.params)
    assert match_to_tol(s2.value, s1.value)
    assert match_to_tol(s2.error, s1.error)

# end def
