#!/usr/bin/env python

import io
import sys
from numpy import inf
from pytest import raises

from stalk.params import Parameter
from stalk.params import BondLength
from stalk.params.parameter import ParameterLimitException
from stalk.util.util import match_to_tol

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test Parameter class
def test_Parameter():

    # Cannot construct without a value
    with raises(TypeError):
        Parameter()
    # end with

    # Cannot construct with non-scalar
    with raises(ValueError):
        Parameter([])
    # end with

    # test empty
    p = Parameter(0.0)
    assert p.value == 0.0
    assert p.error == 0.0
    assert p.unit == ''
    assert p.label == 'p'
    assert p.limits[0] == -inf
    assert p.limits[1] == inf
    
    with raises(TypeError):
        p.limits = 2.0
    # end with

    # test nominal
    value = 1.0
    error = 2.0
    limits = (-3.0, 5.0)
    label = 'label'
    unit = 'unit'
    p = Parameter(value, error, label=label, unit=unit, limits=limits)
    assert p.value == value
    assert p.error == error
    assert p.label == label
    assert p.unit == unit
    assert p.limits == limits

    # test shifting
    shift = 3.0
    p.shift(shift)
    assert p.value == value + shift
    assert p.error == 0.0
    
    # Trying to shift further violates limits
    with raises(ParameterLimitException):
        p.shift(limits[1])
    # end with
    with raises(ParameterLimitException):
        p.shift(-100.0)
    # end with
    # Original value is preserved
    assert p.value == value + shift

    # test printing
    test_stdout = io.StringIO()
    sys.stdout = test_stdout
    print(p)
    sys.stdout = sys.__stdout__
    param_str = test_stdout.getvalue()
    param_str_ref = 'label         4.0000 unit       \n'
    assert param_str == param_str_ref
# end def


# Test BondLength class
def test_BondLength():

    # Cannot construct without a value
    with raises(TypeError):
        BondLength()
    # end with

    # Cannot construct with non-scalar
    with raises(ValueError):
        BondLength([])
    # end with

    # test defaults
    pairs = ([0.0, 1.0], [1.0, 0.0])
    p = BondLength(pairs)
    assert match_to_tol(p.value, 2.0**0.5)
    assert p.error == 0.0
    assert p.unit == 'A'
    assert p.label == 'd'
    assert p.limits[0] == 0.0
    assert p.limits[1] == inf

    # test multiple pairs
    pairs = [([0.0, 1.0], [1.0, 0.0]), ([3.0, -5.0], [4.0, -4.0])]
    limits = (1.0, 5.0)
    label = 'label'
    unit = 'unit'
    p = BondLength(pairs, label=label, unit=unit, limits=limits)
    assert match_to_tol(p.value, 2.0**0.5)
    assert p.error == 0.0
    assert p.label == label
    assert p.unit == unit
    assert p.limits == limits

# end def