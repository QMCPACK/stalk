#!/usr/bin/env python

import io
import sys
from numpy import inf, pi
from pytest import raises

from stalk.params import Parameter
from stalk.params import BondLength
from stalk.params.parameter import BondAngle, ParameterLimitException, PhaseAngle
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

    # Cannot construct with empty non-scalar
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


# Test BondAngle class
def test_BondAngle():

    # Cannot construct without a value
    with raises(TypeError):
        BondAngle()
    # end with

    # Cannot construct with non-scalar
    with raises(ValueError):
        BondAngle([])
    # end with

    # test defaults
    triplets = ([0.0, 1.0], [1.0, 0.0], [0.0, -1.0])
    p = BondAngle(triplets)
    assert match_to_tol(p.value, 90.0)
    assert p.error == 0.0
    assert p.unit == 'ang'
    assert p.label == 'a'
    assert p.limits[0] == 0.0
    assert p.limits[1] == 180.0

    # test multiple angles
    triplets = [([0.0, 1.0], [1.0, 0.0], [0.0, -1.0]), ([0.0, 1.0], [-1.0, 0.0], [0.0, -1.0])]
    limits = (1.0, 2.0)
    label = 'label'
    unit = 'rad'
    p = BondAngle(triplets, label=label, unit=unit, limits=limits)
    assert match_to_tol(p.value, pi / 2)
    assert p.error == 0.0
    assert p.label == label
    assert p.unit == unit
    assert p.limits == limits

# end def


# Test PhaseAngle class
def test_PhaseAngle():

    # Cannot construct without a value
    with raises(TypeError):
        PhaseAngle()
    # end with

    # Cannot construct with non-scalar
    with raises(ValueError):
        PhaseAngle([])
    # end with

    # test defaults
    value = 1.0
    p = PhaseAngle(value)
    assert p.value == value
    assert p.error == 0.0
    assert p.unit == 'rad'
    assert p.label == 't'
    assert p.limits[0] == -pi
    assert p.limits[1] == pi

    # test alternative units
    label = 'label'
    unit = 'ang'
    value = 512.0
    p = PhaseAngle(value, label=label, unit=unit)
    assert p.value == value % 360.0
    assert p.error == 0.0
    assert p.label == label
    assert p.unit == unit
    assert p.limits == (-180.0, 180.0)
    # test negative angles
    p.value = 200.0
    assert p.value == -160.0

# end def
