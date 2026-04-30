#!/usr/bin/env python

from pytest import raises

from stalk.params.linesearch_point import LineSearchPoint

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test LineSearchPoint class
def test_LineSearchPoint():

    # Empty init leads to error
    with raises(TypeError):
        LineSearchPoint()
    # end with

    # Test faulty format
    with raises(TypeError):
        LineSearchPoint([0])
    # end with

    # Test nominal grid init
    offset = 1.0
    point = LineSearchPoint(offset)
    assert point.offset == offset
    assert point.value is None
    assert point.error == 0.0
    assert point.enabled
    assert not point.valid

    # Test grid init with faulty value
    with raises(ValueError):
        LineSearchPoint(offset, value=[0.0])
    # end with

    # Test nominal grid init with value
    value = 2.0
    point = LineSearchPoint(offset, value=value)
    assert point.offset == offset
    assert point.value == value
    assert point.error == 0.0
    assert point.enabled
    assert point.valid

    # Test grid init with faulty error
    with raises(ValueError):
        LineSearchPoint(offset, error=[0.0])
    # end with
    with raises(ValueError):
        LineSearchPoint(offset, error=-1.0e-9)
    # end with

    # Test nominal grid init with value and error
    value = 2.0
    error = 3.0
    point = LineSearchPoint(offset, value=value, error=error)
    assert point.offset == offset
    assert point.value == value
    assert point.error == error
    assert point.enabled
    assert point.valid

    # Test comparison operator
    assert LineSearchPoint(1.0) == LineSearchPoint(1.0 + 1e-10)
    assert not LineSearchPoint(1.0) == LineSearchPoint(1.0 - 1e-8)

    # Test ordering operator
    assert LineSearchPoint(1.0) < LineSearchPoint(1.0 + 1e-10)
    assert not LineSearchPoint(1.0) < LineSearchPoint(1.0 - 1e-10)

# end def
