#!/usr/bin/env python3

from pytest import raises
from stalk.params.PesResult import PesResult

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


def test_PesResult():
    # Test degraded
    # Cannot init empty
    with raises(TypeError):
        PesResult()
    # end with

    val = 1.0
    err = 2.0
    err_default = 0.0
    sigma = 3.0

    # test nominal (no error)
    res0 = PesResult(val)
    assert res0.value == val
    assert res0.error == err_default
    res0.add_sigma(sigma)
    assert res0.error == sigma
    # Only test that value has been changed but not by how much
    assert res0.value != val

    # Test nominal (with error)
    res1 = PesResult(val, err)
    assert res1.value == val
    assert res1.error == err
    # Add zero sigma and expect no effect
    res1.add_sigma(0.0)
    assert res1.value == val
    assert res1.error == err
    res1.add_sigma(sigma)
    assert res1.error == (err**2 + sigma**2)**0.5
    # Only test that value has been changed but not by how much
    assert res1.value != val

    # Test rescale
    res2 = PesResult(val, err)
    scale = 1.234
    res2.rescale(scale)
    assert res2.value == val / scale
    assert res2.error == err / scale

    with raises(ValueError):
        res1.add_sigma([])
    # end with
    with raises(ValueError):
        res1.add_sigma(-1e-9)
    # end with

    with raises(ValueError):
        PesResult(0.0, [])
    # end with
    with raises(ValueError):
        PesResult(0.0, -1e-9)
    # end with

# end def
