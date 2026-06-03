#!/usr/bin/env python

from pytest import raises
from numpy import nan, ndarray, isnan

from stalk.io.txt_data import TxtData

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test PesLoader class
def test_TxtData():

    # Test with arbitrary value
    pl = TxtData('suffix')
    assert pl.suffix == 'suffix'
    assert pl.scale == 1.0

    with raises(ValueError):
        pl = TxtData(123)
    # end with

    with raises(ValueError):
        pl.scale = 0.0
    # end with

    # test loading with missing file w/ default
    res = pl.load_result('tests/unit_tests/assets/missing_file.dat', default=[nan])
    assert isinstance(res, ndarray)
    assert isnan(res[0])

    # test loading with missing file w/o default
    with raises(FileNotFoundError):
        res = pl.load_result('tests/unit_tests/assets/missing_file.dat')
    # end with

    # test loading with existing file
    pl = TxtData(suffix='energy.dat', scale=2.0)
    res = pl.load_result('tests/unit_tests/assets')
    E_ref, err_ref = 15.0, 0.1
    assert isinstance(res, ndarray)
    assert res[0] == E_ref / 2
    assert res[1] == err_ref / 2

    # Test pointing to a file directly
    res = pl.load_result('tests/unit_tests/assets/energy.dat', rescale=False)
    assert isinstance(res, ndarray)
    assert res[0] == E_ref
    assert res[1] == err_ref

# end def
