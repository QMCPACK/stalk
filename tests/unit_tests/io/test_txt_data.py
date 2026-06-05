#!/usr/bin/env python

from pytest import raises
from numpy import nan, ndarray, isnan

from stalk.io.txt_data import TxtData
from stalk.util.util import match_to_tol

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test PesLoader class
def test_TxtData(tmp_path):

    # Test default values
    pl = TxtData()
    assert pl.suffix == 'data.dat'
    assert pl.scale == 1.0
    # Test custom values
    pl = TxtData('suffix.dat', 2.0)
    assert pl.suffix == 'suffix.dat'
    assert pl.scale == 2.0

    with raises(ValueError):
        pl = TxtData(123)
    # end with

    with raises(ValueError):
        pl.scale = 0.0
    # end with

    # test treating a missing file
    pl = TxtData()
    path = tmp_path
    # Looks for 'data.dat' in tmp_path
    assert not pl.exists(path)
    assert pl.get_filename(path) == tmp_path / 'data.dat'
    assert pl.get_filename(path / 'missing_dir') == tmp_path / 'missing_dir' / 'data.dat'
    # Loading with a default value
    res = pl.load_result(tmp_path / 'missing', default=nan)
    assert isnan(res)

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

    # test saving to file
    pl = TxtData(suffix='normal.dat', scale=2.0)
    data = [1.0, 0.1]
    pl.save_result(tmp_path, data=data)
    data_rescale = pl.load_result(tmp_path, rescale=True)
    assert match_to_tol(data_rescale * 2, data)
    # Testing overwrite and treating None data
    pl.save_result(tmp_path, data=None, overwrite=False)
    # The old data is recovered without overwrite
    data_load = pl.load_result(tmp_path, rescale=False)
    assert match_to_tol(data_load, data)
    pl.save_result(tmp_path, data=None, overwrite=True)
    data_load = pl.load_result(tmp_path)
    assert isnan(data_load).all()

# end def
