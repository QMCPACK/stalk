#!/usr/bin/env python

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises

from stalk.params.pes_function import NotEvaluatedException
from stalk.params.pes_result import PesResult
from stalk.io.pes_loader import PesLoader


# Test PesLoader class
def test_PesLoader():

    # Test default settings
    pl = PesLoader()
    assert pl.scale == 1.0
    assert pl.suffix == 'energy.dat'

    # Try out alternative values
    pl = PesLoader(suffix='e.dat', scale=2.0, arg='test')
    assert pl.scale == 2.0
    assert pl.args.get('arg') == 'test'

    path = 'tests/unit_tests/assets'
    # Missing a file always raises exception
    with raises(NotEvaluatedException):
        pl.load(path)
    # end with
    # The file is found by the default name energy.dat
    pl.suffix = 'energy.dat'
    # But loading fails with the extra argument to np.loadtxt
    with raises(TypeError):
        pl.load(path)
    # end with
    # Clean it up, and the loading should succeed
    pl.args = {}
    res = pl.load(path)
    assert isinstance(res, PesResult)
    # See tests/unit_tests/assets/energy.dat: scaling by 1/2.0 is applied
    E_ref, err_ref = 7.5, 0.05
    assert res.value == E_ref
    assert res.error == err_ref

# end def
