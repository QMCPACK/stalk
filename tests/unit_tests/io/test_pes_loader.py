#!/usr/bin/env python

from pytest import warns
from numpy import isnan

from stalk.params.parameter_set import ParameterSet
from stalk.params.pes_result import PesResult
from stalk.io.pes_loader import PesLoader

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test PesLoader class
def test_PesLoader():

    # Test return by with a custom loader function
    def test_loader(structure: ParameterSet, add=0):
        return PesResult(float(len(structure.file_path)) + add)
    # end def

    # Manually replace load method for testing
    pl = PesLoader(suffix='energy.dat')

    path = 'tests/unit_tests/assets'
    res = pl.load(path)
    assert isinstance(res, PesResult)
    # See tests/unit_tests/assets/energy.dat
    E_ref, err_ref = 15.0, 0.1
    assert res.value == E_ref
    assert res.error == err_ref

    # Test overriding arg
    scale = 1.4
    res2 = pl.load(path, scale=scale)
    assert isinstance(res2, PesResult)
    assert res2.value == E_ref / scale
    assert res2.error == err_ref / scale

    # Test loading add_sigma > 0
    sigma = 1.23
    res_sigma = pl.load(path, sigma=sigma)
    assert isinstance(res_sigma, PesResult)
    assert res_sigma.value != E_ref
    assert res_sigma.error == (err_ref**2 + sigma**2)**0.5

    # Test not finding the file
    with warns(UserWarning):
        res_missing = pl.load('missing')
        assert isnan(res_missing.value)
    # end with

# end def
