#!/usr/bin/env python

from pytest import raises

from stalk.ls.polynomial_fit import PolynomialFit
from stalk.ls.fitting_function import FittingFunction
from stalk.ls.ls_settings import LsSettings

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test LsSettings class
def test_LsSettings():

    # Test empty init / defaults
    with raises(TypeError):
        # Must always provide fitting function/kind
        LsSettings()
    # end with
    ls = LsSettings(fit_kind='pf3')
    assert ls.N == 200
    assert ls.sgn == 1
    assert ls.fraction == 0.025
    assert ls.fit_func == PolynomialFit(3)

    # Test wrong values
    with raises(ValueError):
        LsSettings(N=0)
    # end with
    with raises(ValueError):
        LsSettings(N=None)
    # end with
    with raises(ValueError):
        LsSettings(sgn=0.9)
    # end with
    with raises(ValueError):
        LsSettings(sgn=None)
    # end with
    with raises(ValueError):
        LsSettings(fraction=0.0)
    # end with
    with raises(ValueError):
        LsSettings(fraction=0.5)
    # end with
    with raises(ValueError):
        LsSettings(fraction=None)
    # end with

    # Dummy fitting function handle
    def dummy_fit():
        pass
    # end def

    # test initializations of fitting functions
    N = 50
    sgn = -1
    fraction = 0.4
    func = FittingFunction(dummy_fit, {'test': 0})
    ls1 = LsSettings(
        fit_func=func,
        N=N,
        sgn=sgn,
        fraction=fraction
    )
    assert ls1.N == N
    assert ls1.fraction == fraction
    assert ls1.sgn == sgn
    assert isinstance(ls1.fit_func, FittingFunction)
    assert ls1.fit_func.func is dummy_fit
    assert ls1.fit_func.args['test'] == 0
    # This should initialize to pf4 and ignore fit_args
    ls2 = LsSettings(
        fit_kind='pf4',
        fit_args={'pfn': 5}
    )
    assert ls2.fit_func == PolynomialFit(4)
    with raises(TypeError):
        # Cannot set unrecognized format
        LsSettings(fit_kind='error')
    # end with

    # Test copy and __eq__
    # Plain copy shoud always be equal
    assert ls1 == ls1.copy()
    # Override with same values should also be equal
    assert ls1 == ls1.copy(N=N, sgn=sgn, fraction=fraction, fit_func=func)
    # Override with different values is not equal
    assert ls1 != ls1.copy(N=2 * N)
    assert ls1 != ls1.copy(sgn=-sgn)
    assert ls1 != ls1.copy(fraction=0.5 * fraction)
    assert ls1 != ls1.copy(fit_kind='pf2')

# end def
