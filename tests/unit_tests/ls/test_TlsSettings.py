#!/usr/bin/env python

from pytest import raises, warns
from numpy import linspace
import warnings
from scipy.interpolate import CubicSpline, PchipInterpolator

from stalk.ls import PolynomialFit
from stalk.ls.FittingFunction import FittingFunction
from stalk.ls.TlsSettings import TlsSettings
from stalk.util.util import match_to_tol

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test TlsSettings class
def test_TlsSettings():

    # Test empty init / defaults
    with raises(TypeError):
        # Must always provide fitting function/kind
        TlsSettings()
    # end with

    # Test defaults
    ls = TlsSettings(fit_kind='pf3')
    assert ls.sgn == 1
    assert ls.fraction == 0.025
    assert ls.fit_func == PolynomialFit(3)
    assert ls.bias_mix == 0.0
    assert ls.bias_order == 1
    assert ls.Gs is None
    assert ls.M == 0
    assert ls.N == 0
    assert ls.interp is None
    assert ls.interp_kind is None
    assert ls.target.x0 == 0.0
    assert ls.target.y0 == 0.0

    # Test wrong values
    with raises(ValueError):
        TlsSettings(fit_kind='pf2', bias_order=0)
    # end with
    with raises(ValueError):
        TlsSettings(fit_kind='pf2', bias_mix=-0.1)
    # end with
    with raises(ValueError):
        TlsSettings(fit_kind='pf2', Gs=[0.1, 0.2, 0.3])
    # end with
    with warns(UserWarning):
        offsets = [0.0]
        ls.get_safe_offsets(offsets)
        assert offsets is offsets
    # end with

    # Dummy fitting function handle
    def dummy_fit():
        pass
    # end def

    # test initializations of fitting functions
    N = 50
    sgn = -1
    fraction = 0.4
    bias_mix = 0.1
    bias_order = 2
    x0 = 0.1
    y0 = 0.2
    M = 4
    N = 5
    func = FittingFunction(dummy_fit, {'test': 0})
    ls1 = TlsSettings(
        fit_func=func,
        sgn=sgn,
        fraction=fraction,
        bias_mix=bias_mix,
        bias_order=bias_order,
        target_x0=x0,
        target_y0=y0,
        M=M,
        N=N,
        Gs=None
    )
    assert ls1.N == N
    assert ls1.fraction == fraction
    assert ls1.sgn == sgn
    assert isinstance(ls1.fit_func, FittingFunction)
    assert ls1.fit_func.func is dummy_fit
    assert ls1.fit_func.args['test'] == 0
    assert ls1.bias_mix == bias_mix
    assert ls1.bias_order == bias_order
    assert ls1.Gs is not None
    assert ls1.M == M
    assert ls1.N == N
    assert ls.interp is None
    assert ls1.interp_kind is None
    assert ls1.target.x0 == x0
    assert ls1.target.y0 == y0

    # Test direct supply of Gs
    ls2 = TlsSettings(
        fit_func=func,
        Gs=ls1.Gs
    )
    assert ls2.M == M
    assert ls2.N == N

    # Test copy and __eq__
    # Plain copy should always be equal
    assert ls1 == ls1.copy()
    # Override with same values should also be equal
    assert ls1 == ls1.copy(
        sgn=sgn,
        fraction=fraction,
        fit_func=func,
        bias_mix=bias_mix,
        bias_order=bias_order,
        target_x0=x0,
        target_y0=y0,
        M=M,
        N=N,
    )
    # Same number of samples results in same Gs
    assert ls1 == ls1.copy(Gs=None, M=M, N=N)
    # Supplying Gs takes precedence
    assert ls1 == ls1.copy(Gs=ls1.Gs)
    # Different M or N are not equal
    assert ls1 != ls1.copy(Gs=None, M=(M + 1), N=N)
    assert ls1 != ls1.copy(Gs=None, M=M, N=(N + 1))
    # Override with different values is not equal
    assert ls1 != ls1.copy(fraction=fraction + 0.01)  # LsSettings unequal
    assert ls1 != ls1.copy(bias_mix=2 * bias_mix)
    assert ls1 != ls1.copy(bias_order=2 * bias_order)
    assert ls1 != ls1.copy(bias_order=2 * bias_order)
    # Cannot override target
    assert ls1 == ls1.copy(target_x0=1.0)

    # interpolation tests
    offsets = linspace(-2, 2, 5)
    values = 0.2 * offsets**2
    interp_spl_ext = CubicSpline(
        offsets,
        values,
        extrapolate=True
    )
    interp_spl_noext = CubicSpline(
        offsets,
        values,
        extrapolate=False
    )
    interp_pch_ext = PchipInterpolator(
        offsets,
        values,
        extrapolate=True
    )
    ls3 = TlsSettings(
        fit_func=func,
        sgn=sgn,
        fraction=fraction,
        bias_mix=bias_mix,
        bias_order=bias_order,
        target_x0=x0,
        target_y0=y0,
    )
    ls3.interp = interp_spl_ext

    assert ls3.interp_kind == 'cubic'
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        offsets_out = ls3.get_safe_offsets(offsets.copy())
        assert match_to_tol(offsets_out, offsets)
    # end with
    with warns(UserWarning, match='Extrapolating'):
        dp = 0.02
        offsets_out = ls3.get_safe_offsets(offsets + dp)
        assert match_to_tol(offsets_out, offsets + dp)
    # end with
    with warns(UserWarning, match='Extrapolating'):
        dp = -0.02
        offsets_out = ls3.get_safe_offsets(offsets + dp)
        assert match_to_tol(offsets_out, offsets + dp)
    # end with
    ls3.interp = interp_spl_noext
    with warns(UserWarning, match='Reset offset='):
        dp = 0.02
        offsets_out = ls3.get_safe_offsets(offsets + dp)
        assert match_to_tol(offsets_out[:-1], offsets[:-1] + dp)
        assert match_to_tol(offsets_out[-1], offsets[-1])
    # end with
    with warns(UserWarning, match='Reset offset'):
        dp = -0.02
        offsets_out = ls3.get_safe_offsets(offsets + dp)
        assert match_to_tol(offsets_out[1:], offsets[1:] + dp)
        assert match_to_tol(offsets_out[0], offsets[0])
    # end with

    ls3.interp = interp_pch_ext
    assert ls3.interp_kind == 'pchip'
    ls3.interp = None
    assert ls3.interp_kind is None

# end def
