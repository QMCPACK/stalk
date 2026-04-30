#!/usr/bin/env python

from numpy import isnan, linspace, random
from scipy.interpolate import PchipInterpolator, CubicSpline
from pytest import raises

from stalk.ls.linesearch_grid import LineSearchGrid
from stalk.ls.polynomial_fit import PolynomialFit
from stalk.params.linesearch_point import LineSearchPoint
from stalk.util import match_to_tol
from stalk.ls import TargetLineSearchBase

from ..assets.fitting_pf2 import generate_exact_pf2, generate_exact_pf3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# test TargetLineSearchBase class
def test_TargetLineSearchBase():

    with raises(TypeError):
        # Cannot init without fit_func/kind
        TargetLineSearchBase()
    # end with

    # Test init of defaults
    tls = TargetLineSearchBase(fit_kind='pf3')
    assert tls.settings.N == 200
    assert tls.settings.fit_func == PolynomialFit(3)
    assert tls.settings.sgn == 1
    assert tls.settings.fraction == 0.025
    assert tls.target_settings.Gs is None
    assert tls.target_settings.fit_func == PolynomialFit(3)
    assert tls.target_settings.sgn == 1
    assert tls.target_settings.fraction == 0.025
    assert tls.target_settings.M == 0
    assert tls.target_settings.N == 0
    assert tls.target_settings.bias_mix == 0.0
    assert tls.target_settings.bias_order == 1
    assert tls.target_settings.target.x0 == 0.0
    assert tls.target_settings.target.y0 == 0.0
    assert tls.target_settings.target.x0_err == 0.0
    assert tls.target_settings.target.y0_err == 0.0
    assert tls.target_settings.interp is None
    assert tls.target_settings.interp_kind is None
    assert tls.target_interp is None
    assert not tls.valid_target

    # Test init of non-defaults
    N = 100
    sgn = -1
    fraction = 0.2
    bias_mix = 0.1
    bias_order = 2
    x0 = 0.1
    y0 = 0.2
    tls = TargetLineSearchBase(
        fit_kind='pf3',
        N=N,
        sgn=sgn,
        fraction=fraction,
        bias_mix=bias_mix,
        bias_order=bias_order,
        target_x0=x0,
        target_y0=y0
    )
    assert tls.settings.N == N
    assert tls.settings.fit_func == PolynomialFit(3)
    assert tls.settings.sgn == sgn
    assert tls.settings.fraction == fraction
    assert tls.target_settings.fit_func == PolynomialFit(3)
    assert tls.target_settings.sgn == sgn
    assert tls.target_settings.fraction == fraction
    assert tls.target_settings.bias_mix == bias_mix
    assert tls.target_settings.bias_order == bias_order
    assert tls.target_settings.target.x0 == x0
    assert tls.target_settings.target.y0 == y0

    # Test init with only offsets input
    offsets = linspace(-0.5, 0.5, 21)
    tls = TargetLineSearchBase(
        fit_kind='pf3',
        offsets=offsets
    )
    assert not tls.valid
    # Cannot reset interpolation while empty
    with raises(AssertionError):
        tls.reset_interpolation()
    # end with
    # Target evaluation results in None
    assert isnan(tls.evaluate_target([0.0]))
    assert isnan(tls.evaluate_target(0.0))
    # Bias assessment results in nan
    assert isnan(tls.compute_bias([LineSearchGrid()]))
    # Errorbar assessment results in nan, nan
    assert isnan(tls.compute_errorbar(LineSearchGrid())[0])
    assert isnan(tls.compute_errorbar(LineSearchGrid())[1])
    # Total error assessment results in nan
    assert isnan(tls.compute_error(LineSearchGrid()))
    # Cannot extrapolate error
    with raises(AssertionError):
        tls.bracket_target_bias()
    # end with

    # Test nominal init with reference potential data
    grid, ref = generate_exact_pf2(1.23, 2.34, N=21, error=0.1)
    bias_mix = 0.1
    bias_order = 1
    interpolate_kind = 'pchip'
    fraction = 0.05
    fit_kind = 'pf2'
    N = 20
    tls = TargetLineSearchBase(
        offsets=grid.offsets,
        values=grid.values,
        bias_mix=bias_mix,
        bias_order=bias_order,
        interpolate_kind=interpolate_kind,
        fraction=fraction,
        fit_kind=fit_kind,
        N=N
    )
    assert tls.target_settings.bias_mix == bias_mix
    assert tls.valid_target
    # Test interpolant
    assert isinstance(tls.target_interp, PchipInterpolator)
    # Test Fitting function
    assert tls.fit_res.fraction == fraction
    assert tls.fit_res.analyzed

    # Test evaluation
    assert match_to_tol(tls.evaluate_target(grid.offsets), grid.values)
    assert tls.evaluate_target(grid.offsets[4]) == grid.values[4]
    assert isnan(tls.evaluate_target(grid.offsets[0] - 1e-6))
    assert isnan(tls.evaluate_target(grid.offsets[-1] + 1e-6))
    # Test assessment of bias
    bias = tls.compute_bias(grid)
    assert match_to_tol(bias, ref.x0 + bias_mix * ref.y0)
    # Test assessment of errorbars
    Gs = random.randn(30, len(tls))
    errorbar_x, errorbar_y = tls.compute_errorbar(grid, Gs=Gs)
    assert errorbar_x > 0.0
    assert errorbar_y > 0.0
    # Test assessment of total error
    error = tls.compute_error(grid, Gs=Gs)
    assert match_to_tol(error, bias + errorbar_x)

    # Test reset interpolation
    tls.reset_interpolation(interpolate_kind='cubic')
    assert isinstance(tls.target_settings.interp, CubicSpline)
    with raises(ValueError):
        tls.reset_interpolation('error')
    # end with
    assert match_to_tol(tls.evaluate_target(grid.offsets), grid.values)
    assert tls.evaluate_target(grid.offsets[4]) == grid.values[4]
    assert isnan(tls.evaluate_target(grid.offsets[0] - 1e-6))
    assert isnan(tls.evaluate_target(grid.offsets[-1] + 1e-6))

    # Test bracket target bias
    biased_grid, ref = generate_exact_pf3(1.23, 262.3, N=11)
    tls = TargetLineSearchBase(
        offsets=biased_grid.offsets,
        values=biased_grid.values,
        fit_kind='pf2'
    )
    with raises(ValueError):
        tls.bracket_target_bias(bracket_fraction=0.0)
    # end with
    with raises(ValueError):
        tls.bracket_target_bias(bracket_fraction=1.0)
    # end with
    # At first the target biases are nonzero
    assert tls.target_settings.target.x0 != ref.x0
    assert tls.target_settings.target.y0 != ref.y0
    tls.bracket_target_bias(fit_kind='pf2')
    # After bracketing, the biases should be close to zero
    assert match_to_tol(tls.target_settings.target.x0, ref.x0, tol=1e-4)
    assert match_to_tol(tls.target_settings.target.y0, ref.y0, tol=1e-4)

    # Test addition of a point
    assert len(tls) == 11
    assert len(tls.target_interp.x) == 11
    offset = 3.45
    val = 7.0
    old_val = tls.evaluate_target(offset)
    # The target is updated if both offset and value are provided
    tls.add_point(LineSearchPoint(offset, val))
    new_val = tls.evaluate_target(offset)
    assert old_val != new_val
    assert new_val == val
    assert len(tls) == 12
    assert len(tls.target_interp.x) == 12

# end def
