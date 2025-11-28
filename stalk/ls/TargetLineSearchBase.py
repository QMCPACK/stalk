#!/usr/bin/env python3
'''TargetLineSearch classes for the assessment and evaluation of fitting errors'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array, isscalar, nan
from scipy.interpolate import CubicSpline, PchipInterpolator

from stalk.ls.LineSearchBase import LineSearchBase
from stalk.ls.LineSearchGrid import LineSearchGrid
from stalk.ls.TlsSettings import TlsSettings


# Class for line-search with resampling and bias assessment against target
class TargetLineSearchBase(LineSearchBase):
    _target_settings: TlsSettings

    def __init__(
        self,
        offsets=None,
        values=None,
        errors=None,
        interpolate_kind='cubic',
        bias_mix=0.0,
        bias_order=1,
        target_x0=0.0,
        target_y0=0.0,
        fit_kind=None,
        fit_func=None,
        fit_args={},
        **ls_args,
        # fraction=0.025, sgn=1, fit_kind='pf3', fit_func=None, fit_args={}, N=200
    ):
        LineSearchBase.__init__(
            self,
            offsets=offsets,
            values=values,
            errors=errors,
            fit_kind=fit_kind,
            fit_func=fit_func,
            fit_args=fit_args,
            **ls_args
        )
        self._target_settings = TlsSettings(
            bias_mix=bias_mix,
            bias_order=bias_order,
            target_x0=target_x0,
            target_y0=target_y0,
            # sgn, fraction and fit_func are inherited from LsSettings
            sgn=self.settings.sgn,
            fraction=self.settings.fraction,
            fit_func=self.settings.fit_func,
            # M, N, Gs, and interp are set later upon optimization
        )
        if self.valid:
            self.reset_interpolation(interpolate_kind=interpolate_kind)
        # end if
    # end def

    @property
    def target_settings(self):
        return self._target_settings
    # end def

    @target_settings.setter
    def target_settings(self, target_settings):
        # Only update if the new settings differ from the old ones
        if not isinstance(target_settings, TlsSettings):
            raise TypeError("The settings must be derived from TlsSettings!")
        elif target_settings != self.target_settings:
            self._target_settings = target_settings
        # end if
    # end def

    @property
    def target_interp(self):
        if hasattr(self, 'target_settings'):
            return self.target_settings.interp
        # end if
    # end def

    @property
    def valid_target(self):
        return self.target_interp is not None
    # end def

    def add_point(self, point):
        super().add_point(point)
        if self.valid_target:
            self.reset_interpolation(
                interpolate_kind=self.target_settings.interp_kind
            )
        # end if
    # end def

    def reset_interpolation(
        self,
        interpolate_kind='cubic'
    ):
        if not self.valid:
            raise AssertionError("Must provide values before interpolation")
        # end if
        if interpolate_kind == 'pchip':
            self.target_settings.interp = PchipInterpolator(
                self.valid_offsets,
                self.valid_values,
                extrapolate=False
            )
        elif interpolate_kind == 'cubic':
            self.target_settings.interp = CubicSpline(
                self.valid_offsets,
                self.valid_values,
                extrapolate=False
            )
        else:
            raise ValueError("Could not recognize interpolate kind" + str(interpolate_kind))
        # end if
    # end def

    def evaluate_target(self, offsets):
        if isscalar(offsets):
            return self._evaluate_target_point(offsets)
        else:
            return array([self._evaluate_target_point(offset) for offset in offsets])
        # end if
    # end def

    def _evaluate_target_point(self, offset):
        if not self.valid_target:
            return nan
        elif not self.target_interp.extrapolate and (
                offset < self.target_interp.x.min() or offset > self.target_interp.x.max()):
            return nan
        else:
            return self.target_interp(offset)
        # end if
    # end def

    def compute_error(
        self,
        grid: LineSearchGrid,
        **ls_overrides
        # sgn=1, fraction=0.025, fit_func=None, bias_mix=None, bias_order=None,
        # N=200, Gs=None, fraction=None
    ):
        if not self.valid_target:
            return nan
        # end if
        bias = self.compute_bias(grid, **ls_overrides)
        errorbar_x, errorbar_y = self.compute_errorbar(
            grid,
            **ls_overrides
        )
        return bias + errorbar_x
    # end def

    def compute_bias(
        self,
        grid: LineSearchGrid,
        **ls_overrides
        # sgn=1, fraction=0.025, fit_func=None, bias_mix=None, bias_order=None,
        # N=200, Gs=None, fraction=None
    ):
        # Use stored settings and override where requested. The overrides are checked
        # in the copy constructor
        if not self.valid_target:
            return nan
        # end if
        settings = self.target_settings.copy(**ls_overrides)
        try:
            bias_x, bias_y = self._compute_xy_bias(grid, settings)
        except AssertionError:
            # If fitting altogether fails, assign maximum bias
            bias_x, bias_y = 1e100, 1e100
        # end try
        bias_tot = abs(bias_x) + settings.bias_mix * abs(bias_y)
        return bias_tot
    # end def

    def _compute_xy_bias(
        self,
        grid: LineSearchGrid,
        settings: TlsSettings,
    ):
        if settings.interp is None:
            return nan, nan
        # end if
        # Begin from target
        offsets0 = grid.offsets.copy()
        x0 = 0.0
        # Repeat search 'bias_order' times to simulate how bias is self-induced
        for i in range(settings.bias_order):
            offsets = self.target_settings.get_safe_offsets(offsets0 + x0)
            grid = LineSearchGrid(offsets)
            grid.values = self.evaluate_target(grid.offsets)
            res = settings.fit_func.find_minimum(
                grid,
                sgn=settings.sgn
            )
            x0 = res.x0
            y0 = res.y0
        # end for
        # Offset bias
        bias_x = x0 - settings.target.x0
        # Value bias
        bias_y = y0 - settings.target.y0
        return bias_x, bias_y
    # end def

    def compute_errorbar(
        self,
        grid: LineSearchGrid,
        **ls_overrides
        # sgn=1, fraction=0.025, fit_func=None, bias_mix=None, bias_order=None,
        # N=200, Gs=None, fraction=None
    ):
        # Use stored settings and override where requested. The overrides are checked
        # in the copy constructor
        if not self.valid_target:
            return nan, nan
        # end if
        settings = self.target_settings.copy(**ls_overrides)
        offsets = grid.valid_offsets
        values = self.evaluate_target(offsets)
        if values is not None:
            res = settings.fit_func.find_noisy_minimum(
                grid=LineSearchGrid(
                    offsets,
                    values=values,
                    errors=grid.valid_errors
                ),
                sgn=settings.sgn,
                N=None,
                Gs=settings.Gs[:, :len(offsets)],
                fraction=settings.fraction
            )
            return res.x0_err, res.y0_err
        # end if
    # end def

    def bracket_target_bias(
        self,
        bracket_fraction=0.5,
        M=7,
        max_iter=10,
        bias_order=1,  # Bias order 1 is preferred
        **ls_overrides
        # sgn=1, fraction=0.025, fit_func=None, bias_mix=None
        # N=200, Gs=None, fraction=None
    ):
        if bracket_fraction <= 0 or bracket_fraction >= 1.0:
            raise ValueError("Must be 0 < bracket_fraction < 1")
        # end if
        if not self.valid_target:
            raise AssertionError("Target data is not valid yet.")
        # end if
        settings = self.target_settings.copy(
            bias_order=bias_order,
            **ls_overrides
        )
        x0 = settings.target.x0
        y0 = settings.target.y0

        R_this = self.valid_R_max
        # Begin from initial x0 result
        bias_x = self.fit_res.x0
        for i in range(max_iter):
            offsets = self._make_offsets_R(R_this, M) + bias_x
            bias_x, bias_y = self._compute_xy_bias(
                LineSearchGrid(offsets),
                settings
            )
            R_this *= bracket_fraction
        # end for
        self.target_settings.target.x0 = x0 + bias_x
        self.target_settings.target.y0 = y0 + bias_y
    # end def

# end class
