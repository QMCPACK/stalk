#!/usr/bin/env python3
'''TargetLineSearch classes for the assessment and evaluation of fitting errors'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from matplotlib import pyplot as plt
from numpy import array, isscalar, linspace, nan
from numpy import ndarray

from stalk.ls.ErrorSurface import ErrorSurface
from stalk.ls.LineSearchGrid import LineSearchGrid
from stalk.ls.TlsSettings import TlsSettings
from stalk.util.util import FF
from .LineSearch import LineSearch
from .TargetLineSearchBase import TargetLineSearchBase


# Class for error scan line-search
class TargetLineSearch(TargetLineSearchBase, LineSearch):
    _Gs: ndarray = None  # N x M set of correlated random fluctuations for the grid
    _epsilon = None  # optimized target error
    _W_opt = None  # W to meet epsilon
    _sigma_opt = None  # sigma to meet epsilon
    _error_surface: ErrorSurface = None  # Error surface

    def __init__(
        self,
        # kwargs related to LineSearch
        structure=None,
        hessian=None,
        d=None,
        path='',
        interactive=False,
        # sigma=0.0
        offsets=None,
        M=7,
        W=None,
        R=None,
        pes=None,
        bracket=True,
        # kwargs related to TargetLineSearchBase
        bias_order=1,
        bias_mix=0.0,
        interpolate_kind='cubic',
        fit_kind=None,
        fit_func=None,
        fit_args={},
        # Line-search kwargs
        **ls_args
        # offsets=None, values=None, errors=None, fraction=0.025, sgn=1, N=200, Gs=None
    ):
        # Set bias_mix and target_fit in TargetLineSearchBase
        TargetLineSearchBase.__init__(
            self,
            bias_mix=bias_mix,
            bias_order=bias_order,
            fit_kind=fit_kind,
            fit_func=fit_func,
            fit_args=fit_args,
        )
        # The necessary init is done in LineSearch class
        LineSearch.__init__(
            self,
            structure=structure,
            hessian=hessian,
            d=d,
            M=M,
            W=W,
            R=R,
            offsets=offsets,
            pes=pes,
            interactive=interactive,
            path=path,
            **ls_args,
        )
        # Finally, attempt to reset interpolation
        if self.valid:
            self.reset_interpolation(interpolate_kind=interpolate_kind)
            if bracket:
                self.bracket_target_bias()
            # end if
        # end if
    # end def

    @property
    def target_settings(self):
        return self._target_settings
    # end def

    @property
    def setup(self):
        return self.target_settings is not None and self.M > 0
    # end def

    @target_settings.setter
    def target_settings(self, target_settings):
        if not isinstance(target_settings, TlsSettings):
            raise TypeError("The settings must be derived from TlsSettings!")
        else:
            self._target_settings = target_settings
            # Clear the error surface when new settings are introduced
            self._error_surface = None
        # end if
    # end def

    @property
    def sigma_opt(self):
        return self._sigma_opt
    # end def

    @sigma_opt.setter
    def sigma_opt(self, sigma_opt):
        if isscalar(sigma_opt) and sigma_opt > 0.0 or sigma_opt is None:
            self._sigma_opt = sigma_opt
        else:
            raise ValueError("sigma_opt must be > 0.0 but " + str(sigma_opt) + " was given.")
        # end if
    # end def

    @property
    def W_opt(self):
        return self._W_opt
    # end def

    @W_opt.setter
    def W_opt(self, W_opt):
        if isscalar(W_opt) and W_opt > 0.0 or W_opt is None:
            self._W_opt = W_opt
        else:
            raise ValueError("W_opt must be > 0.0 but " + str(W_opt) + " was given.")
        # end if
    # end def

    @property
    def epsilon(self):
        return self._epsilon
    # end def

    @epsilon.setter
    def epsilon(self, epsilon):
        if isscalar(epsilon) and epsilon > 0.0 or epsilon is None:
            self._epsilon = epsilon
        else:
            raise ValueError("Epsilon must be > 0.0 but " + str(epsilon) + " was given.")
        # end if
    # end def

    @property
    def grid_opt(self):
        if not self.optimized:
            return None
        # end if
        offsets = self.figure_out_offsets(M=self.M, W=self.W_opt)
        errors = self.target_settings.M * [self.sigma_opt]
        grid = LineSearchGrid(
            offsets=offsets,
            errors=errors
        )
        return grid
    # end def

    @property
    def R_opt(self):
        return self._W_to_R(self.W_opt)
    # end def

    @property
    def resampled(self):
        return isinstance(self._error_surface, ErrorSurface)
    # end def

    @property
    def optimized(self):
        return self.W_opt is not None and self.sigma_opt is not None
    # end def

    @property
    def Gs(self):
        return self.target_settings.Gs
    # end def

    @Gs.setter
    def Gs(self, Gs):
        self.target_settings.Gs = Gs
        self._error_surface = None
    # end def

    @property
    def M(self):
        return self.target_settings.M
    # end def

    @property
    def error_surface(self):
        return self._error_surface
    # end def

    def compute_bias_of(
        self,
        M=None,
        R=None,
        W=None,
        num_W=10,  # number of W points if no grids are provided
        **ls_args  # sgn, fraction, fit_func, N, Gs, bias_order, bias_mix
    ):
        M = M if M is not None else len(self)

        Rs = []
        Ws = []
        if R is None:
            if W is None:
                # By default, compute bias for a range of W values
                for W in linspace(0.0, self.valid_W_max, num_W):
                    Ws.append(W)
                    Rs.append(self._W_to_R(W))
                # end for
            elif isscalar(W):
                Ws.append(W)
                Rs.append(self._W_to_R(W))
            else:
                Ws = W
                Rs = [self._W_to_R(W) for W in Ws]
            # end if
        else:
            if isscalar(R):
                Ws.append(self._R_to_W(R))
                Rs.append(R)
            else:
                Rs = R
                Ws = [self._R_to_W(R) for R in Rs]
            # end if
        # end if

        biases = []
        for W in Ws:
            offsets = self.figure_out_adjusted_offsets(M=M, W=W)
            values = self.evaluate_target(offsets)
            bias = self.compute_bias(
                grid=LineSearchGrid(offsets=offsets, values=values),
                **ls_args
            )
            biases.append(bias)
        # end for
        return array(Ws), array(Rs), array(biases)
    # end def

    def figure_out_adjusted_offsets(
        self,
        **grid_args  # M, R, W, offsets
    ):
        adjusted_offsets = self.figure_out_offsets(**grid_args) + self.target_settings.target.x0
        return self.target_settings.get_safe_offsets(adjusted_offsets)
    # end def

    def optimize(
        self,
        epsilon,
        max_rounds=10,  # maximum number of rounds
        skip_setup=False,  # If confident
        **kwargs
        # W_resolution=0.1, S_resolution=0.1, verbosity=1
        # fit_kind=None, fit_func=None, fit_args={}, fraction=None,
        # generate_args: W_num, W_max, sigma_num, sigma_max, noise_frac, M, N, Gs
        # bias_mix, bias_order
    ):
        """Optimize W and sigma to a given target error epsilon > 0."""
        if not self.valid_target:
            raise AssertionError("Must have valid target data before optimization.")
        # end if
        if not isscalar(epsilon) or epsilon <= 0.0:
            raise ValueError("Must provide epsilon > 0.")
        # end if
        if max_rounds <= 0:
            raise ValueError('Must provide max_rounds > 0')
        # end if

        # Skip if already optimized and allowed to skip setup (overlooks all overrides)
        if not self.optimized or not skip_setup:
            self.setup_optimization(**kwargs)
        # end if
        # Find W and sigma that maximize sigma
        for round in range(max_rounds):
            # Probe accuracy requests from ErrorSurface
            W_vals, sigma_vals = self.error_surface.request_points(epsilon)
            if len(W_vals) == 0 and len(sigma_vals) == 0:
                # Break when no more accuracy can be requested
                break
            else:
                # Insert x-cols and y-rows
                for W_val in W_vals:
                    self.insert_W_data(W_val)
                # end for
                for sigma_val in sigma_vals:
                    self.insert_sigma_data(sigma_val)
                # end for
            # end if
        # end while

        W_opt, sigma_opt = self.error_surface.evaluate_target(epsilon)
        if W_opt == 0.0 or sigma_opt == 0.0:
            msg = f"Warning! Optimization to epsilon={epsilon} resulted in "
            msg += f"W_opt={W_opt} and sigma_opt={sigma_opt} which is numerically "
            msg += "unfeasible. Check the error surface or epsilon!"
            print(msg)
            self.W_opt = None
            self.sigma_opt = None
            self.epsilon = None
        else:
            self.W_opt = W_opt
            self.sigma_opt = sigma_opt
            self.epsilon = epsilon
        # end if
    # end def

    def setup_optimization(
        self,
        W_num=3,
        W_max=None,
        sigma_num=3,
        sigma_max=None,
        noise_frac=0.05,
        W_resolution=0.1,
        S_resolution=0.1,
        verbosity=1,
        **ls_overrides
        # fit_kind=None, fit_func=None, fit_args={}, fraction=0.025, Gs=None,
        # M=None, N=None, bias_mix=0.0, bias_order=1
    ):
        if not self.valid_target:
            raise AssertionError("Must have valid target data before setup.")
        # end if
        # Create new settings. If they do not match the previous ones (checked in setter):
        # -> clear E_mat and regenerate Gs and regenerate the error surface
        target_settings = self.target_settings.copy(**ls_overrides)
        if target_settings != self.target_settings:
            self.target_settings = target_settings
        # end if
        # Make sure error surface is generated
        if not self.resampled:
            self.generate_error_surface(
                W_num=W_num,
                W_max=W_max,
                sigma_num=sigma_num,
                sigma_max=sigma_max,
                noise_frac=noise_frac,
                W_resolution=W_resolution,
                S_resolution=S_resolution,
                verbosity=verbosity
            )
        # end if
    # end def

    # X: grid of W values; Y: grid of sigma values; E: grid of total errors
    #   if Gs is not provided, use M and N
    def generate_error_surface(
        self,
        W_num=3,
        W_max=None,
        sigma_num=3,
        sigma_max=None,
        noise_frac=0.05,
        W_resolution=0.1,
        S_resolution=0.1,
        verbosity=1
    ):
        if not self.valid_target:
            raise AssertionError("Must have valid target data before generating error.")
        elif not self.setup:
            raise AssertionError("Must setup target line-search with M > 2")
        # end if
        if W_resolution >= 0.5 or W_resolution <= 0.0:
            raise ValueError('W resolution must be 0.0 < W_resolution < 0.5')
        # end if
        if S_resolution >= 0.5 or S_resolution <= 0.0:
            raise ValueError('S resolution must be 0.0 < S_resolution < 0.5')
        # end if

        W_max = W_max if W_max is not None else self.W_max
        sigma_max = sigma_max if sigma_max is not None else W_max * noise_frac

        if W_max <= 0.0:
            raise ValueError('Must provide W_max > 0')
        # end if
        if sigma_max <= 0.0:
            raise ValueError('Must provide sigma_max > 0')
        # end if

        # Initial W and sigma grids
        self._error_surface = ErrorSurface(
            X_res=W_resolution,
            Y_res=S_resolution,
            verbosity=verbosity
        )
        # Start from adding the first row: sigma=0 -> plain bias
        Ws = linspace(0.0, W_max, W_num)
        for W in Ws:
            if W > 0.0:
                self.insert_W_data(W)
            # end if
        # end for

        # Then, append the noisy rows
        sigmas = linspace(0.0, sigma_max, sigma_num)
        for sigma in sigmas[1:]:
            self.insert_sigma_data(sigma)
        # end for
    # end def

    def insert_sigma_data(self, sigma):
        if not (self.resampled and isscalar(sigma) and sigma > 0):
            raise AssertionError('Must have resampled data and scalar sigma > 0')
        # end if
        E_row = [self._compute_target_error(W, sigma) for W in self.error_surface.Xs]
        self.error_surface.insert_row(sigma, E_row)
    # end def

    def insert_W_data(self, W):
        if not (self.resampled and isscalar(W) and W > 0 and W <= self.W_max):
            raise AssertionError(f'Must have resampled data and scalar 0 < W <= W_max, W_max={self.W_max}')
        # end if
        E_col = [self._compute_target_error(W, sigma) for sigma in self.error_surface.Ys]
        self.error_surface.insert_col(W, E_col)
    # end def

    # Compute fitting bias and error using consistent parameters
    def _compute_target_error(self, W, sigma):
        offsets = self.figure_out_adjusted_offsets(W=W, M=self.M)
        values = self.evaluate_target(offsets)
        errors = self.M * [sigma]
        grid = LineSearchGrid(offsets=offsets, values=values, errors=errors)
        error = self.compute_error(grid)
        return error
    # end def

    def to_settings(self):
        result = {
            'd': self.d,
            'fraction': self.target_settings.fraction,
            'fit_kind': self.target_settings.fit_func.kind,
            'sgn': self.target_settings.sgn,
            'M': self.target_settings.M,
            'W': self.W_opt,
            'sigma': self.sigma_opt,
            'hessian': self.hessian,
            'structure': self.structure,
        }
        return result
    # end def

    def statistical_cost(self):
        """Return statistical cost based on sigma and M"""
        if not self.optimized:
            return -1.0
        # end if
        return self.M * self.sigma_opt**-2
    # end def

    def plot(
        self,
        ax=None,
        **kwargs
    ):
        if ax is None:
            f, ax = self._create_plot()
        # end if
        target = None
        if self.optimized:
            offsets = self.figure_out_adjusted_offsets(W=self.W_opt, M=self.M)
            values = self.evaluate_target(offsets)
            grid = LineSearchGrid(offsets, values)
            target = self.target_settings.fit_func.find_minimum(grid)
            errors = self.M * [self.sigma_opt]
            offsets -= self.target_settings.target.x0
            ax.errorbar(
                offsets,
                values,
                errors,
                linestyle='None',
                marker='*',
                color='r',
                label='Optimized grid'
            )
        # end if
        LineSearch.plot(self, ax=ax, target=target)
    # end def

    def plot_error_surface(
        self,
        ax=None
    ):
        if not self.optimized:
            print('Must optimize before plotting error surface')
            return
        # end if
        if ax is None:
            f, ax = self._create_plot(xlabel='Grid extent W', ylabel='Input noise')
            ax.set_title(f'Error surface: {repr(self)}, epsilon={FF.format(self.epsilon)}')
        # end if
        T = self.error_surface.T_mat
        X = self.error_surface.X_mat
        Y = self.error_surface.Y_mat
        Z = self.error_surface.E_mat
        Z[~T] = nan
        ax.contourf(X, Y, Z)
        ax.contour(X, Y, Z, [self.epsilon], colors='k')
        ax.plot(X.flatten(), Y.flatten(), 'k.', alpha=0.3)
        ax.plot(self.W_opt, self.sigma_opt, 'ko')
        plt.tight_layout()
    # end def

    def _reset_resampling(self):
        self.error_surface = None
    # end def

    def __str__(self):
        string = LineSearch.__str__(self)
        if self.optimized:
            string += '\n  W_opt: ' + FF.format(self.W_opt)
            string += '\n  sigma_opt: ' + FF.format(self.sigma_opt)
            string += '\n  epsilon: ' + FF.format(self.epsilon)
        return string
    # end def

# end class
