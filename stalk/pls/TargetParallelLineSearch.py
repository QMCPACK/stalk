#!/usr/bin/env python3
'''TargetParallelLinesearch class for assessment of parallel mixed errors

This is the surrogate model used to inform and optimize a parallel line-search.
'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import argmin, array, isscalar, mean, linspace

from stalk.ls.LineSearchGrid import LineSearchGrid
from stalk.ls.LsSettings import LsSettings
from stalk.util import get_fraction_error
from stalk.ls import TargetLineSearch
from stalk.pls import ParallelLineSearch
from stalk.util.util import FF, FFS, FI, FIS, FP, FPS


class TargetParallelLineSearch(ParallelLineSearch):
    ls_type = TargetLineSearch
    _ls_list: list[TargetLineSearch] = []
    _temperature = None
    epsilon_p = None
    error_p = None
    error_d = None

    # Return a list of enabled line-searches
    @property
    def ls_list(self):
        return self._ls_list
    # end def

    @property
    def epsilon_d(self):
        epsilon_d = [tls.epsilon for tls in self.ls_list]
        if all(eps is None for eps in epsilon_d):
            return None
        else:
            return epsilon_d
        # end if
    # end def

    @property
    def optimized(self):
        return len(self) > 0 and all([tls.optimized for tls in self.ls_list])
    # end def

    @property
    def x_targets(self):
        return [tls.target_settings.target.x0 for tls in self.ls_list]
    # end def

    @x_targets.setter
    def x_targets(self, targets):
        for x0, ls in zip(targets, self.ls_list):
            ls.target_settings.target.x0 = x0
        # end for
    # end def

    @property
    def temperature(self):
        return self._temperature
    # end def

    @temperature.setter
    def temperature(self, temperature):
        if temperature is None:
            self._temperature = None
        elif isscalar(temperature) and temperature > 0:
            self._temperature = temperature
        else:
            raise ValueError("Temperature must be positive.")
        # end if
    # end def

    @property
    def statistical_cost(self):
        return sum([tls.statistical_cost() for tls in self.ls_list])
    # end def

    @property
    def M(self):
        return [tls.M for tls in self.ls_list]
    # end def

    @property
    def Gs(self):
        return [tls.Gs for tls in self.ls_list]
    # end def

    @property
    def W_opt(self):
        if self.optimized:
            return [tls.W_opt for tls in self.ls_list]
        else:
            return None
        # end if
    # end def

    @property
    def sigma_opt(self):
        if self.optimized:
            return [tls.sigma_opt for tls in self.ls_list]
        else:
            return None
        # end if
    # end def

    def __init__(
        self,
        path='surrogate',
        structure=None,
        hessian=None,
        targets=None,
        interpolate_kind='cubic',
        **pls_args
        # windows=None, window_frac=0.25, noises=None, add_sigma=False, no_eval=False
        # pes=None, pes_func=None, pes_args={}, loader=None, interactive=False,
        # M=7, fit_kind='pf3', fit_func=None, fit_args={}, N=200, Gs=None, fraction=0.025
    ):
        ParallelLineSearch.__init__(
            self,
            path=path,
            structure=structure,
            hessian=hessian,
            **pls_args
        )
        if targets is not None:
            self.x_targets = targets
        # end if
        for tls in self.ls_list:
            if tls.valid:
                tls._search_and_store()
                tls.reset_interpolation(interpolate_kind=interpolate_kind)
            # end if
        # end for
    # end def

    def optimize(
        self,
        reoptimize=True,
        windows=None,
        noises=None,
        epsilon_p=None,
        epsilon_d=None,
        temperature=None,
        noise_frac=0.1,
        resolution=0.01,
        starting_mix=0.5,
        write=None,
        overwrite=False,
        **ls_args
        # M=7, fit_kind=None, fit_func=None, fit_args={}, Gs=None, fraction=0.025
        # bias_mix=0.0, bias_order=1, noise_frac=0.05,
        # W_resolution=0.05, S_resolution=0.05, verbosity=1, max_rounds=10
        # W_num=3, W_max=None, sigma_num=3, sigma_max=None
    ):
        if self.optimized and not reoptimize:
            warnings.warn('Already optimized, use reoptimize = True to reoptimize.')
            return
        # end if
        if not self.evaluated:
            raise AssertionError('Cannot optimize before data has been evaluated.')
        # end if
        if windows is not None and noises is not None:
            self.optimize_windows_noises(
                windows,
                noises,
                **ls_args
            )
        elif temperature is not None:
            self.optimize_temperature(
                temperature,
                noise_frac=noise_frac,
                **ls_args
            )
        elif epsilon_d is not None:
            self.optimize_epsilon_d(
                epsilon_d,
                noise_frac=noise_frac,
                **ls_args
            )
        elif epsilon_p is not None:
            self.optimize_epsilon_p(
                epsilon_p,
                noise_frac=noise_frac,
                resolution=resolution,
                starting_mix=starting_mix,
                **ls_args
            )
        else:
            raise AssertionError('Optimizer constraint not identified')
        # end if
        # Finalize and store the result, write to disk if requested
        self._finalize_optimization(write=write, overwrite=overwrite)
    # end def

    def optimize_windows_noises(
        self,
        windows,
        noises,
        Gs=None,
        **ls_args
        # fit_kind=None, fit_func=None, fit_args={},
        # N=500, M=7, fraction=0.025,
        # bias_mix=0.0, bias_order=1, noise_frac=0.05,
        # W_resolution=0.05, S_resolution=0.05, max_rounds=10
        # W_num=3, W_max=None, sigma_num=3, sigma_max=None
    ):
        # If provided, distribute Gs per line-search; if not, provide None
        print("Optimizing to windows + noises:")
        fmt = ' tls=' + FI + ' (window=' + FF + ', noise=' + FF + ')'
        if Gs is None:
            Gs = len(windows) * [None]
        # end if
        for window, sigma, tls, Gs_this in zip(windows, noises, self.ls_list, Gs):
            # No optimization necessary if the windows, noises are readily provided but
            # generating error surface to store all required settings
            print(fmt.format(tls.d, window, sigma))
            tls.setup_optimization(
                Gs=Gs_this,
                **ls_args
            )
            tls.W_opt = window
            tls.sigma_opt = sigma
        # end for
    # end def

    def optimize_epsilon_d(
        self,
        epsilon_d,
        Gs=None,
        skip_setup=False,
        **ls_args
        # N=500, M=7, fraction=0.025, fit_kind=None, fit_func=None,
        # fit_args={}, bias_mix=0.0, bias_order=1, noise_frac=0.05,
        # W_resolution=0.05, S_resolution=0.05, max_rounds=10
        # W_num=3, W_max=None, sigma_num=3, sigma_max=None
    ):
        # If provided, distribute Gs per line-search; if not, provide None
        print("Optimizing to epsilon_d")
        fmt = ' tls=' + FI + ' (epsilon=' + FF + ')'
        if Gs is None:
            Gs = len(epsilon_d) * [None]
        # end if
        for epsilon, tls, Gs_this in zip(epsilon_d, self.ls_list, Gs):
            print(fmt.format(tls.d, epsilon))
            tls.optimize(epsilon, Gs=Gs_this, skip_setup=skip_setup, **ls_args)
        # end for
        # These will be reset now and updated later if applicable
        self.epsilon_p = None
        self.temperature = None
    # end def

    def optimize_temperature(
        self,
        temperature,
        **ls_args
        # Gs=None, N=500, M=7, fraction=0.025, fit_kind=None, fit_func=None,
        # fit_args={}, bias_mix=0.0, bias_order=1, noise_frac=0.05,
        # W_resolution=0.05, S_resolution=0.05, max_rounds=10
        # W_num=3, W_max=None, sigma_num=3, sigma_max=None
    ):
        print(('Calculating epsilon_d for temperature=' + FF).format(temperature))
        epsilon_d = self._get_thermal_epsilon_d(temperature)
        self.optimize_epsilon_d(epsilon_d, **ls_args)
        self.temperature = temperature
    # end def

    def optimize_epsilon_p(
        self,
        epsilon_p,
        starting_mix=0.5,
        thermal=False,
        Gs=None,
        resolution=0.01,
        verbosity=1,
        **ls_args,
        # N=500, M=7, fraction=0.025, fit_kind=None, fit_func=None,
        # fit_args={}, bias_mix=0.0, bias_order=1, noise_frac=0.05,
        # W_resolution=0.05, S_resolution=0.05, max_rounds=10
        # W_num=3, W_max=None, sigma_num=3, sigma_max=None
    ):
        if Gs is None:
            Gs = len(epsilon_p) * [None]
        # end if
        for tls, Gs_this in zip(self.ls_list, Gs):
            tls.setup_optimization(
                Gs=Gs_this,
                verbosity=verbosity,
                **ls_args
            )
        # end for
        epsilon_p = array(epsilon_p, dtype=float)
        if thermal:
            if verbosity >= 1:
                print("Optimizing to epsilon_p with the thermal constraint")
                self._print_epsilon(epsilon_p)
            # end if
            epsilon_d_opt, T = self._optimize_epsilon_p_thermal(
                epsilon_p,
                resolution=resolution,
                verbosity=verbosity,
            )
            self.optimize_epsilon_d(epsilon_d_opt, skip_setup=True, **ls_args)
            self.temperature = T
        else:
            if verbosity >= 1:
                print("Optimizing to epsilon_p using line-search")
                self._print_epsilon(epsilon_p)
            # end if
            epsilon_d_opt = self._optimize_epsilon_p_ls(
                epsilon_p,
                starting_mix=starting_mix,
                resolution=resolution,
                verbosity=verbosity,
            )
            self.optimize_epsilon_d(epsilon_d_opt, skip_setup=True, **ls_args)
        # end if
        self.epsilon_p = epsilon_p
    # end def

    def _optimize_epsilon_p_thermal(
        self,
        epsilon_p,
        resolution=0.01,  # Relative temperature resolution
        verbosity=1,
    ):
        # initial temperature
        T = self._get_epsilon_p_temperature(epsilon_p) * resolution
        # Initial guess to enter first round
        error_p = epsilon_p * 0
        # First loop: increase T until the errors are no longer capped
        while all(error_p - epsilon_p < 0.0):
            epsilon_d = self._get_thermal_epsilon_d(T)
            if verbosity >= 2:
                max_diff = (error_p - epsilon_p).max()
                print(('  T=' + FF + ', max(error_p-epsilon_p)=' + FF).format(T, max_diff))
            # end if
            error_p = self._resample_errors_p_of_d(epsilon_d, max_rounds=4)
            T *= 1.5
        # end while
        # Second loop: decrease T until the errors are capped
        while not all(error_p - epsilon_p < 0.0):
            T *= 1 - resolution
            epsilon_d = self._get_thermal_epsilon_d(T)
            error_p = self._resample_errors_p_of_d(epsilon_d)
        # end while
        return self._get_thermal_epsilon_d(T), T
    # end def

    def _optimize_epsilon_p_ls(
        self,
        epsilon_p,
        resolution=0.01,
        it_max=10,
        starting_mix=0.5,
        cost_factor=0.5,
        verbosity=1,
    ):
        U = self.hessian.directions
        # Starting mixure of bias + noise
        epsilon_d0 = starting_mix * (
            abs(U @ epsilon_p) + (1 - starting_mix) * U.T @ U @ epsilon_p
        )

        # Internal line-search cost-function: mean-squared relative error + cost
        # constribution. The cost function is inclined to make error_p slightly
        def cost(error_p, rel_stat_cost):
            rel_error = (error_p - epsilon_p) / epsilon_p
            return mean(rel_error**2)**0.5 + rel_stat_cost * cost_factor
        # end def

        epsilon_d_opt = array(epsilon_d0)
        # Initial optimization to get statistical cost
        if verbosity >= 2:
            print("Initial optimization:")
            self._print_epsilon(epsilon_d0, 'd')
        # end if
        for tls, epsilon in zip(self.ls_list, epsilon_d0):
            tls.optimize(epsilon, skip_setup=True)
        # end for

        for it in range(it_max):
            # Squeeze the epsilon_d offset
            old_stat_cost = self.statistical_cost
            R_factor = 0.9**(it + 1)
            epsilon_d_old = epsilon_d_opt.copy()
            # sequential line-search from d0...dD
            for d in range(len(epsilon_d_opt)):
                epsilon_d = epsilon_d_opt.copy()
                R_min = epsilon_d[d] * (1 - R_factor)
                R_max = epsilon_d[d] * (1 + R_factor)
                epsilons = linspace(R_min, R_max, 5)
                costs = []
                for eps_this in epsilons:
                    epsilon_d[d] = eps_this
                    error_p = self._resample_errors_p_of_d(epsilon_d)
                    rel_stat_cost = self.statistical_cost / old_stat_cost
                    costs.append(cost(error_p, rel_stat_cost))
                # end for
                # Argmin is more robust to the uneven shape
                new_epsilon = epsilons[argmin(costs)]
                epsilon_d_opt[d] = new_epsilon
            # end for
            error_p = self._resample_errors_p_of_d(epsilon_d_opt)
            cost_it = cost(error_p, 0.0)
            if verbosity >= 2:
                print(('  iter=' + FI + 'cost=' + FF).format(it, cost_it))
                self._print_epsilon(epsilon_d_opt, 'd')
            # end if
            diff_epsilon_d = mean(abs(epsilon_d_old - epsilon_d_opt))
            if cost_it < resolution or diff_epsilon_d < resolution / 10:
                break
            # end if
        # end for
        # scale down
        while any((error_p - epsilon_p) > 0.0):
            if mean(error_p / epsilon_p) > 5 * resolution:
                epsilon_d_opt *= 1 - 5 * resolution
            else:
                epsilon_d_opt *= 1 - resolution
            # end if
            error_p = self._resample_errors_p_of_d(epsilon_d_opt)
        # end for
        return epsilon_d_opt
    # end def

    def validate(self, N=500, thr=1.1):
        """Validate optimization by independent random resampling
        """
        if not self.optimized:
            raise AssertionError('Must be optimized first')
        # end if
        ref_error_p, ref_error_d = self._resample_errors(N=N)
        valid = not (
            any(ref_error_p > self.error_p * thr) or any(ref_error_d > self.error_d * thr)
        )
        return valid
    # end def

    def bracket_target_biases(
        self,
        **kwargs
        # bias_order=1, mx_iter=10, M=7, bracket_fraction=0.5
        # sgn=1, fraction=0.025, fit_func=None, bias_mix=None
        # N=200, Gs=None, fraction=None
    ):
        for tls in self.ls_list:
            tls.bracket_target_bias(**kwargs)
        # end for
    # end def

    def copy(
        self,
        path='',
        **kwargs
        # pes=None, pes_func=None, pes_args={}
    ):
        if not self.optimized:
            raise AssertionError("Must optimize surrogate before copying")
        # end if
        pls = ParallelLineSearch.copy(
            self,
            path=path,
            # Copy optimized windows, noises
            windows=self.W_opt,
            noises=self.sigma_opt,
            **kwargs
            # pes=None, pes_func=None, pes_args={}
        )
        # Copy each optimized line-search settings
        for ls_new, tls in zip(pls.ls_list, self.ls_list):
            ls_new._settings = LsSettings.copy(tls.target_settings)
        # end for
        return pls
    # end def

    # Finalize optimization by computing error estimates
    def _finalize_optimization(self, write=None, overwrite=False):
        # The errors are calculated strictly based on settings stored in line-searches
        errors = self._resample_errors()
        self.error_d, self.error_p = errors

        print('  Optimization complete:')
        self._print_optimization('d', self.D_list, self.epsilon_d, self.error_d)
        if self.epsilon_p is not None:
            self._print_optimization('p', range(len(self.epsilon_p)), self.epsilon_p, self.error_p)
        # end if

        if isinstance(write, str):
            self.write_to_disk(fname=write, overwrite=overwrite)
        # end if
    # end def

    def _print_optimization(self, label, i_list, epsilon, error):
        # Vertical print
        if epsilon is None:
            return
        # end if
        print(('    ' + FIS + FFS + FFS + FPS).format(label, 'target', 'error', 'rel. '))
        for i, eps, err in zip(i_list, epsilon, error):
            rel_err = err / eps
            print(('    ' + FI + FF + FF + FP).format(i, eps, err, rel_err * 100))
        # end for
    # end def

    def _print_epsilon(self, epsilon, label='p'):
        for p, eps in enumerate(epsilon):
            print(('  {}=' + FI + FF).format(label, p, eps))
        # end for
    # end def

    def _resample_errors(self):
        biases_d, biases_p = self._compute_bias()
        errorbar_d, errorbar_p = self._resample_errorbars()
        error_d = biases_d + errorbar_d
        error_p = biases_p + errorbar_p
        return error_d, error_p
    # end def

    def _compute_bias(self):
        biases_d = []
        for tls in self.ls_list:
            offsets = tls.figure_out_adjusted_offsets(W=tls.W_opt, M=tls.M)
            values = tls.evaluate_target(offsets)
            grid = LineSearchGrid(offsets=offsets, values=values)
            biases_d.append(tls.compute_bias(grid))
        # end for
        biases_d = array(biases_d)
        biases_p = self._calculate_params_next(self.params, biases_d)
        biases_p -= self.params
        return biases_d, biases_p
    # end def

    # based on windows and noises
    def _resample_errorbars(self, N=None):
        # provide correlated sampling
        fraction = self.ls_list[0].target_settings.fraction
        # list of distributions of minima per direction, parameter
        x0s_d, x0s_p = [], []
        # list of statistical errorbars per direction, parameter
        biases_d, errorbar_d, errorbar_p = [], [], []
        for tls in self.ls_list:
            # Get and store x0 distribution based on the optimized fit
            offsets = tls.figure_out_adjusted_offsets(W=tls.W_opt, M=tls.M)
            values = tls.evaluate_target(offsets)
            errors = tls.M * [tls.sigma_opt]
            grid = LineSearchGrid(offsets=offsets, values=values, errors=errors)
            bias_d = tls._compute_xy_bias(
                tls.grid_opt,
                tls.target_settings
            )[0]
            biases_d.append(bias_d)
            if N is None:
                # Correlated sampling
                x0s = tls.target_settings.fit_func.get_x0_distribution(grid, Gs=tls.Gs)
            else:
                # Validation
                x0s = tls.target_settings.fit_func.get_x0_distribution(grid, N=N)
            # end if
            x0s_d.append(x0s)
            errorbar_d.append(get_fraction_error(x0s - bias_d, fraction)[1])
        # end for
        biases_p = self._calculate_params_next(self.params, biases_d)
        # Use directional x0 distribution to linearly map to parameter errorbars
        for x0 in array(x0s_d).T:
            x0_p = self._calculate_params_next(0, x0)
            x0s_p.append(x0_p)
        # end for
        x0s_p = array(x0s_p) - biases_p
        errorbar_p = [get_fraction_error(x0s, fraction)[1] for x0s in array(x0s_p).T]
        return array(errorbar_d), array(errorbar_p)
    # end def

    def _get_thermal_epsilon_d(self, temperature):
        return [(temperature / abs(Lambda))**0.5 for Lambda in self.Lambdas]
    # end def

    def _get_epsilon_d_temperature(self, epsilon_d):
        # return mean temperature given epsilon_d
        return mean(epsilon_d**2 * self.Lambdas)
    # end def

    def _get_thermal_epsilon_p(self, temperature):
        return [(temperature / abs(Lambda))**0.5 for Lambda in self.hessian.hessian.diagonal]
    # end def

    def _get_epsilon_p_temperature(self, epsilon_p):
        # return mean temperature given epsilon_p
        return mean(epsilon_p**2 * self.hessian.hessian.diagonal())
    # end def

    def _resample_errors_p_of_d(
        self,
        epsilon_d,
        **kwargs  # max_rounds=10
    ):
        for epsilon, tls, in zip(epsilon_d, self.ls_list):
            tls.optimize(epsilon, skip_setup=True, **kwargs)
        # end for
        return self._resample_errors()[1]
    # end def

    def plot_error_surfaces(self, **kwargs):
        for tls in self.ls_list:
            tls.plot_error_surface(**kwargs)
        # end for
    # end def

# end class
