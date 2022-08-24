#!/usr/bin/env python3

from numpy import array, mean, linspace, argmin
from functools import partial
from scipy.optimize import broyden1

from lib.util import get_fraction_error
from lib.targetlinesearch import TargetLineSearch
from lib.parallellinesearch import ParallelLineSearch


class TargetParallelLineSearch(ParallelLineSearch):
    epsilon_p = None
    epsilon_d = None
    error_p = None
    error_d = None
    temperature = None
    window_frac = None
    targets = None
    # FLAGS
    optimized = False

    def __init__(
        self,
        structure = None,
        hessian = None,
        targets = None,
        **kwargs
    ):
        ParallelLineSearch.__init__(self, structure = structure, hessian = hessian, **kwargs)
        self.set_targets(targets)
        self.ls_list = self._generate_tls_list(**kwargs)
    #end def

    def set_targets(self, targets):
        # whether target or not
        if targets is None:
            targets = self.D * [0.0]
        else:
            assert(len(targets) == self.D)
        #end if
        self.targets = array(targets)
    #end def

    def optimize(
        self,
        windows = None,
        noises = None,
        epsilon_p = None,
        epsilon_d = None,
        temperature = None,
        **kwargs,
    ):
        """Optimize parallel line-search for noise using different constraints. The optimizer modes are called in the following order of priority based on input parameters provided:
  1) windows, noises (list, list)
     when windows and noises per DIRECTION are provided, they are allocated to each DIRECTION and parameter errors resampled
  2) temperature (float > 0)
    temperature: windows and noises per DIRECTION are optimized to meet DIRECTIONAL tolerances based on thermal equipartition, parameter errors resampled
  3) epsilon_d (list)
    windows and noises per DIRECTION are optimized to meet DIRECTIONAL tolerances
  4) epsilon_p (list)
    windows and noises per DIRECTION are optimized for total predicted cost while meeting PARAMETER tolerances
    Guided by a keyword:
      kind = 'ls': [DEFAULT] use an auxiliary line-search method for optimization
      kind = 'thermal': find maximum temperature to maintain PARAMETER errors
      kind = 'broyden1': use broyden1 minimization to find the optimal solution (very unstable)
useful keyword arguments:
  M = number of points
  N = number of points for resampling
  bias_mix = mixing of energy bias to parameter bias
  fit_kind = fitting function
        """
        if windows is not None and noises is not None:
            self.optimize_windows_noises(windows, noises, **kwargs)
        elif temperature is not None:
            self.optimize_thermal(temperature, **kwargs)
        elif epsilon_d is not None:
            self.optimize_epsilon_d(epsilon_d, **kwargs)
        elif epsilon_p is not None:
            self.optimize_epsilon_p(epsilon_p, **kwargs)
        else:
            raise AssertionError('Optimizer constraint not identified')
        #end if
    #end def

    # All optimizer methods conclude here, where the final errors and key parameters are also stored to the object
    def optimize_windows_noises(self, windows, noises, M = None, fit_kind = None, bias_mix = None, **kwargs):
        M = M if M is not None else self.M
        fit_kind = fit_kind if fit_kind is not None else self.fit_kind
        self.error_d, self.error_p = self._errors_windows_noises(windows, noises, M = M, fit_kind = fit_kind, bias_mix = bias_mix, **kwargs)
        self.M = M
        self.fit_kind = fit_kind
        self.windows = windows
        self.noises = noises
        self.optimized = True
    #end def

    def _errors_windows_noises(self, windows, noises, Gs = None, fit_kind = None, M = None, **kwargs):
        Gs_d = Gs if Gs is not None else self.D * [None]
        return self._resample_errors(windows, noises, Gs = Gs_d, M = M, fit_kind = fit_kind, **kwargs)
    #end def

    def optimize_thermal(self, temperature, **kwargs):
        assert temperature > 0, 'Temperature must be positive'
        self.optimize_epsilon_d(self._get_thermal_epsilon_d(temperature), **kwargs)
    #end def

    # TODO: fixed-point method, also need to init error matrices
    def optimize_epsilon_p(
        self,
        epsilon_p,
        kind = 'ls',  # try line-search by default
        **kwargs,
    ):
        epsilon_d0 = epsilon_p.copy()  # TODO: fix
        if kind == 'ls':
            epsilon_d_opt = self._optimize_epsilon_p_ls(epsilon_p, epsilon_d0, **kwargs)
        elif kind == 'thermal':
            epsilon_d_opt = self._optimize_epsilon_p_thermal(epsilon_p, **kwargs)
        elif kind == 'broyden1':
            # Current: broyden1, probably fails to converge
            validate_epsilon_d = partial(self._resample_errors_p_of_d, target = array(epsilon_p), **kwargs)
            epsilon_d_opt = broyden1(validate_epsilon_d, epsilon_d0, f_tol = 1e-3, verbose = True)
        else:
            raise AssertionError('Fixed-point kind not recognized')
        #end if
        kwargs_d = kwargs.copy()
        kwargs_d.update(fix_res = False)
        self.optimize_epsilon_d(epsilon_d_opt, **kwargs_d)
        self.epsilon_p = epsilon_p
    #end def

    # TODO: check for the first step
    def _optimize_epsilon_p_thermal(self, epsilon_p, T0 = 0.00001, dT = 0.000005, verbose = False, **kwargs):
        T = T0
        error_p = array([-1, -1])  # init
        #first = True
        while all(error_p < 0.0):
            try:
                epsilon_d = self._get_thermal_epsilon_d(T)
                error_p = self._resample_errors_p_of_d(epsilon_d, target = epsilon_p, verbose = verbose)
                if verbose:
                    print('T = {} highest error {} %'.format(T, (error_p + epsilon_p) / epsilon_p * 100))
                #end if
                T += dT
            except AssertionError:
                if verbose:
                    print('T = {} skipped'.format(T))
                #end if
                T += dT
            #end try
        #end while
        T -= 2 * dT
        return self._get_thermal_epsilon_d(T)
    #end def

    def _optimize_epsilon_p_ls(
        self,
        epsilon_p,
        epsilon_d0,
        thr = None,
        it_max = 10,
        **kwargs
    ):
        thr = thr if thr is not None else mean(epsilon_p) / 20

        def cost(derror_p):
            return sum(derror_p**2)**0.5
        #end def
        epsilon_d_opt = array(epsilon_d0)
        for it in range(it_max):
            coeff = 0.5**(it + 1)
            epsilon_d_old = epsilon_d_opt.copy()
            # sequential line-search from d0...dD
            for d in range(len(epsilon_d_opt)):
                epsilon_d = epsilon_d_opt.copy()
                epsilons = linspace(epsilon_d[d] * (1 - coeff), (1 + coeff) * epsilon_d[d], 10)
                costs = []
                for s in epsilons:
                    epsilon_d[d] = s
                    derror_p = self._resample_errors_p_of_d(epsilon_d, target = epsilon_p, fix_res = False, **kwargs)
                    costs.append(cost(derror_p))
                #end for
                epsilon_d_opt[d] = epsilons[argmin(costs)]
            #end for
            derror_p = self._resample_errors_p_of_d(epsilon_d_opt, target = epsilon_p, **kwargs)
            cost_it = cost(derror_p)
            # scale down
            if cost_it < thr or sum(abs(epsilon_d_old - epsilon_d_opt)) < thr / 100:
                break
            #end if
        #end for
        for c in range(100):
            if any(derror_p > 0.0):
                epsilon_d_opt = [e * 0.99 for e in epsilon_d_opt]
                derror_p = self._resample_errors_p_of_d(epsilon_d_opt, target = epsilon_p, **kwargs)
            else:
                break
            #end if
        #end for
        return epsilon_d_opt
    #end def

    def optimize_epsilon_d(
        self,
        epsilon_d,
        Gs = None,
        **kwargs,
    ):
        Gs_d = Gs if Gs is not None else self.D * [None]
        assert len(Gs_d) == self.D, 'Must provide list of Gs equal to the number of directions'
        windows, noises = [], []
        for epsilon, ls, Gs in zip(epsilon_d, self.ls_list, Gs_d):
            ls.optimize(epsilon, Gs = Gs, **kwargs)
            windows.append(ls.W_opt)
            noises.append(ls.sigma_opt)
        #end for
        self.optimize_windows_noises(windows, noises, Gs = Gs_d, **kwargs)
        self.epsilon_d = epsilon_d
    #end def

    def get_Gs(self):
        return [ls.Gs for ls in self.ls_list]
    #end def

    def validate(self, N = 500, verbose = False, thr = 1.1):
        """Validate optimization by independent random resampling"""
        assert self.optimized, 'Must be optimized first'
        ref_error_p, ref_error_d = self._resample_errors(self.windows, self.noises, Gs = None, N = N)
        valid = True
        if verbose:
            print('Parameter errors:')
            print('  {:<2s}  {:<10s}  {:<10s}  {:<5s}'.format('#', 'validation', 'corr. samp.', 'ratio'))
        #end if
        for p, ref, corr in zip(range(len(ref_error_p)), ref_error_p, self.error_p):
            ratio = ref / corr
            valid_this = ratio < thr
            valid = valid and valid_this
            if verbose:
                print('  {:<2d}  {:<10f}  {:<10f}  {:<5f}'.format(p, ref, corr, ratio))
            #end if
        #end for
        if verbose:
            print('Direction errors:')
            print('  {:<2s}  {:<10s}  {:<10s}  {:<5s}'.format('#', 'validation', 'corr. samp.', 'ratio'))
        #end if
        for d, ref, corr in zip(range(len(ref_error_d)), ref_error_d, self.error_d):
            ratio = ref / corr
            valid_this = ratio < thr
            valid = valid and valid_this
            if verbose:
                print('  {:<2d}  {:<10f}  {:<10f}  {:<5f}'.format(d, ref, corr, ratio))
            #end if
        #end for
        return valid
    #end def

    def _get_thermal_epsilon_d(self, temperature):
        return [(temperature / Lambda)**0.5 for Lambda in self.Lambdas]
    #end def

    def _get_thermal_epsilon(self, temperature):
        return [(temperature / Lambda)**0.5 for Lambda in self.hessian.diagonal]
    #end def

    # Important: overrides ls_list generator
    def _generate_ls_list(self, **kwargs):
        return None
    #end def

    def _generate_tls_list(self, **kwargs):
        noises = self.noises if self.noisy else self.D * [None]
        ls_list = []
        for d, window, noise, target in zip(range(self.D), self.windows, noises, self.targets):
            tls = TargetLineSearch(
                structure = self.structure,
                hessian = self.hessian,
                d = d,
                W = window,
                sigma = noise,
                target = target,
                **kwargs)
            ls_list.append(tls)
        #end for
        return ls_list
    #end def

    def compute_bias_p(self, **kwargs):
        return self.compute_bias(**kwargs)[1]
    #end def

    def compute_bias_d(self, **kwargs):
        return self.compute_bias(**kwargs)[0]
    #end def

    def compute_bias(self, windows = None, **kwargs):
        windows = windows if windows is not None else self.windows
        return self._compute_bias(windows, **kwargs)
    #end def

    def _compute_bias(self, windows, **kwargs):
        bias_d = []
        for W, tls in zip(windows, self.ls_list):
            assert W <= tls.W_max, 'window is larger than W_max'
            grid = tls._figure_out_grid(W = W)[0]
            bias_d.append(tls.compute_bias(grid, **kwargs)[0])
        #end for
        bias_d = array(bias_d)
        bias_p = self._calculate_params_next(self.get_params(), self.get_directions(), bias_d) - self.get_params()
        return bias_d, bias_p
    #end def

    # based on windows, noises
    def _resample_errorbars(self, windows, noises, Gs = None, N = None, M = None, fit_kind = None, **kwargs):
        Gs_d = self.D * [None] if Gs is None else Gs  # provide correlated sampling
        biases_d, biases_p = self._compute_bias(windows, fit_kind = fit_kind, **kwargs)  # biases per direction, parameter
        x0s_d, x0s_p = [], []  # list of distributions of minima per direction, parameter
        errorbar_d, errorbar_p = [], []  # list of statistical errorbars per direction, parameter
        for tls, W, noise, bias_d, Gs in zip(self.ls_list, windows, noises, biases_d, Gs_d):
            assert W <= tls.W_max, 'window is larger than W_max'
            grid, M = tls._figure_out_grid(W = W, M = M)
            values = tls.evaluate_target(grid)
            errors = M * [noise]
            Gs = Gs if Gs is not None else tls.Gs
            x0s = tls.get_x0_distribution(grid = grid, values = values, errors = errors, Gs = Gs, N = N, fit_kind = fit_kind)
            x0s_d.append(x0s)
            errorbar_d.append(get_fraction_error(x0s - bias_d, self.fraction)[1])
        #end for
        # parameter errorbars
        for x0 in array(x0s_d).T:
            x0s_p.append(self._calculate_params_next(-biases_p, self.directions, x0))  # TODO: could this be vectorized?
        #end for
        errorbar_p = [get_fraction_error(x0s, self.fraction)[1] for x0s in array(x0s_p).T]
        return array(errorbar_d), array(errorbar_p)
    #end def

    def _resample_errors(self, windows, noises, **kwargs):
        bias_d, bias_p = self._compute_bias(windows, **kwargs)
        errorbar_d, errorbar_p = self._resample_errorbars(windows, noises, **kwargs)
        error_d = abs(bias_d) + errorbar_d
        error_p = abs(bias_p) + errorbar_p
        return error_d, error_p
    #end def

    def _resample_errors_p_of_d(self, epsilon_d, target = 0.0, **kwargs):
        windows, noises = self._windows_noises_of_epsilon_d(epsilon_d, **kwargs)
        return self._resample_errors(windows, noises, **kwargs)[1] - target
    #end def

    def _windows_noises_of_epsilon_d(
        self,
        epsilon_d,
        **kwargs,
    ):
        windows, noises = [], []
        for epsilon, ls, in zip(epsilon_d, self.ls_list):
            W_opt, sigma_opt = ls.maximize_sigma(epsilon, **kwargs)  # no altering the error
            #W_opt, sigma_opt = ls.interpolate_max_sigma(abs(epsilon))
            windows.append(W_opt)
            noises.append(sigma_opt)
        #end for
        return windows, noises
    #end def

    def get_biases_d(self, windows = None):
        windows = windows if windows is not None else self.windows
        biases_d = array([ls.compute_bias_of(W = ls.W_opt)[0][0] for ls in self.ls_list])
        return biases_d
    #end def

    def get_biases_p(self, windows = None):
        biases_p = self._calculate_shifts(self.directions, self.get_biases_d(windows = windows))
        return biases_p
    #end def

#end class