#!/usr/bin/env python3
'''TargetLineSearch classes for the assessment and evaluation of fitting errors
'''

from numpy import array, argsort, isscalar, linspace, append, nan, isnan, where
from numpy import random, argmin, equal

from .LineSearch import LineSearch
from .TargetLineSearchBase import TargetLineSearchBase

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Class for error scan line-search
class TargetLineSearch(TargetLineSearchBase, LineSearch):

    E_mat = None  # resampled W-sigma matrix of errors
    W_mat = None  # resampled W-mesh
    S_mat = None  # resampled sigma-mesh
    T_mat = None  # resampled trust-mesh (whether error is reliable)
    Gs = None  # N x M set of correlated random fluctuations for the grid
    epsilon = None  # optimized target error
    W_opt = None  # W to meet epsilon
    sigma_opt = None  # sigma to meet epsilon
    # FLAGS
    resampled = False
    optimized = False

    def __init__(
        self,
        structure=None,
        hessian=None,
        d=None,
        M=7,
        W=None,  # characteristic window
        R=None,  # max displacement
        grid=None,  # manual set of shifts
        values=None,
        mode='jobs',
        pes=None,
        **kwargs,  # some redundancy in submitting kwargs
    ):
        # provide target_grid, target_values explicitly
        TargetLineSearchBase.__init__(self, **kwargs)
        LineSearch.__init__(
            self,
            structure=structure,
            hessian=hessian,
            d=d,
            M=M,
            W=W,
            R=R,
            grid=grid,
            values=values,
            mode=mode,
            pes=pes,
            **kwargs,
        )
    # end def

    def compute_bias_of(
        self,
        R=None,
        W=None,
        verbose=False,
        bias_mix=None,
        **kwargs
    ):
        bias_mix = bias_mix if bias_mix is not None else self.bias_mix
        biases_x, biases_y, biases_tot = [], [], []
        if verbose:
            print((4 * '{:<10s} ').format('R', 'bias_x', 'bias_y', 'bias_tot'))
        # end if
        if R is not None:
            Rs = array([R]) if isscalar(R) else R
            for R in Rs:
                R = max(R, 1e-5)  # for numerical stability
                grid, M = self._figure_out_grid(R=R, **kwargs)
                bias_x, bias_y, bias_tot = self.compute_bias(
                    grid, bias_mix=bias_mix, **kwargs)
                biases_x.append(bias_x)
                biases_y.append(bias_y)
                biases_tot.append(bias_tot)
                if verbose:
                    print((4 * '{:<10f} ').format(R, bias_x, bias_y, bias_tot))
                # end if
            # end for
        elif W is not None:
            Ws = array([W]) if isscalar(W) else W
            if verbose:
                print((4 * '{:<10s} ').format('W',
                      'bias_x', 'bias_y', 'bias_tot'))
            # end if
            for W in Ws:
                W = max(W, 1e-5)  # for numerical stability
                grid, M = self._figure_out_grid(W=W, **kwargs)
                bias_x, bias_y, bias_tot = self.compute_bias(
                    grid, bias_mix=bias_mix, **kwargs)
                biases_x.append(bias_x)
                biases_y.append(bias_y)
                biases_tot.append(bias_tot)
                if verbose:
                    print((4 * '{:<10f} ').format(W, bias_x, bias_y, bias_tot))
                # end if
            # end for
        # end if
        return array(biases_x), array(biases_y), array(biases_tot)
    # end def

    # override TargetLineSearchBase class
    def set_target(self, grid, values, **kwargs):
        if values is None or all(equal(array(values), None)):
            return
        # end if
        TargetLineSearchBase.set_target(self, grid, values, **kwargs)
        if self.structure_list is None:
            return
        # end if
        for s, v in zip(self.structure_list, values):
            s.value = v
            s.error = 0.0
        # end for
    # end def

    def evaluate_target(self, grid):
        try:
            return TargetLineSearchBase.evaluate_target(self, grid)
        except AssertionError:
            print('  W_max and R_max are respectively {} and {}'.format(
                self.W_max, self.R_max))
            return None
        # end try
    # end def

    @property
    def R_max(self):
        return min([-self.grid.min(), self.grid.max()])
    # end def

    @property
    def W_max(self):
        return self._R_to_W(self.R_max)
    # end def

    def _W_sigma_of_epsilon(self, epsilon, **kwargs):
        W, sigma = self._contour_max_sigma(epsilon)
        return W, sigma
    # end def

    # X: grid of W values; Y: grid of sigma values; E: grid of total errors
    #   if Gs is not provided, use M and N
    def generate_W_sigma_data(
        self,
        W_num=10,
        W_max=None,
        sigma_num=10,
        sigma_max=None,
        noise_frac=0.05,
        Gs=None,
        M=None,
        N=None,
        fit_kind=None,
        verbose=False,
        **kwargs
    ):
        if Gs is None:
            M = M if M is not None else self.M
            self.regenerate_Gs(M, N)  # set (or reset) Gs and M appropriately
        else:  # Gs is not None
            N, M = Gs.shape
            assert N > 1, 'Must provide N > 1'
            assert M > 1, 'Must provide M > 1'
            self.Gs = Gs
            self.M = M
        # end if
        if verbose:
            print(
                '  tls{}: generating {}-by-{} error surface with N={}'.format(self.d, W_num, sigma_num, N))
        # end if
        W_max = W_max if W_max is not None else self.W_max
        fit_kind = fit_kind if fit_kind is not None else self.fit_kind
        sigma_max = sigma_max if sigma_max is not None else W_max * noise_frac
        # starting window array: sigma = 0, so only bias
        Ws = linspace(0.0, W_max, W_num)
        sigmas = linspace(0.0, sigma_max, sigma_num)
        errors = array(self.M * sigmas[0])
        E_this = [self._compute_error(self._make_grid_W(
            W, self.M), errors, Gs=self.Gs, fit_kind=fit_kind, **kwargs) for W in Ws]
        self.fit_kind = fit_kind
        self.E_mat = array([E_this])
        self.W_mat = array([Ws])
        self.S_mat = array([W_num * [sigmas[0]]])
        self.T_mat = self._generate_T_mat()
        # append the rest
        for sigma in sigmas[1:]:
            self._insert_sigma_data(sigma, Gs=self.Gs, fit_kind=fit_kind)
        # end for
        self.resampled = True
    # end def

    def _generate_T_mat(self):
        return self.W_mat >= self.S_mat
    # end def

    def _check_Gs_M_N(self, Gs, M, N):
        """Return True if Gs (and derived quantities) do not need regeneration, otherwise return False"""
        if self.Gs is None:
            return False
        # end if
        Gs = Gs if Gs is not None else self.Gs  # only Gs will matter when provided
        M = M if M is not None else self.M  # check user input
        N = N if N is not None else len(self.Gs)  # check user input
        return (len(Gs) == N and len(Gs[0]) == M)
    # end def

    def regenerate_Gs(self, M, N):
        """Regenerate and save Gs array"""
        assert N is not None and N > 1, 'Must provide N > 1'
        assert M is not None and M > 1, 'Must provide M > 1'
        self.Gs = random.randn(N, M)
        self.M = M
    # end def

    def insert_sigma_data(self, sigma, **kwargs):
        self._insert_sigma_data(
            sigma, Gs=self.Gs, fit_kind=self.fit_kind, **kwargs)
    # end def

    def insert_W_data(self, W, **kwargs):
        assert W < self.W_max, 'Cannot resample past W_max={:<9f}; extend the target data'.format(
            self.W_max)
        self._insert_W_data(W, Gs=self.Gs, fit_kind=self.fit_kind, **kwargs)
    # end def

    def _insert_sigma_data(self, sigma, **kwargs):
        Ws = self.W_mat[0]
        Es = [self._compute_error(self._make_grid_W(
            W, self.M), self.M * [sigma], **kwargs) for W in Ws]
        W_mat = append(self.W_mat, [Ws], axis=0)
        S_mat = append(self.S_mat, [len(Ws) * [sigma]], axis=0)
        E_mat = append(self.E_mat, [Es], axis=0)
        idx = argsort(S_mat[:, 0])
        self.W_mat = W_mat[idx]
        self.S_mat = S_mat[idx]
        self.E_mat = E_mat[idx]
        self.T_mat = self._generate_T_mat()
    # end def

    def _insert_W_data(self, W, **kwargs):
        sigmas = self.S_mat[:, 0]
        grid = self._make_grid_W(W, self.M)
        Es = [self._compute_error(grid, self.M * [sigma], **kwargs)
              for sigma in sigmas]
        W_mat = append(self.W_mat, array([len(sigmas) * [W]]).T, axis=1)
        S_mat = append(self.S_mat, array([sigmas]).T, axis=1)
        E_mat = append(self.E_mat, array([Es]).T, axis=1)
        idx = argsort(W_mat[0])
        self.W_mat = W_mat[:, idx]
        self.S_mat = S_mat[:, idx]
        self.E_mat = E_mat[:, idx]
        self.T_mat = self._generate_T_mat()
    # end def

    def optimize(self, epsilon, **kwargs):
        """Optimize W and sigma to a given target error epsilon > 0."""
        self.W_opt, self.sigma_opt = self.maximize_sigma(epsilon, **kwargs)
        self.epsilon = epsilon
        self.optimized = True
    # end def

    def maximize_sigma(
        self,
        epsilon,
        allow_override=True,  # allow regeneration of errors
        fix_res=True,
        Gs=None,
        M=None,
        N=None,
        verbose=True,
        low_thr=0.9,
        **kwargs
    ):
        """Optimize W and sigma based on maximizing sigma."""
        if self.resampled:
            if not self._check_Gs_M_N(Gs, M, N):
                if allow_override:
                    self.generate_W_sigma_data(
                        Gs=Gs, M=M, N=N, verbose=verbose, **kwargs)
                else:
                    raise AssertionError('Requested inconsistent resampling.')
                # end if
            # end if
        else:
            self.generate_W_sigma_data(
                Gs=Gs, M=M, N=N, verbose=verbose, **kwargs)
        # end if
        W, sigma, E, errs = self._maximize_y(epsilon, low_thr=low_thr)
        if 'not_found' in errs:
            raise AssertionError(
                '  tls{}: W, sigma not found for epsilon = {}. Check minimum bias and raise epsilon.'.format(self.d, epsilon))
        # end if
        if fix_res:
            while 'x_underflow' in errs and self._fix_x_underflow(W, verbose=verbose):
                W, sigma, E, errs = self._maximize_y(epsilon, low_thr=low_thr)
            # end while
            while 'y_underflow' in errs and self._fix_y_underflow(sigma, verbose=verbose):
                W, sigma, E, errs = self._maximize_y(epsilon, low_thr=low_thr)
            # end while
            while 'y_overflow' in errs and self._fix_y_overflow(sigma, verbose=verbose):
                W, sigma, E, errs = self._maximize_y(epsilon, low_thr=low_thr)
            # end while
            while 'low_res' in errs and self._fix_low_res(epsilon, verbose=verbose):
                W, sigma, E, errs = self._maximize_y(epsilon, low_thr=low_thr)
            # end while
        return W, sigma
    # end def

    def _fix_x_underflow(self, W_this, verbose=True):
        W_new = self.W_mat[0, 1] / 2
        if W_new < self.W_max * 1e-3:
            if verbose:
                print('  tls{}: W underflow: did not add W = {}'.format(
                    self.d, W_new.round(7)))
            # end if
            return False
        # end if
        self.insert_W_data(W_new)
        if verbose:
            print('  tls{}: W underflow: added W = {} to resampling grid'.format(
                self.d, (W_new).round(7)))
        # end if
        return True
    # end def

    def _fix_y_underflow(self, S_this, verbose=True):
        S_new = self.S_mat[1, 0] / 2
        if S_new < self.W_max * 1e-8:
            return False
        # end if
        self.insert_sigma_data(S_new)
        if verbose:
            print('  tls{}: Sigma underflow: added sigma = {} to resampling grid'.format(
                self.d, (S_new).round(7)))
        # end if
        return True
    # end def

    def _fix_y_overflow(self, S_this, verbose=True):
        S_new = S_this * 2
        if S_new > self.W_max:
            if verbose:
                print('  tls{}: Sigma overflow: did not add sigma = {} to resampling grid'.format(
                    self.d, S_new))
            # end if
            return False
        # end if
        self.insert_sigma_data(S_new)
        if verbose:
            print('  tls{}: Sigma overlow: added sigma = {} to resampling grid'.format(
                self.d, (S_new).round(7)))
        # end if
        return True
    # end def

    def _fix_low_res(self, epsilon, verbose=True):
        status = False
        Wi, Si = self._argmax_y(self.E_mat, self.T_mat, epsilon)
        W_this, S_this = self.W_mat[0, Wi], self.S_mat[Si, 0]
        if Si < self.S_mat.shape[0] - 1:
            S_new = (self.S_mat[Si, 0] + self.S_mat[Si + 1, 0]) / 2
            if abs(S_new - S_this) > self.S_mat.max() * 5e-3:
                self.insert_sigma_data(S_new)
                status = True
                if verbose:
                    print(
                        '  tls{}: Low-res: added sigma = {} to resampling grid'.format(self.d, S_new.round(7)))
                # end if
            # end if
        # end if
        W_lo = (self.W_mat[0, Wi - 1] + self.W_mat[0, Wi]) / 2
        if Wi < len(self.W_mat[0]) - 1:  # whether to add high W value
            W_hi = (self.W_mat[0, Wi] + self.W_mat[0, Wi + 1]) / 2
            if abs(W_hi - W_this) > self.W_max * 2e-2:
                status = True
                self.insert_W_data(W_hi)
                if verbose:
                    print(
                        '  tls{}: Low-res: added W = {} to resampling grid'.format(self.d, W_hi.round(7)))
                # end if
            # end if
        # end if
        if abs(W_lo - W_this) > self.W_max * 2e-2:
            self.insert_W_data(W_lo)
            status = True
            if verbose:
                print(
                    '  tls{}: Low-res: added W = {} to resampling grid'.format(self.d, W_lo.round(7)))
            # end if
        # end if
        return status
    # end def

    def interpolate_max_sigma(self, epsilon, low_thr=0.9, **kwargs):
        W, sigma, E, errs = self._maximize_y(epsilon, low_thr=low_thr)
        # TODO: bilinear interpolation
        if any([err in errs for err in ['y_overflow', 'x_underflow']]):
            return W, sigma
        else:
            Wi, Si = self._argmax_y(self.E_mat, epsilon)
            a = (self.E_mat[Si + 1, Wi] - epsilon) / \
                (self.E_mat[Si, Wi] - epsilon)
            sigma = a * (self.S_mat[Si, Wi] - self.S_mat[Si + 1, Wi])
            return W, sigma
        # end if
    # end def

    def _maximize_y(self, epsilon, low_thr=0.9):
        """Return X, Y, and E values """
        assert self.resampled, 'Must resample errors first!'
        assert low_thr < 0.99, 'Threshold limit too high'
        assert epsilon > 0, 'epsilon must be positive'
        errs = []
        xi, yi = self._argmax_y(self.E_mat, self.T_mat, epsilon)
        if isnan(xi):
            return nan, nan, nan, ['not_found']
        # end if
        E0 = self.E_mat[yi, xi]
        x0 = self.W_mat[yi, xi]
        y0 = self.S_mat[yi, xi]
        if xi == 0:
            errs.append('x_underflow')
        elif xi == self.E_mat.shape[1] - 1:
            errs.append('x_overflow')
        # end if
        if yi == 0:
            errs.append('y_underflow')
        elif yi == self.E_mat.shape[0] - 1:
            errs.append('y_overflow')
        # end if
        if E0 / epsilon < low_thr:
            errs.append('low_res')
        # end if
        return x0, y0, E0, errs
    # end def

    def _argmax_y(self, E, T, epsilon):
        """Return indices to the highest point in E matrix that is lower than epsilon"""
        xi, yi = nan, nan
        W = self.W_mat[0]
        for i in range(len(E), 0, -1):  # from high to low
            err = where((E[i - 1] < epsilon) & (T[i - 1]))
            if len(err[0]) > 0:
                yi = i - 1
                # xi = err[0][argmax(E[i - 1][err[0]])]
                # take the middle
                xi = err[0][argmin(
                    abs(W[err[0]] - (W[err[0][0]] + W[err[0][-1]]) / 2))]
                break
            # end if
        # end for
        return xi, yi
    # end def

    def statistical_cost(self, sigma=None, M=None):
        """Return statistical cost based on sigma and M"""
        sigma = sigma if sigma is not None else self.sigma_opt
        M = M if M is not None else self.M
        return M * sigma**-2
    # end def

    def plot_error_surface(
        self,
        ax=None,
        **kwargs
    ):
        if not self.optimized:
            print('Must optimize before plotting error surface')
            return
        # end if
        from matplotlib import pyplot as plt
        if ax is None:
            f, ax = plt.subplots(1, 1)
        # end if
        T = self.T_mat
        X = self.W_mat
        Y = self.S_mat
        Z = self.E_mat
        Z[~T] = nan
        ax.contourf(X, Y, Z)
        ax.contour(X, Y, Z, [self.epsilon], colors='k')
        ax.plot(X.flatten(), Y.flatten(), 'k.', alpha=0.3)
        ax.plot(self.W_opt, self.sigma_opt, 'ko')
    # end def

# end class
