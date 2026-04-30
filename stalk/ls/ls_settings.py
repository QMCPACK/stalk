#!/usr/bin/env python3
'''Line-search settings'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import isscalar
from stalk.ls.morse_fit import MorseFit
from stalk.ls.polynomial_fit import PolynomialFit
from stalk.ls.fitting_function import FittingFunction
from stalk.ls.spline_fit import SplineFit


class LsSettings():

    _fit_func: FittingFunction
    _fraction: float
    _sgn: int
    _N: int

    @property
    def fraction(self):
        return self._fraction
    # end def

    @fraction.setter
    def fraction(self, fraction):
        if isinstance(fraction, float) and fraction > 0 and fraction < 0.5:
            self._fraction = fraction
        else:
            raise ValueError("Must provide 0 < fraction < 0.5.")
        # end if
    # end def

    @property
    def N(self):
        return self._N
    # end def

    @N.setter
    def N(self, N):
        if isinstance(N, int) and N > 1:
            self._N = N
        else:
            raise ValueError("Must provide N > 1")
        # end if
    # end def

    @property
    def sgn(self):
        return self._sgn
    # end def

    @sgn.setter
    def sgn(self, sgn):
        if isscalar(sgn) and abs(sgn) == 1:
            self._sgn = int(sgn)
        else:
            raise ValueError('Must provide sgn as 1 or -1')
        # end if
    # end def

    @property
    def fit_func(self):
        return self._fit_func
    # end def

    def __init__(
        self,
        N=200,
        fraction=0.025,
        sgn=1,
        fit_kind=None,
        fit_func=None,
        fit_args={},
    ):
        self.N = N
        self.sgn = sgn
        self.fraction = fraction
        self._set_fit_func(
            fit_func=fit_func,
            fit_kind=fit_kind,
            fit_args=fit_args,
        )
    # end def

    def copy(
        self,
        **ls_overrides
    ):
        ls_args = {
            'fraction': self.fraction,
            'sgn': self.sgn,
            'fit_func': self.fit_func,
            'N': self.N
        }
        ls_args.update(**ls_overrides)
        return LsSettings(**ls_args)
    # end def

    def _set_fit_func(self, fit_func, fit_args, fit_kind):
        # Fit kind (str) takes precedence
        if fit_kind == 'morse':
            fit_func = MorseFit()
        elif fit_kind == 'spline':
            fit_func = SplineFit()
        elif hasattr(fit_kind, "__iter__") and 'pf' in fit_kind:
            fit_func = PolynomialFit(int(fit_kind[2:]))
        elif callable(fit_func):
            # Next, check if fit_func is provided as a function
            fit_func = FittingFunction(fit_func, args=fit_args)
        elif isinstance(fit_func, FittingFunction):
            # Fitting function is good as is
            pass
        else:
            raise TypeError('Fit kind {} not recognized'.format(fit_kind))
        # end fi
        self._fit_func = fit_func
    # end def

    def __str__(self):
        string = 'Fit: {}, N: {}, fraction: {}, sgn: {}'.format(
            self.fit_func.kind,
            self.N,
            self.fraction,
            self.sgn
        )
        return string
    # end def

    # Return true if the settings of self and other are consistent
    def __eq__(
        self,
        other
    ):
        if not isinstance(other, LsSettings):
            return False
        # end if
        result = self.N == other.N
        result &= self.fit_func == other.fit_func
        result &= self.fraction == other.fraction
        result &= self.sgn == other.sgn
        return result
    # end def

# end class
