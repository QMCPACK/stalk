#!/usr/bin/env python3
'''Class for polynomial fitting for curve minimum'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import argmin, polyder, polyfit, polyval, roots, where
from stalk.ls import FittingFunction
from stalk.ls.polynomial_result import PolynomialResult


class PolynomialFit(FittingFunction):
    _n = None

    def __init__(self, n):
        self.n = n
    # end def

    @property
    def n(self):
        return self._n
    # end def

    @property
    def kind(self):
        return f'pf{self.n}'
    # end def

    @n.setter
    def n(self, n):
        if not isinstance(n, int):
            raise TypeError(f'Polynomial degree n must be integer, provided: {n}')
        # end if
        if n < 2:
            raise ValueError(f'Polynomial degree n must be higher than 1, provided: {n}')
        # end if
        self._n = n
    # end def

    def _eval_function(self, offsets, values) -> PolynomialResult:
        if len(offsets) <= self.n:
            raise ValueError(f'Fitting of {self.n} degree polynomial to {len(offsets)} points is under-determined. Aborting.')
        # end if
        pf = polyfit(offsets, values, self.n)
        pfd = polyder(pf)
        r = roots(pfd)
        d = polyval(polyder(pfd), r)
        # filter real minima (maxima with sgn < 0)
        x_mins = r[where((r.imag == 0) & (d > 0))].real
        if len(x_mins) > 0:
            y_mins = polyval(pf, x_mins)
            imin = argmin(abs(x_mins))
        else:
            if abs(offsets).max() > 1e-4:
                # Only warn if using nonzero grid
                warnings.warn('The fit minimum not found inside grid but at the boundary!')
            # end if
            x_mins = [min(offsets), max(offsets)]
            y_mins = polyval(pf, x_mins)
            imin = argmin(y_mins)  # pick the lowest/highest energy
        # end if
        y0 = y_mins[imin]
        x0 = x_mins[imin]
        res = PolynomialResult(x0, y0, fit=pf)
        return res
    # end def

    def __eq__(self, other):
        if not isinstance(other, PolynomialFit):
            return False
        # end if
        result = self.n == other.n
        return result
    # end def

# end class
