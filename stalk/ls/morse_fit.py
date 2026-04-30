#!/usr/bin/env python3
'''Class for polynomial fitting for curve minimum'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import dot
from scipy.optimize import minimize

from stalk.ls import FittingFunction
from stalk.ls.morse_result import MorseResult
from stalk.ls.polynomial_fit import PolynomialFit
from stalk.util.util import match_to_tol, morse


class MorseFit(FittingFunction):
    p0 = None

    def __init__(self, p0=None):
        self.p0 = p0
    # end def

    @property
    def kind(self):
        return 'morse'
    # end def

    def _eval_function(self, offsets, values) -> MorseResult:
        if len(offsets) <= 3:
            raise ValueError(f'Fitting of Morse potential to {len(offsets)} points is under-determined. Aborting.')
        # end if
        # Set x0 so that the whole grid is on the positive side
        x_shift = -min(offsets) + (max(offsets) - min(offsets)) * 0.5
        if self.p0 is None:
            pf2 = PolynomialFit(2)._eval_function(offsets + x_shift, values)
            Einf = values[-1]
            a_arg = 2 * (Einf - pf2.y0) / pf2.fit[2]
            if a_arg > 0.0:
                a = a_arg**0.5
            else:
                a = 1.0
            # end if
            p0 = [pf2.x0, a, Einf - pf2.y0, Einf]
        else:
            p0 = self.p0
        # end if

        # Least-squares error
        def cost(p):
            diff = values - morse(p, offsets + x_shift)
            return dot(diff, diff)
        # end def

        fit_res = minimize(cost, x0=p0).x
        x0 = fit_res[0] - x_shift
        y0 = fit_res[3] - fit_res[2]
        fit_res[0] -= x_shift
        res = MorseResult(x0, y0, fit=fit_res)
        return res
    # end def

    def __eq__(self, other):
        if not isinstance(other, MorseFit):
            return False
        # end if
        if self.p0 is None:
            result = other.p0 is None
        else:
            result = match_to_tol(self.p0, other.p0) and match_to_tol(self.x0, other.x0)
        # end if
        return result
    # end def

# end class
