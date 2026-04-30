#!/usr/bin/env python3
'''Class for polynomial fitting for curve minimum'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import argmin, where
from scipy.interpolate import CubicSpline

from stalk.ls import FittingFunction
from stalk.ls.spline_result import SplineResult


class SplineFit(FittingFunction):

    def __init__(self):
        pass
    # end def

    @property
    def kind(self):
        return 'spline'
    # end def

    def _eval_function(self, offsets, values) -> SplineResult:
        if len(offsets) <= 3:
            raise ValueError(f'Fitting of a spline to {len(offsets)} points is under-determined. Aborting.')
        # end if
        pf = CubicSpline(offsets, values)
        pfd = pf.derivative()
        pfdd = pfd.derivative()
        r = pfd.roots()
        x_mins = r[where(pfdd(r) > 0)]
        if len(x_mins) > 0:
            y_mins = pf(x_mins)
            imin = argmin(abs(x_mins))
        else:
            warnings.warn('The fit minimum not found inside grid but at the boundary!')
            x_mins = [min(offsets), max(offsets)]
            y_mins = pf(x_mins)
            imin = argmin(y_mins)  # pick the lowest/highest energy
        # end if
        y0 = y_mins[imin]
        x0 = x_mins[imin]
        res = SplineResult(x0, y0, fit=pf)
        return res
    # end def

    def __eq__(self, other):
        return isinstance(other, SplineFit)
    # end def

# end class
