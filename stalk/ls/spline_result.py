#!/usr/bin/env python3
'''Class for spline curve minimum and error bars'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from scipy.interpolate import CubicSpline

from stalk.ls import FittingResult


class SplineResult(FittingResult):
    fit: CubicSpline

    def get_hessian(self, x):
        h = self.fit.derivative(nu=2)(x)
        return h
    # end def

    def get_force(self, x):
        return -self.fit.derivative(nu=1)(x)
    # end def

    def get_values(self, x):
        return self.fit(x)
    # end def

# end class
