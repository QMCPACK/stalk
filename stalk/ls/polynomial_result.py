#!/usr/bin/env python3
'''Class for polynomial curve minimum and error bars'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import polyder, polyval

from stalk.ls import FittingResult


class PolynomialResult(FittingResult):

    def get_hessian(self, x):
        pfh = polyder(polyder(self.fit))
        return polyval(pfh, x)
    # end def

    def get_force(self, x):
        pfd = -polyder(self.fit)
        return polyval(pfd, x)
    # end def

    def get_values(self, x):
        return polyval(self.fit, x)
    # end def

# end class
