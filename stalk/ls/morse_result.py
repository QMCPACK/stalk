#!/usr/bin/env python3
'''Class for Morse curve minimum and error bars'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from stalk.ls import FittingResult
from stalk.util.util import morse


class MorseResult(FittingResult):

    def get_hessian(self, x, dx=1e-4):
        # For the moment, just do finite-difference
        ddx = (self.get_force(x + dx, dx=dx) - self.get_force(x, dx=dx)) / dx
        return -ddx
    # end def

    def get_force(self, x, dx=1e-4):
        # For the moment, just do finite-difference
        dx = (self.get_values(x + dx) - self.get_values(x)) / dx
        return -dx
    # end def

    def get_values(self, x):
        return morse(self.fit, x)
    # end def

# end class
