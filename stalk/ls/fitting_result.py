#!/usr/bin/env python3
'''Generic class for curve minimum and error bars'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


class FittingResult():
    # TODO: this implementation is a stub

    fraction = None
    x0 = None
    x0_err = 0.0
    y0 = None
    y0_err = 0.0
    fit = None

    def __init__(
        self,
        x0,
        y0,
        x0_err=0.0,
        y0_err=0.0,
        fit=None,
        fraction=0.025
    ):
        self.fraction = fraction
        self.x0 = x0
        self.y0 = y0
        self.x0_err = x0_err
        self.y0_err = y0_err
        self.fit = fit
    # end def

    @property
    def analyzed(self):
        return self.x0 is not None
    # end def

    def get_hessian(self, x):
        raise NotImplementedError("The Hessian is not implemented for generic class")
    # end def

    def get_force(self, x):
        raise NotImplementedError("The force is not implemented for generic class")
    # end def

    def get_values(self, x):
        raise NotImplementedError("The evaluation is not implemented for generic class")
    # end def

# end class
