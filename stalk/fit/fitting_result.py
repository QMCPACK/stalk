#!/usr/bin/env python3
'''Generic class for curve minimum and error bars'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


from numpy import ndarray
from pathlib import Path

from stalk.io.cacheable import Cacheable
from stalk.io.txt_data import TxtData


class FittingResult(Cacheable):
    fraction = None
    x0 = None
    x0_err = 0.0
    y0 = None
    y0_err = 0.0
    fit = None
    # Flag indicating if the result is at the boundary
    boundary = False

    def __init__(
        self,
        x0: float = None,
        y0: float = None,
        x0_err: float = 0.0,
        y0_err: float = 0.0,
        fit: ndarray = None,
        fraction: float = 0.025,
    ):
        self.fraction = fraction
        self.x0 = x0
        self.y0 = y0
        self.x0_err = x0_err
        self.y0_err = y0_err
        self.fit = fit
        Cacheable.__init__(
            self,
            x0=TxtData('x0.out'),
            y0=TxtData('y0.out'),
            fit=TxtData('fit.out'),
        )
    # end def

    @property
    def analyzed(self):
        return self.x0 is not None
    # end def

    def save_result(self, path: str | Path, overwrite: bool = True) -> None:
        """Save the fitting result to disk."""
        if self.analyzed:
            self.save(
                path=path,
                overwrite=overwrite,
                x0=[self.x0, self.x0_err],
                y0=[self.y0, self.y0_err],
                fit=self.fit,
            )
        else:
            raise AssertionError("Fitting result is not analyzed, cannot save.")
        # end if
    # end def

    def load_result(self, path: str | Path) -> None:
        """Load the fitting result from disk."""
        x0res = self.load(path, 'x0')
        self.x0 = x0res[0]
        self.x0_err = x0res[1]
        y0res = self.load(path, 'y0')
        self.y0 = y0res[0]
        self.y0_err = y0res[1]
        self.fit = self.load(path, 'fit')
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
