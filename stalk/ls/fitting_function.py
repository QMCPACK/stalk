#!/usr/bin/env python3
'''Generic class for fitting for curve minimum'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import ndarray, random, array

from stalk.ls.fitting_result import FittingResult
from stalk.ls.linesearch_grid import LineSearchGrid
from stalk.util import get_fraction_error


class FittingFunction():
    func = None
    args = {}

    def __init__(
        self,
        func=None,
        args={}
    ):
        if callable(func):
            self.func = func
            self.args = args
        else:
            raise ValueError("The fitting function must be callable.")
        # end if
    # end def

    @property
    def kind(self):
        # NOTE: using this argument does cannot create another FittingFunction
        return str(self.func)
    # end def

    def find_minimum(
        self,
        grid,
        sgn=1,
    ):
        if not isinstance(grid, LineSearchGrid):
            raise ValueError("Fitting function input must be inherited from LineSearchGrid.")
        # end if
        return self._eval_function(grid.valid_offsets, grid.valid_values * sgn)
    # end def

    def find_noisy_minimum(
        self,
        grid,
        sgn=1,
        N=200,
        Gs=None,
        fraction=0.025
    ):
        if not isinstance(grid, LineSearchGrid):
            raise ValueError("Fitting function input must be inherited from LineSearchGrid.")
        # end if
        result = self._eval_function(grid.valid_offsets, grid.valid_values * sgn)
        # If errors present, resample errorbars; if not, errors default to 0
        if grid.noisy:
            x0s, y0s = self.get_distribution(
                grid,
                sgn=sgn,
                N=N,
                Gs=Gs
            )
            result.x0_err = get_fraction_error(x0s - result.x0, fraction=fraction)[1]
            result.y0_err = get_fraction_error(y0s - result.y0, fraction=fraction)[1]
        # end if
        result.fraction = fraction
        return result
    # end def

    # Return a random resampled distribution of x0, y0 results based on the input grid
    # (offsets, values, errors) realized on the fitting function
    def get_distribution(
        self,
        grid,
        sgn=1,
        N=200,
        Gs: ndarray = None,
    ):
        if not isinstance(grid, LineSearchGrid):
            raise ValueError("Fitting function input must be inherited from LineSearchGrid.")
        # end if
        if Gs is None:
            if isinstance(N, int) and N > 0:
                Gs = random.randn(N, len(grid.valid_errors))
            else:
                raise ValueError("Must provide either N > 0 or an array of G displacements")
            # end if
        elif Gs.shape[1] != len(grid.valid_errors):
            raise AssertionError("Must provide Gs that are consistent with valid data.")
        # end if
        x0_distribution = []
        y0_distribution = []
        fit_distribution = []
        for G in Gs:
            values = sgn * grid.valid_values + grid.valid_errors * G
            result_this = self._eval_function(grid.valid_offsets, values)
            x0_distribution.append(result_this.x0)
            y0_distribution.append(result_this.y0)
            fit_distribution.append(result_this.fit)
        # end for
        return array(x0_distribution), array(y0_distribution)
    # end def

    def get_x0_distribution(
        self,
        grid,
        **kwargs,  # sng=1, N=200, Gs=None
    ):
        return self.get_distribution(grid, **kwargs)[0]
    # end def

    def get_y0_distribution(
        self,
        grid,
        **kwargs  # sng=1, N=200, Gs=None
    ):
        return self.get_distribution(grid, **kwargs)[1]
    # end def

    def _eval_function(self, offsets, values) -> FittingResult:
        x0, y0, fit = self.func(offsets, values, **self.args)
        return FittingResult(x0, y0, fit=fit)
    # end def

    def __eq__(self, other):
        if not isinstance(other, FittingFunction):
            return False
        # end if
        result = self.func is other.func
        for key, val in self.args.items():
            result &= key in other.args and val == other.args[key]
        # end for
        return result
    # end def

# end class
