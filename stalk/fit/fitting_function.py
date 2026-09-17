#!/usr/bin/env python3
'''Generic class for fitting for curve minimum'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import ndarray, array

from stalk.fit.fitting_result import FittingResult
from stalk.util import get_fraction_error
from stalk.util.function_caller import FunctionCaller
from stalk.util.noise import Noise, NoiseFactory


class FittingFunction(FunctionCaller):
    _result_class = FittingResult

    @property
    def kind(self):
        # NOTE: using this argument does cannot create another FittingFunction
        return str(self.func)
    # end def

    def find_minimum(
        self,
        grid: ndarray,  # grid or offsets
        values: ndarray | None = None,
        errors: ndarray | None = None,
        N: int = 200,
        Gs: ndarray | None = None,
        sgn: int = 1,
        fraction: float = 0.025,
        noise: str | Noise = 'std',
    ):
        offsets, values, errors = self._sanitize_inputs(grid, values, errors)
        result = self._eval_function(offsets * sgn, values)
        # If errors present, resample errorbars; if not, errors default to 0
        if errors is not None:
            x0s, y0s = self.get_distribution(
                offsets,
                values,
                errors,
                N=N,
                Gs=Gs,
                sgn=sgn,
                noise=noise,
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
        grid: ndarray,  # grid or offsets
        values: ndarray | None = None,
        errors: ndarray | None = None,
        N: int = 200,
        Gs: ndarray | None = None,
        sgn: int = 1,
        noise: str | Noise = 'std',
    ):
        offsets, values, errors = self._sanitize_inputs(grid, values, errors)
        if Gs is None:
            if isinstance(N, int) and N > 0:
                Gs = NoiseFactory.create(noise).generate(N, len(errors))
            else:
                raise ValueError("Must provide either N > 0 or an array of G displacements")
            # end if
        elif Gs.shape[1] != len(errors):
            raise AssertionError("Must provide Gs that are consistent with valid data.")
        # end if
        x0_distribution = []
        y0_distribution = []
        fit_distribution = []
        for G in Gs:
            values_this = sgn * values + errors * G
            result_this = self._eval_function(offsets, values_this)
            x0_distribution.append(result_this.x0)
            y0_distribution.append(result_this.y0)
            fit_distribution.append(result_this.fit)
        # end for
        return array(x0_distribution), array(y0_distribution)
    # end def

    def get_x0_distribution(
        self,
        *args,  # grid, values, errors
        **kwargs,  # N=200, Gs=None, sgn=1, noise='std'
    ):
        return self.get_distribution(*args, **kwargs)[0]
    # end def

    def get_y0_distribution(
        self,
        *args,  # grid, values, errors
        **kwargs,  # N=200, Gs=None, sgn=1, noise='std'
    ):
        return self.get_distribution(*args, **kwargs)[1]
    # end def

    def _sanitize_inputs(
        self,
        grid: ndarray,  # grid or offsets
        values: ndarray | None = None,
        errors: ndarray | None = None,
    ):
        if hasattr(grid, 'valid_args'):
            # If grid is a LineSearchGrid, use its valid arguments
            offsets, values, errors = grid.valid_args
        elif values is None:
            raise TypeError("Must provide either a grid or offsets and values.")
        else:
            offsets = grid
            if len(offsets) != len(values):
                raise ValueError("Offsets and values must be of the same length.")
            # end if
        # end if
        if errors is not None and len(offsets) != len(errors):
            raise ValueError("Offsets, values and errors must be of the same length.")
        # end if
        return offsets, values, errors
    # end def

    def _eval_function(self, offsets, values) -> FittingResult:
        x0, y0, fit = self.func(offsets, values, **self.args)
        return self._result_class(x0, y0, fit=fit)
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
