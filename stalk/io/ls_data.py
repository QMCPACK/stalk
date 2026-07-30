#!/usr/bin/env python3
'''Class for saving/loading line-search data and fitting results'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array

from stalk.io.txt_data import TxtData
from stalk.ls.fitting_result import FittingResult
from stalk.ls.linesearch_base import LineSearchBase


class LineSearchData():
    # Line-search data file for saving/loading grid, values and errors
    ls_file: TxtData = None
    # Data file for saving/loading the result coordinates
    result_file: TxtData = None
    # Data file for saving/loading the fitted parameters
    fit_file: TxtData = None

    def __init__(
        self,
        label: str = 'ls',
    ):
        self.label = label
        self.ls_file = TxtData(f'{label}.dat')
        self.result_file = TxtData(f'{label}_res.dat')
        self.fit_file = TxtData(f'{label}_fit.dat')
    # end def

    def save(
        self,
        grid: LineSearchBase,
        path: str,
        overwrite: bool = True
    ) -> None:
        data = array([grid.offsets, grid.values, grid.errors]).T
        self.ls_file.save_result(path, data, overwrite=overwrite)
        if grid.fit_res is not None:
            # Save fitting result if available
            fit_res = grid.fit_res
            x0 = array([[fit_res.x0, fit_res.y0], [fit_res.x0_err, fit_res.y0_err]])
            self.result_file.save_result(path, x0, overwrite=overwrite)
            self.fit_file.save_result(path, fit_res.fit, overwrite=overwrite)
        # end if
    # end def

    def load(
        self,
        path: str,
    ) -> LineSearchBase:
        '''Load line-search data from file and return a LineSearchBase instance of desired type.'''
        data = self.ls_file.load_result(path, None)
        if data is None:
            return None
        # end if
        result_grid = LineSearchBase(
            offsets=data[:, 0],
            values=data[:, 1],
            errors=data[:, 2]
        )
        x0 = self.result_file.load_result(path, None)
        fit = self.fit_file.load_result(path, None)
        if x0 is not None and fit is not None:
            fit_res = FittingResult(
                x0=x0[0, 0],  # x0
                y0=x0[0, 1],  # y0
                x0_err=x0[1, 0],  # x0_err
                y0_err=x0[1, 1],  # y0_err
                fit=fit,
            )
            result_grid.fit_res = fit_res
        # end if
        return result_grid
    # end def

# end class
