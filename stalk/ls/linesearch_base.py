#!/usr/bin/env python3
'''Class for line-search along direction in abstract context'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import linspace
from matplotlib import pyplot as plt
from stalk.ls.fitting_result import FittingResult
from stalk.ls.linesearch_grid import LineSearchGrid
from stalk.ls.ls_settings import LsSettings
from stalk.util.util import FF, FU


class LineSearchBase(LineSearchGrid):
    _settings: LsSettings
    fit_res: FittingResult

    def __init__(
        self,
        offsets=None,
        values=None,
        errors=None,
        fraction=0.025,
        sgn=1,
        fit_kind='pf3',
        fit_func=None,
        fit_args={},
        N=200
    ):
        LineSearchGrid.__init__(self, offsets)
        self._settings = LsSettings(
            fraction=fraction,
            sgn=sgn,
            N=N,
            fit_func=fit_func,
            fit_args=fit_args,
            fit_kind=fit_kind,
        )
        self.fit_res = None
        if values is not None:
            self.values = values
            if errors is not None:
                self.errors = errors
            # end if
            self._search_and_store()
        # end if
    # end def

    @property
    def settings(self):
        return self._settings
    # end def

    @property
    def x0(self):
        return None if self.fit_res is None else self.fit_res.x0
    # end def

    @property
    def x0_err(self):
        return None if self.fit_res is None else self.fit_res.x0_err
    # end def

    @property
    def y0(self):
        return None if self.fit_res is None else self.fit_res.y0
    # end def

    @property
    def y0_err(self):
        return None if self.fit_res is None else self.fit_res.y0_err
    # end def

    def search_with_error(
        self,
        **ls_overrides
    ):
        settings = self.settings.copy(**ls_overrides)
        res = settings.fit_func.find_noisy_minimum(
            self,
            sgn=settings.sgn,
            fraction=settings.fraction,
            N=settings.N
        )
        return res
    # end def

    def search(
        self,
        **ls_overrides
    ):
        settings = self.settings.copy(**ls_overrides)
        res = settings.fit_func.find_minimum(
            self,
            sgn=settings.sgn
        )
        return res
    # end def

    def _search_and_store(self):
        """Perform line-search with the preset values and settings, saving the result to self."""
        self.fit_res = self.search_with_error()
    # end def

    def reset_search(self, x0=0.0, y0=0.0):
        self.fit_res.x0 = x0
        self.fit_res.x0_err = 0.0
        self.fit_res.y0 = y0
        self.fit_res.y0_err = 0.0
        return self
    # end def

    def _make_offsets_R(self, R, M):
        if R < 1e-6:
            raise ValueError("R must be larger than 1e-6")
        # end if
        offsets = linspace(-R, R, M)
        return offsets
    # end def

    def plot(
        self,
        ax=None,
        color='tab:blue',
        target=None,
        **kwargs
    ):
        if not self.valid:
            warnings.warn("Cannot plot without valid data.")
            return
        # end if
        if ax is None:
            f, ax = self._create_plot(**kwargs)
        # end if
        if target is None:
            target = self.fit_res
        # end if
        LineSearchGrid.plot(self, ax=ax, color=color, **kwargs)
        if self.fit_res is not None:
            ax.errorbar(
                target.x0,
                target.y0,
                target.y0_err,
                xerr=target.x0_err,
                linestyle='none',
                marker='x',
                color=color,
                label='Fitted minimum'
            )
            xgrid = self._get_plot_grid()
            ygrid = target.get_values(xgrid)
            ax.plot(
                xgrid,
                ygrid,
                linestyle='--',
                color=color,
                label='Fitted curve'
            )
        # end if
        plt.tight_layout()
    # end def

    def _get_plot_grid(self, fraction=0.1):
        w = (self.offsets.max() - self.offsets.min()) * fraction
        grid = linspace(self.offsets.min() - w, self.offsets.max() + w, 201)
        return grid
    # end def

    def __str__(self):
        string = LineSearchGrid.__str__(self)
        string += '\n  ' + str(self.settings)
        if self.x0 is None:
            string += '\n  x0: not set'
        else:
            string += '\n  x0: ' + FF.format(self.x0)
            if self.fit_res.x0_err > 0:
                string += FU.format(self.x0_err)
            # end if
        # end if
        if self.y0 is None:
            string += '\n  y0: not set'
        else:
            string += '\n  y0: ' + FF.format(self.y0)
            if self.y0_err > 0:
                string += FU.format(self.y0_err)
            # end if
        # end if
        return string
    # end def


# end class
