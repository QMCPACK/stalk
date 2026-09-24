#!/usr/bin/env python3
'''Class for line-search along direction in abstract context'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path
import warnings
from numpy import linspace, isscalar, array
from matplotlib import pyplot as plt

from stalk.io.cacheable import Cacheable
from stalk.fit.fitting_result import FittingResult
from stalk.io.txt_data import TxtData
from stalk.ls.linesearch_grid import LineSearchGrid
from stalk.ls.ls_settings import LsSettings
from stalk.params.parameter_set import ParameterSet
from stalk.pes.pes_function import PesFunction
from stalk.util.util import FF, FU


class LineSearchBase(Cacheable, LineSearchGrid[ParameterSet]):
    _settings: LsSettings = None
    _sigma = 0.0  # Target errorbar
    _fit_res: FittingResult = None

    def __init__(
        self,
        offsets=None,
        values=None,
        errors=None,
        fraction=0.025,
        sigma=0.0,
        sgn=1,
        fit_kind='pf3',
        fit_func=None,
        fit_args={},
        N=200,
        store=True,
        noisy=True,
        # try to load after init
        path: str | None = None,
        check_offsets=False,
    ):
        Cacheable.__init__(
            self,
            data=TxtData('ls.out'),
        )
        LineSearchGrid.__init__(self, offsets)
        self.sigma = sigma
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
            self.search(store=store, noisy=noisy)
        # end if
        # Try to load from disk
        self.try_load_result(path, check_offsets=check_offsets)
    # end def

    @property
    def fit_res(self):
        return self._fit_res
    # end def

    @fit_res.setter
    def fit_res(self, res: FittingResult):
        if res is not None and not isinstance(res, FittingResult):
            raise TypeError("fit_res must be a FittingResult instance or None.")
        # end if
        self._fit_res = res
    # end def

    @property
    def settings(self) -> LsSettings:
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

    @property
    def sigma(self):
        return self._sigma
    # end def

    @sigma.setter
    def sigma(self, sigma):
        if isscalar(sigma) and sigma >= 0.0:
            self._sigma = sigma
        else:
            raise ValueError("Sigma must be >= 0.0")
        # end if
    # end def

    def evaluate(
        self,
        pes: PesFunction,
        path=None,
        var_eff_map=None,
        interactive=False,
        dep_jobs=[],
        add_sigma=False,
        warn_limit=2.0,
        **kwargs,  # etc.
    ):
        '''Evaluate the PES on the line-search grid using an evaluation function.'''
        if not self.shifted:
            raise AssertionError('The line-search grid must be generated before evaluation!')
        # end if
        # Collect enabled structures, assign sigmas
        structures = self.collect_enabled()
        pes.evaluate_all(
            structures,
            path=path,
            var_eff_map=var_eff_map,
            interactive=interactive,
            dep_jobs=dep_jobs,
            add_sigma=add_sigma,
            warn_limit=warn_limit,
            **kwargs
        )
        if self.evaluated:
            self.search()
        else:
            print(f'{repr(self)} missing results for the following structures:')
            for point in self.grid:
                if point.valid:
                    continue
                else:
                    print(f'  {point.offset}')
                # end if
            # end for
        # end if
    # end def

    def search(
        self,
        store=True,
        noisy=True,
        **ls_overrides
    ):
        settings = self.settings.copy(**ls_overrides)
        if noisy:
            res = settings.fit_func.find_minimum(
                self,
                sgn=settings.sgn,
                fraction=settings.fraction,
                N=settings.N
            )
        else:
            res = settings.fit_func.find_minimum(
                self.valid_offsets,
                self.valid_values,
                # no errors
                sgn=settings.sgn
            )
        # end if
        if store:
            self.fit_res = res
        # end if
        return res
    # end def

    def reset_search(self, fit_res: FittingResult = None) -> None:
        self.fit_res = fit_res
    # end def

    def finalize(self) -> None:
        if self.evaluated:
            self.search(store=True, noisy=True)
        else:
            warnings.warn("Cannot finalize without valid data.")
        # end if
    # end def

    def save_result(self, path: str | Path, overwrite: bool = True) -> None:
        if path is None:
            return
        # end if
        path = Path(path)
        data = array([self.offsets, self.values, self.errors]).T
        self.save(path=path, data=data, overwrite=overwrite)
        if self.fit_res is not None:
            # Save fitting result if available
            self.fit_res.save_result(path, overwrite=overwrite)
        # end if
    # end def

    def try_load_result(self, path: str | Path | None, check_offsets: bool = False) -> None:
        if path is None:
            return
        # end if
        path = Path(path)
        data = self.load(path, 'data')
        if data is not None:
            if not check_offsets or len(self) == len(data[0]) and all(self.offsets == data[0]):
                self.grid = data[:, 0]
                self.values = data[:, 1]
                self.errors = data[:, 2]
                print(f'Loaded line-search data from {path}.')
                fit_res = self.settings.fit_func._result_class()
                fit_res.load_result(path)
                if fit_res.analyzed:
                    self.fit_res = fit_res
                    print(f'Loaded fitting result from {path}.')
                # end if
            else:
                msg = 'Offsets of the loaded line-search do not match the current offsets. '
                msg += f'Loaded offsets: {data[0]}, current offsets: {self.offsets}'
                raise ValueError(msg)
            # end if
        # end if
    # end def

    def _make_offsets_R(self, R: float, M: int):
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

    # Override to associate sigma with the structures before evaluation
    def collect_enabled(self) -> list[ParameterSet]:
        structures = super().collect_enabled()
        for structure in structures:
            structure.sigma = self.sigma
        # end for
        return structures
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
