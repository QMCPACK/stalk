#!/usr/bin/env python3
'''ParallelLineSearch class for simultaneous linesearches along conjugate directions'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import ndarray, array
from textwrap import indent

from stalk.io.ls_data import LineSearchData
from stalk.params.pes_function import NotEvaluatedException, PesFunction
from stalk.util import get_fraction_error
from stalk.params import ParameterSet
from stalk.params import ParameterHessian
from stalk.ls import LineSearch


class ParallelLineSearch():
    ls_type = LineSearch
    _ls_list: list[LineSearch] = []  # list of line-search objects
    _hessian = None  # hessian object
    _structure = None  # eqm structure
    _structure_next = None  # next structure
    _path = None

    def __init__(
        self,
        # PLS arguments
        path='pls',
        hessian=None,
        structure=None,
        windows=None,
        window_frac=0.25,
        noises=None,
        # LineSearch args
        **ls_args
        # M=7, fit_kind='pf3', fit_func=None, fit_args={}, N=200, Gs=None, fraction=0.025
    ):
        self.path = path
        if structure is not None:
            self.structure = structure
        # end if
        if hessian is not None:
            self.hessian = hessian
        # end if
        if self.setup:
            self.initialize(
                windows,
                noises,
                window_frac,
                **ls_args
            )
        # end if
        # Make sure that LS is solved after successful loading
        if self.evaluated and self.structure_next is None:
            self._solve_ls()
        # end if
    # end def

    @property
    def path(self):
        return self._path
    # end def

    @path.setter
    def path(self, path):
        if isinstance(path, str):
            self._path = path
        else:
            raise ValueError('path must be str')
        # end if
    # end def

    # Return True if the parallel line-search has starting structure and Hessian
    @property
    def setup(self):
        return self.structure is not None and self.hessian is not None
    # end def

    # Return True if the line-search structures have been shifted
    @property
    def shifted(self):
        return len(self) > 0 and all([ls.shifted for ls in self.ls_list])
    # end def

    # Return a list of all line-searches
    @property
    def ls_list(self):
        return self._ls_list
    # end def

    @property
    def D(self):
        if self.hessian is None:
            return 0
        else:
            return len(self.hessian)
        # end if
    # ed def

    @property
    def evaluated(self):
        return len(self) > 0 and all([ls.evaluated for ls in self.ls_list])
    # end def

    @property
    def hessian(self):
        return self._hessian
    # end def

    @hessian.setter
    def hessian(self, hessian):
        if isinstance(hessian, ndarray):
            hessian = ParameterHessian(hessian=hessian)
        elif not isinstance(hessian, ParameterHessian):
            raise ValueError('Hessian matrix is not supported')
        # end if
        if self._hessian is not None:
            pass  # TODO: check for constistency
        # end if
        self._hessian = hessian
        if self.structure is None:
            self.structure = hessian.structure
        # end if
        # TODO: propagate Hessian information to ls_list?
    # end def

    @property
    def structure(self):
        return self._structure
    # end def

    @structure.setter
    def structure(self, structure):
        if not isinstance(structure, ParameterSet):
            raise TypeError("Structure must be inherited from ParameterSet clas")
        # end if
        self._structure = structure.copy(label='eqm')
        # Upon change, reset line-searches according to old windows/noises, if present
        if self.shifted:
            windows = self.windows
            noises = self.noises
            self._reset_ls_list(windows, noises)
        # end if
    # end def

    @property
    def structure_next(self) -> ParameterSet:
        return self._structure_next
    # end def

    @property
    def Lambdas(self) -> ndarray:
        if self.hessian is None:
            return array([])
        else:
            return array(self.hessian.enabled_lambdas)
        # end if
    # end def

    @property
    def windows(self) -> list[float]:
        result = []
        for ls in self.ls_list:
            if isinstance(ls, LineSearch):
                window = ls.W_max
            else:
                window = None
            # end if
            result.append(window)
        # end if
        return result
    # end def

    @property
    def noises(self) -> list[float]:
        return [ls.sigma for ls in self.ls_list]
    # end def

    @property
    def noises_min(self) -> float:
        return array([ls.sigma for ls in self.ls_list]).min()
    # end def

    @property
    def D_list(self) -> list[int]:
        return [d for d in range(len(self.hessian)) if self.hessian.enabled[d]]
    # end def

    def initialize(
        self,
        windows=None,
        noises=None,
        window_frac=None,
        **ls_args
        # M=7, fit_kind='pf3', fit_func=None, fit_args={}, N=200, Gs=None, fraction=0.025
    ) -> None:
        if windows is None:
            windows = abs(self.Lambdas)**0.5 * window_frac
        # end if
        noises = self._for_all_ls(noises, default=0.0)
        self._reset_ls_list(windows, noises, **ls_args)
    # end def

    def _reset_ls_list(
        self,
        windows,
        noises,
        M=7,
        **ls_args,
        # fit_kind='pf3', fit_func=None, fit_args={}, N=200, Gs=None, fraction=0.025
    ) -> None:
        M = self._for_all_ls(M)
        ls_list = []
        for d, window, noise in zip(self.D_list, windows, noises):
            # Only add if enabled by the Hessian
            if self.hessian.enabled[d]:
                # Try to load from disk
                ls_load = LineSearchData(label=f'ls{d}').load(path=self.path)
                # Create new line-search object
                ls = self.ls_type(
                    structure=self.structure,
                    hessian=self.hessian,
                    d=d,
                    sigma=noise,
                    W=window,
                    M=M[d],
                    **ls_args
                )
                if ls_load is not None and len(ls) == len(ls_load):
                    if not all(ls.offsets == ls_load.offsets):
                        raise ValueError('Offsets of the loaded line-search do not match the current offsets')
                    # end if
                    print(f'{self.path}ls{d}: Line-search data loaded from disk.')
                    ls.values = ls_load.values
                    ls.errors = ls_load.errors
                    ls.fit_res = ls_load.fit_res
                # end if
                ls_list.append(ls)
            # end if
        # end for
        self._ls_list = ls_list
        # Reset next structure if re-initialized
        self._structure_next = None
    # end def

    def evaluate(
        self,
        pes: PesFunction,
        add_sigma=False,
        interactive=False,
        dep_jobs=[],
        warn_limit=2.0,
        var_eff_map=None,
    ) -> None:
        if not self.shifted:
            raise AssertionError("Must have shifted structures first!")
        # end if
        structures, sigmas = self._collect_enabled()
        pes.evaluate_all(
            structures,
            sigmas=sigmas,
            path=self.path,
            add_sigma=add_sigma,
            interactive=interactive,
            dep_jobs=dep_jobs,
            var_eff_map=var_eff_map,
            warn_limit=warn_limit,
        )
        if not all([s.value is not None and s.enabled for s in structures]):
            print('Cannot solve the line-searches, as not all structures were successfully evaluated.')
            return
        # end if
        self._solve_ls()
    # end def

    def evaluate_eqm(
        self,
        pes: PesFunction,
        add_sigma=False,
        interactive=False,
        dep_jobs=[],
        var_eff_map=None,
    ):
        pes.evaluate(
            self.structure,
            sigma=self.noises_min,
            path=self.path,
            add_sigma=add_sigma,
            interactive=interactive,
            dep_jobs=dep_jobs,
            var_eff_map=var_eff_map,
        )
    # end def

    def _collect_enabled(self) -> tuple[list[ParameterSet], list[float]]:
        structures = []
        sigmas = []
        sigma_eqm = self.noises_min
        for ls in self.ls_list:
            for structure in ls.grid:
                structures += [structure]
                if structure.is_eqm:
                    sigmas += [sigma_eqm]
                else:
                    sigmas += [ls.sigma]
                # end if
            # end for
        # end for
        return structures, sigmas
    # end def

    def _solve_ls(self):
        # Set the eqm energy and solve the line-searches
        for ls in self.ls_list:
            ls._search_and_store()
            eqm = ls.get(0.0)
            if eqm is not None:
                self.structure.value = eqm.value
                self.structure.error = eqm.error
            # end if
        # end for
        # Calculate next params
        params_next, params_next_err = self.calculate_next_params()  # **kwargs
        self._structure_next = self.structure.copy(
            params=params_next,
            params_err=params_next_err
        )
    # end def

    @property
    def noisy(self):
        return any([ls.noisy for ls in self.ls_list])
    # end def

    @property
    def params(self):
        if self.structure is not None:
            return self.structure.params
        # end if
    # end def

    @property
    def params_err(self):
        if self.structure is not None:
            return self.structure.params_err
        # end if
    # end def

    def calculate_next_params(
        self,
        N=200,
        Gs=None,
        fraction=0.025
    ):
        # deterministic
        params = self.params
        shifts = self.shifts
        params_next = self._calculate_params_next(params, shifts)
        # stochastic
        if self.noisy:
            x0s = []
            for ls in self.ls_list:
                x0s.append(ls.settings.fit_func.get_x0_distribution(ls, N=N, Gs=Gs))
            # end if
            x0s = array(x0s).T
            dparams = []
            for shifts_this in x0s:
                dparams.append(
                    self._calculate_params_next(
                        params,
                        shifts_this
                    ) - params_next
                )
            # end for
            dparams = array(dparams).T
            params_next_err = array(
                [get_fraction_error(p, fraction=fraction)[1] for p in dparams]
            )
        else:
            params_next_err = array(self.D * [0.0])
        # end if
        return params_next, params_next_err
    # end def

    def ls(self, i) -> LineSearch:
        if i < 0 or i >= len(self.ls_list):
            raise ValueError("Must choose line-search between 0 and " + str(len(self.ls_list)))
        # end if
        return self.ls_list[i]
    # end def

    def _calculate_params_next(self, params, shifts):
        return params + shifts @ self.hessian.enabled_directions
    # end def

    @property
    def shifts(self):
        return array([ls.x0 for ls in self.ls_list])
    # end def

    def copy(
        self,
        path,
        structure=None,
        hessian=None,
        windows=None,
        noises=None,
        M=None,
    ):
        structure = structure if structure is not None else self.structure
        hessian = hessian if hessian is not None else self.hessian
        windows = windows if windows is not None else self.windows
        noises = noises if noises is not None else self.noises
        pls_args = {}
        if M is not None:
            pls_args['M'] = M
        # end if
        copy_pls = ParallelLineSearch(
            path=path,
            structure=structure,
            hessian=hessian,
            windows=windows,
            noises=noises,
            **pls_args,
        )
        for ls, ls_new in zip(self.ls_list, copy_pls.ls_list):
            ls_new._settings = ls._settings
        # end for
        return copy_pls
    # end def

    def propagate(
        self,
        pes: PesFunction,
        next_path=None,
        overwrite=True,
        add_sigma=False,
        interactive=False,
        **kwargs  # dep_jobs=[], var_eff_map=None
    ):
        if not self.evaluated:
            self.evaluate(pes=pes, add_sigma=add_sigma, interactive=interactive, **kwargs)
        # end if
        if not self.evaluated:
            raise NotEvaluatedException("Cannot propagate, as not all line-searches were successfully evaluated.")
        # end if
        next_path = next_path if next_path is not None else self.path + '_next/'
        # Write to disk
        for ls in self.ls_list:
            LineSearchData(label=f'ls{ls.d}').save(ls, path=self.path, overwrite=overwrite)
        # end if
        pls_next = self.copy(
            next_path,
            structure=self.structure_next
        )
        return pls_next
    # end def

    def plot(
        self,
        **kwargs  # TODO: list kwargs
    ):
        for ls in self.ls_list:
            ls.plot(**kwargs)
        # end for
    # end def

    def _for_all_ls(self, data, default=None):
        if isinstance(data, list):
            if len(data) == len(self.hessian):
                return data
            elif len(data) == len(self.D_list):
                result = []
                for d in range(len(self.hessian)):
                    if d in self.D_list:
                        result.append(data.pop(0))
                    else:
                        result.append(default)
                    # end if
                # end for
                return result
            else:
                raise ValueError("Data must be a list of length equal to the number of line-searches or enabled directions")
            # end if
        elif data is None:
            return len(self.hessian) * [default]
        else:
            return len(self.hessian) * [data]
        # end if
    # end def

    def __str__(self):
        string = self.__class__.__name__
        if self.ls_list is None:
            string += '\n  Line-searches: None'
        else:
            string += '\n  Line-searches:\n'
            string += indent('\n'.join([str(ls) for ls in self.ls_list]), '    ')
        # end if
        # TODO
        return string
    # end def

    def __len__(self):
        return len(self.ls_list)
    # end def

# end class
