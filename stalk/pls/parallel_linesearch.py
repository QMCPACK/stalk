#!/usr/bin/env python3
'''ParallelLineSearch class for simultaneous linesearches along conjugate directions'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import ndarray, array
from textwrap import indent
from dill import dumps, loads
from os import makedirs, path

from stalk.params.pes_function import PesFunction
from stalk.util import get_fraction_error
from stalk.params import ParameterSet
from stalk.params import ParameterHessian
from stalk.ls import LineSearch
from stalk.util.util import directorize


class ParallelLineSearch():
    ls_type = LineSearch
    _ls_list: list[LineSearch] = []  # list of line-search objects
    _hessian = None  # hessian object
    _structure = None  # eqm structure
    _structure_next = None  # next structure
    _path = ''
    _pes: PesFunction = None

    # Try to load the instance from file before ordinary init
    def __new__(cls, path='', load=None, *args, **kwargs):
        if load is None:
            return super().__new__(cls)
        else:
            # Try to load a pickle file from disk.
            try:
                fname = directorize(path) + load
                with open(fname, mode='rb') as f:
                    data = loads(f.read(), ignore=False)
                # end with
                if isinstance(data, cls):
                    return data
                else:
                    raise TypeError("The loaded file is not the same kind!")
                # end if
            except FileNotFoundError:
                return super().__new__(cls)
            # end try
        # end if
    # end def

    def __init__(
        self,
        # PLS arguments
        path='pls',
        hessian=None,
        structure=None,
        windows=None,
        window_frac=0.25,
        noises=None,
        add_sigma=False,
        no_eval=False,
        pes=None,
        pes_func=None,
        pes_args={},
        interactive=False,
        load=None,  # eliminate loading arg
        # LineSearch args
        **ls_args
        # M=7, fit_kind='pf3', fit_func=None, fit_args={}, N=200, Gs=None, fraction=0.025
    ):
        if load is not None and self.pes is not None:
            # Proxies of successful loading from disk
            return
        # end if
        if isinstance(pes, PesFunction):
            self.pes = pes
        else:
            # If none are provided, raises TypeError
            self.pes = PesFunction(pes_func, pes_args)
        # end if
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
            if self.shifted and not no_eval:
                # Successful evaluation leads to estimation of next structure
                self.evaluate(add_sigma=add_sigma, interactive=interactive)
            # end if
        # end if
    # end def

    @property
    def pes(self):
        return self._pes
    # end def

    @pes.setter
    def pes(self, pes):
        if isinstance(pes, PesFunction):
            self._pes = pes
        else:
            raise TypeError("Must provide PES that is inherited from PesFunction.")
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
    def structure_next(self):
        return self._structure_next
    # end def

    @property
    def Lambdas(self):
        if self.hessian is None:
            return array([])
        else:
            return array(self.hessian.enabled_lambdas)
        # end if
    # end def

    @property
    def windows(self):
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
    def noises(self):
        return [ls.sigma for ls in self.ls_list]
    # end def

    @property
    def noises_min(self):
        return array([ls.sigma for ls in self.ls_list]).min()
    # end def

    @property
    def D_list(self):
        return [d for d in range(len(self.hessian)) if self.hessian.enabled[d]]
    # end def

    def initialize(
        self,
        windows=None,
        noises=None,
        window_frac=None,
        **ls_args
        # M=7, fit_kind='pf3', fit_func=None, fit_args={}, N=200, Gs=None, fraction=0.025
    ):
        if windows is None:
            windows = abs(self.Lambdas)**0.5 * window_frac
        # end if
        if noises is None:
            noises = self.D * [0.0]
        # end if
        self._reset_ls_list(windows, noises, **ls_args)
    # end def

    def _reset_ls_list(
        self,
        windows,
        noises,
        **ls_args,
        # M=7, fit_kind='pf3', fit_func=None, fit_args={}, N=200, Gs=None, fraction=0.025
    ):
        ls_list = []
        for d, window, noise in zip(self.D_list, windows, noises):
            # Only add if enabled by the Hessian
            if self.hessian.enabled[d]:
                ls = self.ls_type(
                    structure=self.structure,
                    hessian=self.hessian,
                    d=d,
                    sigma=noise,
                    W=window,
                    **ls_args
                )
                ls_list.append(ls)
            # end if
        # end for
        self._ls_list = ls_list
        # Reset next structure if re-initialized
        self._structure_next = None
    # end def

    def evaluate(
        self,
        add_sigma=False,
        interactive=False,
        dep_jobs=[],
        var_eff_map=None,
    ):
        if not self.shifted:
            raise AssertionError("Must have shifted structures first!")
        # end if
        structures, sigmas = self._collect_enabled()
        self.pes.evaluate_all(
            structures,
            sigmas,
            add_sigma=add_sigma,
            path=self.path,
            interactive=interactive,
            dep_jobs=dep_jobs,
            var_eff_map=var_eff_map,
        )
        # Set the eqm energy
        for ls in self.ls_list:
            eqm = ls.find_point(0.0)
            if eqm is not None:
                self.structure.value = eqm.value
                self.structure.error = eqm.error
                break
            # end if
        # end for
        self._solve_ls()
        # Calculate next params
        params_next, params_next_err = self.calculate_next_params()  # **kwargs
        self._structure_next = self.structure.copy(
            params=params_next,
            params_err=params_next_err
        )
    # end def

    def evaluate_eqm(
        self,
        add_sigma=False,
        interactive=False,
        dep_jobs=[],
        var_eff_map=None,
    ):
        self.pes.evaluate(
            self.structure,
            sigma=array(self.noises).min(),
            path=self.path,
            add_sigma=add_sigma,
            interactive=interactive,
            dep_jobs=dep_jobs,
            var_eff_map=var_eff_map,
        )
    # end def

    def _collect_enabled(self):
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
        for ls in self.ls_list:
            ls._search_and_store()
        # end for
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
        pes=None
    ):
        structure = structure if structure is not None else self.structure
        hessian = hessian if hessian is not None else self.hessian
        windows = windows if windows is not None else self.windows
        noises = noises if noises is not None else self.noises
        pes = pes if pes is not None else self.pes
        copy_pls = ParallelLineSearch(
            path=path,
            structure=structure,
            hessian=hessian,
            windows=windows,
            noises=noises,
            no_eval=True,
            pes=pes
        )
        for ls, ls_new in zip(self.ls_list, copy_pls.ls_list):
            ls_new._settings = ls._settings
        # end for
        # If no new pes was supplied, use the old one
        if copy_pls.pes is None:
            copy_pls.pes = self.pes
        # end if
        return copy_pls
    # end def

    def propagate(
        self,
        path=None,
        write=True,
        overwrite=True,
        add_sigma=False,
        fname='pls.p',
        interactive=False,
        **kwargs  # dep_jobs=[], var_eff_map=None
    ):
        if not self.evaluated:
            self.evaluate(add_sigma=add_sigma, interactive=interactive, **kwargs)
        # end if
        path = path if path is not None else self.path + '_next/'
        # Write to disk
        if write:
            self.write_to_disk(fname=fname, overwrite=overwrite)
        # end if
        # check if manually providing structure
        pls_next = self.copy(
            path,
            structure=self.structure_next
        )
        return pls_next
    # end def

    def write_to_disk(self, fname='data.p', overwrite=False):
        fpath = directorize(self.path) + fname
        if path.exists(fpath) and not overwrite:
            print(f'File {fpath} exists. To overwrite, run with overwrite = True')
            return
        # end if
        makedirs(self.path, exist_ok=True)
        with open(fpath, mode='wb') as f:
            f.write(dumps(self, byref=True))
        # end with
    # end def

    def plot(
        self,
        **kwargs  # TODO: list kwargs
    ):
        for ls in self.ls_list:
            ls.plot(**kwargs)
        # end for
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
