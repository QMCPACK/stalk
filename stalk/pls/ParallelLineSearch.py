#!/usr/bin/env python3
'''ParallelLineSearch class for simultaneous linesearches along conjugate directions'''

from numpy import ndarray, array, matmul
from os import makedirs, path
from dill import dumps, loads
from textwrap import indent

from stalk.util import get_fraction_error, directorize
from stalk.params import ParameterSet
from stalk.params import ParameterHessian
from stalk.ls import LineSearch, LineSearchDummy
from .PesSampler import PesSampler

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Class for a bundle of parallel line-searches
class ParallelLineSearch(PesSampler):

    skip_init = False
    ls_type = LineSearch
    ls_list: list[LineSearch] = []  # list of line-search objects
    hessian = None  # hessian object
    Lambdas = None
    directions = None
    D = None
    M = None  # number of grid points
    path = None
    windows = None
    noises = None
    fraction = None
    structure = None  # eqm structure
    structure_next = None  # next structure
    x_unit = None
    E_unit = None
    noisy = False  # flag whether deterministic or noisy
    msg_setup = 'Setup: first set_structure() and set_hessian() with valid input'
    msg_shifted = 'Shifted: first set_windows() to define displacements'
    msg_loaded = 'Loaded: first load_results() with valid input'

    # Try to load the instance from file before ordinary init
    def __new__(cls, load=None, mode=None, *args, **kwargs):
        if load is None:
            return super().__new__(cls)
        else:
            # Try to load a pickle file from disk.
            try:
                with open(load, mode='rb') as f:
                    data = loads(f.read())
                # end with
                data.skip_init = True
                return data
            except FileNotFoundError:
                # If only trying to load a file, stop right here
                if mode == 'load':
                    return None
                else:
                    return super().__new__(cls)
                # end if
            # end try
        # end if
    # end def

    def __init__(
        self,
        # Remove disk loading argument
        load=None,
        # PesSampler arguments
        mode=None,
        pes=None,
        pes_func=None,
        pes_args={},
        loader=None,
        load_func=None,
        load_args={},
        # PLS arguments
        hessian=None,
        structure=None,
        windows=None,
        window_frac=0.25,
        noises=None,
        path='pls',
        M=7,
        fit_kind='pf3',
        x_unit='A',
        E_unit='Ry',
        fraction=0.025,
        shift_params=None,
        D=None
    ):
        if self.skip_init:
            return
        # end if
        PesSampler.__init__(self, mode, pes=pes, pes_func=pes_func, pes_args=pes_args,
                            loader=loader, load_func=load_func, load_args=load_args)
        self.x_unit = x_unit
        self.E_unit = E_unit
        self.path = directorize(path)
        self.fraction = fraction
        self.M = M
        self.fit_kind = fit_kind
        if structure is not None:
            self.set_structure(structure, shift_params)
        # end if
        if hessian is not None:
            self.set_hessian(hessian)
        # end if
        if self.status.setup:
            self.guess_windows(windows, window_frac, no_reset=True, D=D)
            self.set_noises(noises, D=D)
        # end if
        if self.mode == 'pes':
            if self.status.shifted:
                self.generate_pls_jobs()
            # end if
            if self.status.generated:
                self.load_results()
            # end if
        # end if
    # end def

    def _setup(self):
        # test for proper hessian and structure
        if isinstance(self.hessian, ParameterHessian) and isinstance(self.structure, ParameterSet):
            # test for compatibility
            if not len(self.hessian.Lambda) == len(self.structure.params):
                return False
            # end if
        else:
            return False
        # end if
        return PesSampler._setup(self)
    # end def

    # ie displacement of positions
    def _shifted(self):
        # check windows
        if self.windows is None or not len(self.windows) == self.D or self.M is None:
            return False
        # end if
        # check noises
        if self.noisy and self.noises is None:
            return False
        # end if
        return PesSampler._shifted(self)
    # end def

    # def _generated(self): not overridden
    def _loaded(self):
        if not all([ls.loaded for ls in self.ls_list]):
            return False
        # end if
        return PesSampler._loaded(self)
    # end def

    def _analyzed(self):
        if not all([ls.analyzed for ls in self.ls_list]):
            return False
        # end if
        return PesSampler._analyzed(self)
    # end def

    def set_hessian(self, hessian):
        self._avoid_protected()
        if hessian is None:
            return
        elif isinstance(hessian, ndarray):
            hessian = ParameterHessian(
                hessian=hessian, x_unit=self.x_unit, E_unit=self.E_unit)
        elif isinstance(hessian, ParameterHessian):
            pass
        else:
            raise ValueError('Hessian matrix is not supported')
        # end if
        self.hessian = hessian
        self.Lambdas = hessian.Lambda
        self.directions = hessian.get_directions()
        self.D = hessian.D
        self.cascade()
    # end def

    def set_structure(self, structure, shift_params=None):
        self._avoid_protected()
        if structure is None:
            return
        # end if
        assert isinstance(
            structure, ParameterSet), 'Structure must be ParameterSet object'
        self.structure = structure.copy(label='eqm')
        if shift_params is not None:
            self.structure.shift_params(shift_params)
        # end if
        self.cascade()
    # end def

    def guess_windows(self, windows, window_frac, **kwargs):
        self._avoid_protected()
        self._require_setup()
        if windows is None:
            windows = abs(self.Lambdas)**0.5 * window_frac
            self.windows_frac = window_frac
        # end if
        self.set_windows(windows, **kwargs)
    # end def

    def set_windows(self, windows, no_reset=False, **kwargs):
        self._avoid_protected()
        self._require_setup()
        if windows is not None:
            assert windows is not None or len(
                windows) == self.D, 'length of windows differs from the number of directions'
            self.windows = array(windows)
        # end if
        self.cascade()
        if not no_reset:
            self.reset_ls_list(**kwargs)  # always reset ls_list
        # end if
    # end def

    def set_noises(self, noises, **kwargs):
        self._avoid_protected()
        self._require_setup()
        if noises is None:
            self.noisy = False
            self.noises = None
        else:
            assert (len(noises) == self.D)
            self.noisy = True
            self.noises = array(noises)
        # end if
        self.cascade()
        self.reset_ls_list(**kwargs)  # always reset ls_list
    # end def

    def reset_ls_list(self, D=None, **kwargs):
        self._avoid_protected()
        self._require_shifted()
        if D is None:
            D = range(self.D)
        else:
            D = [d if d in D else None for d in range(self.D)]
        # end if
        assert len(D) == self.D, 'len(D) must match D'
        noises = self.noises if self.noisy else self.D * [None]
        ls_list = []
        for d, window, noise in zip(D, self.windows, noises):
            if d is None:
                ls = LineSearchDummy(d=d)
            else:
                ls = self.ls_type(
                    structure=self.structure,
                    hessian=self.hessian,
                    d=d,
                    W=window,
                    M=self.M,
                    fit_kind=self.fit_kind,
                    sigma=noise,
                    **kwargs)
            ls_list.append(ls)
        # end for
        self.ls_list = ls_list
        self.cascade()
    # end def

    def copy(self, path='', c_noises=1.0, update_hessian=False, squeeze=1.0, **kwargs):
        if self.noises is None:
            noises = None
        else:
            noises = array(
                [noise * c_noises for noise in self.noises]) * squeeze
        # end if
        if self.windows is None:
            windows = None
        else:
            windows = array(self.windows) * squeeze
        # end if
        ls_args = {
            'path': path,
            'structure': self.structure,
            'hessian': self.hessian,
            'windows': windows,
            'noises': noises,
            'M': self.M,
            'fit_kind': self.fit_kind,
            'pes': self.pes,
            'loader': self.loader,
            'mode': self.mode,
        }
        if update_hessian:
            ls_args['hessian'] = self.updated_hessian()
        # end if
        ls_args.update(**kwargs)
        pls_next = ParallelLineSearch(**ls_args)
        return pls_next
    # end def

    def propagate(self, path=None, protect=True, write=True, **kwargs):
        if not self.status.analyzed:
            return
        # end if
        path = path if path is not None else self.path + '_next/'
        # check if manually providing structure
        if 'structure' in kwargs.keys():
            # TODO assert
            pls_next = self.copy(path=path, **kwargs)
        else:
            pls_next = self.copy(
                path=path, structure=self.structure_next, **kwargs)
        # end if
        self.status.protected = protect
        if write:
            self.write_to_disk()
        # end if
        return pls_next
    # end def

    def generate_pls_jobs(self, **kwargs):
        self._require_shifted()

        if (self.mode == 'pes'):
            self.status.generated = True
            return
        # end if

        # Find the smallest noise target for the eqm job
        sigma_min = None if not self.noisy else self.noises.min()

        # Generate eqm jobs
        eqm_jobs = self.ls_list[0].generate_eqm_jobs(
            self.pes, path=self.path, sigma=sigma_min)
        jobs = eqm_jobs

        # Generate line-search jobs
        for ls in self.ls_list:
            jobs += ls.generate_ls_jobs(self.pes, path=self.path,
                                        eqm_jobs=eqm_jobs, **kwargs)
        # end for
        self.status.generated = True
        return jobs
    # end def

    def generate_eqm_jobs(self, **kwargs):
        sigma_min = None if not self.noisy else self.noises.min()
        eqm_jobs = self.ls_list[0].generate_eqm_jobs(
            self.pes, path=self.path, sigma=sigma_min, **kwargs)
        return eqm_jobs
    # end def

    def run_jobs(self, interactive=True, eqm_only=False, **kwargs):
        if eqm_only:
            jobs = self.generate_eqm_jobs(**kwargs)
        else:
            jobs = self.generate_pls_jobs(**kwargs)
        # end if
        if jobs is None or jobs == []:
            return
        # end if
        if interactive:
            print('About to submit the following new jobs:')
            any_new = False
            for job in jobs:
                if not job.submitted:
                    print('  {}'.format(job.path))
                    any_new = True
                # end if
            # end for
            if any_new:
                if input('proceed? (Y/n) ') in ['n', 'N']:
                    exit()
                # end if
            # end if
        # end if
        from nexus import run_project
        run_project(jobs)
    # end def

    # can either load based on analyze_func or by providing values/errors
    def load_results(self, loader=None, values=None, errors=None, add_sigma=False, **kwargs):
        if self.status.protected:
            return
        # end if

        if self.mode == 'pes':
            values_ls, errors_ls = [], []
            for ls in self.ls_list:
                ls.evaluate_pes(pes_eval=self.pes, add_sigma=add_sigma)
                values_ls.append(ls.valid_values)
                errors_ls.append(ls.valid_errors)
            # end for
        else:
            # Unless a list/array of per-direction values/errors is provided, init to None
            values_ls = values if values is not None else self.D * [None]
            errors_ls = errors if errors is not None else self.D * [None]
        # end if

        loaded = True
        loader = loader if loader is not None else self.loader
        for ls, values, errors in zip(self.ls_list, values_ls, errors_ls):
            if not isinstance(ls, LineSearch):
                continue
            # end if
            loaded_this = ls.set_results(
                ls.grid,
                values=values,
                errors=errors,
                add_sigma=add_sigma,
                **kwargs)
            loaded = loaded and loaded_this
        # end for
        if loaded:
            self.find_eqm_value()
            self.status.generated = True
            self.status.loaded = True
            self.calculate_next()
        # end if
        self.cascade()
    # end def

    def load_eqm_results(self, add_sigma=False, **kwargs):
        if self.status.protected:
            return
        # end if
        E, err = self.ls_list[0].load_eqm_results(
            self.loader,
            path=self.path,
            add_sigma=add_sigma,
            **kwargs
        )
        self.structure.value = E
        self.structure.error = err
    # end def

    def find_eqm_value(self):
        E, err = None, None
        for ls in filter(lambda x: not isinstance(x, LineSearchDummy), self.ls_list):
            eqm = ls.find_point(0.0)
            E, err = eqm.value, eqm.error
        # end for
        self.structure.value = E
        self.structure.error = err
    # end def

    def calculate_next(self, **kwargs):
        self._avoid_protected()
        self._require_loaded()
        params_next, params_next_err = self.calculate_next_params(**kwargs)
        self.structure_next = self.structure.copy(
            params=params_next, params_err=params_next_err)
        self.cascade()
    # end def

    def ls(self, i) -> LineSearch:
        return self.ls_list[i]
    # end def

    def get_next_params(self):
        assert self.structure_next is not None, 'Next structure has not been computed yet'
        return self.structure_next.params
    # end def

    def calculate_next_params(self, **kwargs):
        # deterministic
        params_next = self._calculate_params_next(
            self.get_params(), self.get_directions(), self.get_shifts())
        # stochastic
        if self.noisy:
            params_next_err = self._calculate_params_next_error(
                self.get_params(), self.get_directions(), params_next, **kwargs)
        else:
            params_next_err = array(self.D * [0.0])
        # end if
        return params_next, params_next_err
    # end def

    def get_params(self):
        return self.structure.params
    # end def

    def get_params_err(self):
        err = self.structure.params_err
        if err is None:
            return self.structure.params * 0.0
        else:
            return err
        # end if
    # end def

    def get_shifted_params(self, i=None):
        if i is None:
            return [ls.get_shifted_params() for ls in self.ls_list]
        else:
            return self.ls(i).get_shifted_params()
        # end if
    # end def

    def get_directions(self, d=None):
        return self.hessian.get_directions(d)
    # end def

    def get_shifts(self):
        self._require_loaded()
        return self._get_shifts()
    # end def

    def _get_shifts(self):
        return array([ls.get_x0(err=False) for ls in self.ls_list])
    # end def

    def _calculate_params_next(self, params, directions, shifts):
        return params + self._calculate_shifts(directions, shifts)
    # end def

    def _calculate_shifts(self, directions, shifts):
        return shifts @ directions
    # end def

    def _calculate_params_next_error(self, params, directions, params_next, N=200, fraction=0.025, **kwargs):
        x0s_d = self._get_x0_distributions(N=N)
        params_d = []
        for x0s in x0s_d:
            params_d.append(self._calculate_params_next(
                params, directions, x0s) - params_next)
        # end for
        params_next_err = [get_fraction_error(p, fraction=fraction)[
            1] for p in array(params_d).T]
        return array(params_next_err)
    # end def

    def _get_x0_distributions(self, N=200, **kwargs):
        return array([ls.get_x0_distribution(errors=ls.errors, N=N, **kwargs) for ls in self.ls_list]).T
    # end def

    def write_to_disk(self, fname='data.p', overwrite=False):
        if path.exists(self.path + fname) and not overwrite:
            print('File {} exists. To overwrite, run with overwrite = True'.format(
                self.path + fname))
            return
        # end if
        makedirs(self.path, exist_ok=True)
        with open(self.path + fname, mode='wb') as f:
            f.write(dumps(self, byref=True))
        # end with
    # end def

    def __str__(self):
        string = self.__class__.__name__
        if self.ls_list is None:
            string += '\n  Line-searches: None'
        else:
            string += '\n  Line-searches:\n'
            string += indent('\n'.join(['#{:<2d} {}'.format(ls.d, str(ls)) for ls in filter(
                lambda x: not isinstance(x, LineSearchDummy), self.ls_list)]), '    ')
        # end if
        # TODO
        return string
    # end def

    def plot(self, **kwargs):
        for ls in self.ls_list:
            ls.plot(**kwargs)
        # end for
    # end def

    def updated_hessian(self, method='powell', **kwargs):
        H = self.hessian.hessian
        U = self.hessian.U
        if method == 'powell':
            d = matmul(U, array([[ls.x0 for ls in self.ls_list]]).T)
            y = matmul(U, array([[ls.get_force() for ls in self.ls_list]]).T)
            u = d * (d.T @ d)**-0.5
            j = y - H @ d
            dH = j @ u.T + u @ j.T - d.T @ j * u @ u.T
        else:
            print('Method {} not implemented'. format(method))
            dH = 0.0
        # end if
        return ParameterHessian(hessian=H + dH, structure=self.structure_next)
    # end def

# end class
