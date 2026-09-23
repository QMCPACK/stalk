#!/usr/bin/env python3
"""ParameterHessian class to consider Hessians according to a ParameterSet mapping."""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path
import warnings
from numpy import array, linalg, diag, isscalar, ndarray, zeros, ones, where, mean, polyfit

from stalk.io.txt_data import TxtData
from stalk.params.parameter_structure import ParameterStructure
from stalk.pes.pes_function import NotEvaluatedException, PesFunction
from stalk.util import bipolyfit
from stalk.params.parameter_set import ParameterSet
from stalk.io.cacheable import Cacheable


class ParameterHessian(Cacheable):
    _structure: ParameterSet = None
    _hessian: ndarray = None
    _Lambda: ndarray = None
    _U: ndarray = None
    _enabled: list[bool] = None
    _disable_limit: float = None

    def __init__(
        self,
        structure: ParameterSet,
        hessian: ndarray | None = None,
        disable_limit: float = 0.0,
    ):
        # The ParameterHessian must have a ParameterSet structure
        self.structure = structure
        self.disable_limit = disable_limit
        self.hessian = hessian
        Cacheable.__init__(
            self,
            hessian=TxtData('hessian.out'),
        )
    # end def

    @property
    def structure(self) -> ParameterSet:
        return self._structure
    # end def

    @structure.setter
    def structure(self, structure: ParameterSet) -> None:
        if isinstance(structure, ParameterStructure):
            if structure.require_consistent and not structure.check_consistency():
                raise AssertionError('The structure is not consistent! Aborting.')
            # end if
            self._structure = structure
        elif isinstance(structure, ParameterSet):
            self._structure = structure
        else:
            raise TypeError(f'Structure must be a ParameterSet, provided {type(structure)}')
        # end if
    # end def

    @property
    def hessian(self) -> ndarray | None:
        if self._hessian is None:
            # Default Hessian
            d = len(self.structure)
            return diag(d * [1.0])
        else:
            return self._hessian
        # end if
    # end def

    @hessian.setter
    def hessian(self, hessian: ndarray) -> None:
        if hessian is None:
            self._hessian = None
            self._enabled = None
            return
        # end if
        d = len(self.structure)
        hessian = array(hessian)
        if len(hessian.shape) != 2 or hessian.shape[0] != d or hessian.shape[1] != d:
            raise AssertionError(f'The Hessian must be {d}x{d} array, provided: {hessian.shape}')
        # end if
        self._hessian = hessian
        Lambda, U = linalg.eig(self.hessian)
        self._Lambda = Lambda
        self._U = U
        # List of disabled directions
        self._enabled = len(self) * [True]
        # Disable directions softer than the disable limit
        lmax = self.lambdas.max()
        for h in range(len(self)):
            condition = abs(self.lambdas[h]) / lmax
            if condition < self.disable_limit:
                self.enabled[h] = False
                print(f'Disabled direction {h} with Lambda/max(Lambda)={condition:.4} < {self.disable_limit}')
            # end if
        # end for
    # end def

    @property
    def disable_limit(self) -> float:
        return self._disable_limit
    # end def

    @disable_limit.setter
    def disable_limit(self, disable_limit):
        if isscalar(disable_limit) and disable_limit >= 0.0:
            self._disable_limit = disable_limit
        else:
            raise TypeError(f'Disable limit must be nonnegative, provided {disable_limit}')
        # end if
    # end def

    @property
    def enabled(self):
        return self._enabled
    # end def

    @property
    def U(self) -> ndarray | None:
        return self._U
    # end def

    @property
    def directions(self) -> ndarray | None:
        if self.U is not None:
            return self.U.T
        # end if
    # end def

    @property
    def enabled_directions(self) -> ndarray | None:
        return self.U.T[where(self.enabled)]
    # end def

    @property
    def lambdas(self) -> ndarray | None:
        return self._Lambda
    # end def

    @property
    def enabled_lambdas(self) -> ndarray | None:
        return self.lambdas[where(self.enabled)]
    # end def

    # Only preserved for backward compatibility
    def init_hessian_array(self, hessian):
        self.hessian = hessian
    # end def

    def compute_fdiff(
        self,
        pes: PesFunction,
        path: Path | str | None = None,
        dp=0.01,
        dpos_mode=False,
        **kwargs,
    ):
        if self.load_hessian(path):
            return
        # end if
        P = len(self)

        # Figure out finite differences
        if isscalar(dp):
            dps = array(P * [dp])
        elif len(dp) == P:
            dps = array(dp)
        else:
            raise ValueError(f'Error: Provided {len(dp)} dps for {P} directions! Aborting.')
        # end if

        # Get list of displacements and structures
        dp_list, structure_list = self._get_fdiff_data(dps, dpos_mode=dpos_mode)
        pes.evaluate_all(
            structure_list,
            path=path,
            **kwargs
        )
        if not all([s.value is not None and s.enabled for s in structure_list]):
            print('Did not update the Hessian, as not all structures were evaluated successfully!')
            return
        # end if
        # Issue warning when eqm energy is not the apparent minimum
        self._warn_energy(structure_list)

        # Pick those displacements and energies that were successfully computed
        energies = []
        pdiffs = []
        for dp, s in zip(dp_list, structure_list):
            if s.value is not None and s.enabled:
                pdiffs.append(dp)
                energies.append(s.value)
            # end if
        # end for
        pdiffs = array(pdiffs)
        energies = array(energies)

        params = self.structure.params
        if P == 1:  # for 1-dimensional problems
            pf = polyfit(pdiffs[:, 0], energies, 2)
            hessian = array([[pf[0]]])
        else:
            hessian = zeros((P, P))
            pfs = [[] for p in range(P)]
            for p0, param0 in enumerate(params):
                for p1, param1 in enumerate(params):
                    if p1 <= p0:
                        continue
                    # end if
                    # filter out the values where other parameters were altered
                    ids = ones(len(pdiffs), dtype=bool)
                    for p in range(P):
                        if p == p0 or p == p1:
                            continue
                        # end if
                        ids = ids & (abs(pdiffs[:, p]) < 1e-10)
                    # end for
                    XY = pdiffs[where(ids)]
                    E = energies[where(ids)]
                    X = XY[:, p0]
                    Y = XY[:, p1]
                    pf = bipolyfit(X, Y, E, 2, 2)
                    hessian[p0, p1] = pf[4]
                    hessian[p1, p0] = pf[4]
                    pfs[p0].append(2 * pf[6])
                    pfs[p1].append(2 * pf[2])
                # end for
            # end for
            for p0 in range(P):
                hessian[p0, p0] = mean(pfs[p0])
            # end for
        # end if
        self.hessian = hessian
        self.save_hessian(path)
    # end def

    def _get_fdiff_data(self, dps, dpos_mode=False):
        dp_list = [0.0 * dps]
        structure_list = [self.structure.copy(label='eqm')]

        def shift_params(id_ls, dp_ls):
            dparams = array(len(dps) * [0.0])
            label = 'eqm'
            for p, dp in zip(id_ls, dp_ls):
                dparams[p] += dp
                label += f'_{self.structure.params_list[p].label}'
                if dp > 0:
                    label += '+'
                # end if
                label += f'{dp}'
            # end for
            structure_new = self.structure.copy(label=label)
            if isinstance(structure_new, ParameterStructure):
                structure_new.shift_params(dparams, dpos_mode=dpos_mode)
            else:
                structure_new.shift_params(dparams)
            # end if
            structure_list.append(structure_new)
            dp_list.append(dparams)
        # end def

        for p0, dp0 in enumerate(dps):
            shift_params([p0], [+dp0])
            shift_params([p0], [-dp0])
            for p1, dp1 in enumerate(dps):
                if p1 <= p0:
                    continue
                # end if
                shift_params([p0, p1], [+dp0, +dp1])
                shift_params([p0, p1], [+dp0, -dp1])
                shift_params([p0, p1], [-dp0, +dp1])
                shift_params([p0, p1], [-dp0, -dp1])
            # end for
        # end for
        return dp_list, structure_list
    # end def

    def load_hessian(self, path: str) -> bool:
        if path is None:
            return False
        # end if
        hessian = self.load(path, key='hessian')
        params = self.structure.load(key='params', path=path)
        if isinstance(params, ndarray) and isinstance(self.structure, ParameterSet):
            self.structure.params = params
            print(f'Loaded Hessian parameters from {path}.')
            if isinstance(hessian, ndarray):
                self.hessian = hessian
                print(f'Loaded Hessian from {path}.')
                return True
            # end if
        # end if
        return False
    # end def

    def save_hessian(self, path: Path | str) -> None:
        if self.hessian is not None:
            self.save(path, hessian=self.hessian)
            self.structure.save(path=path, params=self.structure.params)
            print(f'Saved Hessian to {path}.')
        # end if
    # end def

    def _warn_energy(self, structure_list: list[ParameterSet]):
        eqm_value = structure_list[0].value
        if eqm_value is None:
            raise NotEvaluatedException('Cannot issue energy warning, eqm energy is not available!')
        # end if
        self.structure.value = eqm_value
        for structure in structure_list[0:]:
            if structure.value < eqm_value:
                warnings.warn(f'E({structure.label})={structure.value} < E(eqm)={eqm_value}!')
            # end if
        # end for
    # end def

    def __len__(self) -> int:
        if self.structure is None:
            return 0
        else:
            return len(self.structure)
        # end if
    # end def

    def __str__(self) -> str:
        string = self.__class__.__name__
        if self.hessian is not None:
            string += '\n  hessian:'
            for h in self.hessian:
                string += ('\n    ' + len(h) * '{:<+1.6f} ').format(*tuple(h))
            # end for
            string += '\n  Conjugate directions:'
            string += '\n    Lambda     Direction'
            for Lambda, direction, enabled in zip(self.lambdas, self.directions, self.enabled):
                string += '\n    {:<8f}   '.format(Lambda)
                string += (len(direction) * '{:<+1.6f} ').format(*tuple(direction))
                if not enabled:
                    string += ' <- disabled'
                # end if
            # end for
        else:
            string += '\n  hessian: not set'
        # end if
        return string
    # end def

# end class
