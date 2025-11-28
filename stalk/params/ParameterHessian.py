#!/usr/bin/env python3
"""ParameterHessian class to consider Hessians according to a ParameterSet mapping."""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import array, linalg, diag, isscalar, ndarray, zeros, ones, where, mean, polyfit

from stalk.params.ParameterStructure import ParameterStructure
from stalk.params.PesFunction import PesFunction
from stalk.util import bipolyfit
from stalk.params.ParameterSet import ParameterSet


class ParameterHessian():
    _structure = None
    _hessian = None
    _Lambda = None
    _U = None
    require_consistent = True

    def __init__(
        self,
        hessian=None,
        structure=None,
        require_consistent=True,
    ):
        self.require_consistent = require_consistent
        self.structure = structure
        # Hessian must be set after structure
        if self.structure is not None:
            self.hessian = hessian
        # end if
    # end def

    @property
    def structure(self):
        return self._structure
    # end def

    @structure.setter
    def structure(self, structure):
        if structure is None:
            self.reset()
        elif isinstance(structure, ParameterSet):
            if not structure.check_consistency():
                if self.require_consistent:
                    raise AssertionError('The structure is not consistent! Aborting.')
                else:
                    warnings.warn("The structure is not consistent!")
                # end if
            # end if
            self._structure = structure
        else:
            raise TypeError('Structure must be inherited from ParameterSet class.')
        # end if
    # end def

    @property
    def hessian(self):
        return self._hessian
    # end def

    @hessian.setter
    def hessian(self, hessian):
        if self.structure is None:
            raise AssertionError('Cannot set Hessian without setting structure first!')
        # end if
        d = len(self.structure)
        if hessian is None:
            # Default Hessian
            self._hessian = diag(d * [1.0])
        else:
            hessian = array(hessian)
            if len(hessian.shape) != 2 or hessian.shape[0] != d or hessian.shape[1] != d:
                raise AssertionError(f'The Hessian must be {d}x{d} array, provided: {hessian.shape}')
            # end if
            self._hessian = hessian
        # end if
        Lambda, U = linalg.eig(self.hessian)
        self._Lambda = Lambda
        self._U = U
    # end def

    @property
    def U(self) -> ndarray:
        return self._U
    # end def

    @property
    def directions(self) -> ndarray:
        if self.U is not None:
            return self.U.T
        # end if
    # end def

    @property
    def lambdas(self) -> ndarray:
        return self._Lambda
    # end def

    def reset(self):
        self._structure = None
        self._hessian = None
        self._Lambda = None
        self._U = None
    # end def

    # Only preserved for backward compatibility
    def init_hessian_array(self, hessian):
        self.hessian = hessian
    # end def

    def compute_fdiff(
        self,
        pes: PesFunction,
        structure=None,
        dp=0.01,
        dpos_mode=False,
        **kwargs,
    ):
        if structure is not None:
            self.structure = structure
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
        pes.evaluate_all(structure_list, **kwargs)
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
    # end def

    def _get_fdiff_data(self, dps, dpos_mode=False):
        dp_list = [0.0 * dps]
        structure_list = [self.structure.copy(label='eqm')]

        def shift_params(id_ls, dp_ls):
            dparams = array(len(dps) * [0.0])
            label = 'eqm'
            for p, dp in zip(id_ls, dp_ls):
                dparams[p] += dp
                label += '_p{}'.format(p)
                if dp > 0:
                    label += '+'
                # end if
                label += '{}'.format(dp)
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

    def _warn_energy(self, structure_list: list[ParameterSet]):
        eqm_value = structure_list[0].value
        self.structure.value = eqm_value
        for structure in structure_list[0:]:
            if structure.value < eqm_value:
                warnings.warn(f'Offset energy lower than eqm: E({structure.label})={structure.value} < E(eqm)={eqm_value}!')
            # end if
        # end for
    # end def

    def __len__(self):
        if self.structure is None:
            return 0
        else:
            return len(self.structure)
        # end if
    # end def

    def __str__(self):
        string = self.__class__.__name__
        if self.hessian is not None:
            string += '\n  hessian:'
            for h in self.hessian:
                string += ('\n    ' + len(h) * '{:<8f} ').format(*tuple(h))
            # end for
            string += '\n  Conjugate directions:'
            string += '\n    Lambda     Direction'
            for Lambda, direction in zip(self.lambdas, self.directions):
                string += ('\n    {:<8f}   ' + len(direction) * '{:<+1.6f} ').format(Lambda, *tuple(direction))
            # end for
        else:
            string += '\n  hessian: not set'
        # end if
        return string
    # end def

# end class
