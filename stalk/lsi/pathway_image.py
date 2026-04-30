#!/usr/bin/env python3
'''LineSearchIteration class for treating iteration of subsequent parallel linesearches'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from copy import copy
from os import makedirs
from numpy import savetxt, loadtxt, array
from functools import partial

from numpy import ndarray, zeros
from stalk.lsi.linesearch_iteration import LineSearchIteration
from stalk.params.parameter_hessian import ParameterHessian
from stalk.params.parameter_set import ParameterSet
from stalk.params.pes_function import PesFunction
from stalk.pls.surrogate import Surrogate
from stalk.util.util import get_fraction_error, orthogonal_subspace_basis


class PathwayImage():
    _path = None  # Path is set upon Hessian calculation
    _lsi: LineSearchIteration = None  # Line-search iteration
    _structure: ParameterSet = None
    _hessian: ParameterHessian = None
    _tangent = None
    _subspace = None
    _reaction_coordinate = None
    _surrogate = None

    def __init__(
        self,
        structure: ParameterSet,
        reaction_coordinate=0
    ):
        self._structure = structure
        self._reaction_coordinate = reaction_coordinate
    # end def

    @property
    def lsi(self):
        return self._lsi
    # end def

    @property
    def structure(self):
        return self._structure
    # end def

    @property
    def hessian(self):
        return self._hessian
    # end def

    @property
    def surrogate(self):
        return self._surrogate
    # end def

    @property
    def reaction_coordinate(self):
        return self._reaction_coordinate
    # end def

    @property
    def structure_init(self):
        if self._tangent is None:
            return self.lsi.structure_init
        else:
            return extend_structure(self.structure, self.lsi.structure_init, self._subspace)
        # end if
    # end def

    @property
    def structure_final(self):
        if self._tangent is None:
            return self.lsi.structure_final
        else:
            return extend_structure_errors(
                self.structure,
                self.lsi.structure_final,
                self._subspace
            )
        # end if
    # end def

    def calculate_hessian(
        self,
        tangent,
        pes: PesFunction,
        path='',
        **hessian_args  # dp=0.01, dpos_mode=False, structure=None
    ):
        hessian_file = f'{path}/hessian.dat'
        makedirs(path, exist_ok=True)
        if tangent is None:
            # points A and B are calculated in full parametric space
            hessian = ParameterHessian(structure=self.structure)
            subspace = None
            pes_comp = pes
        else:
            # Calculate orthogonal subspace exluding the tangent direction
            subspace = orthogonal_subspace_basis(tangent)
            # Make a copy of the PesFunction wrapper, then replace the func
            pes_comp = copy(pes)
            pes_comp.func = partial(extended_pes, self.structure, subspace, pes)
            # Subspace parameter set is a zero-centered p-1 vector
            structure_sub = ParameterSet(zeros(len(subspace)))
            hessian = ParameterHessian(structure=structure_sub)
        # end if
        try:
            hessian_array = loadtxt(hessian_file, ndmin=2)
            hessian.init_hessian_array(hessian_array)
        except FileNotFoundError:
            hessian.compute_fdiff(
                pes=pes_comp,
                **hessian_args
            )
            savetxt(hessian_file, hessian.hessian)
        # end try
        self._path = path
        self._subspace = subspace
        self._tangent = tangent
        self._hessian = hessian
    # end def

    def generate_surrogate(
        self,
        pes: PesFunction = None,
        overwrite=False,
        **surrogate_args
    ):
        if self._tangent is None:
            pes_sub = pes
        else:
            # Make a copy of the PesFunction wrapper, then replace the func
            pes_sub = copy(pes)
            pes_sub.func = partial(extended_pes, self.structure, self._subspace, pes)
        # end if
        path = '{}surrogate'.format(self._path)
        surrogate = Surrogate(
            load='data.p',
            path=path,
            hessian=self.hessian,
            pes=pes_sub,
            **surrogate_args
        )
        surrogate.write_to_disk(overwrite=overwrite)
        self._surrogate = surrogate
    # end def

    def optimize_surrogate(
        self,
        overwrite=True,
        **optimize_args
    ):
        self.surrogate.bracket_target_biases()
        self.surrogate.optimize(**optimize_args)
        self.surrogate.write_to_disk(overwrite=overwrite)
    # end def

    def run_linesearch(
        self,
        num_iter=3,
        path='lsi',
        pes: PesFunction = None,
        add_sigma=False,
        **lsi_args
    ):
        if self._tangent is None:
            pes_comp = pes
        else:
            # Make a copy of the PesFunction wrapper, then replace the func
            pes_comp = copy(pes)
            pes_comp.func = partial(extended_pes, self.structure, self._subspace, pes)
        # end if
        lsi = LineSearchIteration(
            path=self._path + path,
            surrogate=self.surrogate,
            pes=pes_comp,
            **lsi_args
        )
        for i in range(num_iter):
            lsi.propagate(i, add_sigma=add_sigma)
        # end for
        self._lsi = lsi
    # end def

    def __lt__(self, other):
        return (hasattr(other, 'reaction_coordinate') and
                self.reaction_coordinate > other.reaction_coordinate)
    # end def

# end class


def extend_structure(structure0: ParameterSet, structure_sub: ParameterSet, subspace):
    structure = structure0.copy(label=structure_sub.label)
    structure.shift_params(structure_sub.params @ subspace)
    return structure
# end def


def extend_structure_errors(
    structure0: ParameterSet,
    structure_sub: ParameterSet,
    subspace,
    N=200,
    fraction=0.025
):
    ps_sub = structure_sub.get_params_distribution(N=N)
    ps = []
    for p_sub in ps_sub:
        ps.append(structure0.params + p_sub @ subspace)
    # end for
    ps = array(ps).T
    params_err = [get_fraction_error(p, fraction=fraction)[1] for p in ps]
    structure = structure0.copy(label=structure_sub.label)
    structure.set_params(
        structure0.params + structure_sub.params @ subspace,
        params_err=params_err
    )
    structure.value = structure_sub.value
    structure.error = structure_sub.error
    return structure
# end def


def extended_pes(
    structure: ParameterSet,
    subspace: ndarray,
    pes: PesFunction,
    structure_sub: ParameterSet,
    **kwargs
):
    new_structure = extend_structure(structure, structure_sub, subspace)
    return pes.func(new_structure, **kwargs)
# end def
