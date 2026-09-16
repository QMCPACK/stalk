#!/usr/bin/env python3

from numpy import pi

from stalk import ParameterStructure
from stalk import BondLength
from stalk import BondAngle

from params import forward, backward, relax_pyscf, pes_dict, relax_dict


# Let us initiate a ParameterStructure object that implements the parametric mappings
params_init = [
    BondLength(0.97, label='r_OH'),
    BondAngle(104.0 / 180 * pi, unit='rad', label='a_HOH')
]
elem = ['O'] + 2 * ['H']
structure = ParameterStructure(
    forward=forward,
    backward=backward,
    params=params_init,
    elem=elem,
    units='A'
)

# Treat a collection relaxed geometries based on alternative XC functionals
structure_relax = {}
for xc, pes in pes_dict.items():
    structure_relax[xc] = relax_dict[xc](
        path=f'{xc}',
        relax_func=relax_pyscf,
        structure=structure.copy(label='relax'),
    )
    print(f'Relaxed parameters ({xc}):')
    print(structure_relax[xc].params)
# end for
