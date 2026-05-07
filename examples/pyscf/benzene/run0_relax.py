#!/usr/bin/env python3

from stalk import ParameterStructure
from stalk import XyzGeometry
from stalk import BondLength
from stalk.util import Bohr

from params import forward, backward, relax_pyscf


# Let us initiate a ParameterStructure object that implements the parametric mappings
units = 'B'
params_init = [
    BondLength(2.651, label='r_CC', unit=units),
    BondLength(2.055, label='r_CH', unit=units),
]
elem = 6 * ['C'] + 6 * ['H']
structure_init = ParameterStructure(
    forward=forward,
    backward=backward,
    params=params_init,
    elem=elem,
    units=units,
)

xyz = XyzGeometry(suffix='relax.xyz', scale=Bohr)
structure_relax = xyz.load_or_relax(
    path='./',
    relax_func=relax_pyscf,
    structure=structure_init
)

if __name__ == '__main__':
    print('Initial parameters:')
    print(structure_init.params)
    print('Relaxed parameters:')
    print(structure_relax.params)
# end if
