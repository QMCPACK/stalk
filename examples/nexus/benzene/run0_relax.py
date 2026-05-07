#!/usr/bin/env python3

from stalk.nexus import NexusStructure
from stalk import BondLength

from params import forward, backward, relax_pyscf

interactive = __name__ == "__main__"

# Let us initiate a NexusStructure object that implements the parametric mappings
units = 'B'
params_init = [
    BondLength(2.651, label='r_CC', unit=units),
    BondLength(2.055, label='r_CH', unit=units),
]
elem = 6 * ['C'] + 6 * ['H']
structure_init = NexusStructure(
    forward=forward,
    backward=backward,
    params=params_init,
    elem=elem,
    units=units,
)

structure_relax = structure_init.copy()

relax_pyscf.relax(
    structure_relax,
    path='relax/',
    interactive=interactive,
)

if interactive:
    print('Initial params:')
    print(structure_init.params)
    print('Relaxed params:')
    print(structure_relax.params)
# end if
