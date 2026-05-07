#!/usr/bin/env python3

from stalk import BondLength
from stalk.nexus import NexusStructure

from params import forward, backward, relax_pyscf

interactive = __name__ == "__main__"

# Let us initiate a NexusStructure object that implements the parametric mappings
elem = 2 * ['C']
structure_init = NexusStructure(
    forward=forward,
    backward=backward,
    params=[BondLength(1.54, label='r_CC', unit='A')],
    elem=elem,
    units='A'
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
