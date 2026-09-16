#!/usr/bin/env python3

from numpy import array

from stalk.nexus import NexusStructure
from stalk import Parameter

from params import forward, backward, relax_pwscf

interactive = __name__ == "__main__"

# Let us initiate a NexusStructure object that implements the parametric mappings
params_init = array([1.7])
elem = 2 * ['C']
structure_init = NexusStructure(
    forward=forward,
    backward=backward,
    params=[Parameter(1.7, label='a')],
    elem=elem,
    units='A'
)

structure_relax = structure_init.copy()

structure_relax = relax_pwscf(
    structure_init.copy(label='relax'),
    path='./',
    interactive=interactive,
)

if interactive:
    print('Initial params:')
    print(structure_init.params)
    print('Relaxed params:')
    print(structure_relax.params)
# end if
