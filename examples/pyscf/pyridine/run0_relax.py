#!/usr/bin/env python3

from numpy import pi, sin, cos

from stalk import ParameterStructure
from stalk import XyzGeometry
from stalk import Parameter

from params import forward, backward, relax_pyscf


# Let us initiate a ParameterStructure object that implements the parametric mappings
r_CH = 1.08
pi_3 = pi / 3
params_init = [
    Parameter(1.335 * sin(pi_3), label='x_C04'),
    Parameter(1.335 * cos(pi_3), label='y_C04'),
    Parameter(1.335 * sin(pi_3), label='x_C13'),
    Parameter(1.395 + 1.335 * cos(pi_3), label='y_C13'),
    Parameter(1.395 + 1.335 * cos(pi_3) + 1.390 * cos(pi_3), label='y_C2'),
    Parameter((1.335 + r_CH) * sin(pi_3), label='x_H04'),
    Parameter((1.335 + r_CH) * cos(pi_3), label='y_H04'),
    Parameter((1.335 + r_CH) * sin(pi_3), label='x_H13'),
    Parameter(1.395 + (1.335 + r_CH) * cos(pi_3), label='y_H13'),
    Parameter(r_CH + 1.395 + 1.335 * cos(pi_3) + 1.390 * cos(pi_3), label='y_H2'),
]
elem = ['N'] + 5 * ['C'] + 5 * ['H']
structure_init = ParameterStructure(
    forward=forward,
    backward=backward,
    params=params_init,
    elem=elem,
    units='A'
)

xyz = XyzGeometry(suffix='relax.xyz')
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
