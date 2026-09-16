#!/usr/bin/env python3

from numpy import pi

from stalk import ParameterStructure
from stalk import XyzGeometry
from stalk import BondLength
from stalk import Parameter

from params import forward, backward, relax_pbe

# Let us initiate a ParameterStructure object that implements the parametric mappings
elem = ['N'] + 3 * ['H']
structure_a_init = ParameterStructure(
    forward=forward,
    backward=backward,
    params=[
        BondLength(1.04, label='r_NH'),
        Parameter(1.2, label='a_Hx', unit='rad')
    ],
    elem=elem,
    units='A',
    label='pointA',
)
structure_b_init = ParameterStructure(
    forward=forward,
    backward=backward,
    params=[
        BondLength(1.04, label='r_NH'),
        Parameter(pi - 1.2, label='a_Hx', unit='rad')
    ],
    elem=elem,
    units='A',
    label='pointB',
)

# The suffix 'relax.xyz' is hardcoded to the relaxation function
xyz = XyzGeometry(suffix='relax.xyz')
structure_a = relax_pbe(
    structure=structure_a_init,
    path='relax',
)
structure_b = relax_pbe(
    structure=structure_b_init,
    path='relax',
)
