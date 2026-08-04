#!/usr/bin/env python3

from numpy import pi

from stalk import ParameterStructure
from stalk import XyzGeometry
from stalk import BondLength
from stalk import Parameter

from params import forward, backward, relax_pbe

# Let us initiate a ParameterStructure object that implements the parametric mappings
elem = ['N'] + 3 * ['H']
structure_A_init = ParameterStructure(
    forward=forward,
    backward=backward,
    params=[
        BondLength(1.04, label='r_NH'),
        Parameter(1.2, label='a_Hx', unit='rad')
    ],
    elem=elem,
    units='A',
)
structure_B_init = ParameterStructure(
    forward=forward,
    backward=backward,
    params=[
        BondLength(1.04, label='r_NH'),
        Parameter(pi - 1.2, label='a_Hx', unit='rad')
    ],
    elem=elem,
    units='A',
)

# The suffix 'relax.xyz' is hardcoded to the relaxation function
xyz = XyzGeometry(suffix='relax.xyz')
structure_a = xyz.load_or_relax(
    path='pointA',
    relax_func=relax_pbe,
    structure=structure_A_init,
)
structure_b = xyz.load_or_relax(
    path='pointB',
    relax_func=relax_pbe,
    structure=structure_B_init,
)
