#!/usr/bin/env python3

from os import makedirs

from stalk import ParameterStructure
from stalk import XyzGeometry
from stalk import BondLength
from stalk import Parameter

from params import forward, backward, relax_pyscf

# Generate base directory for point A
basedir = 'pointA'
makedirs(basedir, exist_ok=True)

# Let us initiate a ParameterStructure object that implements the parametric mappings
elem = ['N'] + 3 * ['H']
structure_init = ParameterStructure(
    forward=forward,
    backward=backward,
    params=[
        BondLength(1.04, label='r_NH'),
        Parameter(1.2, label='a_Hx', unit='rad')
    ],
    elem=elem,
    units='A',
)

outfile = f'{basedir}/relax.xyz'
xyz = XyzGeometry(suffix=outfile)
structure_relax = xyz.load_or_relax(
    path='./',
    relax_func=relax_pyscf,
    structure=structure_init,
    xc='pbe',
    outfile=outfile,
)
