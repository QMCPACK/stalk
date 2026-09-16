#!/usr/bin/env python3

from numpy import array

from stalk import ParameterStructure
from stalk import XyzGeometry
from stalk import BondLength
from stalk import Parameter

from params import forward, backward


# Let us initiate a ParameterStructure object that implements the parametric mappings
params_init = array([
    BondLength(1.8, label='r_CCl'),
    Parameter(0.2, label='z_CH'),
    Parameter(1.0, label='xy_CH'),
])
elem = ['C'] + ['Cl'] + 3 * ['H']
structure = ParameterStructure(
    forward=forward,
    backward=backward,
    params=params_init,
    elem=elem,
    units='A',
    label='init'
)

xyz_r = XyzGeometry(suffix='relax.xyz')
xyz_i = XyzGeometry(suffix='structure.xyz')
relaxdir = 'relax'
try:
    # Try to load the file from 'relax/relax.xyz'
    geom = xyz_r.load(relaxdir)
except FileNotFoundError:
    # If not found, create the initial structure and instruct to relax it
    xyz_i.write(structure, relaxdir)
    print(f'Wrote initial structure to {relaxdir}/structure.xyz')
    print(f'Next, run "python3 relax_geometry.py {relaxdir}/structure.xyz" to relax.')
    exit(0)
# end try
print('Initial params:')
print(structure.params)
print('Relaxed structure:')
new_params = structure.map_forward(geom.get_pos())
structure_relax = structure.copy(params=new_params)
print(structure_relax)
