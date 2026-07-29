#!/usr/bin/env python3

from numpy import array
from pyscf.geomopt.geometric_solver import optimize
from pyscf.gto.mole import tofile

from stalk import ParameterStructure
from stalk import XyzGeometry
from stalk import BondLength
from stalk import Parameter

from params import forward, backward, kernel_pyscf


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

outfile = 'relax.xyz'
try:
    geom = XyzGeometry().load(outfile)
except FileNotFoundError:
    mf = kernel_pyscf(structure=structure, xc='pbe')
    mf.kernel()
    mol_eq = optimize(mf, maxsteps=100, constraints='chloromethane_constraints.txt')
    # Write to external file
    tofile(mol_eq, outfile, format='xyz')
    geom = XyzGeometry().load(outfile)
# end try
new_params = structure.map_forward(geom.get_pos())
print('Initial params:')
print(structure.params)
print('Relaxed structure:')
structure_relax = structure.copy(params=new_params)
print(structure_relax)
