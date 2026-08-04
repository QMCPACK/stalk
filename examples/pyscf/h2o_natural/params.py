#!/usr/bin/env python3

from numpy import array, sin, cos, ndarray, pi

from pyscf import dft
from pyscf import gto
from pyscf.geomopt.geometric_solver import optimize
from pyscf.gto.mole import tofile

from stalk import ParameterStructure
from stalk import PesFunction
from stalk import BondLength
from stalk import BondAngle


# Natural forward mapping using bond lengths and angles
def forward(pos: ndarray) -> ndarray:
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    O0 = pos[0]
    H0 = pos[1]
    H1 = pos[2]

    # for redundancy, calculate mean bond lengths
    r = BondLength([
        (O0, H0),
        (O0, H1)
    ])
    a = BondAngle((H0, O0, H1), unit='rad')
    params = [r, a]
    return params
# end def


# Backward mapping: produce array of atomic positions from parameters
def backward(params: ndarray) -> ndarray:
    r = params[0]
    # Transform bond angle to triangular angle
    a = (pi - params[1]) / 2
    # place O at origin, and Hs in xy plane
    O0 = [0.0, 0.0, 0.0]
    H1 = [r * sin(a), r * cos(a), 0.0]
    H2 = [r * sin(a), -r * cos(a), 0.0]
    pos = array([O0, H1, H2])
    return pos
# end def


def kernel_pyscf(structure: ParameterStructure):
    atom = []
    for el, pos in zip(structure.elem, structure.pos):
        atom.append([el, tuple(pos)])
    # end for
    mol = gto.Mole()
    mol.atom = atom
    mol.verbose = 1
    mol.basis = 'ccpvdz'
    mol.unit = 'A'
    mol.ecp = 'ccecp'
    mol.charge = 0
    mol.spin = 0
    mol.symmetry = False
    mol.cart = True
    mol.build()

    mf = dft.RKS(mol)
    return mf
# end def


def relax_pyscf(structure: ParameterStructure, outfile='relax.xyz', xc='pbe'):
    mf = kernel_pyscf(structure=structure)
    mf.xc = xc
    mf.kernel()
    mol_eq = optimize(mf, maxsteps=100, constraints='h2o_constraints.txt')
    # Write to external file
    tofile(mol_eq, outfile, format='xyz')
# end def


def pes_pyscf(structure: ParameterStructure, xc='pbe', **kwargs):
    mf = kernel_pyscf(structure=structure)
    mf.xc = xc
    e_scf = mf.kernel()
    return e_scf, 0.0
# end def


# Treat a collection PESs based on alternative XC functionals
xcs = ['pbe', 'b3lyp', 'lda']
# Colors to go with different functionals for plotting
colors = ['tab:blue', 'tab:orange', 'tab:green']
pes_dict = {}
co_dict = {}
for xc, color in zip(xcs, colors):
    co_dict[xc] = color
    pes_dict[xc] = PesFunction(pes_pyscf, xc=xc, create_files=True)
# end for
