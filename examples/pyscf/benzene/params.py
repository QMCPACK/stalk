#!/usr/bin/env python3

from numpy import array, sin, pi, cos, ndarray

from pyscf import dft
from pyscf import gto
from pyscf.geomopt.geometric_solver import optimize
from pyscf.gto.mole import tofile

from stalk import ParameterStructure
from stalk import PesFunction
from stalk import BondLength


# Forward mapping: produce parameter values from an array of atomic positions
def forward(pos: ndarray):
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    C0 = pos[0]
    C1 = pos[1]
    C2 = pos[2]
    C3 = pos[3]
    C4 = pos[4]
    C5 = pos[5]
    H0 = pos[6]
    H1 = pos[7]
    H2 = pos[8]
    H3 = pos[9]
    H4 = pos[10]
    H5 = pos[11]

    # for redundancy, calculate mean bond lengths
    # 0) from neighboring C-atoms
    r_CC = BondLength([
        (C0, C1),
        (C1, C2),
        (C2, C3),
        (C3, C4),
        (C4, C5),
        (C5, C0)
    ], label='r_CC', unit='B', tol=1e-4)
    # 1) from corresponding H-atoms
    r_CH = BondLength([
        (C0, H0),
        (C1, H1),
        (C2, H2),
        (C3, H3),
        (C4, H4),
        (C5, H5)
    ], label='r_CH', unit='B', tol=1e-4)
    params = [r_CC, r_CH]
    return params
# end def


# Backward mapping: produce array of atomic positions from parameters
def backward(params: ndarray):
    r_CC = params[0]
    r_CH = params[1]
    # place atoms on a hexagon in the xy-directions
    hex_xy = array([[cos(3 * pi / 6), sin(3 * pi / 6), 0.],
                    [cos(5 * pi / 6), sin(5 * pi / 6), 0.],
                    [cos(7 * pi / 6), sin(7 * pi / 6), 0.],
                    [cos(9 * pi / 6), sin(9 * pi / 6), 0.],
                    [cos(11 * pi / 6), sin(11 * pi / 6), 0.],
                    [cos(13 * pi / 6), sin(13 * pi / 6), 0.]])
    # C-atoms are one C-C length apart from origin
    pos_C = hex_xy * r_CC
    # H-atoms one C-H length apart from C-atoms
    pos_H = hex_xy * (r_CC + r_CH)
    pos = array([pos_C, pos_H])
    return pos
# end def


def kernel_pyscf(structure: ParameterStructure, xc='pbe'):
    atom = []
    for el, pos in zip(structure.elem, structure.pos):
        atom.append([el, tuple(pos)])
    # end for
    mol = gto.Mole()
    mol.atom = atom
    mol.verbose = 3
    mol.basis = 'ccpvdz'
    mol.unit = 'B'
    mol.ecp = 'ccecp'
    mol.charge = 0
    mol.spin = 0
    mol.symmetry = False
    mol.cart = True
    mol.build()

    mf = dft.RKS(mol)
    mf.xc = xc
    return mf
# end def


def relax_pyscf(structure: ParameterStructure, outfile='relax.xyz'):
    mf = kernel_pyscf(structure=structure)
    mf.kernel()
    mol_eq = optimize(mf, maxsteps=100, constraints='benzene_constraints.txt')
    # Write to external file
    tofile(mol_eq, outfile, format='xyz')
# end def


def pes_pyscf(structure: ParameterStructure, **kwargs):
    print(f'Computing: {structure.label}')
    mf = kernel_pyscf(structure=structure)
    e_scf = mf.kernel()
    return e_scf, 0.0
# end def


pes = PesFunction(pes_pyscf)
