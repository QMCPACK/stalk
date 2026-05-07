#!/usr/bin/env python3

from numpy import array, sin, cos, ndarray, pi
from pyscf import dft, gto

from stalk import FilesPes
from stalk import ParameterStructure
from stalk import BondLength
from stalk import BondAngle
from stalk import Parameter


# Natural forward mapping using bond lengths and angles
def forward_natural(pos: ndarray):
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    C = pos[0]
    Cl = pos[1]
    H0 = pos[2]
    H1 = pos[3]
    H2 = pos[4]

    r_CCl = BondLength((C, Cl), label='r_CCl')
    r_CH = BondLength([
        (C, H0),
        (C, H1),
        (C, H2),
    ], label='r_CH', tol=1e-4)
    a = BondAngle([
        (H0, C, H1),
        (H1, C, H2),
        (H2, C, H0),
    ], unit='rad', tol=1e-4)
    params = [r_CCl, r_CH, a]
    return params
# end def


# Auxiliary parameter mapping using z, xy distances to set H atoms
def forward(pos: ndarray):
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    C = pos[0]
    Cl = pos[1]
    H0 = pos[2]
    H1 = pos[3]
    H2 = pos[4]

    r_CCl = BondLength((C, Cl), label='r_CCl')
    z_CH = Parameter([
        (H0 - C)[2],
        (H1 - C)[2],
        (H2 - C)[2],
    ], label='z_CH', tol=1e-4)
    xy_CH = Parameter([
        ((C[0] - H0[0])**2 + (C[1] - H0[1])**2)**0.5,
        ((C[0] - H1[0])**2 + (C[1] - H1[1])**2)**0.5,
        ((C[0] - H2[0])**2 + (C[1] - H2[1])**2)**0.5,
    ], label='xy_CH', tol=1e-3)
    params = [r_CCl, z_CH, xy_CH]
    return params
# end def


# Backward mapping: produce array of atomic positions from parameters
def backward(params: ndarray):
    r_CCl = params[0]
    z = params[1]
    xy = params[2]
    aux_ang = 2 * pi / 3

    # place atoms on an equilateral triangle in the xy-directions
    C = [0.0, 0.0, 0.0]
    Cl = [0.0, 0.0, -r_CCl]
    H0 = [xy, 0.0, z]
    H1 = [xy * cos(aux_ang), xy * sin(aux_ang), z]
    H2 = [xy * cos(-aux_ang), xy * sin(-aux_ang), z]
    pos = array([C, Cl, H0, H1, H2])
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
    mol.unit = 'A'
    mol.ecp = 'ccecp'
    mol.charge = 0
    mol.spin = 0
    mol.symmetry = False
    mol.build()

    mf = dft.RKS(mol)
    mf.xc = xc
    return mf
# end def


pes_xyz = FilesPes()
