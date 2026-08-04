#!/usr/bin/env python3

from numpy import array, ndarray
from os import makedirs

from pyscf import dft
from pyscf import gto
from pyscf.geomopt.geometric_solver import optimize
from pyscf.gto.mole import tofile

from stalk import ParameterStructure
from stalk import PesFunction
from stalk import Parameter


# Forward mapping: produce parameter values from an array of atomic positions
def forward(pos: ndarray) -> ndarray:
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    N0 = pos[0]
    C0 = pos[1]
    C1 = pos[2]
    C2 = pos[3]
    C3 = pos[4]
    C4 = pos[5]
    H0 = pos[6]
    H1 = pos[7]
    H2 = pos[8]
    H3 = pos[9]
    H4 = pos[10]

    x_C04 = Parameter([
        (C0 - N0)[0],
        -(C4 - N0)[0]
    ], tol=1e-3)
    y_C04 = Parameter([
        (C0 - N0)[1],
        (C4 - N0)[1]
    ], tol=1e-3)
    x_C13 = Parameter([
        (C1 - N0)[0],
        -(C3 - N0)[0]
    ], tol=1e-3)
    y_C13 = Parameter([
        (C1 - N0)[1],
        (C3 - N0)[1]
    ], tol=1e-3)
    y_C2 = (C2 - N0)[1]
    x_H04 = Parameter([
        (H0 - N0)[0],
        -(H4 - N0)[0]
    ], tol=1e-3)
    y_H04 = Parameter([
        (H0 - N0)[1],
        (H4 - N0)[1]
    ], tol=1e-3)
    x_H13 = Parameter([
        (H1 - N0)[0],
        -(H3 - N0)[0]
    ], tol=1e-3)
    y_H13 = Parameter([
        (H1 - N0)[1],
        (H3 - N0)[1]
    ], tol=1e-3)
    y_H2 = Parameter((H2 - N0)[1])
    params = [
        x_C04,
        y_C04,
        x_C13,
        y_C13,
        y_C2,
        x_H04,
        y_H04,
        x_H13,
        y_H13,
        y_H2
    ]
    return params
# end def


# Backward mapping: produce array of atomic positions from parameters
def backward(params: ndarray) -> ndarray:
    x_C04, y_C04, x_C13, y_C13, y_C2, x_H04, y_H04, x_H13, y_H13, y_H2 = tuple(params)
    N0 = [0.0, 0.0, 0.0]
    C0 = [x_C04, y_C04, 0.0]
    C1 = [x_C13, y_C13, 0.0]
    C2 = [0.0, y_C2, 0.0]
    C3 = [-x_C13, y_C13, 0.0]
    C4 = [-x_C04, y_C04, 0.0]
    H0 = [x_H04, y_H04, 0.0]
    H1 = [x_H13, y_H13, 0.0]
    H2 = [0.0, y_H2, 0.0]
    H3 = [-x_H13, y_H13, 0.0]
    H4 = [-x_H04, y_H04, 0.0]
    pos = array([
        N0,
        C0,
        C1,
        C2,
        C3,
        C4,
        H0,
        H1,
        H2,
        H3,
        H4
    ])
    return pos
# end def


def kernel_pyscf(structure: ParameterStructure, **kwargs):
    atom = []
    for el, pos in zip(structure.elem, structure.pos):
        atom.append([el, tuple(pos)])
    # end for
    mol = gto.Mole()
    mol.atom = atom
    mol.verbose = 2
    mol.basis = 'ccecp-ccpvdz'
    mol.unit = 'A'
    mol.ecp = 'ccecp'
    mol.charge = 0
    mol.spin = 0
    mol.symmetry = False
    mol.cart = True
    mol.build()

    mf = dft.RKS(mol)
    mf.xc = 'pbe'
    return mf
# end def


def relax_pyscf(structure: ParameterStructure, **kwargs):
    mf = kernel_pyscf(structure=structure, **kwargs)
    mf.kernel()
    mol_eq = optimize(mf, maxsteps=100)
    # Write to external file
    makedirs(structure.path, exist_ok=True)
    tofile(mol_eq, f'{structure.path}/relax.xyz', format='xyz')
    return mf.e_tot, 0.0
# end def


def pes_pyscf(structure: ParameterStructure, **kwargs):
    mf = kernel_pyscf(structure=structure, **kwargs)
    e_scf = mf.kernel()
    return e_scf, 0.0
# end def


pes = PesFunction(pes_pyscf, create_files=True)
