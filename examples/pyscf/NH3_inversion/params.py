#!/usr/bin/env python3

from numpy import sin, cos, ndarray
from os import makedirs

from pyscf import dft
from pyscf import gto
from pyscf import grad
from pyscf.geomopt.geometric_solver import optimize
from pyscf.gto.mole import tofile

from ase.calculators.calculator import Calculator, all_changes
from ase import Atoms
from ase.constraints import FixAtoms

from stalk import angle
from stalk import rotate_2d
from stalk import ParameterStructure
from stalk import PesFunction
from stalk import BondLength
from stalk import Parameter


# Forward mapping: produce parameter values from an array of atomic positions
def forward(pos: ndarray):
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    N0 = pos[0]
    H0 = pos[1]
    H1 = pos[2]
    H2 = pos[3]

    # for redundancy, calculate mean bond lengths
    r = BondLength([
        (N0, H0),
        (N0, H1),
        (N0, H2),
    ], label='r_NH', tol=1e-4)
    # Calculate angle between x axis and H -> should be around 90 degrees
    x = [1.0, 0.0, 0.0]
    d0 = H0 - N0
    d1 = H1 - N0
    d2 = H2 - N0
    a = Parameter([
        angle(d0, x, units='rad'),
        angle(d1, x, units='rad'),
        angle(d2, x, units='rad'),
    ], label='a_Hx', tol=1e-4)
    params = [r, a]
    return params
# end def


# Backward mapping: produce array of atomic positions from parameters
def backward(params: ndarray):
    r = params[0]
    a = params[1]
    N0 = [0.0, 0.0, 0.0]
    x = r * cos(a)
    xy = [r * sin(a), 0.0]
    H0 = [x, *rotate_2d(xy, 0, units='ang')]
    H1 = [x, *rotate_2d(xy, 120, units='ang')]
    H2 = [x, *rotate_2d(xy, -120, units='ang')]
    pos = [N0, H0, H1, H2]
    return pos
# end def


def kernel_pyscf(
    positions,
    elem,
    xc,
    basis,
    ecp,
):
    atom = []
    for el, pos in zip(elem, positions):
        atom.append([el, tuple(pos)])
    # end for
    mol = gto.Mole()
    mol.atom = atom
    mol.verbose = 3
    mol.basis = basis
    mol.unit = 'A'
    mol.ecp = ecp
    mol.charge = 0
    mol.spin = 0
    mol.symmetry = False
    mol.cart = True
    mol.build()

    mf = dft.RKS(mol)
    mf.xc = xc
    return mf
# end def


def relax_pyscf(
    structure: ParameterStructure,
    xc='pbe',
    basis='ccecpccpvdz',
    ecp='ccecp',
    **kwargs
):
    mf = kernel_pyscf(
        structure.pos,
        structure.elem,
        xc=xc,
        basis=basis,
        ecp=ecp,
    )
    mf.kernel()
    mol_eq = optimize(mf, maxsteps=100, constraints='nh3_constraints.txt')
    # Write to external file
    makedirs(structure.path, exist_ok=True)
    tofile(mol_eq, f'{structure.path}/relax.xyz', format='xyz')
    e_scf = mf.e_tot
    return e_scf, 0.0
# end def


def pes_pyscf(
    structure: ParameterStructure,
    xc='pbe',
    basis='ccecpccpvdz',
    ecp='ccecp',
    **kwargs
):
    mf = kernel_pyscf(
        structure.pos,
        structure.elem,
        xc=xc,
        basis=basis,
        ecp=ecp,
    )
    e_scf = mf.kernel()
    return e_scf, 0.0
# end def


# PBE PES
pes_pbe = PesFunction(pes_pyscf, create_files=True, xc='pbe', basis='ccecpccpvdz', ecp='ccecp')
relax_pbe = PesFunction(relax_pyscf, create_files=True, xc='pbe', basis='ccecpccpvdz', ecp='ccecp')
# LDA PES
pes_lda = PesFunction(pes_pyscf, create_files=True, xc='lda', basis='ccecpccpvdz', ecp='ccecp')
relax_lda = PesFunction(relax_pyscf, create_files=True, xc='lda', basis='ccecpccpvdz', ecp='ccecp')


# To run a NEB calculation with ASE, we need to define a calculator that can compute
# energies and forces. Here we define a simple wrapper around the PySCF calculations.
class PySCF(Calculator):
    implemented_properties = ['energy', 'forces']
    pes_args = {}

    def __init__(
        self,
        pes_args={},
        **kwargs
    ):
        Calculator.__init__(self, **kwargs)
        self.pes_args = pes_args
    # end def

    def calculate(
        self,
        atoms=None,
        properties=['energy'],
        system_changes=all_changes
    ):
        Calculator.calculate(self, atoms, properties, system_changes)

        # Get positions and symbols from the Atoms object
        positions = self.atoms.get_positions()
        elem = self.atoms.get_chemical_symbols()

        # Calculate energy and forces here using your custom method
        energy = self.compute_energy(positions, elem)
        forces = self.compute_forces(positions, elem)

        # Set results to be used by ASE
        self.results = {'energy': energy, 'forces': forces}
    # end def

    def compute_energy(self, positions, elem):
        mf = kernel_pyscf(positions, elem, **self.pes_args)
        return mf.kernel()
    # end def

    def compute_forces(self, positions, elem):
        mf = kernel_pyscf(positions, elem, **self.pes_args)
        mf.kernel()
        mf_grad = grad.RKS(mf)
        forces = -mf_grad.kernel()
        return forces
    # end def

# end class


def neb_image(structure: ParameterStructure, pes_args={}):
    atoms = Atoms(structure.elem, positions=structure.pos)
    atoms.calc = PySCF(pes_args)
    # Fix N atom
    atoms.set_constraint(FixAtoms(indices=[0]))
    return atoms
# end def
