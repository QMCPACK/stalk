#!/usr/bin/env python3

import argparse
from pathlib import Path
from pyscf.geomopt.geometric_solver import optimize
from pyscf.gto.mole import tofile

from stalk import XyzGeometry
from stalk import ParameterStructure

from params import kernel_pyscf

parser = argparse.ArgumentParser(description='Relax an initial XYZ geometry')
parser.add_argument('filename', help='Initial file')
parser.add_argument('--xc', default='pbe', help='XC functional')


relaxfile = 'relax.xyz'

if __name__ == '__main__':
    args = parser.parse_args()
    xyzfile = Path(args.filename)
    if xyzfile.name == relaxfile:
        print('Cannot treat files named "relax.xyz" as they would be overwritten.')
        exit(0)
    # end if
    xyzpath = xyzfile.parent

    # Load the initial geometry
    geom = XyzGeometry(suffix=xyzfile.name).load(xyzpath)
    structure = ParameterStructure(
        pos=geom.get_pos(),
        elem=geom.get_elem(),
    )

    # Calculate the energy using PySCF kernel
    xc = args.xc
    print(f'Relaxing: {xyzfile} ({xc})')
    mf = kernel_pyscf(structure=structure, xc=xc)
    mf.kernel()
    mol_eq = optimize(mf, maxsteps=100, constraints='chloromethane_constraints.txt')

    # Write to external file
    relaxfile = xyzpath / 'relax.xyz'
    tofile(mol_eq, relaxfile, format='xyz')
# end if
