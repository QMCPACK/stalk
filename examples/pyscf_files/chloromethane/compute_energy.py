#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np

from stalk import XyzGeometry
from stalk import PesResult

from params import kernel_pyscf
from run0_relax import structure_relax

parser = argparse.ArgumentParser(description='Compute energies from XYZ structures')
parser.add_argument('filename', nargs='+', help='Structure files')
parser.add_argument('--xc', default='pbe', help='XC functional')


if __name__ == '__main__':
    args = parser.parse_args()
    for fname in args.filename:
        if not fname.endswith('structure.xyz'):
            print(f'Skipping {fname}')
            continue
        # end if
        efile = fname.replace('structure.xyz', 'energy.dat')
        sfile = fname.replace('structure.xyz', 'sigma.dat')

        geom = XyzGeometry(suffix='').load(fname)
        structure = structure_relax.copy(pos=geom.get_pos())

        xc = args.xc
        print(f'Computing: {fname} ({xc})')
        mf = kernel_pyscf(structure=structure, xc=xc)
        e_scf = mf.kernel()
        energy = PesResult(e_scf)

        if Path(sfile).exists():
            sigma = float(np.loadtxt(sfile))
            energy.add_sigma(sigma)
        # end if
        np.savetxt(efile, [energy.value, energy.error])
    # end for
# end if
