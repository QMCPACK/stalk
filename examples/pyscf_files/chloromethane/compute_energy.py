#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np

from stalk import XyzGeometry
from stalk import PesResult
from stalk import ParameterStructure

from params import kernel_pyscf

parser = argparse.ArgumentParser(description='Compute energies from XYZ structures')
parser.add_argument('filename', nargs='+', help='Structure files')
parser.add_argument('--xc', default='pbe', help='XC functional')


if __name__ == '__main__':
    args = parser.parse_args()
    for xyzname in args.filename:
        # Only treat .xyz files
        if not xyzname.endswith('.xyz'):
            print(f'Skipping {xyzname}')
            continue
        # end if
        xyzfile = Path(xyzname)
        xyzpath = xyzfile.parent

        # Load the geometry
        geom = XyzGeometry(suffix=xyzfile.name).load(xyzpath)
        structure = ParameterStructure(
            pos=geom.get_pos(),
            elem=geom.get_elem(),
        )

        # Calculate the energy using PySCF kernel
        xc = args.xc
        print(f'Computing: {xyzname} ({xc})')
        mf = kernel_pyscf(structure=structure, xc=xc)
        e_scf = mf.kernel()
        energy = PesResult(e_scf)

        # Add sigma to energy if it exists
        sfile = xyzpath / 'sigma.in'
        if Path(sfile).exists():
            sigma = float(np.loadtxt(sfile))
            energy.add_sigma(sigma)
        # end if

        # Write the energy to file
        efile = xyzpath / 'value.in'
        np.savetxt(efile, [energy.value, energy.error])
    # end for
# end if
