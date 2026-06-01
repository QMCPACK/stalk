#!/usr/bin/env python3

import numpy as np

from stalk import ParameterHessian

from params import vqe_pes
from run0_relax import s_relax, directory


interactive = __name__ == '__main__'

# Step 1: obtain parametric subspace Hessian in the LDA PES (or load)
hessian = ParameterHessian(structure=s_relax)
# Disable soft directions with low condition number
hessian_file = f'{directory}hessian.dat'
try:
    hessian_array = np.loadtxt(hessian_file, ndmin=2)
    hessian.init_hessian_array(hessian_array)
    print(f'Loaded Hessian from: {hessian_file}')
except FileNotFoundError:
    print('Computing Hessian with finite-difference method:')
    hessian.compute_fdiff(pes=vqe_pes)
    np.savetxt(hessian_file, hessian.hessian)
# end try
if interactive:
    print(hessian)
# end if
