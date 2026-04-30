#!/usr/bin/env python

from numpy import loadtxt, savetxt

from stalk import ParameterHessian

from params import pes
from run0_relax import structure_relax


# Disable line-search directions whose condition number is lower than this:
hessian = ParameterHessian(structure=structure_relax, disable_limit=0.05)
hessian_file = 'hessian.dat'
try:
    hessian_array = loadtxt(hessian_file)
    hessian.init_hessian_array(hessian_array)
    print(f'Loaded Hessian from: {hessian_file}')
except FileNotFoundError:
    print('Computing Hessian with finite-difference method:')
    hessian.compute_fdiff(pes=pes)
    savetxt(hessian_file, hessian.hessian)
# end try

print('Hessian:')
print(hessian)
