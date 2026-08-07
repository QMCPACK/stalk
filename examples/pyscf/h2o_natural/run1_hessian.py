#!/usr/bin/env python

from stalk import ParameterHessian

from params import pes_dict
from run0_relax import structure_relax

# Treat a collection parameter hessians based on alternative XC functionals
hessians = {}
for xc, pes in pes_dict.items():
    hessian = ParameterHessian(structure=structure_relax[xc])
    hessian_dir = f'hessian/{xc}'
    if not hessian.load(hessian_dir):
        hessian.compute_fdiff(pes=pes, path=hessian_dir, dp=0.01)
    # end if
    print(hessian)
    hessians[xc] = hessian
# end for
