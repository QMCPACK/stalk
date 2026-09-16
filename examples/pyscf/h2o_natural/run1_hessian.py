#!/usr/bin/env python3

from stalk import ParameterHessian

from params import pes_dict
from run0_relax import structure_relax

# Treat a collection parameter hessians based on alternative XC functionals
hessians = {}
for xc, pes in pes_dict.items():
    hessian = ParameterHessian(structure=structure_relax[xc])
    hessian.compute_fdiff(
        pes=pes,
        path=f'{xc}/hessian',
        dp=0.01
    )
    hessians[xc] = hessian
    print(hessian)
# end for
