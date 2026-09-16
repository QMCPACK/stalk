#!/usr/bin/env python3

from stalk import ParameterHessian

from params import pes_xyz
from run0_relax import structure_relax


hessian = ParameterHessian(structure=structure_relax)
hessian.compute_fdiff(
    pes=pes_xyz,
    path='hessian',
    dp=0.01
)
print(hessian)
