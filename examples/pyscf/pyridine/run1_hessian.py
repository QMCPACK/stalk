#!/usr/bin/env python3

from stalk import ParameterHessian

from params import pes_pbe
from run0_relax import structure_relax


hessian = ParameterHessian(structure=structure_relax)
hessian.compute_fdiff(
    pes=pes_pbe,
    path='hessian',
    dp=0.01
)
print(hessian)
