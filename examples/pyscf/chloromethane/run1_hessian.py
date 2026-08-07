#!/usr/bin/env python

from stalk import ParameterHessian

from params import pes_pbe
from run0_relax import structure_relax


hessian_dir = 'hessian/'
hessian = ParameterHessian(structure=structure_relax)
if not hessian.load(hessian_dir):
    hessian.compute_fdiff(pes=pes_pbe, path=hessian_dir, dp=0.01)
# end if
print(hessian)
