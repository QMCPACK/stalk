#!/usr/bin/env python

from stalk import ParameterHessian

from params import pes_pwscf
from run0_relax import structure_relax

interactive = __name__ == "__main__"
hessian_dir = 'hessian/'
hessian = ParameterHessian(structure=structure_relax)
if not hessian.load(hessian_dir):
    hessian.compute_fdiff(
        pes=pes_pwscf,
        path=hessian_dir,
        dp=0.001,
        interactive=interactive,
    )
# end if
if interactive:
    print(hessian)
# end if
