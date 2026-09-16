#!/usr/bin/env python3

from stalk import ParameterHessian

from params import pes_pyscf
from run0_relax import structure_relax

interactive = __name__ == "__main__"
hessian = ParameterHessian(structure=structure_relax)
hessian.compute_fdiff(
    pes=pes_pyscf,
    path='hessian',
    dp=0.01,
    interactive=interactive,
)
if interactive:
    print(hessian)
# end if
