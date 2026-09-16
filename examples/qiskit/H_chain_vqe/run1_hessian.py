#!/usr/bin/env python3

from stalk import ParameterHessian

from params import vqe_pes
from run0_relax import s_relax, directory


interactive = __name__ == '__main__'
hessian = ParameterHessian(structure=s_relax)
hessian.compute_fdiff(
    pes=vqe_pes,
    path=f'{directory}hessian',
    dp=0.01
)
if interactive:
    print(hessian)
# end if
