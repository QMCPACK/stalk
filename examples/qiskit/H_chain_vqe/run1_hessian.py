#!/usr/bin/env python3

from stalk import ParameterHessian

from params import vqe_pes
from run0_relax import s_relax, directory


interactive = __name__ == '__main__'
hessian_dir = f'{directory}hessian/'
hessian = ParameterHessian(structure=s_relax)
if not hessian.load_hessian(hessian_dir):
    hessian.compute_fdiff(pes=vqe_pes, path=hessian_dir, dp=0.01)
# end if
if interactive:
    print(hessian)
# end if
