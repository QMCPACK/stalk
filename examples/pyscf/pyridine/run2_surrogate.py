#!/usr/bin/env python

from stalk import Surrogate

from params import pes
from run1_hessian import hessian


surrogate_file = 'surrogate.p'
surrogate = Surrogate(
    path='surrogate/',
    fit_kind='pf3',
    load=surrogate_file,
    structure=hessian.structure,
    hessian=hessian,
    pes=pes,
    window_frac=0.5,  # maximum displacement relative to Lambda of each direction
    M=15  # number of points per direction to sample
)
surrogate.write_to_disk('pes.p')

epsilon_p = 10 * [0.01]
surrogate.optimize(
    epsilon_p=epsilon_p,
    fit_kind='pf3',
    M=7,
    N=400,
    reoptimize=False,
    write=surrogate_file,
)

if __name__ == '__main__':
    print(surrogate)
    print(surrogate.sigma_opt)
    surrogate.plot()
# end if
