#!/usr/bin/env python

from matplotlib import pyplot as plt

from stalk import Surrogate

from params import pes_pyscf
from run1_hessian import hessian

interactive = __name__ == "__main__"

surrogate_file = 'surrogate.p'
surrogate = Surrogate(
    path='surrogate/',
    fit_kind='pf3',
    load=surrogate_file,
    structure=hessian.structure,
    hessian=hessian,
    pes=pes_pyscf,
    window_frac=0.5,  # maximum displacement relative to Lambda of each direction
    M=15,  # number of points per direction to sample
    interactive=interactive
)

epsilon_p = [0.02, 0.02]
surrogate.optimize(
    epsilon_p=epsilon_p,
    fit_kind='pf3',
    M=7,
    N=400,
    noise_frac=0.05,
    reoptimize=False,
    write=surrogate_file,
)

if interactive:
    print(surrogate)
    surrogate.plot_error_surfaces()
    plt.show()
# end if
