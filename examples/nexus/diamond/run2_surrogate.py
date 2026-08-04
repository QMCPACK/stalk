#!/usr/bin/env python

from matplotlib import pyplot as plt

from stalk import Surrogate

from params import pes_pwscf
from run1_hessian import hessian

interactive = __name__ == "__main__"

surrogate = Surrogate(
    path='surrogate/',
    fit_kind='pf3',
    structure=hessian.structure,
    hessian=hessian,
    window_frac=0.3,  # maximum displacement relative to Lambda of each direction
    M=15,  # number of points per direction to sample
)
surrogate.evaluate(pes=pes_pwscf, interactive=interactive)

epsilon_p = [0.01]
surrogate.optimize(
    epsilon_p=epsilon_p,
    fit_kind='pf3',
    M=7,
    N=500,
    reoptimize=False,
)

if interactive:
    print(surrogate)
    surrogate.plot_error_surfaces()
    plt.show()
# end if
