#!/usr/bin/env python

from matplotlib import pyplot as plt

from stalk import Surrogate

from params import pes_xyz
from run1_hessian import hessian


surrogate = Surrogate(
    path='surrogate/',
    fit_kind='pf3',
    structure=hessian.structure,
    hessian=hessian,
    window_frac=0.2,
    M=15
)
surrogate.evaluate(pes=pes_xyz)

epsilon_p = [0.02, 0.02, 0.02]
surrogate.optimize(
    epsilon_p=epsilon_p,
    fit_kind='pf3',
    M=7,
    N=400,
    reoptimize=False,
)

if __name__ == "__main__":
    print(surrogate)
    surrogate.plot()
    plt.legend()
    plt.show()
# end if
