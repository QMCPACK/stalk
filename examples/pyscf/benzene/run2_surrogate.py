#!/usr/bin/env python

from stalk import Surrogate

from params import pes
from run1_hessian import hessian


surrogate = Surrogate(
    path='surrogate/',
    fit_kind='pf3',
    structure=hessian.structure,
    hessian=hessian,
    window_frac=0.5,  # maximum displacement relative to Lambda of each direction
    M=15,  # number of points per direction to sample
)
surrogate.evaluate(pes=pes)

epsilon_p = [0.02, 0.02]
surrogate.optimize(
    epsilon_p=epsilon_p,
    fit_kind='pf3',
    M=7,
    N=400,
    reoptimize=False,
    logger=3,
)
print(surrogate)
