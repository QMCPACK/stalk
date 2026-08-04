#!/usr/bin/env python

import numpy as np

from stalk import Surrogate

from params import pes_dict
from run1_hessian import hessians


# Treat a collection surrogates based on alternative XC functionals
surrogates = {}
epsilon_p = np.array([0.01, 0.01])

for xc, pes in pes_dict.items():
    surrogate_dir = f'surrogate/{xc}'
    # Characterize PES
    surrogate = Surrogate(
        path=surrogate_dir,
        fit_kind='pf3',
        structure=hessians[xc].structure,
        hessian=hessians[xc],
        window_frac=0.3,
        M=15
    )
    surrogate.evaluate(pes)

    # Optimize to tolerances
    surrogate.optimize(
        epsilon_p=epsilon_p / 2,
        fit_kind='pf3',
        M=7,
        N=400,
        reoptimize=False,
    )
    surrogate.optimize(epsilon_p=epsilon_p)
    print(f'Surrogate model ({xc})')
    print(surrogate)
    surrogates[xc] = surrogate
# end for
