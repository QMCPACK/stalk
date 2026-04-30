#!/usr/bin/env python

from os import makedirs
import numpy as np

from stalk import Surrogate

from params import pes_dict
from run1_hessian import hessians


surrogate_dir = 'surrogate/'
makedirs(surrogate_dir, exist_ok=True)

# Treat a collection surrogates based on alternative XC functionals
surrogates = {}
epsilon_p = np.array([0.01, 0.01])

for xc, pes in pes_dict.items():
    surrogate_file = f'{xc}.p'
    # Characterize PES
    surrogate = Surrogate(
        path=surrogate_dir,
        fit_kind='pf3',
        load=surrogate_file,
        structure=hessians[xc].structure,
        hessian=hessians[xc],
        pes=pes,
        window_frac=0.3,
        M=15
    )
    surrogate.bracket_target_biases()

    # Optimize to tolerances
    surrogate.optimize(
        epsilon_p=epsilon_p / 2,
        fit_kind='pf3',
        M=7,
        N=400,
        reoptimize=False,
        write=surrogate_file,
    )
    surrogate.optimize(epsilon_p=epsilon_p)
    surrogate.write_to_disk(surrogate_file, overwrite=True)
    print(f'Surrogate model ({xc})')
    print(surrogate)
    surrogates[xc] = surrogate
# end for
