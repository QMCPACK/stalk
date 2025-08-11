#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt

from stalk import LineSearchIteration

from params import pes
from run2_surrogate import surrogate


shifted_structure = surrogate.structure.copy()
shifted_structure.shift_params(0.1 * np.random.randn(10))
# Then generate line-search iteration object based on the shifted surrogate
srg_ls = LineSearchIteration(
    surrogate=surrogate,
    structure=shifted_structure,
    path='srg_ls',
    pes=pes,
)
# Propagate the parallel line-search (compute values, analyze, then move on) 4 times
#   add_sigma = True means that target errorbars are used to simulate random noise
for i in range(4):
    srg_ls.propagate(i, add_sigma=True)
# end for
# Evaluate the latest eqm structure
srg_ls.pls().evaluate_eqm(add_sigma=True)

if __name__ == '__main__':
    # Print the line-search performance
    print(srg_ls)
    srg_ls.plot(target=surrogate)
    plt.show()
# end if
