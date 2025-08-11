#!/usr/bin/env python3

from matplotlib import pyplot as plt

from stalk import LineSearchIteration

from params import forward_natural, pes_xyz
from run2_surrogate import surrogate


shifted_structure = surrogate.structure.copy()
shifted_structure.shift_params([0.1, -0.1, 0.1])
# Then generate line-search iteration object based on the shifted surrogate
srg_ls = LineSearchIteration(
    surrogate=surrogate,
    structure=shifted_structure,
    path='srg_ls',
    pes=pes_xyz,
)
# Propagate the parallel line-search (compute values, analyze, then move on) 4 times
#   add_sigma = True means that target errorbars are used to simulate random noise
for i in range(4):
    srg_ls.propagate(i, add_sigma=True)
# end for
# Evaluate the latest eqm structure
srg_ls.pls().evaluate_eqm(add_sigma=True)

if __name__ == "__main__":
    # Print the line-search performance
    print(srg_ls)
    print('Original energy and params:')
    print(surrogate.ls(0).target_settings.target.y0, surrogate.structure.params)
    # Remap to natural parameters
    params_natural = surrogate.structure.remap_forward(forward_natural)
    print(f'C-Cl bond length: {params_natural[0]}')
    print(f'C-H bond length: {params_natural[1]}')
    print(f'C-H bond angle: {params_natural[2]}')
    srg_ls.plot()
    plt.show()
# end if
