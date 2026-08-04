#!/usr/bin/env python3

from matplotlib import pyplot as plt

from stalk import LineSearchIteration

from params import forward_natural, pes_pbe
from run2_surrogate import surrogate


shifted_structure = surrogate.structure.copy()
shifted_structure.shift_params([0.1, -0.1, 0.1])
# Then generate line-search iteration object based on the shifted surrogate
srg_ls = LineSearchIteration(
    surrogate=surrogate,
    structure=shifted_structure,
    path='lsi',
)
# Propagate the parallel line-search (compute values, analyze, then move on) 4 times
#   add_sigma = True means that target errorbars are used to simulate random noise
for i in range(4):
    srg_ls.propagate(pes_pbe, i, add_sigma=True)
# end for
# Evaluate the latest eqm structure
srg_ls.pls().evaluate_eqm(add_sigma=True)

if __name__ == "__main__":
    # Print the line-search performance
    print(srg_ls)
    print('Original energy and params:')
    print(surrogate.ls(0).target_settings.target.y0, surrogate.structure.params)
    # Remap to natural parameters
    p_natural_init = srg_ls.structure_init.remap_forward(forward_natural)
    p_natural_final = srg_ls.structure_final.remap_forward(forward_natural)
    print(f'{p_natural_init[0]} -> {p_natural_final[0]}')
    print(f'{p_natural_init[1]} -> {p_natural_final[1]}')
    print(f'{p_natural_init[2]} -> {p_natural_final[2]}')
    srg_ls.plot(target=surrogate.structure)
    plt.show()
# end if
