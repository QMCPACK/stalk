#!/usr/bin/env python3

from matplotlib import pyplot as plt

from stalk import LineSearchIteration

from params import pes_pwscf
from run2_surrogate import surrogate

interactive = __name__ == "__main__"

shifted_structure = surrogate.structure.copy()
shifted_structure.shift_params([0.1])
# Then generate line-search iteration object based on the shifted surrogate
srg_ls = LineSearchIteration(
    surrogate=surrogate,
    structure=shifted_structure,
    path='srg_ls',
    pes=pes_pwscf,
)
# Propagate the parallel line-search (compute values, analyze, then move on) 4 times
#   add_sigma = True means that target errorbars are used to simulate random noise
for i in range(4):
    srg_ls.propagate(i, add_sigma=True, interactive=interactive)
    if interactive:
        print(srg_ls)
        srg_ls.pls(i).plot()
        plt.show()
    # end if
# end for
# Evaluate the latest eqm structure
srg_ls.pls().evaluate_eqm(add_sigma=True, interactive=interactive)

# Print the line-search performance
if interactive:
    print(srg_ls)
    srg_ls.plot(target=surrogate.structure)
    plt.show()
# end if
