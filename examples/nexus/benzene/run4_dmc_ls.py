#!/usr/bin/env python3

from matplotlib import pyplot as plt

from stalk import LineSearchIteration

from params import pes_dmc
from run2_surrogate import surrogate

interactive = __name__ == "__main__"

# Run a snapshot job to sample effective variance w.r.t relative DMC samples
var_eff_map = pes_dmc.get_var_eff_map(
    structure=surrogate.structure,
    path='dmc/var_eff',
    samples=10,
    interactive=interactive,
)
# Add job dependencies to recycle Jastrow
dep_jobs = surrogate.structure.jobs

# Then generate line-search iteration object based on the shifted surrogate
dmc_ls = LineSearchIteration(
    surrogate=surrogate,
    path='dmc/lsi',
    var_eff_map=var_eff_map,
)
# Propagate the parallel line-search (compute values, analyze, then move on) 4 times
#   add_sigma = True means that target errorbars are used to simulate random noise
for i in range(3):
    dmc_ls.propagate(pes_dmc, i, interactive=interactive, dep_jobs=dep_jobs)
    if interactive:
        print(dmc_ls)
        dmc_ls.pls(i).plot()
        plt.show()
    # end if
# end for

# Evaluate all
eqms = [pls.structure for pls in dmc_ls.pls_list]
pes_dmc.evaluate_all(
    eqms,
    sigmas=len(eqms) * [0.002],
    path='dmc/eqms/',
    interactive=interactive,
    dep_jobs=dep_jobs
)

# Print the line-search performance
if interactive:
    print(dmc_ls)
# end if
