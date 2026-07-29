#!/usr/bin/env python3

from matplotlib import pyplot as plt

from stalk import LineSearchIteration

from params import pes_dmc
from run2_surrogate import surrogate

interactive = __name__ == "__main__"

# Run a snapshot job to sample effective variance w.r.t relative DMC samples
# Choosing higher tiling makes for higher accuracy and higher cost
pes_dmc.args['tile_opt'] = 4
pes_dmc.args['twist_grid'] = (2, 2, 2)
# Rescale the result appropriately
pes_dmc.loader.scale = 4
var_eff = pes_dmc.get_var_eff(
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
    pes=pes_dmc,
    var_eff=var_eff,
)
# Propagate the parallel line-search (compute values, analyze, then move on) 4 times
for i in range(2):
    dmc_ls.propagate(
        pes_dmc,
        i,
        interactive=interactive,
        dep_jobs=dep_jobs
    )
    if interactive:
        print(dmc_ls)
        dmc_ls.pls(i).plot()
        plt.show()
    # end if
# end for
# Evaluate the latest eqm structure
dmc_ls.pls().evaluate_eqm(
    pes_dmc,
    interactive=interactive,
    dep_jobs=dep_jobs
)

# Print the line-search performance
if interactive:
    print(dmc_ls)
    dmc_ls.plot(target=surrogate.structure)
    plt.show()
# end if
