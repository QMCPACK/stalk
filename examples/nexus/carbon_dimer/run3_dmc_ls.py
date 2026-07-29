#!/usr/bin/env python3

from matplotlib import pyplot as plt

from stalk import LineSearch

from params import pes_dmc
from run2_surrogate import surrogate
from stalk.params.effective_variance_map import EffectiveVarianceMap

interactive = __name__ == "__main__"

structure_qmc = surrogate.structure.copy(label='eqm')
# Run a snapshot job to sample effective variance w.r.t relative DMC samples
var_eff5 = pes_dmc.get_var_eff(
    structure_qmc,
    path='dmc/var_eff/5',
    samples=5,
    interactive=interactive
)
var_eff10 = pes_dmc.get_var_eff(
    structure_qmc.copy(),
    path='dmc/var_eff/10',
    samples=10,
    interactive=interactive
)
var_eff_map = EffectiveVarianceMap(structure_qmc, var_eff5 + var_eff10)

# Finally, perform line-search iteration based on surrogate settings and DMC PES
dmc_ls = LineSearch(**surrogate.to_settings())
dmc_ls.evaluate(
    pes_dmc,
    path='dmc_ls',
    dep_jobs=structure_qmc.jobs,
    interactive=interactive,
    var_eff_map=var_eff_map,
)

# Print the line-search performance
if interactive:
    print(dmc_ls)
    dmc_ls.plot()
    plt.show()
# end if
