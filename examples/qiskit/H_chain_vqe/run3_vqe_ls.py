#!/usr/bin/env python3

from matplotlib import pyplot as plt
from pathlib import Path
import numpy as np

from stalk import LineSearchIteration

from params import vqe_pes_noisy
from run2_surrogate import surrogate, directory


interactive = __name__ == '__main__'

# Step 3: Optimize VQE by line-search from a shifted position
# Shift the starting position for some challenge
structures, lsis, ntot = [], [], 5
sigma = 0.1
for n in range(ntot):
    structure = surrogate.structure.copy()
    sfile = f'{directory}/{n}_shifts_init.dat'
    if Path(sfile).exists():
        shifts = np.loadtxt(sfile)
        print(f'Loaded initial shifts from {sfile}')
    else:
        shifts = sigma * np.random.randn(len(structure))
        np.savetxt(sfile, shifts)
    # end if
    structure.shift_params(shifts)
    structures.append(structure)
    lsi = LineSearchIteration(
        surrogate=surrogate,
        structure=structure,
        path=f'{directory}STALK/lsi{n}',
    )
    lsis.append(lsi)
    # Propagate the parallel line-search (compute values, analyze, then move on) 5 times
    for i in range(5):
        lsi.propagate(vqe_pes_noisy, i)
        # Plot each line-search
        if interactive:
            lsi.pls(-2).plot()
            plt.title(f'Line-search {i} for structure {n}')
            plt.show()
        # end if
    # end for
    # Evaluate the latest eqm structure
    lsi.pls().evaluate_eqm(vqe_pes_noisy)

    if interactive:
        print(lsi)
        # Plot the convergence of parameters and energy with parallel iteration count
        lsi.plot(bundle=False, target=surrogate.structure)
        # Consider converged after 2 iterations
        lsi.transient = 2
        s_final = lsi.structure_final
        for p, (param, param_ref) in enumerate(zip(surrogate.structure.params, s_final.params)):
            diff = abs(param - param_ref)
            tolerance = surrogate.epsilon_p[p]
            print(f'  |theta{p} - ref| = {diff} < {tolerance}? {diff < tolerance}')
        # end for
        print(f'  E_ref={surrogate.structure.value} vs E_vqe={s_final.value} +/- {s_final.error}')
        plt.show()
    # end if
# end for
