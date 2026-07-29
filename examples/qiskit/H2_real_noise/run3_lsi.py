#!/usr/bin/env python3

from matplotlib import pyplot as plt
from pathlib import Path
import numpy as np

from stalk import LineSearchIteration
from stalk import Surrogate
from stalk import PesFunction
from stalk import ParameterStructure

from params import pes_backend
from run2_surrogate import surrogate, directory


interactive = __name__ == '__main__'


# Step 3: Optimize VQE by line-search from a shifted position
# Shift the starting position for some challenge
def run_stalk_lsi(
    surrogate: Surrogate,
    structure,
    n: int,
    temperature=0.001,
    niter: int = 5,
    backend_kind='generic',
):
    surrogate.optimize(temperature=temperature, reoptimize=True)
    path = f'{directory}STALK/{backend_kind}/T{temperature}/lsi{n}'
    pes = PesFunction(func=pes_backend, kernel_args={}, backend=backend_kind)
    lsi = LineSearchIteration(
        surrogate=surrogate,
        structure=structure,
        path=path,
    )
    # Propagate the parallel line-search (compute values, analyze, then move on) 5 times
    for i in range(niter):
        lsi.propagate(pes, i)
    # end for
    # Evaluate the latest eqm structure
    lsi.pls().evaluate_eqm(pes)
    return lsi
# end for


def get_structure(structure_orig: ParameterStructure, directory, n, sigma=0.4):
    sfile = f'{directory}{n}_shifts_init.dat'
    if Path(sfile).exists():
        shifts = np.loadtxt(sfile)
        print(f'Loaded initial shifts from {sfile}')
    else:
        shifts = sigma * np.random.randn(len(surrogate.structure))
        np.savetxt(sfile, shifts)
    # end if
    structure = structure_orig.copy()
    structure.shift_params(shifts)
    return structure
# end def


# Run the algorithm a few times with default settings to test out
lsis, ntot = [], 5
backend_kind = 'generic'
temperatures = [0.0005, 0.001, 0.002, 0.004]
for n in range(ntot):
    lsi_row = []
    structure = get_structure(surrogate.structure, directory, n, sigma=0.4)
    for temperature in temperatures:
        print(f'Running STALK line-search {n} with temperature = {temperature}')
        lsi = run_stalk_lsi(
            surrogate=surrogate,
            structure=structure,
            n=n,
            niter=5,
            temperature=temperature,
            backend_kind=backend_kind,
        )
        lsi_row.append(lsi)
    # end for
    lsis.append(lsi_row)
# end for

if interactive:
    # Plot the data of each line-search
    for lsi_row in lsis:
        for lsi in lsi_row:
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
        # end for
    # end for
    plt.show()
# end if
