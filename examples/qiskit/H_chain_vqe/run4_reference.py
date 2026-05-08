#!/usr/bin/env python3

from os import makedirs
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from qiskit_algorithms.optimizers import CG, COBYLA, POWELL

from params import kernel_vqe
from run0_relax import s_init, directory
from run3_vqe_ls import lsis, structures


params_history = []
energy_history = []


def callback(eval_count, params, value, estimator_result):
    print(f"Energy = {value:.6f} Ha ({eval_count})")
    energy_history.append(value)
    params_history.append(params)
# end def


# Exact solver energy for reference
ansatz, operator, estimator, exact_solver, problem = kernel_vqe(s_init, exact=True)
exact_result = exact_solver.solve(problem)
exact_energy = exact_result.groundenergy

# Resample alternative optimizers so many times
optimizers = []
# using the same starting point as STALK?
same_init = False
# COBYLA
optimizers.append(('COBYLA', COBYLA(maxiter=400), [], [], []))
# CG
optimizers.append(('CG', CG(maxiter=20), [], [], []))
# POWELL
optimizers.append(('POWELL', POWELL(maxiter=20), [], [], []))

# Loop through all optimizers
for label, optimizer, edata, pdata, podata in optimizers:
    for n, structure in enumerate(structures):
        makedirs(f'{directory}/{label}/', exist_ok=True)
        if same_init:
            efile = f'{directory}{label}/{n}_energy.dat'
            pfile = f'{directory}{label}/{n}_params.dat'
            pofile = f'{directory}{label}/{n}_params_opt.dat'
            initial_point = structure.params
        else:
            efile = f'{directory}{label}/{n}_energy_random.dat'
            pfile = f'{directory}{label}/{n}_params_random.dat'
            pofile = f'{directory}{label}/{n}_params_opt_random.dat'
            initial_point = None
        # end if
        if Path(efile).exists() and Path(pfile).exists() and Path(pofile).exists():
            energies = np.loadtxt(efile)
            params_raw = np.loadtxt(pfile)
            params_opt = np.loadtxt(pofile)
            print(f'Loaded {label}, round {n} from disk')
        else:
            params_history.clear()
            energy_history.clear()
            ansatz, operator, estimator, vqe_solver, problem = kernel_vqe(
                s_init,
                callback=callback,
                optimizer=optimizer,
                # To be fair, use the same initial point
                initial_point=initial_point,
            )
            print(f'Starting {label}, round {n}:')
            vqe_result = vqe_solver.solve(problem)
            print(f'Finished {label}, round {n}')
            energies = np.array(energy_history)
            params_raw = np.array(params_history)
            params_opt = vqe_result.raw_result.optimal_point
            np.savetxt(efile, energies)
            np.savetxt(pfile, params_raw)
            np.savetxt(pofile, params_opt)
        # end if
        edata.append(energies)
        pdata.append(params_raw)
        podata.append(params_opt)
    # end for
# end for


# Plot energies
plt.figure()
plt.xlabel('VQE iteration')
plt.ylabel('Energy (Ha)')
if same_init:
    plt.title('Energy convergence (same init)')
else:
    plt.title('Energy convergence (random init)')
# end if
maxlen = 0
for label, optimizer, edata, pdata, podata in optimizers:
    p = plt.plot(np.nan, label=label)
    for row in edata:
        maxlen = max(maxlen, len(row))
        p = plt.plot(row, alpha=1.0 / len(edata), color=p[0].get_color())
    # end for
# end for

p = plt.plot(np.nan, marker='o', label='STALK')
for lsi in lsis:
    ls_energies, ls_errs, ls_evals, tot_evals = [], [], [], 0
    for pls in lsi.pls_list:
        ls_evals.append(tot_evals)
        tot_evals += sum([len(ls) for ls in pls.ls_list])
        ls_energies.append(pls.structure.value)
        ls_errs.append(pls.structure.error)
    # end for
    plt.errorbar(
        ls_evals,
        ls_energies,
        ls_errs,
        marker='o',
        color=p[0].get_color(),
    )
    maxlen = max(maxlen, tot_evals)
# end for

# Exact reference line
plt.plot([0, maxlen], 2 * [exact_energy], 'k-', label='Exact solver')

plt.legend()
plt.show()
