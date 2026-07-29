#!/usr/bin/env python3

from os import makedirs
from qiskit import transpile
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import CG, COBYLA, POWELL
from qiskit.primitives import BackendEstimatorV2

from params import get_backend, kernel_vqe, pes_exact
from run0_relax import s_init, directory
from run3_lsi import get_structure, run_stalk_lsi, surrogate


def run_classical_optimizer(
    optimizer,
    structure,
    label,
    n,
    precision=0.01,
    same_init=True,
    backend_kind='generic',
):
    # Record the history of parameters and energies for each optimizer
    params_history = []
    energy_history = []

    def callback(eval_count, params, value, estimator_result):
        print(f"Energy = {value:.6f} Ha ({eval_count})")
        energy_history.append(value)
        params_history.append(params)
    # end def

    path = f'{directory}/{label}/{backend_kind}/prec_{precision}/'
    makedirs(path, exist_ok=True)
    if same_init:
        efile = f'{path}/energy_{n}.dat'
        pfile = f'{path}/params_{n}.dat'
        pofile = f'{path}/params_opt_{n}.dat'
        initial_point = structure.params
    else:
        efile = f'{path}/energy_{n}_random.dat'
        pfile = f'{path}/params_{n}_random.dat'
        pofile = f'{path}/params_opt_{n}_random.dat'
        initial_point = None
    # end if
    if Path(efile).exists() and Path(pfile).exists() and Path(pofile).exists():
        energies = np.loadtxt(efile)
        params_raw = np.loadtxt(pfile)
        params_opt = np.loadtxt(pofile)
        print(f'Loaded {label}, round {n} from disk')
    else:
        ansatz, mapper, q_hamiltonian = kernel_vqe(
            s_init,
            callback=callback,
            optimizer=optimizer,
            # To be fair, use the same initial point
            initial_point=initial_point,
        )
        print(f'Starting {label}, round {n}:')
        backend, initial_layout = get_backend(backend_kind, q_hamiltonian)
        transpiled_ansatz = transpile(
            ansatz,
            backend=backend,
            initial_layout=initial_layout,
            optimization_level=3
        )
        estimator = BackendEstimatorV2(
            backend=backend,
            options={"default_precision": precision}
        )
        vqe = VQE(
            estimator,
            transpiled_ansatz,
            optimizer,
            callback=callback,
            initial_point=initial_point
        )
        vqe_result = vqe.compute_minimum_eigenvalue(q_hamiltonian)
        print(f'Finished {label}, round {n}')
        energies = np.array(energy_history)
        params_raw = np.array(params_history)
        params_opt = vqe_result.optimal_point
        np.savetxt(efile, energies)
        np.savetxt(pfile, params_raw)
        np.savetxt(pofile, params_opt)
    # end if
    return energies, params_raw, params_opt
# end def


# Get the exact energy for reference
exact_energy, sigma = pes_exact(s_init)

# Use the following initial structures for all methods
structures = []
for n in range(5):
    structure = get_structure(surrogate.structure, directory, n, sigma=0.4)
    structures.append(structure)
# end for

# First run the STALK line-searches to get the data for comparison
temperatures = [0.0005, 0.001, 0.002, 0.004]
lsis_fake, lsis_back = [], []
for temperature in temperatures:
    lsi_fake_row, lsi_back_row = [], []
    for n, structure in enumerate(structures):
        print(f'Running STALK line-search {n} on fake_aphrodite at T = {temperature}')
        lsi_fake = run_stalk_lsi(
            surrogate=surrogate,
            structure=structure,
            n=n,
            niter=5,
            temperature=temperature,
            backend_kind='fake_aphrodite',
        )
        lsi_fake_row.append(lsi_fake)
        print(f'Running STALK line-search {n} on generic backend at T = {temperature}')
        lsi_back = run_stalk_lsi(
            surrogate=surrogate,
            structure=structure,
            n=n,
            niter=5,
            temperature=temperature,
            backend_kind='generic',
        )
        lsi_back_row.append(lsi_back)
    # end for
    lsis_fake.append(lsi_fake_row)
    lsis_back.append(lsi_back_row)
# end for


# Compare selected optimizers with STALK, using either the same initial point or
# random initial points
same_init = True
precisions = [0.005, 0.01, 0.02]
refdata_generic = []
# COBYLA
refdata_generic.append(('COBYLA', COBYLA(maxiter=400), [], [], []))
# CG
refdata_generic.append(('CG', CG(maxiter=20), [], [], []))
# POWELL
refdata_generic.append(('POWELL', POWELL(maxiter=20), [], [], []))

# Calculate data for the classical optimizers with generic backend
for label, optimizer, edata, pdata, podata in refdata_generic:
    for precision in precisions:
        edatarow, pdatarow, podatarow = [], [], []
        for n, structure in enumerate(structures):
            energies, params_raw, params_opt = run_classical_optimizer(
                optimizer=optimizer,
                structure=structure,
                label=label,
                n=n,
                precision=precision,
                same_init=same_init,
                backend_kind='generic',
            )
            edatarow.append(energies)
            pdatarow.append(params_raw)
            podatarow.append(params_opt)
        # end for
        edata.append(edatarow)
        pdata.append(pdatarow)
        podata.append(podatarow)
    # end for
# end for


refdata_fake_aphrodite = []
# COBYLA
refdata_fake_aphrodite.append(('COBYLA', COBYLA(maxiter=400), [], [], []))
# CG
refdata_fake_aphrodite.append(('CG', CG(maxiter=20), [], [], []))
# POWELL
refdata_fake_aphrodite.append(('POWELL', POWELL(maxiter=20), [], [], []))

# Calculate data for the classical optimizers with fake Aphrodite backend
for label, optimizer, edata, pdata, podata in refdata_fake_aphrodite:
    for precision in precisions:
        edatarow, pdatarow, podatarow = [], [], []
        for n, structure in enumerate(structures):
            energies, params_raw, params_opt = run_classical_optimizer(
                optimizer=optimizer,
                structure=structure,
                label=label,
                n=n,
                precision=precision,
                same_init=same_init,
                backend_kind='fake_aphrodite',
            )
            edatarow.append(energies)
            pdatarow.append(params_raw)
            podatarow.append(params_opt)
        # end for
        edata.append(edatarow)
        pdata.append(pdatarow)
        podata.append(podatarow)
    # end for
# end for


# Plot energies vs evaluations
plt.figure()
i_prec = 0
plt.title(f'Energy convergence vs evals, prec={precisions[i_prec]}')
plt.xlabel('Energy evaluations')
plt.ylabel('Energy (Ha)')
maxlen = 0
for label, optimizer, edata, pdata, podata in refdata_generic:
    p = plt.plot(np.nan, label=label)
    for row in edata[i_prec]:
        maxlen = max(maxlen, len(row))
        p = plt.plot(row, alpha=1.0 / len(edata), color=p[0].get_color())
    # end for
# end for

i_temp = 0
p = plt.plot(np.nan, marker='o', label=f'STALK (tol={temperatures[i_temp]})')
for lsi in lsis_back[i_temp]:
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


# Plot energies vs projected shots
def to_shots(sigma):
    return sigma**-2
# end def


plt.figure()
plt.title(f'Energy convergence vs shots, prec={precisions[i_prec]}')
plt.xlabel('Relative shots')
plt.ylabel('Energy (Ha)')
maxlen = 0
for label, optimizer, edata, pdata, podata in refdata_generic:
    p = plt.plot(np.nan, label=label)
    for row in edata[i_prec]:
        shots = to_shots(precisions[i_prec]) * np.arange(1, len(row) + 1)
        p = plt.plot(shots, row, alpha=1.0 / len(edata), color=p[0].get_color())
        maxlen = max(maxlen, shots[-1])
    # end for
# end for

p = plt.plot(np.nan, marker='o', label=f'STALK (tol={temperatures[i_temp]})')
for lsi in lsis_back[i_temp]:
    tot_shots = 0
    ls_energies, ls_errs, ls_evals, tot_evals = [], [], [], 0
    for pls in lsi.pls_list:
        ls_evals.append(tot_shots)
        tot_shots = tot_shots + sum([to_shots(ls.sigma) * len(ls) for ls in pls.ls_list])
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
    maxlen = max(maxlen, tot_shots)
# end for
# Exact reference line
plt.plot([0, maxlen], 2 * [exact_energy], 'k-', label='Exact solver')
plt.legend()

plt.show()
