#!/usr/bin/env python3

from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_algorithms.minimum_eigensolvers import NumPyMinimumEigensolver
from qiskit.primitives import StatevectorEstimator

from stalk import ParameterStructure
from stalk import PesFunction


# Generic VQE kernel for the H chain using PySCFDriver
def kernel_vqe(
    structure: ParameterStructure,
    charge=0,
    spin=0,
    basis="sto3g",
    **kwargs
):
    atom = []
    for elem, pos in zip(structure.elem, structure.pos):
        atom.append(f'{elem} {pos[0]} {pos[1]} {pos[2]}')
    # end for
    # Define molecule (linear chain, equal H-H spacing)
    driver = PySCFDriver(
        atom=atom,
        basis=basis,
        charge=charge,
        spin=spin,
    )
    # Get the problem
    problem = driver.run()
    mapper = JordanWignerMapper()
    fermionic_hamiltonian = problem.hamiltonian.second_q_op()
    q_hamiltonian = mapper.map(fermionic_hamiltonian).simplify()

    # Prepare the ansatz
    hf_state = HartreeFock(
        num_spatial_orbitals=problem.num_spatial_orbitals,
        num_particles=problem.num_particles,
        qubit_mapper=mapper,
    )
    ansatz = UCCSD(
        num_spatial_orbitals=problem.num_spatial_orbitals,
        num_particles=problem.num_particles,
        qubit_mapper=mapper,
        initial_state=hf_state,
        reps=1,
    )
    return ansatz, mapper, q_hamiltonian
# end def


def pes_vqe(
    structure: ParameterStructure,
    kernel_args={},  # dict to hold ansatz, operator, estimator
    **kwargs  # charge=0, spin=0, basis="sto3g"
):
    # Cache the VQE kernel to cut redundant operations
    if all([k in kernel_args for k in ['ansatz', 'hamiltonian', 'mapper']]):
        ansatz = kernel_args['ansatz']
        mapper = kernel_args['mapper']
        q_hamiltonian = kernel_args['hamiltonian']
    else:
        ansatz, mapper, q_hamiltonian = kernel_vqe(structure, **kwargs)
        kernel_args['ansatz'] = ansatz
        kernel_args['hamiltonian'] = q_hamiltonian
        kernel_args['mapper'] = mapper
    # end if
    # Count evaluations
    if 'evals' not in kernel_args:
        kernel_args['evals'] = 0
    # end if

    estimator = StatevectorEstimator()
    sigma = structure.sigma

    # Evaluate VQE with given parameters
    params = structure.params.reshape(-1, len(structure.params))
    if sigma > 0:
        # Finite precision
        job = estimator.run([(ansatz, q_hamiltonian, params)], precision=sigma)
    else:
        # Exact precision
        job = estimator.run([(ansatz, q_hamiltonian, params)])
    # end if
    energy = job.result()[0].data.evs[0]
    kernel_args['evals'] += 1

    # In this ideal VQE, sigma=target_precision is also the apparent uncertainty
    printout = 'Run theta ='
    for p in params[0]:
        printout += '  %+5.4f' % p
    # end for
    printout += f", E = {'%+5.4f' % energy} +/- {'%+5.4f' % sigma}"
    printout += f'  ({kernel_args['evals']})'
    print(printout)
    return energy, sigma
# end def


# This is a PES that solves the exact ground state numerically, provided an atomic structure
def pes_exact(
    structure: ParameterStructure,
    **kwargs  # charge=0, spin=0, basis="sto3g", callback=None, exact=False
):
    ansatz, mapper, q_hamiltonian = kernel_vqe(structure, **kwargs)
    solver = NumPyMinimumEigensolver()
    result = solver.compute_minimum_eigenvalue(q_hamiltonian)
    energy = result.eigenvalue.real

    printout = 'Run NumpyMinimumSolver:'
    printout += f" E = {'%+5.4f' % energy}"
    print(printout)
    return energy, 0.0
# end def


# VQE surrogate PES
vqe_pes = PesFunction(func=pes_vqe, kernel_args={})
# Exact PES
exact_pes = PesFunction(func=pes_exact, kernel_args={})
# Create another instance to allow (optionally) a different PES and to reset eval count
#   NB: Using here the same PES, only this time it is noisy
vqe_pes_noisy = PesFunction(func=pes_vqe, kernel_args={})
