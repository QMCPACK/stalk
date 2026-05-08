#!/usr/bin/env python3

from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_algorithms.minimum_eigensolvers import NumPyMinimumEigensolver
from qiskit_algorithms.optimizers import COBYLA
from qiskit_algorithms import VQE
from qiskit.primitives import StatevectorEstimator

from stalk import ParameterStructure
from stalk import PesFunction


# Generic VQE kernel for the H chain using PySCFDriver
def kernel_vqe(
    structure: ParameterStructure,
    charge=0,
    spin=0,
    basis="sto3g",
    callback=None,
    optimizer=COBYLA(maxiter=200),
    initial_point=None,
    exact=False,
    **kwargs,
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
    # Set up the mappers and transformers
    mapper = JordanWignerMapper()
    # Hydrogen simplifies to 1 el per ion
    num_electrons = len(structure.elem)
    transformer = ActiveSpaceTransformer(num_electrons=num_electrons, num_spatial_orbitals=num_electrons)

    # Transform the problem
    problem = transformer.transform(problem)

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
    )

    est = StatevectorEstimator()
    if exact:
        # Set up exact solver
        solver = GroundStateEigensolver(mapper, NumPyMinimumEigensolver())
    else:
        # Set up VQE
        vqe = VQE(est, ansatz, optimizer, callback=callback, initial_point=initial_point)
        solver = GroundStateEigensolver(mapper, vqe)
    # end if
    operator = solver.get_qubit_operators(problem, None)[0]
    return ansatz, operator, est, solver, problem
# end def


def pes_vqe(
    structure: ParameterStructure,
    sigma=0.0,
    path='',
    kernel_args={},  # dict to hold ansatz, operator, estimator
    **kwargs  # charge=0, spin=0, basis="sto3g", callback=None, exact=False
):
    # Cache the VQE kernel to cut redundant operations
    if all([k in kernel_args for k in ['ansatz', 'operator', 'estimator']]):
        ansatz = kernel_args['ansatz']
        operator = kernel_args['operator']
        estimator = kernel_args['estimator']
    else:
        ansatz, operator, estimator, solver, problem = kernel_vqe(structure)
        kernel_args['ansatz'] = ansatz
        kernel_args['operator'] = operator
        kernel_args['estimator'] = estimator
    # end if
    # Count evaluations
    if 'evals' not in kernel_args:
        kernel_args['evals'] = 0
    # end

    # Evaluate VQE with given parameters
    params = structure.params.reshape(-1, len(structure.params))
    if sigma > 0:
        # Finite precision
        job = estimator.run([(ansatz, operator, params)], precision=sigma)
    else:
        # Exact precision
        job = estimator.run([(ansatz, operator, params)])
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
    sigma=0.0,
    path='',
    **kwargs  # charge=0, spin=0, basis="sto3g", callback=None, exact=False
):
    ansatz, operator, estimator, exact_solver, problem = kernel_vqe(structure, exact=True)
    exact_result = exact_solver.solve(problem)
    exact_energy = exact_result.groundenergy

    printout = 'Run NumpyMinimumSolver:'
    printout += f" E = {'%+5.4f' % exact_energy} +/- {'%+5.4f' % sigma}"
    print(printout)
    return exact_energy, sigma
# end def


# VQE surrogate PES
vqe_pes = PesFunction(func=pes_vqe, kernel_args={})
# Exact PES
exact_pes = PesFunction(func=pes_exact, kernel_args={})
# Create another instance to allow (optionally) a different PES and to reset eval count
#   NB: Using here the same PES, only this time it is noisy
vqe_pes_noisy = PesFunction(func=pes_vqe, kernel_args={})
