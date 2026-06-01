#!/usr/bin/env python3

from qiskit import transpile
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit.primitives import BackendEstimatorV2
from qiskit_algorithms.minimum_eigensolvers import NumPyMinimumEigensolver
from qiskit.primitives import StatevectorEstimator
from qiskit.providers.fake_provider import GenericBackendV2
from iqm.qiskit_iqm import IQMFakeAphrodite  # 54-qubit Aphrodite architecture

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
    mapper = ParityMapper(problem.num_particles)
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


# Idealized VQE PES that only considers phenomenological precision
def pes_ideal(
    structure: ParameterStructure,
    sigma=0.0,
    path='',
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


def get_backend(backend, q_hamiltonian):
    if backend is None or backend == 'generic':
        backend = GenericBackendV2(q_hamiltonian.num_qubits)
        initial_layout = None
    elif backend == 'fake_aphrodite':
        # fake IQM Aphrodite PES
        backend = IQMFakeAphrodite()
        initial_layout = [
            backend.qubit_name_to_index("QB5"),
            backend.qubit_name_to_index("QB6")
        ]
    else:
        raise ValueError(f"Unsupported backend: {backend}")
    # end if
    return backend, initial_layout
# end def


# More realostic VQE PES that derives precision from a backend and transpiled circuit
#  NB: cannot handle zero noise (sigma=0)
def pes_backend(
    structure: ParameterStructure,
    sigma=0.0,
    path='',
    backend=None,
    kernel_args={},  # dict to hold ansatz, operator, estimator
    **kwargs  # charge=0, spin=0, basis="sto3g"
):
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

    backend, initial_layout = get_backend(backend, q_hamiltonian)
    transpiled_ansatz = transpile(
        ansatz,
        backend=backend,
        initial_layout=initial_layout,
        optimization_level=3
    )
    # Apply the same layout to the qubit Hamiltonian
    transpiled_hamiltonian = q_hamiltonian.apply_layout(transpiled_ansatz.layout)

    estimator = BackendEstimatorV2(backend=backend)

    # Evaluate VQE with given parameters
    params = structure.params.reshape(-1, len(structure.params))
    job = estimator.run([(transpiled_ansatz, transpiled_hamiltonian, params)], precision=sigma)
    energy = job.result()[0].data.evs[0]
    kernel_args['evals'] += 1

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
    ansatz, mapper, q_hamiltonian = kernel_vqe(structure, **kwargs)
    solver = NumPyMinimumEigensolver()
    result = solver.compute_minimum_eigenvalue(q_hamiltonian)
    energy = result.eigenvalue.real

    printout = 'Run NumpyMinimumSolver:'
    printout += f" E = {'%+5.4f' % energy} +/- {'%+5.4f' % sigma}"
    print(printout)
    return energy, sigma
# end def


# Ideal VQE PES
vqe_pes = PesFunction(func=pes_ideal, kernel_args={})
# Generic face backend PES
backend_pes = PesFunction(
    func=pes_backend,
    kernel_args={},
    backend='generic',
)
aphrodite_pes = PesFunction(
    func=pes_backend,
    kernel_args={},
    backend='fake_aphrodite',
)
