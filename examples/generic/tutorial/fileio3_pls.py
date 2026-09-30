#!/usr/bin/env python3

# This example is meant to accompany the file I/O tutorial available in the STALK
# documentation. It demonstrates some basic principles of file I/O in STALK.

from stalk import ParameterSet
from stalk import PesFunction
from stalk import ParameterHessian
from stalk import ParallelLineSearch
from stalk import morse


# A simple 2D PES function based on two uncoupled Morse potentials
def pes_func(params: ParameterSet, r0=1.0, a=3.0, D=0.5, E_inf=0.0):
    E0 = morse([r0, a, D, E_inf], params[0])
    E1 = morse([r0, a, D, E_inf], params[1])
    return E0 + E1
# end def


# Define a 2-parameter problem and a conforming Hessian.
p = ParameterSet([0.9, 1.1])
pes = PesFunction(pes_func, r0=1.0, a=3.0, D=0.5, E_inf=0.0)

hessian = ParameterHessian(structure=p)
hessian.hessian = [[1.0, 0.0], [0.0, 1.0]]

# Creating a ParallelLineSearch composes two line-searches in orthogonal directions, whose
# evaluation is performed in parallel.
pls = ParallelLineSearch(
    # This will create the following line-search directiories upon evaluation:
    #  io/pls/ls0
    #  io/pls/ls1
    path='io/pls',
    hessian=hessian,
    Rs=[0.2, 0.15],
    sigmas=[0.01, 0.02],
    M=5,
)
print('Parallel-linesearch before evaluation (unless loaded from io/pls):')
print(pls)
input('Press Enter to continue...\n')
# 'Propagation' means evaluation of all line-searches, and then solving for the next
# minimum-energy structure estimate
pls.propagate(pes, add_sigma=True)
print('Parallel-linesearch after evaluation:')
print(pls)

print('\nInitial parameters: p =', pls.structure.params)
print('Final parameters: p =', pls.structure_next.params)
input('Press Enter to continue...\n')

# As before, noisy evaluation will yield different results each time, unless the results
# are serialized to disk. This will happen automatically when the 'path' argument is used.
print('Resampling noisy ParallelLineSearch, with and without I/O serialization...')
for i in range(3):
    pls_no_io = ParallelLineSearch(
        hessian=hessian,
        Rs=[0.2, 0.15],
        sigmas=[0.01, 0.02],
        M=5,
    )
    pls_no_io.propagate(pes, add_sigma=True)
    pls_io = ParallelLineSearch(
        path='io/pls',
        hessian=hessian,
        Rs=[0.2, 0.15],
        sigmas=[0.01, 0.02],
        M=5,
    )
    pls_io.propagate(pes, add_sigma=True)
    print(f'  {i + 1}/3 no I/O params: {pls_no_io.structure_next.params}')
    print(f'  {i + 1}/3 w/ I/O params: {pls_io.structure_next.params}')
