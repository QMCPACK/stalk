#!/usr/bin/env python3

# This example is meant to accompany the file I/O tutorial available in the STALK
# documentation. It demonstrates some basic principles of file I/O in STALK.

from stalk import ParameterSet
from stalk import PesFunction
from stalk import ParameterHessian
from stalk import morse
from stalk import LineSearchIteration


# A simple PES function based on the Morse potential
def pes_func(params: ParameterSet, r0=1.0, a=3.0, D=0.5, E_inf=0.0):
    E0 = morse([r0, a, D, E_inf], params[0])
    E1 = morse([r0, a, D, E_inf], params[1])
    return E0 + E1
# end def


p = ParameterSet([0.9, 1.1])
pes = PesFunction(pes_func, r0=1.0, a=3.0, D=0.5, E_inf=0.0)

hessian = ParameterHessian(structure=p)
hessian.hessian = [[1.0, 0.0], [0.0, 1.0]]

# Create an iteration of noisy parallel line-searches with and without I/O.
for k in range(3):
    lsi_no_io = LineSearchIteration(
        hessian=hessian,
        sigmas=[0.01, 0.01],
    )
    lsi_io = LineSearchIteration(
        hessian=hessian,
        sigmas=[0.01, 0.01],
        path='io/lsi',
    )
    for i in range(3):
        lsi_no_io.propagate(pes, i=i, add_sigma=True)
        lsi_io.propagate(pes, i=i, add_sigma=True)
    # end for
    print(f'\n{k + 1}/3 final parameters after Line-search iteration:')
    # The final params withou serialization will be different because of the random noise
    print(f'  no I/O: {lsi_no_io.structure_final.params} (different each time)')
    # The inal params with I/O will be the same, because the first iteration will be
    # loaded from the disk on subsequent runs.
    print(f'  w/ I/O: {lsi_io.structure_final.params} (same each time)')
# end for
