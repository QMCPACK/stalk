#!/usr/bin/env python3

# This example is meant to accompany the file I/O tutorial available in the STALK
# documentation. It demonstrates some basic principles of file I/O in STALK.

from stalk import ParameterSet
from stalk import PesFunction
from stalk import morse


# A simple PES function based on the Morse potential
def pes_func(params: ParameterSet, r0=1.0, a=3.0, D=0.5, E_inf=0.0):
    E = morse([r0, a, D, E_inf], params[0])
    return E
# end def


p = ParameterSet([1.0])
pes = PesFunction(pes_func, r0=1.0, a=3.0, D=0.5, E_inf=0.0)
sigma = 0.04

# Single evaluation, no I/O (default behavior)
p0 = p.copy(label='no_io')
result0 = pes(p0)
print(f'no I/O energy: {result0.value}')
input('Press Enter to continue...\n')

# Single evaluation, I/O enabled
p1 = p.copy(label='with_io')
result1 = pes(p1, path='io/pes')
print(f'with I/O energy: {result1.value}')
print('See the "io/pes" directory for the generated files.')

# Load the value from disk
p_load = p.copy(label='with_io')
p_load.try_load_value('io/pes/with_io')
print(f'loaded energy: {p_load.value}')
