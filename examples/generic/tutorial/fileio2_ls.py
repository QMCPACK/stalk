#!/usr/bin/env python3

# This example is meant to accompany the file I/O tutorial available in the STALK
# documentation. It demonstrates some basic principles of file I/O in STALK.

from stalk import ParameterSet
from stalk import PesFunction
from stalk import LineSearch
from stalk import morse


# A simple PES function based on the Morse potential
def pes_func(params: ParameterSet, r0=1.0, a=3.0, D=0.5, E_inf=0.0):
    E = morse([r0, a, D, E_inf], params[0])
    return E
# end def


p = ParameterSet([1.0])
pes = PesFunction(pes_func, r0=1.0, a=3.0, D=0.5, E_inf=0.0)
sigma = 0.04

# Resample noiseless line-search, no I/O
print('Resampling noiseless line-search, no I/O...')
for i in range(3):
    ls = LineSearch(structure=p, d=0, M=5, R=4.0, sigma=sigma)
    pes(ls)
    # The result will be same each time, because the PES is deterministic
    print(f'  {i + 1}/3 found minimum: {ls.x0:+.3f} ± {ls.x0_err:.3f}')
# end for
input('Press Enter to continue...\n')

# Resample noiseless line-search, I/O enabled
print('Resampling noiseless line-search, with I/O...')
for i in range(3):
    ls = LineSearch(structure=p, d=0, M=5, R=4.0, sigma=sigma)
    # I/O is enabled by providing a 'path' for the evaluation.
    pes(ls, path='io/ls')
    # At first, the evaluated points will be created. Later, the evaluated points will
    # be loaded from disk. Regardless, the result will remain the same each time.
    print(f'  {i + 1}/3 found minimum: {ls.x0:+.3f} ± {ls.x0_err:.3f}')
# end for
input('Press Enter to continue...\n')

# Now, let us resample a noisy line-search with I/O enabled
print('Resampling noisy line-search, with I/O...')
for i in range(3):
    ls = LineSearch(structure=p, d=0, M=5, R=4.0, sigma=sigma)
    # We will actually use the old results, only this time adding random noise.
    # The noise from 'add_sigma=True' is added _after_ loading, so the results will be
    # different each time. This is an important reminder that the solution of a noisy
    # line-search is also a random variable.
    pes(ls, path='io/ls', add_sigma=True)
    print(f'  {i + 1}/3 found minimum: {ls.x0:+.3f} ± {ls.x0_err:.3f}')
# end for
input('Press Enter to continue...\n')

# Finally, let us serialize the line-search result. Then, we can decide to make evaluation
# of the noisy line-search conditional: only evaluate if the result is not available yet.

print('Resampling noisy line-search, with I/O and result serialization...')
for i in range(3):
    # In addition to before, let us attempt to load the line-search from 'io/ls'
    # At first, there is no data and the loading will fail silently.
    ls = LineSearch(structure=p, d=0, M=5, R=4.0, sigma=sigma, path='io/ls')
    # If fitting result is not available, evaluate the line-search and save the result.
    if ls.fit_res is None:
        pes(ls, add_sigma=True)
        ls.save_result('io/ls')
    # end if
    print(f'  {i + 1}/3 found minimum: {ls.x0:+.3f} ± {ls.x0_err:.3f}')
# end for
