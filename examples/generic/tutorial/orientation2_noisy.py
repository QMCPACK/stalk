#!/usr/bin/env python3

# This example demonstrates the noisy evaluation of a PES in STALK orientation.

import numpy as np
from matplotlib import pyplot as plt

from orientation1_pes import pes_func

from stalk import ParameterSet
from stalk import PesFunction


sigma = 0.1  # noise level for the PES evaluations
p = ParameterSet([1.0], sigma=sigma)  # initial parameter set
pes = PesFunction(pes_func, c=[1.0, 2.0, 3.0])
n = 1000

energies = []
for _ in range(n):
    pes(p, add_sigma=True)
    energies.append(p.value)
# end for
# Calculate the apparent standard deviation of the noisy PES evaluations
sigma_out = np.std(energies)
# Calculate the aggregate mean of the noisy PES evaluations
E_mean = np.mean(energies)
# Calculate the exact PES for reference
E_exact = pes(p, add_sigma=False).value

plt.hist(energies, bins=20)
plt.axvline(E_exact, color='r', linestyle='dashed', linewidth=1, label='Exact PES value')
plt.axvline(E_mean, color='g', linestyle='dashed', linewidth=1, label=f'Mean of {n} PES evaluations')
plt.title(f'PES evaluations with sigma={sigma}, measured sigma_N={sigma_out:.4f}')
plt.legend()
plt.show()
