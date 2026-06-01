#!/usr/bin/env python3

import numpy as np
from os import makedirs

from stalk import ParameterStructure
from stalk import PhaseAngle

from params import vqe_pes, kernel_vqe


interactive = __name__ == '__main__'

# H2 dimer
num_H = 2
# Spacing
d = 0.735
directory = f'H{num_H}/d_{d}/'
makedirs(directory, exist_ok=True)

# Let's start by defining the physical problem and placing the H-chain along z axis
pos = []
for i in range(num_H):
    pos.append([0.0, 0.0, i * d])
# end for
# ParameterStructure will be used here to contain the 'theta' parameters (to be optimized)
# and also the atomic structure (to be held fixed), and there is not requirement of consistency
# between the two
s_init = ParameterStructure(
    label='init',
    pos=pos,
    elem=num_H * ['H'],
    require_consistent=False,
)
# Obtain the number of VQE parameters from the kernel
ansatz = kernel_vqe(s_init)[0]
num_theta = len(ansatz.parameters)

# Start optimization from a semi-random initial guess
theta_init = [PhaseAngle(0.1 * v, label=f't{t}') for t, v in enumerate(np.random.randn(num_theta))]
s_init.params = theta_init

# Step 0: optimize the VQE classically (or load)
relaxfile = f'{directory}relax.dat'
try:
    params_relax = np.loadtxt(relaxfile)
    s_relax = s_init.copy(params_relax, label='relax')
except FileNotFoundError:
    s_relax = s_init.copy(label='relax')
    vqe_pes.relax(s_relax)
    np.savetxt(relaxfile, s_relax.params)
# end try
if interactive:
    print(s_relax)
# end if
