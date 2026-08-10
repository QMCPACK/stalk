#!/usr/bin/env python3

# This example is meant to demonstrate the noisy line-search fitting in STALK orientation.

import numpy as np
from matplotlib import pyplot as plt

from orientation1_pes import pes_func

from stalk import ParameterSet
from stalk import PesFunction
from stalk import LineSearch


# Define the PES and parameter set
c = [1.0, 2.0, 3.0]
pes = PesFunction(pes_func, c=c, create_files=False)
p = ParameterSet([1.0])  # initial parameter set

# Let's define 3 different line-searches with variable grids
offsets1 = np.linspace(-2, 2, 7)
offsets2 = np.linspace(-3, 3, 7)
offsets3 = np.linspace(-4, 4, 7)
sigma = 0.1  # noise level for the PES evaluation
ls1 = LineSearch(p, d=0, offsets=offsets1, sigma=sigma)
ls2 = LineSearch(p, d=0, offsets=offsets2, sigma=sigma)
ls3 = LineSearch(p, d=0, offsets=offsets3, sigma=sigma)
# exact minimum of the PES
x0_exact = -c[1] / (2 * c[0]) - p[0]

res1, res2, res3 = [], [], []
if __name__ == "__main__":
    f, ax = plt.subplots()
    for _ in range(100):
        # Evaluate the line-searches (bypass intermediate results)
        ls1.evaluate(pes, reset_value=True, add_sigma=True)
        ls2.evaluate(pes, reset_value=True, add_sigma=True)
        ls3.evaluate(pes, reset_value=True, add_sigma=True)
        res1.append([ls1.fit_res.x0, ls1.fit_res.y0])
        res2.append([ls2.fit_res.x0, ls2.fit_res.y0])
        res3.append([ls3.fit_res.x0, ls3.fit_res.y0])
    # end for
    res1 = np.array(res1)
    res2 = np.array(res2)
    res3 = np.array(res3)
    print('LS results with max-offsets=2:')
    print(f'  x0={ls1.fit_res.x0} +/- {ls1.fit_res.x0_err}, ')
    print(f'  y0={ls1.fit_res.y0} +/- {ls1.fit_res.y0_err}, ')
    print('LS results with max-offsets=3:')
    print(f'  x0={ls2.fit_res.x0} +/- {ls2.fit_res.x0_err}, ')
    print(f'  y0={ls2.fit_res.y0} +/- {ls2.fit_res.y0_err}, ')
    print('LS results with max-offsets=4:')
    print(f'  x0={ls3.fit_res.x0} +/- {ls3.fit_res.x0_err}, ')
    print(f'  y0={ls3.fit_res.y0} +/- {ls3.fit_res.y0_err}, ')
    ax.scatter(res1[:, 0], res1[:, 1], label='ls1', alpha=0.5)
    ax.scatter(res2[:, 0], res2[:, 1], label='ls2', alpha=0.5)
    ax.scatter(res3[:, 0], res3[:, 1], label='ls3', alpha=0.5)
    ax.set_xlabel('x0')
    ax.set_ylabel('y0')
    plt.axvline(x0_exact, color='r', linestyle='dashed', linewidth=1, label='Exact PES value')
    ax.set_title(f'Line-search results with sigma={sigma}')
    ax.legend()
    plt.show()
# end if
