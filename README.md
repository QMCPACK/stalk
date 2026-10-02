# Surrogate Theory Accelerated Line-search Kit (STALK)

Surrogate Theory Accelerated Line-search Kit (STALK) is a Python implementation of The
Surrogate Hessian Accelerated Parallel Line-search method. The method is intended for
local optimization of noisy, multivariable cost functions that can be approximated with a
surrogate model. It has the following features:
- **Robust**: rarely misses the minimum
- **Fast to converge**: requires minimal iterations to solution
- **Controllable accuracy**: can be tuned to the desired accuracy
- **Cost-efficient**: meets the accuracy with minimal statistical sampling
- **Parallelizable**: numerous evaluations may run in parallel

Typical applications of STALK include:
- Relaxation of atomic structures with noisy energies, e.g., quantum Monte Carlo. 
- Optimization of variational quantum eigensolvers

Advanced features include:
- Statistical analyses of black-box cost functions and line-searches
- Treatment of low-dimensional parametric spaces, Hessians, etc
- Gradient-free transition pathway search

## GETTING STARTED

The code can be installed with pip from PyPi
```bash
pip install stalk-qmc
```
See the [installation
instructions](https://stalk.readthedocs.io/en/latest/installation.html) for more details.

Installation comprises various methods in the `stalk` package, for example:
```python
# Using line-search along 45 degree direction, find the minimum of a simple 2D function
import numpy as np
from stalk import ParameterSet, PesFunction, LineSearch

def my_pes(x: ParameterSet) -> float:
    return (x.params[0])**2 + (x.params[1])**2

pes = PesFunction(my_pes)
p = ParameterSet([2.0, 2.0])
direction = np.array([0.5, 0.5])**0.5
sigma = 0.1  # statistical noise
ls = LineSearch(p, direction=direction, R=4.0, M=7, sigma=sigma)
pes(ls, add_sigma=True)
print(ls)
# The true solution (0, 0) is solved relative to the line-search starting point (2, 2)
# ls.x0 ~ -1 ± 0.05
# ls.y0 ~ -2**1.5 ± 0.1
p0 = p.copy()
p0.shift_params(ls.x0 * direction)
print(p0.params)
# p0.params ~ [0.0, 0.0]
```

See [Documentation](https://stalk.readthedocs.io/en/latest/) (WIP) to learn the basic
concepts, algorithms and operating principles of the code.

See [examples](https://github.com/QMCPACK/stalk/tree/master/examples) to study comprehensive
workflows of the code, including geometry relaxation, transition pathway search and
optimization of the variational quantum eigensolver.

## CITING

Upon publishing results based on the method, we kindly ask you to cite
[The original work](https://doi.org/10.1063/5.0079046).

> Juha Tiihonen, Paul R. C. Kent, and Jaron T. Krogel \
The Journal of Chemical Physics \
156, 054104 (2022)

## SUPPORT

The software and its documentation are under development with no warranties. Support may be
inquired by [contacting the authors](mailto:tiihonen@iki.fi).

## ACKNOWLEDGEMENTS

The authors of this method are Juha Tiihonen, Paul R. C. Kent and Jaron T.
Krogel, working in the Center for Predictive Simulation of Functional Materials
(https://cpsfm.ornl.gov/)

The code is written and developed by:
- Juha Tiihonen

Contributions to the algorithms have been made by:
- Jaron T. Krogel
- Gopal Iyer
- Simon Nirenberg

[Further contributions](CONTRIBUTING.md) are always welcome.

This work has been authored in part by UT-Battelle, LLC, under contract
DE-AC05-00OR22725 with the US Department of Energy (DOE).
