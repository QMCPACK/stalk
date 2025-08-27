#!/usr/bin/env python

from matplotlib import pyplot as plt

from stalk import TargetLineSearch

from params import pes_pyscf
from run1_hessian import hessian

interactive = __name__ == "__main__"

surrogate = TargetLineSearch(
    d=0,
    fit_kind='pf3',
    structure=hessian.structure,
    hessian=hessian,
    R=0.3,
    M=11,
    pes=pes_pyscf,
    path='surrogate',
    interactive=interactive,
)
surrogate.bracket_target_bias()
surrogate.reset_interpolation(interpolate_kind='cubic')

# Set target parameter error tolerances (epsilon): 0.01 Angstrom accuracy
# Then, optimize the surrogate line-search to meet the tolerances given the line-search
surrogate.optimize(
    epsilon=0.02,  # Optimize to 0.02 Angstrom
    fit_kind='pf3',
    M=7,
    N=400,
)

if interactive:
    print(surrogate)
    surrogate.plot()
    surrogate.plot_error_surface()
    plt.show()
# end if
