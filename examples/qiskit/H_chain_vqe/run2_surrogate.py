#!/usr/bin/env python3

from matplotlib import pyplot as plt

from stalk import Surrogate

from params import vqe_pes
from run1_hessian import hessian, directory


interactive = __name__ == '__main__'

# Step 2: generate a surrogate model of the parallel line-search by
# supplying the relaxed structure and its Hessian, and characterizing the
# PES along the conjugate directions
surrogate = Surrogate(
    path=f'{directory}surrogate/',
    fit_kind='pf3',
    hessian=hessian,
    window_frac=1.0,  # maximum displacement relative to Lambda of each direction
    M=11  # number of points per direction to sample
)
surrogate.evaluate(vqe_pes, interactive=interactive)

# Optimize the model such that the parameters can be resolved to accuracy
#   epsilon_p < |bias| + uncertainty
# using minimal statistical cost (maximum tolerated noise)
surrogate.optimize(
    temperature=0.001,
    fit_kind='pf2',  # fitting 3rd order polynomial
    M=5,  # 5 points along each line-search
    N=200,  # correlated error resampling population
    bias_order=1,  # treat "bias-induced" bias
    reoptimize=False,
    overwrite=True,
)

if interactive:
    # Plot the surrogate data of each line-search
    surrogate.plot()
    # Plot the error surface maps, showing the contour of requested tolerance
    # and the chosen point that maximizes input noise (sigma)
    surrogate.plot_error_surfaces()
    plt.show()
# end if
