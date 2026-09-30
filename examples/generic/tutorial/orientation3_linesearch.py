#!/usr/bin/env python3

# This example is meant to demonstrate the line-search optimization method in the STALK
# orientation.

import numpy as np
from matplotlib import pyplot as plt

from orientation1_pes import pes_func

from stalk import ParameterSet
from stalk import LineSearchBase
from stalk import PesFunction
from stalk import LineSearch


# Define the PES and parameter set
pes = PesFunction(pes_func, c=[1.0, 2.0, 3.0])
p = ParameterSet([1.0])  # initial parameter set


# For the sake of comprehension, let's first consider the base LineSearch class. This class
# is agnostic of the PES context, and so it operates on abstract points, whose values are
# supplied manually. It is therefore more clumsy to use but also simple to understand.
if __name__ == "__main__":
    f, ax = plt.subplots()
    # Creating a regular grid of points assumed to contain the minimum
    points = np.linspace(-3, 3, 7)
    # Create a line-search grid based on the points, using 3rd order polynomial fit (pf3)
    lsb = LineSearchBase(points, fit_kind='pf3')
    print('Base Line-search grid before evaluation:')
    print(lsb)
    input('Press Enter to continue...\n')

    # Let us collect the PES values
    values = []
    print('Evaluating the PES at the points of the base line-search grid:')
    for point in points:
        # Copy the parameter set and set it to the new point
        p_this = p.copy(params=[point])
        # Evaluate and record the energy value
        E = pes(p_this)
        values.append(E.value)
        print(f'  E({point}) = {E.value}')
    # end for

    # Supply the values and perform the search
    lsb.values = values
    lsb.search()
    print('Base Line-search grid after evaluation and search:')
    print(lsb)
    input('Press Enter to continue...\n')
    lsb.plot(ax=ax, color='tab:blue')
# end if


# Let's do it again but now with the context-aware LineSearch class. Here, we operate
# directly on the parameter sets, and the search will relative to an initial point
if __name__ == "__main__":
    # Let's do it again but now with the context-aware LineSearch class
    offsets = np.linspace(-3, 3, 7)
    # Now, the search is centered around 'p' searches by the offsets around it.
    # d: The direction index is mandatory. Here, d=0 means that the search direction will
    # be along the first parameter of the 1D parameter set.
    ls = LineSearch(p, d=0, offsets=offsets)
    print('Contextual Line-search before evaluation:')
    print(ls)
    input('Press Enter to continue...\n')

    # When operating with parameter sets (instead of abstract points), we can readily
    # evaluate the PES.
    pes(ls)
    print('Contextual Line-search after evaluation and search:')
    print(ls)
    ls.plot(ax=ax, color='tab:orange')
    plt.legend()
    plt.show()
# end if
