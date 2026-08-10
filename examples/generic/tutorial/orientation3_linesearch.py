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
pes = PesFunction(pes_func, c=[1.0, 2.0, 3.0], create_files=False)
p = ParameterSet([1.0])  # initial parameter set

# Creating a regular grid of points assumed to contain the minimum
points = np.linspace(-3, 3, 7)
# Yes, this is a bit hacky, because we start up with a contextless LineSearchBase, which
# does not know about the parameter set or the PES, just a grid of points and values.
values = [pes(p.copy(params=[point]), reset_value=True).value for point in points]

if __name__ == "__main__":
    f, ax = plt.subplots()
    # Create a line-search grid
    lsb = LineSearchBase(points, fit_kind='pf3')
    print('Base Line-search grid before evaluation:')
    print(lsb)
    # Supply values
    lsb.values = values
    results = lsb.search()
    lsb.fit_res = results
    print('Base Line-search grid after evaluation and search:')
    print(lsb)
    lsb.plot(ax=ax, color='tab:blue')
# end if


if __name__ == "__main__":
    # Let's do it again but now with the context-aware LineSearch class
    offsets = np.linspace(-3, 3, 7)
    # Now, the search is centered around 'p' searches by the offsets around it. Clearer, eh?
    ls = LineSearch(p, d=0, offsets=offsets)
    print('Contextual Line-search before evaluation:')
    # When operating with parameter sets (instead of abstract points), we can readily evaluate
    # the PES.
    ls.evaluate(pes)
    print('Contextual Line-search after evaluation and search:')
    print(ls)
    ls.plot(ax=ax, color='tab:orange')
    plt.legend()
    plt.show()
# end if
