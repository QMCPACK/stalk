#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

# Manual test to inspect printing of object properties

from matplotlib import pyplot as plt
from numpy import linspace

from stalk import LineSearch
from stalk import LineSearchBase
from stalk import LineSearchGrid
from stalk import ParameterHessian
from stalk import ParameterSet
from stalk import PesFunction


def pes_2d(structure: ParameterSet, x0=[0.0, 0.0], **kwargs):
    return sum((structure.params - x0)**2), 0.0
# end def


def plot_linesearchgrid():
    # Plot without data
    offsets = linspace(0.0, 1.0, 11)
    ls = LineSearchGrid(offsets=offsets)
    print("Plotting line-search grid offsets (does nothing but a warning):")
    ls.plot()
    input("Proceed? ")

    ls.values = offsets**2
    print("Plotting line-search grid offsets and data (no error):")
    ls.plot()

    ls.errors = offsets * 0.2
    print("Plotting line-search grid offsets, data and errors:")
    ls.plot()
# end def


def plot_linesearchbase():
    # Print without data
    offsets = linspace(0.0, 1.0, 11)
    ls = LineSearchBase(
        offsets=offsets,
    )
    print("Plotting line-search base offsets (does nothing but a warning):")
    ls.plot()
    input("Proceed? ")

    ls = LineSearchBase(
        offsets=offsets,
        values=(offsets - 0.5)**2
    )
    print("Plotting line-search base offsets and data (no error):")
    ls.plot()

    ls = LineSearchBase(
        offsets=offsets,
        values=(offsets - 0.5)**2,
        errors=offsets * 0.1,
        fit_kind='pf2'
    )
    print("Plotting line-search base offsets, data and errors (pf2):")
    ls.plot()
# end def


def plot_linesearch():
    # Print without data
    x0 = [2.0, 3.0]
    structure = ParameterSet(x0)
    hessian = ParameterHessian(
        structure=structure,
        hessian=[[0.5, 0.2], [0.2, 0.5]],
    )
    ls = LineSearch(
        hessian=hessian,
        d=0,
        W=0.4,
        sigma=0.1
    )
    print("Plotting line-search offsets (no eval):")
    ls.plot()
    input("Proceed? ")

    pes = PesFunction(pes_2d, {'x0': x0})
    ls.evaluate(pes=pes, add_sigma=False)
    print("Plotting line-search after evaluated without sigma:")
    ls.plot()

    ls.reset()
    ls.evaluate(pes=pes, add_sigma=True)
    print("Plotting line-search after evaluated with sigma:")
    ls.plot()
# end def


if __name__ == '__main__':
    plot_linesearchgrid()
    plot_linesearchbase()
    plot_linesearch()
    plt.show()
# end if
