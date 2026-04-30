#!/usr/bin/env python3

# Manual test to inspect printing of object properties

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


def print_linesearchgrid():
    # Print without data
    offsets = linspace(0.0, 1.0, 11)
    ls = LineSearchGrid(offsets=offsets)
    print("Printing line-search grid offsets (no eval):")
    print(ls)
    input("Proceed? ")

    ls.values = offsets**2
    print("Printing line-search grid offsets and data (no error):")
    print(ls)
    input("Proceed? ")

    ls.errors = offsets * 0.2
    print("Printing line-search grid offsets, data and errors:")
    print(ls)
    input("Proceed? ")
# end def


def print_linesearchbase():
    # Print without data
    offsets = linspace(0.0, 1.0, 11)
    ls = LineSearchBase(
        offsets=offsets,
    )
    print("Printing line-search base offsets (no eval):")
    print(ls)
    input("Proceed? ")

    ls = LineSearchBase(
        offsets=offsets,
        values=(offsets - 0.5)**2
    )
    print("Printing line-search base offsets and data (no error):")
    print(ls)
    input("Proceed? ")

    ls = LineSearchBase(
        offsets=offsets,
        values=(offsets - 0.5)**2,
        errors=offsets * 0.1,
        fit_kind='pf2'
    )
    print("Printing line-search base offsets, data and errors (pf2):")
    print(ls)
    input("Proceed? ")
# end def


def print_linesearch():
    # Print without data
    x0 = [2.0, 3.0]
    structure = ParameterSet(x0)
    hessian = ParameterHessian(
        [[0.5, 0.2], [0.2, 0.5]],
        structure=structure,
    )
    ls = LineSearch(
        hessian=hessian,
        d=0,
        W=0.4,
        sigma=0.1
    )
    print("Printing line-search offsets (no eval):")
    print(ls)
    input("Proceed? ")

    pes = PesFunction(pes_2d, {'x0': x0})
    ls.evaluate(pes=pes, add_sigma=False)
    print("Printing line-search after evaluated without sigma:")
    print(ls)
    input("Proceed? ")

    ls.evaluate(pes=pes, add_sigma=True)
    print("Printing line-search after evaluated with sigma:")
    print(ls)
    input("Proceed? ")
# end def


if __name__ == '__main__':
    print_linesearchgrid()
    print_linesearchbase()
    print_linesearch()
# end if
