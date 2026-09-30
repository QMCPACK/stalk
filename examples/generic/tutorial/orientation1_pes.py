#!/usr/bin/env python3

# This example is meant to accompany the orientation tutorial available in the STALK
# documentation. It demonstrates some of the basic concepts and building blocks of the STALK
# framework, including the potential energy surface (PES) and parameter sets and basic
# line-searches.

from stalk import ParameterSet
from stalk import PesFunction


# A raw PES function with fixed coefficients
def pes_func_fixed_c(params: ParameterSet):
    c = [1.0, 2.0, 3.0]
    E = c[0] * params[0]**2 + c[1] * params[0] + c[2]
    return E
# end def


# A raw PES function with kwargs coefficients
def pes_func(params: ParameterSet, c: list = [1.0, 2.0, 3.0]):
    E = c[0] * params[0]**2 + c[1] * params[0] + c[2]
    return E
# end def


# Define properly wrapped-up PES functions
pes_fixed = PesFunction(pes_func_fixed_c)
pes1 = PesFunction(pes_func, c=[1.0, 2.0, 3.0])
pes2 = PesFunction(pes_func, c=[1.0, -2.0, 0.0])


if __name__ == "__main__":
    # Define sets of parameters to evaluate
    pf = ParameterSet([1.0])
    Ef = pes_fixed(pf)
    print('Energy of the PES with fixed coefficients.')
    print(f'  E({pf.params}) = {Ef.value}')
    input('Press Enter to continue...\n')

    # Create a new set of parameters and evaluate the first parameterized PES
    p1 = ParameterSet([1.0])
    E1 = pes1(p1)
    print(f'Energy of the PES1 with coefficients c = {pes1.args['c']}:')
    print(f'  E1({p1.params}) = {E1.value}')
    input('Press Enter to continue...\n')

    # Relax the first PES using a classical algorithm
    print('Relaxing PES1 using BFGS algorithm...')
    pes1.relax(p1, method='BFGS', tol=1e-6)
    print(f'Relaxed parameters: p1* = {p1.params}')
    print(f'Relaxed energy: E1(p1*) = {p1.value}')
    input('Press Enter to continue...\n')

    # Create a new set of parameters and evaluate the second PES
    p2 = ParameterSet([1.0])
    E2 = pes2(p2)
    print(f'Energy of the PES2 with coefficients c = {pes2.args['c']}:')
    print(f'  E2({p2.params}) = {E2.value}')
    input('Press Enter to continue...\n')
    # Relax the second PES using a classical algorithm
    print('Relaxing PES2 using BFGS algorithm...')
    pes2.relax(p2, method='BFGS', tol=1e-6)
    print(f'Relaxed parameters: p2* = {p2.params}')
    print(f'Relaxed energy: E2(p2*) = {p2.value}')
# end if
