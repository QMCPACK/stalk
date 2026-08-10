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
    print(f'Energy of the PES with fixed coeff. p = {pf.params}: E = {Ef.value}')

    p = ParameterSet([1.0])
    E1 = pes1(p)
    E2 = pes2(p)
    print(f'Energy of the PES1 at p = {p.params}: E1 = {E1.value}')
    print(f'Energy of the PES2 at p = {p.params}: E2 = {E2.value}')

    # Relaxing
    print('Relaxing PES1 using BFGS algorithm...')
    pes1.relax(p, method='BFGS', tol=1e-6)
    print('Relaxed parameters for pes1: ', p.params)
    print('Relaxed energy for pes1: ', p.value)
    p2 = p.copy()
    print('Relaxing PES2 using BFGS algorithm...')
    pes2.relax(p2, method='BFGS', tol=1e-6)
    print('Relaxed parameters for pes2: ', p2.params)
    print('Relaxed energy for pes2: ', p2.value)
# end if
