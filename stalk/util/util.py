#!/usr/bin/env python3
"""Various utility functions and constants commonly needed in line-search workflows"""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import exp, median, array, isnan
from numpy import meshgrid, linalg, linspace, dot, eye


Bohr = 0.5291772105638411  # A
Ry = 13.605693012183622  # eV
Hartree = 27.211386024367243  # eV

# Print formats

# Floating-point format (signed)
FF = '{:> 9.4f} '
# Floating-point format (for errors)
FU = '+/- {:<5.4f} '
# String format for float fields (right-aligned)
FFS = '{:>9s} '
# String format for float fields (left-aligned)
FFSL = '{:<9s} '
# Integer format
FI = '{:>4d} '
# String format for integer fields
FIS = '{:>4s} '
# Percentage format (unsigned)
FP = '{:>5.3f}% '
# String format for percentage fields
FPS = '{:>6s} '
# Path+label format
PL = '{}/{}'
# Label format (right-aligned)
FL = '{:>10s} '
# Label format (left-aligned)
FLL = '{:<10s} '
# structure label format
SL = 'd{}_{:+5.4f}'


def get_fraction_error(data, fraction, both=False):
    """Estimate uncertainty from a distribution based on a percentile fraction"""
    if fraction < 0.0 or fraction >= 0.5:
        raise ValueError('Invalid fraction')
    # end if
    data = array(data, dtype=float)
    data = data[~isnan(data)]  # remove nan
    ave = median(data)
    data = data[data.argsort()] - ave  # sort and center
    pleft = abs(data[int((len(data) - 1) * fraction)])
    pright = abs(data[int((len(data) - 1) * (1 - fraction))])
    if both:
        err = [pleft, pright]
    else:
        err = max(pleft, pright)
    # end if
    return ave, err
# end def


def match_to_tol(val1, val2, tol=1e-8):
    """Match the values of two vectors. True if all match, False if not."""
    val1 = array(val1).flatten()
    val2 = array(val2).flatten()
    return abs(val1 - val2).max() < tol
# end def


def bipolynomials(X, Y, nx, ny):
    """Construct a bipolynomial expansion of variables

    XYp = x**0 y**0, x**0 y**1, x**0 y**2, ...
    courtesy of Jaron Krogel"""
    X = X.flatten()
    Y = Y.flatten()
    Xp = [0 * X + 1.0]
    Yp = [0 * Y + 1.0]
    for n in range(1, nx + 1):
        Xp.append(X**n)
    # end for
    for n in range(1, ny + 1):
        Yp.append(Y**n)
    # end for
    XYp = []
    for Xn in Xp:
        for Yn in Yp:
            XYp.append(Xn * Yn)
        # end for
    # end for
    return XYp
# end def bipolynomials


def bipolyfit(X, Y, Z, nx, ny):
    """Fit to a bipolynomial set of variables"""
    XYp = bipolynomials(X, Y, nx, ny)
    p, r, rank, s = linalg.lstsq(array(XYp).T, Z.flatten(), rcond=None)
    return p
# end def bipolyfit


def bipolyval(p, X, Y, nx, ny):
    """Evaluate based on a bipolynomial set of variables"""
    shape = X.shape
    XYp = bipolynomials(X, Y, nx, ny)
    Z = 0 * X.flatten()
    for pn, XYn in zip(p, XYp):
        Z += pn * XYn
    # end for
    Z.shape = shape
    return Z
# end def bipolyval


def bipolymin(p, X, Y, nx, ny, itermax=6, shrink=0.1, npoints=10):
    """Find the minimum of a bipolynomial set of variables"""
    for i in range(itermax):
        Z = bipolyval(p, X, Y, nx, ny)
        X = X.ravel()
        Y = Y.ravel()
        Z = Z.ravel()
        imin = Z.argmin()
        xmin = X[imin]
        ymin = Y[imin]
        zmin = Z[imin]
        dx = shrink * (X.max() - X.min())
        dy = shrink * (Y.max() - Y.min())
        xi = linspace(xmin - dx / 2, xmin + dx / 2, npoints)
        yi = linspace(ymin - dy / 2, ymin + dy / 2, npoints)
        X, Y = meshgrid(xi, yi)
        X = X.T
        Y = Y.T
    # end for
    return xmin, ymin, zmin
# end def bipolymin


def directorize(path: str):
    """If missing, add '/' to the end of path"""
    path = str(path)
    if len(path) > 0 and not path[-1] == '/':
        path += '/'
    # end if
    path = path.replace("//", "/")
    return path
# end def


# Find an orthonormal basis for the subspace orthogonal to 'vec'.
def orthogonal_subspace_basis(vec, threshold=1e-10):
    # Normalize 'vec'
    vec = vec / linalg.norm(vec)

    # Get the dimension of the space
    dim = len(vec)

    # Create an initial basis for the whole space including 'vec'
    initial_basis = [vec] + [eye(dim)[:, i] for i in range(dim)]

    # Orthonormalize the entire basis
    orthonormalized = gram_schmidt(initial_basis)

    # Exclude the orthonormal vector corresponding to 'vec'
    subspace_basis = array([v for v in orthonormalized if abs(dot(v, vec)) < threshold])

    return subspace_basis
# end def


# Perform the Gram-Schmidt process on a list of vectors.
def gram_schmidt(vectors):
    orthonormal_basis = []
    for v in vectors:
        for u in orthonormal_basis:
            v -= dot(u, v) * u
        # end for
        if linalg.norm(v) > 0:
            orthonormal_basis.append(v / linalg.norm(v))
        # end if
    return orthonormal_basis
# end def


def morse(p, r):
    # p0: eqm value, p1: stiffness, p2: well depth, p3: E_inf
    return p[2] * ((1 - exp(-(r - p[0]) / p[1]))**2 - 1) + p[3]
# end def
