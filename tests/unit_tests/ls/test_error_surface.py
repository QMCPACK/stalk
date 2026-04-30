#!/usr/bin/env python3

from numpy import abs

from pytest import raises
from stalk.ls.error_surface import ErrorSurface
from stalk.util.util import match_to_tol

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# test ErrorSurface class
def test_ErrorSurface():

    # Utility functions to simulate error monotonous with w, sigma parameters
    def errors(x, y):
        return abs(y) + x**2
    # end def

    def add_row(y, es: ErrorSurface):
        errs = [errors(x, y) for x in es.Xs]
        es.insert_row(y, errs)
        return errs
    # end def

    def add_col(x, es: ErrorSurface):
        errs = [errors(x, y) for y in es.Ys]
        es.insert_col(x, errs)
        return errs
    # end def

    es = ErrorSurface()
    assert len(es.Xs) == 1
    assert len(es.Ys) == 1
    assert es.X_res == 0.1
    assert es.Y_res == 0.1
    assert es.X_mat.shape[0] == 1
    assert es.X_mat.shape[1] == 1
    assert es.Y_mat.shape[0] == 1
    assert es.Y_mat.shape[1] == 1
    assert es.E_mat.shape[0] == 1
    assert es.E_mat.shape[1] == 1
    assert es.T_mat.shape[0] == 1
    assert es.T_mat.shape[1] == 1
    assert es.X_mat[0, 0] == 0.0
    assert es.Y_mat[0, 0] == 0.0
    assert es.E_mat[0, 0] == 0.0
    assert es.T_mat[0, 0]

    # Add x-point to end
    x0 = 1.0
    errs0 = add_col(x0, es)
    assert match_to_tol(es.E_mat[:, 1], errs0)
    assert match_to_tol(es.Xs, [0.0, x0])
    assert match_to_tol(es.Ys, [0.0])
    # Add x-point to middle
    x1 = 0.5
    errs1 = add_col(0.5, es)
    assert match_to_tol(errs1, es.E_mat[:, 1])
    assert match_to_tol(es.Xs, [0.0, x1, x0])
    assert match_to_tol(es.Ys, [0.0])

    # Add y-point to end
    y2 = 0.4
    errs2 = add_row(y2, es)
    assert match_to_tol(errs2, es.E_mat[-1, :])
    assert match_to_tol(es.Xs, [0.0, x1, x0])
    assert match_to_tol(es.Ys, [0.0, y2])
    # Add y-point to middle
    y3 = 0.2
    errs3 = add_row(y3, es)
    assert match_to_tol(errs3, es.E_mat[1, :])
    assert match_to_tol(es.Xs, [0.0, x1, x0])
    assert match_to_tol(es.Ys, [0.0, y3, y2])

    # Check internal consistency
    for Y_row in es.Y_mat.T:
        assert match_to_tol(Y_row, es.Ys)
    # end for
    for X_col in es.X_mat:
        assert match_to_tol(X_col, es.Xs)
    # end for

    # Cannot add negative or inconsistent data
    with raises(ValueError):
        y = -0.1
        errs = [errors(x, y) for x in es.Xs]
        es.insert_col(y, errs)
    # end with
    with raises(ValueError):
        y = 0.1
        errs = [errors(x, y) for x in es.Xs][:-1]
        es.insert_col(y, errs)
    # end with
    with raises(ValueError):
        x = -0.1
        errs = [errors(x, y) for x in es.Ys]
        es.insert_row(x, errs)
    # end with
    with raises(ValueError):
        x = 0.1
        errs = [errors(x, y) for x in es.Ys][:-1]
        es.insert_row(x, errs)
    # end with

    # Test y-maximization
    # Cannot optimize to sub-zero epsilon
    with raises(ValueError):
        es._argmax_y(-1e-9)
    # end with

    # Small-epsilon limit is in the origin
    xi, yi = es._argmax_y(1e-9)
    assert xi == 0
    assert yi == 0

    # High-epsilon limit is in the opposite corner
    xi, yi = es._argmax_y(1e9)
    assert xi == len(es.Xs) - 1
    assert yi == len(es.Ys) - 1

    # Medium-epsilon is inside E_mat
    xi, yi = es._argmax_y(es.E_mat[1, 1] + 1e-2)
    assert xi == 1
    assert yi == 1

    # TODO: es.request_points
    # es.evaluate_surface
    # es.evaluate_target

# end def
