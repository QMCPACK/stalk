#!/usr/bin/env python3
"""Class for representing a mapping between reducible positions (pos, axes) and irreducible parameters."""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array, ndarray

from stalk.util.function_caller import FunctionCaller
from stalk.util.util import match_to_tol


class ParameterMapping():
    _forward: FunctionCaller = None
    _backward: FunctionCaller = None
    dim = None

    def __init__(
        self,
        forward_func=None,  # pos to params
        forward_args={},
        backward_func=None,  # params to pos
        backward_args={},
        dim=3,
    ):
        self.dim = dim
        if forward_func is not None:
            self.set_forward(forward_func, forward_args)
        # end if
        if backward_func is not None:
            self.set_backward(backward_func, backward_args)
        # end if
    # end def

    @property
    def forward(self):
        return self._forward
    # end def

    @property
    def backward(self):
        return self._backward
    # end def

    @property
    def incomplete(self):
        return self.forward is None or self.backward is None
    # end def

    def set_forward(
        self,
        forward_func,
        forward_args={},
    ):
        if isinstance(forward_func, FunctionCaller):
            self._forward = forward_func
        else:
            self._forward = FunctionCaller(forward_func, **forward_args)
        # end if
    # end def

    def set_backward(
        self,
        backward_func,
        backward_args={},
    ):
        if isinstance(backward_func, FunctionCaller):
            self._backward = backward_func
        else:
            self._backward = FunctionCaller(backward_func, **backward_args)
        # end if
    # end def

    # Perform forward mapping if present and return params; else, return None
    def map_forward(
        self,
        pos: ndarray,
        axes=None,
        **kwargs,
    ):
        if self.forward is None:
            return None
        # end if
        args = self.forward.get_updated(kwargs)
        try:
            params = self.forward.func(pos=pos, axes=axes, **args)
        except TypeError:
            # If axes not supported, try without
            params = self.forward.func(pos=pos, **args)
        # end try
        return array(params, dtype=float).flatten()
    # end def

    # Perform backward mapping if present and return pos, axes; else, return None, None
    def map_backward(
        self,
        params: ndarray,
        **kwargs,  # e.g. axes
    ):
        if self.backward is None:
            return None, None
        # end if
        args = self.backward.get_updated(kwargs)
        result = self.backward.func(params=params, **args)

        if isinstance(result, tuple):
            if len(result) == 2:
                pos = array(result[0], dtype=float).reshape(-1, self.dim)
                axes = array(result[1], dtype=float).reshape(-1, self.dim)
            else:
                raise ValueError(f'The backward func must either return tuple(pos, axes) or pos array, returning {result}')
            # end if
        else:
            pos = array(result, dtype=float).reshape(-1, self.dim)
            axes = None
        # end if
        return pos, axes
    # end def

    def check_pos_consistency(
        self,
        pos: ndarray,
        axes=None,
        tol=1e-6,
    ):
        if self.incomplete:
            # Obviously False if the mapping is incomplete
            return False
        # end if
        # Map pos+args forward
        params = self.map_forward(pos=pos, axes=axes)
        # Then back
        pos_new, axes_new = self.map_backward(params)
        consistent = len(pos_new) == len(pos)
        for p, (po, pn) in enumerate(zip(pos, pos_new)):
            if any(abs(po - pn) > tol):
                consistent = False
                print(f'  Consistency warning, pos[{p}]: {po} - {pn} > {tol}')
            # end if
        # end for
        # Check axes only if both are present
        if axes_new is not None and axes is not None:
            consistent &= match_to_tol(axes_new, axes, tol=tol)
        # end if
        return consistent
    # end def

    def check_params_consistency(
        self,
        params,
        tol=1e-6
    ):
        if self.incomplete:
            # Obviously False if the mapping is incomplete
            return False
        # end if
        # Map params forward
        pos, axes = self.map_backward(params)
        if axes is None:
            params_new = self.map_forward(pos)
        else:
            params_new = self.map_forward(pos, axes=axes)
        # end if
        consistent = len(params) == len(params_new)
        for p, (po, pn) in enumerate(zip(params, params_new)):
            if abs(po - pn) > tol:
                consistent = False
                print(f'  Consistency warning, p[{p}]: {po} - {pn} > {tol}')
            # end if
        # end for
        return consistent
    # end def

# end class
