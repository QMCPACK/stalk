#!/usr/bin/env python3
"""Base class for representing a mapping between reducible positions (pos, axes) and irreducible parameters."""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array, isscalar, diag
from copy import deepcopy

from stalk.util import match_to_tol, get_fraction_error
from stalk.util.util import FF
from stalk.params.ParameterSet import ParameterSet


class ParameterStructure(ParameterSet):
    forward_func = None  # mapping function from pos to params
    backward_func = None  # mapping function from params to pos
    forward_args = None  # kwargs for the forward mapping
    backward_args = None  # kwargs for the backward mapping
    pos = None  # real-space position
    axes = None  # cell axes
    dim = None  # dimensionality
    elem = None  # list of elements
    units = None  # position units
    tol = None  # consistency tolerance
    require_consistent = None  # is consistency required?

    def __init__(
        self,
        forward=None,  # pos to params
        backward=None,  # params to pos
        pos=None,
        axes=None,
        elem=None,
        params=None,
        params_err=None,
        forward_args={},
        backward_args={},
        value=None,
        error=0.0,
        label='',
        units='B',
        dim=3,
        translate=True,  # attempt to translate pos
        tol=1e-7,
        require_consistent=True,
        **kwargs,  # kinds, labels, units
    ):
        self.dim = dim
        self.label = label
        self.tol = tol
        self.units = units
        self.require_consistent = require_consistent
        self.set_forward_func(forward, forward_args)
        self.set_backward_func(backward, backward_args)
        if params is not None:
            self.init_params(params, params_err, **kwargs)
            self.set_params(self.params)
        # end if
        if pos is not None:
            self.set_position(pos, translate=translate)
        # end if
        if axes is not None:
            self.set_axes(axes)
        # end if
        if value is not None:
            self.value = value
            self.error = error
        # end if
        if elem is not None:
            self.set_elem(elem)
        # end if
    # end def

    @property
    def consistent(self):
        return self.check_consistency()
    # end def

    @property
    def periodic(self):
        # By default axes is None, but Nexus may set it to empty list/array
        return self.axes is not None and len(self.axes) > 0
    # end def

    # Set forward mapping function and keyword arguments
    def set_forward_func(self, forward_func, forward_args={}):
        self.forward_func = forward_func
        self.forward_args = forward_args
        # Attempt to update params
        new_params = self.map_forward()
        if new_params is not None:
            self.set_params(new_params)
        # end if
    # end def

    # Set backward mapping function and keyword arguments
    def set_backward_func(self, backward_func, backward_args={}):
        self.backward_func = backward_func
        self.backward_args = backward_args
        # Attempt to update pos+axes
        new_pos, new_axes = self.map_backward()
        if new_pos is not None:
            self.set_position(new_pos, new_axes)
        # end if
    # end def

    # Set a new pos+axes and, if so configured, translate pos through backward mapping.
    def set_position(self, pos, axes=None, translate=True):
        pos = array(pos)
        assert pos.size % self.dim == 0, f'Position vector inconsistent with {self.dim} dimensions!'
        # Set the new pos+axes
        self.pos = array(pos).reshape(-1, self.dim)
        if axes is not None:
            # Set axes without mapping checks or moves
            self.set_axes(axes, check=False)
        # end if

        # If forward_func has been given, update params (but not positions)
        if self.forward_func is not None:
            # Only update params, not positions
            ParameterSet.set_params(self, self.map_forward())
        # end if

        # setting pos will unset value
        self.reset_value()

        # If set up to translate, take another move backward.
        if translate and self.backward_func is not None:
            self.pos, self.axes = self.map_backward()
        # end if
    # end def

    # Set a new pos+axes and, if so configured, translate pos through backward mapping.
    def set_axes(self, axes, check=True):
        if array(axes).size == self.dim:
            axes = diag(axes)
        else:
            axes = array(axes)
            assert axes.size == self.dim**2, f'Axes vector inconsistent with {self.dim} dimensions!'
            axes = array(axes).reshape(self.dim, self.dim)
        # end if
        self.axes = axes

        # Skip parameter updates and checks if they are expected later
        if check:
            # If forward_func has been given, update params; if not, return None
            params = self.map_forward(self.pos, self.axes)
            if params is not None:
                self.set_params(params)
            # end if
            self.reset_value()  # setting axes will reset value
        # end if
    # end def

    def set_elem(self, elem):
        self.elem = array(elem)
    # end def

    def set_params(self, params, params_err=None, dpos_mode=False):
        params_old = self.params
        ParameterSet.set_params(self, params, params_err)
        # After params have been set, attempt to update pos+axes
        if self.backward_func is not None:
            pos_new, self.axes = self.map_backward()
            if dpos_mode:
                self.pos += pos_new - self.map_backward(params_old)[0]
            else:
                self.pos = pos_new
            # end if
        # end if
    # end def

    # Perform forward mapping: if mapping function provided, return new params; else, return None
    def map_forward(self, pos=None, axes=None):
        pos = pos if pos is not None else self.pos
        if self.forward_func is None or pos is None:
            return None
        # end if
        if self.periodic:
            axes = axes if axes is not None else self.axes
            return array(self.forward_func(array(pos), axes, **self.forward_args))
        else:
            return array(self.forward_func(array(pos), **self.forward_args))
        # end if
    # end def

    # Perform backward mapping: if mapping function provided, return new pos+axes; else, return None
    def map_backward(self, params=None):
        params = params if params is not None else self.params
        if self.backward_func is None or params is None:
            return None, None
        # end if
        # Periodic mappings return (pos, axes) and non-periodic only pos
        res = self.backward_func(array(params), **self.backward_args)
        if type(res) is tuple:
            return array(res[0]).reshape(-1, 3), array(res[1]).reshape(-1, 3)
        else:
            return array(res).reshape(-1, 3), None
        # end if
    # end def

    def check_consistency(self, params=None, pos=None, axes=None, tol=None, verbose=False):
        """Check consistency of present forward-backward mapping.
        If params or pos/axes are supplied, check at the corresponding points. If not, check at the present point.
        """
        if not self.require_consistent:
            return True
        # end if
        if self.forward_func is None or self.backward_func is None:
            return False
        # end if
        tol = tol if tol is not None else self.tol
        axes = axes if axes is not None else self.axes
        if pos is None and params is not None:
            return self._check_params_consistency(array(params), tol)
        elif pos is not None and params is None:
            return self._check_pos_consistency(array(pos), array(axes), tol)
        # end if
        # if both params and pos are given, check their internal consistency
        pos = array(pos) if pos is not None else self.pos
        params = array(params) if params is not None else self.params
        if pos is not None and params is not None:
            params_new = array(self.map_forward(pos, axes))
            pos_new, axes_new = self.map_backward(params)
            if self.periodic:
                return match_to_tol(params, params_new, tol) and match_to_tol(pos, pos_new, tol) and match_to_tol(axes, axes_new, tol)
            else:
                return match_to_tol(params, params_new, tol) and match_to_tol(pos, pos_new, tol)
            # end if
        else:
            return False
        # end if
    # end def

    def _check_pos_consistency(self, pos, axes, tol=None):
        tol = tol if tol is not None else self.tol
        if self.periodic:
            params = self.map_forward(pos, axes)
            pos_new, axes_new = self.map_backward(params)
            consistent = match_to_tol(
                pos, pos_new, tol) and match_to_tol(axes, axes_new, tol)
        else:
            params = self.map_forward(pos, axes)
            pos_new, axes_new = self.map_backward(params)
            consistent = match_to_tol(pos, pos_new, tol)
        # end if
        return consistent
    # end def

    def _check_params_consistency(self, params, tol=None):
        tol = tol if tol is not None else self.tol
        pos, axes = self.map_backward(params)
        params_new = self.map_forward(pos, axes)
        return match_to_tol(params, params_new, tol)
    # end def

    def shift_pos(self, dpos, translate=True):
        assert self.pos is not None, 'position has not been set'
        if isscalar(dpos):
            new_pos = self.pos + dpos
        else:
            dpos = array(dpos)
            assert self.pos.size == dpos.size
            new_pos = (self.pos.flatten() + dpos.flatten()).reshape(-1, self.dim)
        # end if
        self.set_position(new_pos, translate=translate)
    # end def

    def shift_params(self, dparams, dpos_mode=False):
        params_old = self.params
        ParameterSet.shift_params(self, dparams)
        # After params have been shifted, attempt to update pos+axes
        if self.backward_func is not None:
            pos_new, self.axes = self.map_backward()
            if dpos_mode:
                self.pos += pos_new - self.map_backward(params_old)[0]
            else:
                self.pos = pos_new
            # end if
        # end if
    # end def

    def copy(
        self,
        params=None,
        params_err=None,
        label=None,
        pos=None,
        axes=None,
        offset=None,
        **kwargs,
    ):
        structure = deepcopy(self)
        if offset is not None:
            structure.offset = offset
        # end if
        if params is not None:
            structure.set_params(params, params_err)
        # end if
        if pos is not None:
            structure.set_position(pos)
        # end if
        if axes is not None:
            structure.set_axes(axes)
        # end if
        if label is not None:
            structure.label = label
        # end if
        return structure
    # end def

    def pos_difference(self, pos_ref):
        dpos = pos_ref.reshape(-1, 3) - self.pos
        return dpos
    # end def

    def jacobian(self, dp=0.001):
        assert self.consistent, 'The mapping must be consistent'
        jacobian = []
        for p in range(len(self.params)):
            params_this = self.params.copy()
            params_this[p] += dp
            pos, axes = self.map_backward(params_this)
            dpos = self.pos_difference(pos)
            jacobian.append(dpos.flatten() / dp)
        # end for
        return array(jacobian).T
    # end def

    def remap_forward(self, forward, N=None, fraction=0.159, **kwargs):
        assert self.consistent, 'The mapping must be consistent'
        pos, axes = self.map_backward()
        if self.periodic:
            params = forward(pos, axes, **kwargs)
        else:
            params = forward(pos, **kwargs)
        # end if
        if N is None:
            return params
        elif sum(self.params_err) > 0:  # resample errorbars
            if self.periodic:
                psdata = [forward(*self.map_backward(p))
                          for p in self.get_params_distribution(N=N)]
            else:
                psdata = [forward(self.map_backward(p)[0])
                          for p in self.get_params_distribution(N=N)]
            # end if
            params_err = array([get_fraction_error(ps, fraction=fraction)[
                               1] for ps in array(psdata).T])
            return params, params_err
        else:  # errors are zero
            return params, 0 * params
        # end if
    # end def

    def __str__(self):
        string = ParameterSet.__str__(self)
        if self.require_consistent:
            string += '\n  consistent: n/a'
        elif self.consistent:
            string += '\n  consistent: yes'
        else:
            string += '\n  consistent: no'
        # end if
        # pos
        if self.pos is None:
            string += '\n  pos: not set'
        else:
            string += '\n  pos ({:d} atoms)'.format(len(self.pos))
            for elem, pos in zip(self.elem, self.pos):
                string += ('\n    {:2s} ' + FF + FF + FF).format(
                    elem, pos[0], pos[1], pos[2])
            # end for
        # end if
        if self.periodic:
            string += '\n  periodic: yes'
            if self.axes is None:
                string += '\n  axes: not set'
            else:
                string += '\n  axes:'
                for axes in self.axes:
                    string += '\n    ' + (FF + FF + FF).format(
                        axes[0], axes[1], axes[2])
                # end for
            # end if
        else:
            string += '\n  periodic: no'
        # end if
        return string
    # end def

# end class
