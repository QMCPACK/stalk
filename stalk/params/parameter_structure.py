#!/usr/bin/env python3
"""Base class for representing a mapping between reducible positions (pos, axes) and irreducible parameters."""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array
from copy import deepcopy

from stalk.params.parameter_mapping import ParameterMapping
from stalk.util import get_fraction_error
from stalk.util.function_caller import FunctionCaller
from stalk.util.util import FF
from stalk.params.parameter_set import ParameterSet


class ParameterStructure(ParameterSet):
    _mapping: ParameterMapping = None
    _pos = None  # real-space position
    _axes = None  # cell axes
    _elem = None  # list of elements
    units = None  # position units
    tol = None  # consistency tolerance
    require_consistent = None  # is consistency required?

    def __init__(
        self,
        pos=None,
        axes=None,
        elem=None,
        params=None,
        params_err=None,
        mapping=None,  # parametric mapping
        forward=None,  # pos to params
        forward_args={},
        backward=None,  # params to pos
        backward_args={},
        dim=3,
        value=None,
        error=0.0,
        label='',
        units='B',
        tol=1e-7,
        require_consistent=True,
    ):
        self.label = label
        self.tol = tol
        self.units = units
        self.require_consistent = require_consistent
        # Initialize parameter mapping
        self._init_mapping(
            mapping,
            forward,
            forward_args,
            backward,
            backward_args,
            dim=dim,
        )
        if params is not None:
            self.params = params
        # end if
        if params_err is not None:
            self.params_err = params_err
        # end if
        if pos is not None:
            # Add here in case axes is needed to update pos
            self._axes = axes
            self.pos = pos
        # end if
        if axes is not None:
            self.axes = axes
        # end if
        if value is not None:
            self.value = value
            self.error = error
        # end if
        if elem is not None:
            self.elem = elem
        # end if
    # end def

    @property
    def mapping(self):
        return self._mapping
    # end def

    @mapping.setter
    def mapping(self, mapping: ParameterMapping):
        if not isinstance(mapping, ParameterMapping):
            raise TypeError(f'The mapping must be ParameterMapping, provided: {mapping}')
        # end if
        self._mapping = mapping
        # Not updating the parameters
    # end def

    @property
    def forward(self):
        return self.mapping.forward
    # end def

    @forward.setter
    def forward(self, forward):
        if isinstance(forward, FunctionCaller) or callable(forward):
            self.mapping.set_forward(forward)
            # Map pos, axes forward if possible
            if self.pos is not None:
                self.params = self.map_forward(self.pos, self.axes)
            # end if
        else:
            raise TypeError(f'Reset forward mapping must be a FunctionCaller, provided: {forward}')
        # end if
        # Not updating the parameters
    # end def

    @property
    def backward(self):
        return self.mapping.backward
    # end def

    @backward.setter
    def backward(self, backward):
        if isinstance(backward, FunctionCaller) or callable(backward):
            self.mapping.set_backward(backward)
            # Map params backward if possible
            if self.params is not None:
                self.pos, self.axes = self.map_backward(self.params)
            # end if
        else:
            raise TypeError(f'Reset backward mapping must be a FunctionCaller, provided: {backward}')
        # end if
        # Not updating the pos, axes here
    # end def

    @property
    def consistent(self):
        if self.elem is None or self.params is None or self.mapping.incomplete:
            # Obviously False if not both pos and params present or mapping incomplete
            return False
        # end if
        consistent = self.mapping.check_params_consistency(self.params, tol=self.tol)
        consistent &= self.mapping.check_pos_consistency(self.pos, self.axes, tol=self.tol)
        return consistent
    # end def

    @property
    def elem(self):
        if self._elem is None:
            if self.pos is None:
                return []
            else:
                return len(self.pos) * [None]
            # end if
        else:
            return self._elem
        # end if
    # end def

    @elem.setter
    def elem(self, elem):
        if elem is None:
            self._elem = None
        elif self.pos is None or len(elem) != len(self.pos):
            raise ValueError('The "elem" list must be the same length as self.pos')
        else:
            # TODO: check actual contents
            self._elem = elem
        # end if
    # end def

    @property
    def dim(self):
        return self.mapping.dim
    # end def

    @property
    def pos(self):
        return self._pos
    # end def

    @pos.setter
    def pos(self, pos):
        if pos is None:
            self._pos = None
        # end if
        pos = array(pos).reshape(-1, self.dim)
        self._pos = pos
        # If the number of positions should change, reset 'elem'
        if len(pos) != len(self.elem):
            self.elem = None
        # end if
        # Map forward if possible
        if self.forward is not None:
            params = self.map_forward(self.pos, self.axes)
            if self.require_consistent:
                self.params = params
            else:
                # Reset params using parent class method
                ParameterSet.params.fset(self, params)
            # end if
        # end if
        self.reset_value()
    # end def

    # Override to add backward synchronization
    @ParameterSet.params.setter
    def params(self, params):
        # First set/init the params
        ParameterSet.params.fset(self, params)
        # Map backward if possible
        if self.backward is not None:
            pos, axes = self.map_backward(self.params)
            self._pos = pos
            self._axes = axes
        # end if
    # end def

    @property
    def axes(self):
        return self._axes
    # end def

    @axes.setter
    def axes(self, axes):
        if axes is None:
            self._axes = None
        else:
            axes = array(axes).reshape(self.dim, self.dim)
            self._axes = axes
            # Map forward if possible
            if self.forward is not None:
                params = self.map_forward(self.pos, axes)
                self.params = params
            # end if
        # end if
        self.reset_value()
    # end def

    @property
    def periodic(self):
        # By default axes is None, but Nexus may set it to empty list/array
        return self.axes is not None and len(self.axes) > 0
    # end def

    def _init_mapping(
        self,
        mapping=None,
        # Keep this constructor for backward compatibility
        forward_func=None,
        forward_args=None,
        backward_func=None,
        backward_args=None,
        dim=3,
    ):
        if not isinstance(mapping, ParameterMapping):
            mapping = ParameterMapping(
                forward_func=forward_func,
                forward_args=forward_args,
                backward_func=backward_func,
                backward_args=backward_args,
                dim=dim,
            )
        # end if
        self.mapping = mapping
    # end def

    def map_forward(self, pos, axes=None, **kwargs):
        return self.mapping.map_forward(pos, axes, **kwargs)
    # end def

    def map_backward(self, params, **kwargs):
        return self.mapping.map_backward(params, **kwargs)
    # end def

    def shift_params(self, shifts, dpos_mode=False):
        params_old = self.params
        ParameterSet.shift_params(self, shifts)
        # Map backward if possible
        if self.backward is not None:
            pos_new, self._axes = self.map_backward(self.params)
            if dpos_mode:
                pos_old = self.map_backward(params_old)[0]
                self._pos += pos_new - pos_old
            else:
                self._pos = pos_new
            # end if
        # end if
    # end def

    # Kept for backward compatibility
    def set_position(self, pos, axes=None, translate=True):
        self.pos = pos
        if axes is not None:
            self.axes = axes
        # end if

        # If set up to translate, take another move backward.
        if translate and self.consistent:
            self.pos, self.axes = self.map_backward(self.params)
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
            structure.params = params
        # end if
        if params_err is not None:
            structure.params_err = params_err
        # end if
        if pos is not None:
            structure.pos = pos
        # end if
        if axes is not None:
            structure.axes = axes
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
        if not self.consistent:
            raise AssertionError('The mapping must be consistent')
        # end if
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
        if not self.consistent:
            raise AssertionError('The mapping must be consistent')
        # end if
        pos, axes = self.map_backward(self.params)
        try:
            params = forward(pos, axes=axes, **kwargs)
        except TypeError:
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
        if not self.require_consistent:
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
