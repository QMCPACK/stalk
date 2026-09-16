#!/usr/bin/env python3
'''Class for a collection of structures.'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from typing import Generic, TypeVar

from stalk.params.linesearch_point import LineSearchPoint

T = TypeVar('T', bound=LineSearchPoint)


class StructureCollection(Generic[T]):
    # List of structures inherited from LineSearchPoint
    _grid: list[T] = []

    def __init__(self, thr=1e-6) -> None:
        self._grid = list[T]()
        self.thr = thr
    # end def

    @property
    def grid(self) -> list[T]:
        '''Return list of points'''
        return [point for point in self._grid]
    # end def

    @grid.setter
    def grid(self, grid: list[T]) -> None:
        self._grid = list[T]()
        for point in grid:
            self.add_point(point)
        # end for
    # end def

    def collect_enabled(self) -> list[T]:
        '''Return a list of enabled structures'''
        return [point for point in self._grid if point.enabled]
    # end def

    def finalize(self) -> None:
        '''Finalization hook post evaluation to be implemented in subclasses.'''
        pass
    # end def

    def add_point(self, point: T) -> bool:
        # Allow addition of floats as LineSearchPoints
        if isinstance(point, float):
            # Will throw exception if T is incompatible
            point = LineSearchPoint(point)
        elif not isinstance(point, LineSearchPoint):
            raise TypeError(f'Point must be a LineSearchPoint or float, got {type(point)}')
        # end if
        # Only add the point if it is not already in the grid
        if point not in self:
            self._grid.append(point)
            # Keep the grid sorted
            self._grid.sort()
            return True
        else:
            return False
        # end if
    # end def

    def __contains__(self, point: T) -> bool:
        # Return True if a point with the same offset is already in the grid
        return any([point == p for p in self._grid])
    # end def

    def __len__(self):
        return len(self.grid)
    # end def

# end class
