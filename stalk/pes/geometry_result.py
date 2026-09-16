#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


from numpy import ndarray


class GeometryResult():
    pos: ndarray = None
    axes: ndarray | None = None
    elem: list[str] = None

    def __init__(self, pos, axes=None, elem=None):
        self.pos = pos
        self.axes = axes
        self.elem = elem
    # end def

    def get_pos(self):
        return self.pos
    # end def

    def get_axes(self) -> ndarray | None:
        return self.axes
    # end def

    def get_elem(self) -> list[str] | None:
        return self.elem
    # end def

    def get_result(self) -> tuple[ndarray, ndarray | None]:
        return self.get_pos(), self.get_axes()
    # end def

    def rescale(self, scale) -> None:
        if self.pos is not None:
            self.pos /= scale
        # end if
        if self.axes is not None:
            self.axes /= scale
        # end if
    # end def

# end class
