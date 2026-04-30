#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings

from nexus import PwscfAnalyzer

from stalk.io.geometry_loader import GeometryLoader
from stalk.params.geometry_result import GeometryResult


class PwscfGeometry(GeometryLoader):

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        suffix='relax.in',
        **kwargs
    ):
        my_args = {'suffix': suffix}
        my_args.update(**args, **kwargs)
        super().__init__(**my_args)
    # end def

    def _load(self, filename: str, **kwargs):
        ai = PwscfAnalyzer(str(filename), **kwargs)
        ai.analyze()

        if not hasattr(ai, "structures") or len(ai.structures) == 0:
            warnings.warn(f"PwscfGeometry loader could not find structures in {str(filename)}")
            return GeometryResult(None, None)
        # end if
        final_structure = ai.structures[len(ai.structures) - 1]
        pos = final_structure.positions
        if 'axes' in final_structure:
            axes = final_structure.axes
        else:
            axes = None
        # end if
        return GeometryResult(pos, axes)
    # end def

# end class
