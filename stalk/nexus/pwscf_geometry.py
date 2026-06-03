#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from nexus import PwscfAnalyzer

from stalk.io.geometry_loader import GeometryLoader
from stalk.params.geometry_result import GeometryResult
from stalk.params.pes_function import NotEvaluatedException


class PwscfGeometry(GeometryLoader):
    _suffix = 'relax.in'

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        scale=1.0,
        c_pos=None,
        **kwargs
    ):
        GeometryLoader.__init__(self, args=args, scale=scale, c_pos=c_pos, **kwargs)
    # end def

    def _load(self, path: str, **kwargs):
        p = self.get_filename(path)
        if p.exists():
            ai = PwscfAnalyzer(str(p), **kwargs)
            ai.analyze()
        else:
            raise NotEvaluatedException(f"PwscfGeometry could not find {p}. Raising exception.")
        # end if

        if not hasattr(ai, "structures") or len(ai.structures) == 0:
            raise NotEvaluatedException(f"PwscfGeometry could not find analyze {p}. Raising exception.")
        # end if
        final_structure = ai.structures[len(ai.structures) - 1]
        pos = final_structure.positions / self.scale
        if 'axes' in final_structure:
            axes = final_structure.axes / self.scale
        else:
            axes = None
        # end if
        return GeometryResult(pos, axes)
    # end def

# end class
