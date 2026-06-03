#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from nexus import PwscfAnalyzer

from stalk.params.pes_function import NotEvaluatedException
from stalk.params.pes_result import PesResult
from stalk.io.pes_loader import PesLoader


class PwscfPes(PesLoader):
    _suffix = 'scf.in'

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        scale=1.0,
        **kwargs,
    ):
        args.update(**kwargs)
        # suffix = None means we'll use class-level default
        suffix = args.pop('suffix', None)
        PesLoader.__init__(self, suffix=suffix, scale=scale, **args)
    # end def

    def _load(self, path: str, **kwargs) -> PesResult:
        p = self.get_filename(path)
        if p.exists():
            ai = PwscfAnalyzer(str(p), **kwargs)
            ai.analyze()
        else:
            raise NotEvaluatedException(f"PwscfPes could not find {p}. Raising exception.")
        # end if

        if not hasattr(ai, "E") or ai.E == 0.0:
            # Analysis has failed
            raise NotEvaluatedException(f"PwscfPes loader could not find energy in {p}. Raising exception.")
        else:
            E = ai.E
        # end if
        Err = 0.0
        result = PesResult(E, Err)
        result.rescale(self.scale)
        return result
    # end def

# end class
