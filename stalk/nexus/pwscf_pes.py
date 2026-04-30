#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import nan

from nexus import PwscfAnalyzer

from stalk.params.pes_result import PesResult
from stalk.io.pes_loader import PesLoader


class PwscfPes(PesLoader):

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        suffix='scf.in',
        **kwargs
    ):
        my_args = {'suffix': suffix}
        my_args.update(**args, **kwargs)
        super().__init__(**my_args)
    # end def

    def _load(self, filename: str, **kwargs) -> PesResult:
        ai = PwscfAnalyzer(filename, **kwargs)
        ai.analyze()
        if not hasattr(ai, "E") or ai.E == 0.0:
            # Analysis has failed
            warnings.warn(f"PwscfPes loader could not find energy in {filename}. Returning None.")
            E = nan
        else:
            E = ai.E
        # end if
        Err = 0.0
        return PesResult(E, Err)
    # end def

# end class
