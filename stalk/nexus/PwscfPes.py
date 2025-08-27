#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import nan
from pathlib import Path

from nexus import PwscfAnalyzer

from stalk.params.ParameterSet import ParameterSet
from stalk.params.PesResult import PesResult
from stalk.util.util import PL
from stalk.io.PesLoader import PesLoader


class PwscfPes(PesLoader):

    def __init__(self, args={}):
        self._func = None
        self.args = args
    # end def

    def _load(self, structure: ParameterSet, suffix='scf.in', **kwargs):
        input_file = Path(PL.format(structure.file_path, suffix))
        # Testing existence here, because Nexus will shut down everything upon failure
        if input_file.exists():
            ai = PwscfAnalyzer(str(input_file), **kwargs)
            ai.analyze()
        else:
            warnings.warn(f"PwscfPes loader could not find {str(input_file)}. Returning None.")
            return PesResult(nan)
        # end if

        if not hasattr(ai, "E") or ai.E == 0.0:
            # Analysis has failed
            warnings.warn(f"PwscfPes loader could not find energy in {str(input_file)}. Returning None.")
            E = nan
        else:
            E = ai.E
        # end if
        Err = 0.0
        return PesResult(E, Err)
    # end def

# end class
