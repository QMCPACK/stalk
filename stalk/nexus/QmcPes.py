#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import nan
from pathlib import Path

from nexus import QmcpackAnalyzer

from stalk.nexus.NexusStructure import NexusStructure
from stalk.params.PesResult import PesResult
from stalk.io.PesLoader import PesLoader
from stalk.util.util import PL


class QmcPes(PesLoader):

    def __init__(self, args={}):
        self._func = None
        self.args = args
    # end def

    def _load(
        self,
        structure: NexusStructure,
        qmc_idx=1,
        suffix='dmc/dmc.in.xml',
        term='LocalEnergy',
        **kwargs  # e.g. equilibration=None
    ):
        input_file = Path(PL.format(structure.file_path, suffix))
        # Testing existence here, because Nexus will shut down everything upon failure
        if input_file.exists():
            ai = QmcpackAnalyzer(str(input_file), **kwargs)
            ai.analyze()
        else:
            warnings.warn(f"QmcPes loader could not find {str(input_file)}. Returning None.")
            return PesResult(nan)
        # end if

        if not hasattr(ai, "qmc") or len(ai.qmc) < qmc_idx or not hasattr(ai.qmc[qmc_idx], "scalars"):
            # Analysis has failed
            warnings.warn(f"QmcPes loader could not find energy in {str(input_file)}. Returning None.")
            return PesResult(nan)
        else:
            return self._analyze_energy_term(ai.qmc[qmc_idx].scalars, term)
        # end if
    # end def

    def _analyze_energy_term(self, scalars, label):
        LE = scalars[label]
        value = LE.mean
        error = LE.error
        return PesResult(value, error)
    # end def

# end class
