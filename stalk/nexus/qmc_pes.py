#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import nan, array
from pathlib import Path

from nexus import QmcpackAnalyzer

from stalk.params.pes_result import PesResult
from stalk.io.pes_loader import PesLoader


class QmcPes(PesLoader):

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        suffix='dmc/dmc.in.xml',
        **kwargs
    ):
        my_args = {'suffix': suffix}
        my_args.update(**args, **kwargs)
        super().__init__(**my_args)
    # end def

    def _load(
        self,
        filename,
        qmc_idx=1,
        term='LocalEnergy',
        twist_averaging=False,
        twist_weights=None,
        **kwargs  # e.g. equilibration=None
    ) -> PesResult:
        # Testing existence here, because Nexus will shut down everything upon failure
        p = Path(filename)
        if p.exists():
            ai = QmcpackAnalyzer(str(p), **kwargs)
            ai.analyze()
        else:
            warnings.warn(f"QmcPes loader could not find {str(p)}. Returning NaN.")
            return PesResult(nan)
        # end if

        if twist_averaging and self._check_bundled(ai):
            return self._perform_twist_averaging(ai, qmc_idx, term, twist_weights)
        else:
            if not hasattr(ai, "qmc") or len(ai.qmc) < qmc_idx or not hasattr(ai.qmc[qmc_idx], "scalars"):
                # Analysis has failed
                warnings.warn(f"QmcPes loader could not find energy in {str(p)}. Returning NaN.")
                return PesResult(nan)
            else:
                return self._analyze_energy_term(ai.qmc[qmc_idx].scalars, term)
            # end if
    # end def

    def _check_bundled(self, ai: QmcpackAnalyzer):
        if not hasattr(ai, "bundled_analyzers") or ai.bundled_analyzers is None:
            warnings.warn("QmcpackAnalyzer could not find twist bundles. Reverting to non-twist energy.")
            return False
        else:
            return True
        # end if
    # end def

    def _analyze_energy_term(self, scalars, label) -> PesResult:
        LE = scalars[label]
        value = LE.mean
        error = LE.error
        return PesResult(value, error)
    # end def

    def _perform_twist_averaging(self, ai: QmcpackAnalyzer, qmc_idx, label, twist_weights):
        if twist_weights is None:
            twist_weights = array(len(ai.bundled_analyzers) * [1])
        # end if
        weighted_sum = 0.0
        weighted_error2 = 0.0
        weight = 0.0
        for analyzer, w in zip(ai.bundled_analyzers, twist_weights):
            res = self._analyze_energy_term(analyzer.qmc[qmc_idx].scalars, label)
            weighted_sum += w * res.value
            weighted_error2 += w * res.error**2
            weight += w
        # end for
        weighted_sum /= weight
        weighted_error = weighted_error2**0.5 / weight
        return PesResult(weighted_sum, weighted_error)
    # end def

# end class
