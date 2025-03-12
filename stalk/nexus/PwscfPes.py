#!/usr/bin/env python3

from nexus import PwscfAnalyzer

from stalk.params.PesResult import PesResult
from stalk.util.util import PL
from stalk.io.PesLoader import PesLoader


class PwscfPes(PesLoader):

    def __init__(self, args={}):
        self._func = None
        self.args = args
    # end def

    def _load(self, path, suffix='scf.in', **kwargs):
        ai = PwscfAnalyzer(PL.format(path, suffix), **kwargs)
        ai.analyze()
        E = ai.E
        Err = 0.0
        return PesResult(E, Err)
    # end def

# end class
