#!/usr/bin/env python3
"""Surrogate Theory Accelerated Line-search Kit"""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"
__version__ = "0.2.0"

from stalk import io
from stalk import ls
from stalk import lsi
from stalk import params
from stalk import pls
from stalk import util
from stalk.params import ParameterHessian
from stalk.params import ParameterSet
from stalk.params import ParameterStructure
from stalk.lsi import LineSearchIteration
from stalk.ls import LineSearch
from stalk.ls import TargetLineSearch
from stalk.pls import ParallelLineSearch
from stalk.pls import TargetParallelLineSearch
try:
    from stalk import nexus
    from stalk.nexus import NexusStructure
    from stalk.nexus import NexusPes
    nexus_enabled = True
except ModuleNotFoundError:
    # Nexus not found
    nexus_enabled = False
    pass
# end try


__all__ = [
    'io',
    'ls',
    'lsi',
    'params',
    'pls',
    'util',
    'ParameterHessian',
    'ParameterSet',
    'ParameterStructure',
    'LineSearchIteration',
    'LineSearch',
    'TargetLineSearch',
    'ParallelLineSearch',
    'TargetParallelLineSearch',
    'nexus',
    'NexusPes',
    'NexusStructure',
]
