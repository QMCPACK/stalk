#!/usr/bin/env python3
"""Surrogate Theory Accelerated Line-search Kit"""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"
__version__ = "0.2.2dev3"

# fit module
from stalk.fit import FittingFunction
from stalk.fit import FittingResult
from stalk.fit import MorseFit
from stalk.fit import MorseResult
from stalk.fit import PolynomialFit
from stalk.fit import PolynomialResult
from stalk.fit import SplineFit
from stalk.fit import SplineResult
# io module
from stalk.io import FilesPes
from stalk.io import GeometryLoader
from stalk.io import GeometryWriter
from stalk.io import XyzGeometry
from stalk.io import TxtData
from stalk.io import write_xyz_sigma
# ls module
from stalk.ls import LineSearch
from stalk.ls import LineSearchBase
from stalk.ls import LineSearchGrid
from stalk.ls import LsSettings
from stalk.ls import TargetLineSearch
from stalk.ls import TargetLineSearchBase
from stalk.ls import TlsSettings
# lsi module
from stalk.lsi import LineSearchIteration
from stalk.lsi import PathwayImage
from stalk.lsi import TransitionPathway
# nexus module
from stalk.nexus import nexus_enabled
if nexus_enabled:
    from stalk.nexus import NexusGeometry
    from stalk.nexus import NexusPes
    from stalk.nexus import NexusStructure
    from stalk.nexus import PwscfGeometry
    from stalk.nexus import PwscfPes
    from stalk.nexus import QmcPes
    from stalk.nexus import XsfGeometry
# end if
# params module
from stalk.params import EffectiveVariance
from stalk.params import EffectiveVarianceMap
from stalk.params import LineSearchPoint
from stalk.params import Parameter
from stalk.params import BondLength
from stalk.params import BondAngle
from stalk.params import PhaseAngle
from stalk.params import ParameterHessian
from stalk.params import ParameterMapping
from stalk.params import ParameterSet
from stalk.params import ParameterStructure
from stalk.params import angle
from stalk.params import bond_angle
from stalk.params import distance
from stalk.params import interpolate_params
from stalk.params import mean_distances
from stalk.params import mean_bond_angles
from stalk.params import mean_param
from stalk.params import periodic_distance
from stalk.params import periodic_bond_angle
from stalk.params import rotate_2d
# pes module
from stalk.pes import GeometryResult
from stalk.pes import PesFunction
from stalk.pes import PesLoader
from stalk.pes import PesResult
from stalk.pes import RelaxFunction
from stalk.pes import StructureCollection
# pls module
from stalk.pls import ParallelLineSearch
from stalk.pls import Surrogate
# util module
from stalk.util import ArgsContainer
from stalk.util import FunctionCaller
from stalk.util import AbsNoise
from stalk.util import Noise
from stalk.util import NoiseFactory
from stalk.util import WhiteNoise
from stalk.util import morse

# Make practical alias
TargetParallelLineSearch = Surrogate

__all__ = [
    # fit module
    'FittingFunction',
    'FittingResult',
    'MorseFit',
    'MorseResult',
    'PolynomialFit',
    'PolynomialResult',
    'SplineFit',
    'SplineResult',
    # io module
    'FilesPes',
    'GeometryLoader',
    'GeometryWriter',
    'XyzGeometry',
    'TxtData',
    'write_xyz_sigma',
    # ls module
    'LineSearch',
    'LineSearchBase',
    'LineSearchGrid',
    'LsSettings',
    'TargetLineSearch',
    'TargetLineSearchBase',
    'TlsSettings',
    # lsi module
    'LineSearchIteration',
    'PathwayImage',
    'TransitionPathway',
    # nexus module
    'NexusGeometry',
    'NexusStructure',
    'NexusPes',
    'PwscfGeometry',
    'PwscfPes',
    'QmcPes',
    'XsfGeometry',
    # params module
    'EffectiveVariance',
    'EffectiveVarianceMap',
    'GeometryResult',
    'LineSearchPoint',
    'Parameter',
    'BondLength',
    'BondAngle',
    'PhaseAngle',
    'ParameterHessian',
    'ParameterMapping',
    'ParameterSet',
    'ParameterStructure',
    'angle',
    'bond_angle',
    'distance',
    'interpolate_params',
    'mean_distances',
    'mean_bond_angles',
    'mean_param',
    'periodic_distance',
    'periodic_bond_angle',
    'rotate_2d',
    # pes module
    'GeometryResult',
    'PesFunction',
    'PesLoader',
    'RelaxFunction',
    'PesResult',
    'StructureCollection',
    # pls module
    'ParallelLineSearch',
    'TargetParallelLineSearch',
    'Surrogate',
    # util module
    'ArgsContainer',
    'FunctionCaller',
    'AbsNoise',
    'Noise',
    'NoiseFactory',
    'WhiteNoise',
    'morse',
]
