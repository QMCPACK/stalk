#!/usr/bin env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path
from numpy import loadtxt, savetxt, nan

from stalk.io.xyz_geometry import XyzGeometry
from stalk.params.parameter_set import ParameterSet
from stalk.params.pes_result import PesResult
from stalk.util.util import directorize


def write_xyz_sigma(
    structure: ParameterSet,
    suffix='structure.xyz',
    sigma=None,
    sigma_suffix='sigma.dat',
    **kwargs
):
    g = XyzGeometry()
    if sigma is not None:
        sigmafilename = directorize(structure.file_path) + sigma_suffix
        savetxt(sigmafilename, [sigma])
    # end if
    g.write(structure=structure, path=structure.file_path, suffix=suffix)
# end def


def load_energy(filename):
    if Path(filename).exists():
        data = loadtxt(filename)
        print(f"Loaded {filename}")
        result = PesResult(data[0], data[1])
    else:
        print(f"Waiting for {filename}")
        # Return null result
        result = PesResult(nan)
    # end if
    return result
# end def
