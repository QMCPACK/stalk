#!/usr/bin/env python3

import sys
import numpy as np

from stalk.io import XyzGeometry
from stalk import ParameterStructure

from unit_tests.assets.test_jobs import efilename, xyzfilename, axesfilename
from unit_tests.assets.h2o import pes_H2O, pos_H2O, elem_H2O
from unit_tests.assets.diamond import pos_diamond, axes_diamond, elem_diamond, pes_diamond

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


if __name__ == '__main__':
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'r') as fhandle:
            # The first line points to structure file
            struct_name = fhandle.readline().replace("\n", "")
            pos = XyzGeometry({'suffix': struct_name}).load('.').get_pos()
            # The second line points to one the hardcoded pes functions
            pes_variable = fhandle.readline().replace("\n", "")
        # end with
        # Evaluate according to preset modes
        if pes_variable == 'h2o':
            # If evaluating PES, write energy and errorbar to disk
            value, error = pes_H2O(pos)
            np.savetxt(efilename, [value, error])
        elif pes_variable == 'relax_h2o':
            value, error = pes_H2O(pos)
            np.savetxt(efilename, [value, error])
            structure = ParameterStructure(pos=pos_H2O, elem=elem_H2O)
            writer = XyzGeometry({'suffix': xyzfilename})
            writer.write(structure, '.')
        elif pes_variable == 'relax_diamond':
            value, error = pes_diamond(axes_diamond[0, 0])
            np.savetxt(efilename, [value, error])
            structure = ParameterStructure(pos=pos_diamond, elem=elem_diamond)
            writer = XyzGeometry({'suffix': xyzfilename})
            writer.write(structure, '.')
            np.savetxt(axesfilename, axes_diamond)
        elif pes_variable == 'evm':
            # If evaluating PES, write energy and errorbar to disk
            value, error = pes_H2O(pos, dvar_eff=[0.1, 0.2])
            np.savetxt(efilename, [value, error])
        else:  # default: dummy
            pass
        # end if
    # end if
# end if
