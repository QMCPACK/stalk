#!/usr/bin/env python3

import os
import numpy as np
from numpy import nan
from pathlib import Path
from nexus import obj, job, settings, generate_physical_system, input_template
from nexus.simulation import GenericSimulation, SimulationAnalyzer, Simulation
from stalk.io.geometry_loader import GeometryLoader
from stalk.io.pes_loader import PesLoader
from stalk.io.xyz_geometry import XyzGeometry
from stalk.nexus.nexus_structure import NexusStructure
from stalk.params.geometry_result import GeometryResult
from stalk.params.pes_result import PesResult

# Add nexus tester to path
app_path = os.path.dirname(__file__) + "/../../nxs_tester.py"

# The hard-coded energy file name
efilename = 'test_energy.dat'
xyzfilename = 'relax.xyz'
axesfilename = 'axes.dat'

# First, the user must set up Nexus according to their computing environment.
testjob = obj(app_command=app_path + ' nxs_test.in', cores=1, ppn=1, serial=True)
nx_settings = obj(
    sleep=1,
    pseudo_dir='./',  # not used in testing
    runs='',  # not used in testing
    results='',  # not used in testing
    status_only=0,
    generate_only=0,
    command_line=False,  # disable for testing
    machine='ws4',
)


# Dummy simulation class to suppress actual identifier check-ups in testing
class DummySimulation(Simulation):
    finished = False

    def __init__(self):
        pass
    # end def
# end class


# Tailor Nexus analyzer for generic testing
class TestAnalyzer(SimulationAnalyzer):
    value = None
    error = 0.0
    pos = None
    axes = None

    def __init__(
        self,
        arg0
    ):
        if isinstance(arg0, Simulation):
            self.path = arg0.path
        elif isinstance(arg0, NexusStructure):
            self.path = arg0.file_path
        else:
            self.path = arg0
        # end if
    # end def

    def analyze(self):
        path = Path(self.path)
        if not path.exists():
            raise AssertionError('The job has not run yet!')
        # end if
        # In case pointing to final file
        if not path.is_dir():
            path = path.parent
        # end if
        e_file = path / efilename
        if e_file.exists():
            res = np.loadtxt(e_file)
            if len(res) == 1:
                self.value = res[0]
            if len(res) > 1:
                self.value = res[0]
                self.error = res[1]
            # end if
        else:
            raise AssertionError('The job has not run yet!')
        # end if
        xyz_file = path / xyzfilename
        axes_file = path / axesfilename
        if xyz_file.exists():
            res = XyzGeometry(suffix=xyzfilename).load(str(path))
            if axes_file.exists():
                res.axes = np.loadtxt(axes_file)
            # end if
            self.pos = res.get_pos()
            self.axes = res.get_axes()
    # end def
# end class


def test_dummy(*args, **kwargs):
    pass
# end def


# Tailor Nexus analyzer for generic testing
class TestLoader(PesLoader):

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        suffix='test_energy.dat',
        **kwargs
    ):
        my_args = {'suffix': suffix}
        my_args.update(**args, **kwargs)
        super().__init__(**my_args)
    # end def

    def _load(self, structure: NexusStructure, produce_fail=False, **kwargs):
        ai = TestAnalyzer(structure, **kwargs)
        ai.analyze()
        if produce_fail:
            return PesResult(nan, 0.0)
        else:
            return PesResult(ai.value, ai.error)
        # end if
    # end def

# end class


# Tailor Nexus analyzer for generic testing
class TestGeometryLoader(GeometryLoader):

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        suffix='relax.xyz',
        **kwargs
    ):
        my_args = {'suffix': suffix}
        my_args.update(**args, **kwargs)
        super().__init__(**my_args)
    # end def

    def _load(self, filename, produce_fail=False, **kwargs):
        ai = TestAnalyzer(filename, **kwargs)
        ai.analyze()
        if produce_fail:
            return GeometryResult(None, None)
        else:
            return GeometryResult(ai.pos, ai.axes)
        # end if
    # end def

# end class


def init_nexus():
    if len(settings) == 0:
        settings(**nx_settings)
    # end if
# end def


def nxs_generic_pes(
    structure,
    path,
    system_args={},
    pes_variable='dummy',
    sigma=None,
    **kwargs
):
    init_nexus()
    system = generate_physical_system(structure=structure, **system_args)
    # TODO: this could be done more elegantly using Nexus templates
    input = "nxs_test.struct.xyz\n" + pes_variable + '\n'
    job_input = input_template(input)
    sim = GenericSimulation(
        system=system,
        job=job(**testjob),
        path=path,
        input=job_input,
        identifier='nxs_test'
    )
    return [sim]
# end def
