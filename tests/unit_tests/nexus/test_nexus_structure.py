#!/usr/bin/env python

from nexus import Structure
from stalk.nexus.nexus_structure import NexusStructure
from stalk.params.parameter_structure import ParameterStructure

from ..assets.h2o import backward_H2O, elem_H2O, forward_H2O, pos_H2O
from ..assets.test_jobs import DummySimulation

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test NexusStructure class
def test_NexusStructure(tmp_path):

    # empty init (matches ParameterStructure)
    s = NexusStructure()
    assert isinstance(s, ParameterStructure)
    assert not s.analyzed
    assert not s.generated
    assert not s.finished
    assert s.jobs is None
    assert isinstance(s.get_nexus_structure(), Structure)

    # Meaningful init
    s = NexusStructure(
        label='label',
        pos=pos_H2O,
        elem=elem_H2O,
        forward=forward_H2O,
        backward=backward_H2O,
        units='A'
    )
    assert isinstance(s.get_nexus_structure(), Structure)

    # Test resetting of jobs and value
    job1 = DummySimulation()
    job2 = DummySimulation()
    s.jobs = [job1, job2]
    assert s.generated
    # Not finished until all jobs are finished
    assert not s.finished
    job1.finished = True
    assert not s.finished
    job2.finished = True
    assert s.finished
    # check analyzed and reset value
    assert not s.analyzed
    s.value = 1.0
    assert s.analyzed
    s.reset_value()
    assert s.jobs is None
    assert s.value is None

    # Test copy
    s.jobs = [job1, job2]
    s_copy = s.copy()
    assert s_copy.jobs is None

# end def
