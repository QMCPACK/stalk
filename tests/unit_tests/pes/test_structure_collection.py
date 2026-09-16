#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises

from stalk.pes.structure_collection import StructureCollection
from stalk.params.linesearch_point import LineSearchPoint
from stalk.params.parameter_set import ParameterSet


# Test StructureCollection class
def test_StructureCollection():
    # Test initialization
    sc = StructureCollection()
    assert sc.grid == []
    assert sc.collect_enabled() == []
    assert sc.thr == 1e-6

    # Test adding points
    p1 = LineSearchPoint(0.0, 1.0, 0.1)
    p2 = LineSearchPoint(1.0, 2.0, 0.2)
    p3 = LineSearchPoint(0.5, 1.5, 0.15)
    assert sc.add_point(p1) is True
    assert sc.add_point(p2) is True
    assert sc.add_point(p3) is True
    assert sc.grid == [p1, p3, p2]  # Should be sorted by offset

    # Test adding duplicate point
    assert sc.add_point(p1) is False
    assert len(sc.grid) == 3

    # Test collect_enabled
    p2.enabled = False
    enabled_points = sc.collect_enabled()
    assert enabled_points == [p1, p3]

    # Test addition of floats
    assert sc.add_point(2.0) is True
    sc2 = StructureCollection[ParameterSet]()
    sc2.add_point(ParameterSet([1.0]))
    # Addition of something else raises exception
    with raises(TypeError):
        sc2.add_point([])
    # end with
# end def
