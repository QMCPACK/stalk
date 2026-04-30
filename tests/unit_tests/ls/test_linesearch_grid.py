#!/usr/bin/env python

from stalk.ls.linesearch_grid import LineSearchGrid
from stalk.util.util import match_to_tol

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test LineSearchGrid class
def test_LineSearchGrid():

    # Empty init
    empty_grid = LineSearchGrid()
    assert len(empty_grid) == 0
    assert len(empty_grid.grid) == 0
    assert len(empty_grid.offsets) == 0
    assert len(empty_grid.values) == 0
    assert len(empty_grid.errors) == 0
    assert not empty_grid.evaluated

    # Init with 1 data
    grid1_in = [1.0]
    grid1 = LineSearchGrid(grid1_in)
    assert len(grid1) == 1
    assert match_to_tol(grid1.offsets, grid1_in)
    assert len(grid1.values) == 1
    assert grid1.values[0] is None
    assert len(grid1.errors) == 1
    assert grid1.errors[0] == 0.0
    assert len(grid1.valid_grid) == 0
    assert len(grid1.valid_offsets) == 0
    assert len(grid1.valid_values) == 0
    assert len(grid1.valid_errors) == 0
    assert not grid1.evaluated
    # Make valid by setting value
    values1_in = [2.0]
    errors1_in = [3.0]
    grid1.values = values1_in
    grid1.errors = errors1_in
    assert match_to_tol(grid1.values, values1_in)
    assert match_to_tol(grid1.errors, errors1_in)
    assert match_to_tol(grid1.valid_offsets, grid1_in)
    assert match_to_tol(grid1.valid_values, values1_in)
    assert match_to_tol(grid1.valid_errors, errors1_in)
    assert grid1.evaluated

    # Init with unsorted data
    grid3_in = [0.0, -1.0, 1.0]
    grid3_sorted = [-1.0, 0.0, 1.0]
    grid3 = LineSearchGrid(grid3_in)
    assert len(grid3) == 3
    assert match_to_tol(grid3.offsets, grid3_sorted)
    assert not grid3.evaluated

# end def
