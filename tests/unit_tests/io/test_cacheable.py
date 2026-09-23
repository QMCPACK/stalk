#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises

from stalk.io.cacheable import Cacheable
from stalk.io.txt_data import TxtData
from stalk.util.util import match_to_tol


# Test Cacheable class
def test_Cacheable(tmp_path):

    # Test empty init
    cache = Cacheable()
    assert isinstance(cache._cache, dict)
    assert len(cache._cache) == 0
    cache.save(tmp_path)  # Should not raise any exceptions
    cache.load(tmp_path)  # Should not raise any exceptions

    # Test init with invalid data type
    with raises(TypeError):
        cache = Cacheable(value=1.0)  # Not a TxtData instance
    # end with

    # Test init with valid TxtData
    cache = Cacheable(
        value=TxtData('value.dat'),
        value2=TxtData('value2.dat'),
    )
    assert len(cache._cache) == 2
    # Loading nonexistent should return dict with None values
    result = cache.load(tmp_path / 'load')
    assert result['value'] is None
    assert result['value2'] is None

    # Writing and reading back values
    value = 1.234
    value2 = [[2, 3, 4, 5], [6, 7, 8, 9]]
    cache.save(tmp_path / 'save', value=value, value2=value2)

    cache2 = Cacheable(
        newvalue=TxtData('value.dat'),
        newvalue2=TxtData('value2.dat'),
    )
    result = cache2.load(tmp_path / 'save')
    assert match_to_tol(result['newvalue'], value)
    # The second field is a scaled TxtData
    assert match_to_tol(result['newvalue2'], value2)
    # Cannot save a nonexisting value
    with raises(KeyError):
        cache.save(tmp_path / 'save', key='nonexisting')
    # end with
    # Cannot load a nonexisting value
    with raises(KeyError):
        cache.load(tmp_path / 'save', key='nonexisting')
    # end with

    # Overwrite test to one field
    new_value = 9.876
    cache.save(tmp_path / 'save', overwrite=False, value=new_value)
    # Without overwrite, matches the old value
    result = cache2.load(tmp_path / 'save', key='newvalue')
    assert match_to_tol(result, value)
    # With overwrite, matches the new value
    cache.save(tmp_path / 'save', overwrite=True, value=new_value)
    result = cache2.load(tmp_path / 'save', key='newvalue')
    assert match_to_tol(result, new_value)

    # Write only one field to another directory
    value = 4
    cache.save(tmp_path / 'save2', value=value)
    # Load all fields
    result = cache.load(tmp_path / 'save2', default=[])
    assert match_to_tol(result['value'], value)
    assert result['value2'] == []  # The other field not written, so should be the default

    # Add a new cache field
    cache.add_cache('value3', TxtData('value3.dat'))
    assert len(cache._cache) == 3
    # Use add_cache to overwrite an existing field
    cache.add_cache('value', TxtData('value_new.dat'))
    assert len(cache._cache) == 3  # Should still be 3, not 4

# end def
