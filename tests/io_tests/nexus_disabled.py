#!/usr/bin/env python

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


from pytest import raises


def disable_nexus():
    from sys import path
    # Temporarily remove nexus from path
    for p in path:
        if '/nexus' in p:
            path.remove(p)
        # end if
    # end for
# end def


# Test Nexus being disabled
if __name__ == '__main__':
    disable_nexus()
    from stalk import nexus
    assert not nexus.nexus_enabled
    # Trying to import specific Nexus classes raises ImportError
    with raises(ImportError):
        from stalk.nexus import NexusPes
        NexusPes()
    # end with
    print("All clear!")
# end if
