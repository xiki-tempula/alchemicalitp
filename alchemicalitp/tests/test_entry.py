"""
Unit and regression test for the alchemicalitp package.
"""

# Import package, test suite, and other packages as needed
import alchemicalitp
import pytest
from pkg_resources import resource_filename

def test_atom_unequal():
    glu = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'GLU.top'))
    assert glu.content_dict['atoms'].content[2] != glu.content_dict['atoms'].content[3]