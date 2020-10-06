"""
Unit and regression test for the alchemicalitp package.
"""

# Import package, test suite, and other packages as needed
import alchemicalitp
import pytest
from pkg_resources import resource_filename

@pytest.fixture
def glu():
    glu = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'GLU.top'))
    return glu

@pytest.fixture
def glh():
    glh = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'GLH.top'))
    return glh

def test_non_index(glu):
    assert glu.content_dict['atoms'][10000] is None

def test_unequal_field(glu, glh):
    assert glu.content_dict['atoms'] != glh.content_dict['atoms']