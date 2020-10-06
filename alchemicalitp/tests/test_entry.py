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

def test_atom_unequal(glu):
    assert glu.content_dict['atoms'].content[2] != glu.content_dict['atoms'].content[3]

def test_pair(glu):
    pair = glu.content_dict['pairs'].content[2]
    assert pair.idx_in([1, 2]) == False
    assert 11 in pair



def test_angle(glu):
    angle = glu.content_dict['angles'].content[2]
    assert angle.idx_in([1, 10]) == False

def test_dihedral(glu):
    assert glu.content_dict['dihedrals'].content[2] != \
           glu.content_dict['dihedrals'].content[3]
    assert 11 in glu.content_dict['dihedrals'].content[2]
    assert glu.content_dict['dihedrals'].content[93] == \
           glu.content_dict['dihedrals'].content[93]
