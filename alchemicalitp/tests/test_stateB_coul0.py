import pytest
import alchemicalitp
from pkg_resources import resource_filename
import os
@pytest.fixture
def urea():
    return alchemicalitp.top.Topology(filename=resource_filename(__name__, 'example/urea.itp'))

def test_stateB_coul0(urea):
    '''Test if add the same state B but only change the charge is fine'''
    topology = urea
    atom = topology.content_dict['atoms'][3]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
               4, 'H', 1, 'URE', 'H11', 4, '0.395055', '1.00800', '', '', '')
    new_top = topology.add_coul0()
    atom = new_top.content_dict['atoms'][3]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
               4, 'H', 1, 'URE', 'H11', 4, '0.395055', '1.00800', 'H', 0, '1.00800')

def test_extract_stateB_coul0(urea):
    topology = urea
    new_top = topology.add_coul0()
    new_top = new_top.to_stateB()
    atom = new_top.content_dict['atoms'][3]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
               4, 'H', 1, 'URE', 'H11', 4, 0, '1.00800', '', '', '')
