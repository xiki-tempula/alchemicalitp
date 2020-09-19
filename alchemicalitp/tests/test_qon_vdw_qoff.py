"""
Unit and regression test for the alchemicalitp package.
"""

# Import package, test suite, and other packages as needed
import alchemicalitp
import pytest
from pkg_resources import resource_filename

@pytest.fixture
def aejn():
    topology = alchemicalitp.top.Topology(filename=resource_filename(__name__, 'AEJN.itp'))
    top_A, top_B = topology.split_coul()
    return topology, top_A, top_B

@pytest.fixture
def ajen():
    topology = alchemicalitp.top.Topology(filename=resource_filename(__name__, 'AJEN.itp'))
    top_A, top_B = topology.split_coul()
    return topology, top_A, top_B

def test_unchanged_atom(aejn):
    topology, top_A, top_B = aejn
    atom = top_A.content_dict['atoms'][1]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (1, 'HC', 1, 'ACE', 'HH31', 1, '0.112300', '1.0080', '', '', '')

def test_changed_atom(aejn):
    # test only charge changes
    topology, top_A, top_B = aejn
    # Test phase 1
    atom = top_A.content_dict['atoms'][7]
    print()
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        7, 'N', 2, 'E2J', 'N', 7, '-0.516300', '14.0100', 'N', '-0.466000', '14.0100')
    # Test phase 2
    atom = top_B.content_dict['atoms'][7]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        7, 'N', 2, 'E2J', 'N', 7, '-0.466000', '14.0100', 'N', '-0.415700', '14.0100')

def test_type_changed_atom(aejn):
    # test types change as well
    topology, top_A, top_B = aejn
    # Test phase 1
    atom = top_A.content_dict['atoms'][9]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        9, 'CT', 2, 'E2J', 'CA', 9, '0.039700', '12.0100', 'CT', '0.027100', '12.0100')
    # Test phase 2
    atom = top_B.content_dict['atoms'][9]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        9, 'CT', 2, 'E2J', 'CA', 9, '0.027100', '12.0100', 'CT', '0.014500', '12.0100')

def test_dum2_atom(aejn):
    # test from dummy to a real atom
    topology, top_A, top_B = aejn
    # Test phase 1
    atom = top_A.content_dict['atoms'][22]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        22, 'DUM_HO', 2, 'E2J', 'DHE2', 22, '0.000000', '1.0000', 'HO', 0, '1.0080')
    # Test phase 2
    atom = top_B.content_dict['atoms'][22]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        22, 'DUM_HO', 2, 'E2J', 'DHE2', 22, '0.000000', '1.0000', 'HO', '0.464100', '1.0080')

def test_2dum_atom(ajen):
    # test from dummy to a real atom
    topology, top_A, top_B = ajen
    # Test phase 1
    atom = top_A.content_dict['atoms'][20]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        20, 'HO', 2, 'J2E', 'HE2', 20, '0.464100', '1.0080', 'DUM_HO', '0.000000', '1.0000')
    # Test phase 2
    atom = top_B.content_dict['atoms'][20]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        20, 'HO', 2, 'J2E', 'HE2', 20, 0, '1.0080', 'DUM_HO', '0.000000', '1.0000')

def test_bonds(aejn):
    # test that bonds are unchanged
    topology, top_A, top_B = aejn
    # Test empty bond
    assert topology.content_dict['bonds'][1] == top_A.content_dict['bonds'][1]
    assert topology.content_dict['bonds'][1] == top_B.content_dict['bonds'][1]
    # Test bond with value
    assert topology.content_dict['bonds'][10] == top_A.content_dict['bonds'][10]
    assert topology.content_dict['bonds'][10] == top_B.content_dict['bonds'][10]

def test_pairs(aejn):
    # test that bonds are unchanged
    topology, top_A, top_B = aejn
    # Test empty pairs
    assert topology.content_dict['pairs'][1] == top_A.content_dict['pairs'][1]
    assert topology.content_dict['pairs'][1] == top_B.content_dict['pairs'][1]

def test_angles(aejn):
    # test that bonds are unchanged
    topology, top_A, top_B = aejn
    # Test empty angles
    assert topology.content_dict['angles'][1] == top_A.content_dict['angles'][1]
    assert topology.content_dict['angles'][1] == top_B.content_dict['angles'][1]
    # Test angles with value
    assert topology.content_dict['angles'][10] == top_A.content_dict['angles'][10]
    assert topology.content_dict['angles'][10] == top_B.content_dict['angles'][10]

def test_dihedrals(aejn):
    # test that bonds are unchanged
    topology, top_A, top_B = aejn
    # Test empty dihedrals
    assert topology.content_dict['dihedrals'][1] == top_A.content_dict['dihedrals'][1]
    assert topology.content_dict['dihedrals'][1] == top_B.content_dict['dihedrals'][1]
    # Test dihedrals with value
    assert topology.content_dict['dihedrals'][10] == top_A.content_dict['dihedrals'][10]
    assert topology.content_dict['dihedrals'][10] == top_B.content_dict['dihedrals'][10]