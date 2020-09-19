"""
Unit and regression test for the alchemicalitp package.
"""

# Import package, test suite, and other packages as needed
import alchemicalitp
import pytest
import os
from pkg_resources import resource_filename

@pytest.fixture
def urea():
    return alchemicalitp.top.Topology(filename=resource_filename(__name__, 'urea.itp'))

def test_defaults(urea):
    topology = urea
    defaults = topology.content_dict['defaults']
    assert (defaults.nbfunc, defaults.comb_rule, defaults.gen_pairs, defaults.fudgeLJ, defaults.fudgeQQ) == (
        '1', '2', 'yes', '0.5000', '0.8333')

def test_atomtypes(urea):
    topology = urea
    atomtype = topology.content_dict['atomtypes'][1]
    assert (atomtype.at_num, atomtype.name, atomtype.mass, atomtype.charge, atomtype.ptype, atomtype.sigma, atomtype.epsilon) == (
        'n4', 'n4', '0.00000', '0.00000', 'A', '3.25000e-01', '7.11280e-01')

def test_moleculetype(urea):
    topology = urea
    moleculetype = topology.content_dict['moleculetype']
    assert (moleculetype.name, moleculetype.nrexcl) == ('URE', '3')

def test_atoms(urea):
    topology = urea
    atom = topology.content_dict['atoms'][0]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
               1, 'C', 1, 'URE', 'C', 1, '0.880229', '12.01000', '', '', '')
    atom = topology.content_dict['atoms'][6]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        7, 'N', 2, 'E2J', 'N', 7, '-0.516300', '14.0100', 'N', '-0.415700', '14.0100')


def test_pairs(urea):
    topology = urea
    pair = topology.content_dict['pairs'][0]
    assert (pair.i, pair.j, pair.func) == (1,7,1)

def test_bonds(urea):
    topology = urea
    bond = topology.content_dict['bonds'][0]
    assert (bond.i, bond.j, bond.func, bond.b0, bond.kb, bond.b0B, bond.kbB) == (1,2,1,'','','','')
    bond = topology.content_dict['bonds'][1]
    assert (bond.i, bond.j, bond.func, bond.b0, bond.kb, bond.b0B, bond.kbB) == (1,3,1,'0.09572','502416.0','','')
    bond = topology.content_dict['bonds'][2]
    assert (bond.i, bond.j, bond.func, bond.b0, bond.kb, bond.b0B, bond.kbB) == (1,4,1,'0.09572','502416.0','0.09572','502416.0')

def test_comment(urea):
    topology = urea
    comment = topology.content_dict['angles'][0]
    assert comment.comment == ' i     j       k       funct   angle   force_constant'

def test_angles(urea):
    topology = urea
    angle = topology.content_dict['angles'][1]
    assert (angle.i, angle.j, angle.k, angle.func, angle.th0, angle.cth, angle.th0B, angle.cthB) == (
        2, 1, 3, 1, '', '', '', '')
    angle = topology.content_dict['angles'][2]
    assert (angle.i, angle.j, angle.k, angle.func, angle.th0, angle.cth, angle.th0B, angle.cthB) == (
        2, 1, 6, 1, '104.52', '628.02', '', '')
    angle = topology.content_dict['angles'][3]
    assert (angle.i, angle.j, angle.k, angle.func, angle.th0, angle.cth, angle.th0B, angle.cthB) == (
        3, 1, 6, 1, '104.52', '628.02', '104.52', '628.02')

def test_dihedrals(urea):
    topology = urea
    dihedral = topology.content_dict['dihedrals'][1]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (2, 1, 3, 4, 9, '', '', '', '', '', '')
    dihedral = topology.content_dict['dihedrals'][2]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (2, 1, 3, 5, 9, 180.0, '0.56484', 4, '', '', '')
    dihedral = topology.content_dict['dihedrals'][3]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (2, 1, 6, 7, 9, 357.2, '1.48473', 3, 180.0, '1.72381', 3)
    dihedral = topology.content_dict['dihedrals'][4]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func,
            dihedral.C0, dihedral.C1, dihedral.C2, dihedral.C3, dihedral.C4, dihedral.C5,
            dihedral.C0B, dihedral.C1B, dihedral.C2B, dihedral.C3B, dihedral.C4B, dihedral.C5B) == \
           (2, 1, 6, 8, 3, 3.68192, -4.35136, 0.00000, 1.33888, 0.00000, 0.00000,
            '', '', '', '', '', '')
    dihedral = topology.content_dict['dihedrals'][5]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (2, 1, 6, 3, 4, '', '', '', '', '', '')
    dihedral = topology.content_dict['dihedrals'][6]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (1, 4, 3, 5, 4, 180.0, '0.56484', 4, '', '', '')
    dihedral = topology.content_dict['dihedrals'][7]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (1, 7, 6, 8, 4, 357.2, '1.48473', 3, 180.0, '1.72381', 3)

def test_write(urea):
    topology = urea
    topology.write('test.itp')
    test_result = alchemicalitp.top.Topology(filename='test.itp')
    assert topology.content_dict['atoms'] == test_result.content_dict['atoms']
    assert topology.content_dict['bonds'] == test_result.content_dict['bonds']
    assert topology.content_dict['pairs'] == test_result.content_dict['pairs']
    assert topology.content_dict['angles'] == test_result.content_dict['angles']
    assert topology.content_dict['dihedrals'] == test_result.content_dict['dihedrals']

def test_writeitp(urea):
    topology = urea
    topology.write('test.itp', format = 'itp')
    test_result = alchemicalitp.top.Topology(filename='test.itp')
    assert topology.content_dict['atoms'] == test_result.content_dict['atoms']
    assert topology.content_dict['bonds'] == test_result.content_dict['bonds']
    assert topology.content_dict['pairs'] == test_result.content_dict['pairs']
    assert topology.content_dict['angles'] == test_result.content_dict['angles']
    assert topology.content_dict['dihedrals'] == test_result.content_dict['dihedrals']
    os.remove('test.itp')

def test_writeitp(urea):
    topology = urea
    topology.write('test.top', format = 'top')
    test_result = alchemicalitp.top.Topology(filename='test.top')
    assert topology.content_dict['defaults'] == test_result.content_dict['defaults']
    assert topology.content_dict['atomtypes'] == test_result.content_dict['atomtypes']
    os.remove('test.top')