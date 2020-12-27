"""
Unit and regression test for the alchemicalitp package.
"""

# Import package, test suite, and other packages as needed
import alchemicalitp
import pytest
from pkg_resources import resource_filename

@pytest.fixture
def dum_to():
    glu = alchemicalitp.top.Topology(filename=resource_filename(__name__, 'GLU.top'))
    glh = alchemicalitp.top.Topology(filename=resource_filename(__name__, 'GLH.top'))
    glu2glh, mapping = glu.add_stateB(glh, [19, None,], [19, 20,])
    return glu, glh, glu2glh, mapping

@pytest.fixture
def to_dum():
    glu = alchemicalitp.top.Topology(filename=resource_filename(__name__, 'GLU.top'))
    glh = alchemicalitp.top.Topology(filename=resource_filename(__name__, 'GLH.top'))
    glh2glu, mapping = glh.add_stateB(glu, [19, 20, ], [19, None, ])
    return glu, glh, glh2glu, mapping

def test_defaults(dum_to):
    glu, glh, glu2glh, mapping = dum_to
    defaults = glu2glh.content_dict['defaults']
    assert (defaults.nbfunc, defaults.comb_rule, defaults.gen_pairs, defaults.fudgeLJ, defaults.fudgeQQ) == (
        '1', '2', 'yes', '0.5', '0.83333333')

def test_mapping_dum_to(dum_to):
    glu, glh, glu2glh, mapping = dum_to
    # Test if full
    mapping1, mapping2 = mapping
    for key in mapping2:
        assert mapping2[key] == key
    for key in mapping1:
        if key < 20:
            assert mapping1[key] == key
        else:
            assert mapping1[key] == key + 1

def test_mapping_to_dum(to_dum):
    glu, glh, glh2glu, mapping = to_dum
    # Test if full
    mapping1, mapping2 = mapping
    for key in mapping1:
        assert mapping1[key] == key
    for key in mapping2:
        if key < 20:
            assert mapping2[key] == key
        else:
            assert mapping2[key] == key + 1

def test_atomtypes(dum_to):
    glu, glh, glu2glh, mapping = dum_to
    # Test total number
    assert len(glu2glh.content_dict['atomtypes']) == 18

def test_unchanged_atom(dum_to):
    glu, glh, glu2glh, mapping = dum_to
    assert glu.content_dict['atoms'][1] == glu2glh.content_dict['atoms'][1]
    assert glh.content_dict['atoms'][1] == glu2glh.content_dict['atoms'][1]

def test_changed_atom(dum_to):
    glu, glh, glu2glh, mapping = dum_to
    atom = glu2glh.content_dict['atoms'][18]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        19, 'O3', 2, 'GLU', 'OE2', 19, '-0.84375000', '16.000000', 'OH', '-0.53240000', '16.000000')

def test_dum2_atom(dum_to):
    glu, glh, glu2glh, mapping = dum_to
    atom = glu2glh.content_dict['atoms'][19]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        20, 'DUM', 2, 'GLH', 'HE2', 20, 0.0, '1.008000', 'HO', '0.40610000', '1.008000')

def test_2dum_atom(to_dum):
    glu, glh, glh2glu, mapping = to_dum
    atom = glh2glu.content_dict['atoms'][19]
    assert (atom.nr, atom.type, atom.resnr, atom.residue, atom.atom, atom.cgnr, atom.charge, atom.mass,
            atom.typeB, atom.chargeB, atom.massB) == (
        20, 'HO', 2, 'GLH', 'HE2', 20, '0.40610000', '1.008000', 'DUM', 0.0, '1.008000')

def test_unchanged_atom(dum_to):
    glu, glh, glu2glh, mapping = dum_to
    assert glu.content_dict['atoms'][1] == glu2glh.content_dict['atoms'][1]
    assert glh.content_dict['atoms'][1] == glu2glh.content_dict['atoms'][1]

def test_unchanged_bonds(to_dum):
    # test that bonds are unchanged
    glu, glh, glh2glu, mapping = to_dum
    # Test empty pairs
    assert glh2glu.content_dict['bonds'][1] == glu.content_dict['bonds'][1]
    assert glh2glu.content_dict['bonds'][1] == glh.content_dict['bonds'][1]

def test_changed_bonds(to_dum):
    # test that bonds are changed
    glu, glh, glh2glu, mapping = to_dum
    bond = glh.content_dict['bonds'][18]
    assert (bond.i, bond.j, bond.func, bond.b0, bond.kb, bond.b0B, bond.kbB) == (
        17, 19, 1, '0.13640', '376560.000000', '', '')
    bond = glu.content_dict['bonds'][18]
    assert (bond.i, bond.j, bond.func, bond.b0, bond.kb, bond.b0B, bond.kbB) == (
        17, 19, 1, '0.12500', '548940.800000', '', '')
    bond = glh2glu.content_dict['bonds'][18]
    assert (bond.i, bond.j, bond.func, bond.b0, bond.kb, bond.b0B, bond.kbB) == (
        17, 19, 1, '0.13640', '376560.000000', '0.12500', '548940.800000')

def test_2dum_bonds(to_dum):
    # test that bonds change to dummy
    glu, glh, glh2glu, mapping = to_dum
    bond = glh.content_dict['bonds'][19]
    assert (bond.i, bond.j, bond.func, bond.b0, bond.kb, bond.b0B, bond.kbB) == (
        19, 20, 1, '0.09600', '462750.400000', '', '')
    bond = glu.content_dict['bonds'][19]
    assert (bond.i, bond.j, bond.func, bond.b0, bond.kb, bond.b0B, bond.kbB) == (
        21, 22, 1, '0.12290', '476976.000000', '', '')
    bond = glh2glu.content_dict['bonds'][19]
    assert (bond.i, bond.j, bond.func, bond.b0, bond.kb, bond.b0B, bond.kbB) == (
        19, 20, 1, '0.09600', '462750.400000', '', '')

def test_dum2_bonds(dum_to):
    # test that bonds change to dummy
    glu, glh, glu2glh, mapping = dum_to
    bond = glu2glh.content_dict['bonds'][19]
    assert (bond.i, bond.j, bond.func, bond.b0, bond.kb, bond.b0B, bond.kbB) == (
        19, 20, 1, '0.09600', '462750.400000', '', '')

def test_unchanged_pairs(to_dum):
    # test that bonds are unchanged
    glu, glh, glh2glu, mapping = dum_to
    # Test empty pairs
    assert glh2glu.content_dict['pairs'][1] == glu.content_dict['pairs'][1]
    assert glh2glu.content_dict['pairs'][1] == glh.content_dict['pairs'][1]

def test_unchanged_pairs(to_dum):
    # test that bonds are unchanged
    glu, glh, glh2glu, mapping = to_dum
    # Test empty pairs
    assert glh2glu.content_dict['pairs'][1] == glu.content_dict['pairs'][1]
    assert glh2glu.content_dict['pairs'][1] == glh.content_dict['pairs'][1]

def test_2dum_pairs(to_dum):
    # test that bonds are unchanged
    glu, glh, glh2glu, mapping = to_dum
    # Test empty pairs
    assert glh2glu.content_dict['pairs'][43] == glh.content_dict['pairs'][43]

def test_dum2_pairs(dum_to):
    # test that bonds are unchanged
    glu, glh, glu2glh, mapping = dum_to
    # Test empty pairs
    assert glu2glh.content_dict['pairs'][43] == glh.content_dict['pairs'][43]

def test_unchanged_angles(to_dum):
    # test that bonds are unchanged
    glu, glh, glh2glu, mapping = to_dum
    # Test empty pairs
    assert glh2glu.content_dict['angles'][1] == glu.content_dict['angles'][1]
    assert glh2glu.content_dict['angles'][1] == glh.content_dict['angles'][1]

def test_changed_angles(to_dum):
    # test that bonds are changed
    glu, glh, glh2glu, mapping = to_dum
    angle = glh.content_dict['angles'][9]
    assert (angle.i, angle.j, angle.k, angle.func, angle.th0, angle.cth, angle.th0B, angle.cthB) == (
        5, 7, 9, 1, '119.6100511', '411.090552', '', '')
    angle = glu.content_dict['angles'][9]
    assert (angle.i, angle.j, angle.k, angle.func, angle.th0, angle.cth, angle.th0B, angle.cthB) == (
        5, 7, 9, 1, '122.5600526', '376.057083', '', '')
    angle = glh2glu.content_dict['angles'][9]
    assert (angle.i, angle.j, angle.k, angle.func, angle.th0, angle.cth, angle.th0B, angle.cthB) == (
        5, 7, 9, 1, '119.6100511', '411.090552', '122.5600526', '376.057083')

def test_2dum_angles(to_dum):
    # test that bonds are changed
    glu, glh, glh2glu, mapping = to_dum
    angle = glh.content_dict['angles'][34]
    assert (angle.i, angle.j, angle.k, angle.func, angle.th0, angle.cth, angle.th0B, angle.cthB) == (
        17, 19, 20, 1, '113.0000484', '418.400000', '', '')
    angle = glu.content_dict['angles'][34]
    assert (angle.i, angle.j, angle.k, angle.func, angle.th0, angle.cth, angle.th0B, angle.cthB) == (
        18, 17, 19, 1, '126.0000540', '669.440000', '', '')
    angle = glh2glu.content_dict['angles'][34]
    assert (angle.i, angle.j, angle.k, angle.func, angle.th0, angle.cth, angle.th0B, angle.cthB) == (
        17, 19, 20, 1, '113.0000484', '418.400000', '', '')

def test_dum2_angles(dum_to):
    # test that bonds are changed
    glu, glh, glu2glh, mapping = dum_to
    angle = glu2glh.content_dict['angles'][34]
    assert (angle.i, angle.j, angle.k, angle.func, angle.th0, angle.cth, angle.th0B, angle.cthB) == (
        17, 19, 20, 1, '113.0000484', '418.400000', '', '')

def test_unchanged_dihedrals(to_dum):
    # test that bonds are unchanged
    glu, glh, glh2glu, mapping = to_dum
    # Test empty pairs
    assert glh2glu.content_dict['dihedrals'][1] == glu.content_dict['dihedrals'][1]
    assert glh2glu.content_dict['dihedrals'][1] == glh.content_dict['dihedrals'][1]

def test_changed_dihedrals(to_dum):
    # test that bonds are changed
    glu, glh, glh2glu, mapping = to_dum
    dihedral = glh.content_dict['dihedrals'][11]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) ==(5, 7, 9, 10, 1, 0.0, '1.6182875', 2, '', '', '')
    dihedral = glu.content_dict['dihedrals'][11]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (5, 7, 9, 10, 1, 0.0, '2.3499018', 2, '', '', '')
    dihedral = glh2glu.content_dict['dihedrals'][11]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (5, 7, 9, 10, 1, 0.0, '1.6182875', 2, 0.0, '2.3499018', 2)

def test_dihedral2dum(to_dum):
    # test the dihedral exits in one top but doesn't exist in the next due to the absence of the dummy atom
    # Thus, to ensure that the dummy atom stay in its place, the dihedral needs to be enforced
    glu, glh, glh2glu, mapping = to_dum
    dihedral = glh.content_dict['dihedrals'][72]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) ==(14, 17, 19, 20, 1, 180.0000771, '10.1273302', 2, '', '', '')
    dihedral = glu.content_dict['dihedrals'][88]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (14, 18, 17, 19, 4, 180.0000771, '43.9320000', 2, '', '', '')
    dihedral = glh2glu.content_dict['dihedrals'][72]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (14, 17, 19, 20, 1, 180.0000771, '10.1273302', 2, '', '', '')

def test_dum2dihedral(dum_to):
    # the reverse of the test_dihedral2dum, where the dihedral is not present in the first top due to the
    # atom not being there
    # Thus, to ensure that the dummy atom stay in its place, the dihedral needs to be enforced
    glu, glh, glu2glh, mapping = dum_to
    dihedral = glu2glh.content_dict['dihedrals'][72]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (14, 17, 19, 20, 1, 180.0000771, '10.1273302', 2, '', '', '')

def test_dihedral2zero(to_dum):
    # The dihedral in the first is in absence in the second top due to topology difference
    # should be set to zero in this case
    glu, glh, glh2glu, mapping = to_dum
    dihedral = glh.content_dict['dihedrals'][74]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) ==(15, 14, 17, 18, 1, 0.0, '-0.5486898', 2, '', '', '')
    dihedral = glu.content_dict['dihedrals'][74]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (16, 14, 17, 18, 1, 0.0, '-0.5294015', 2, '', '', '')
    dihedral = glh2glu.content_dict['dihedrals'][74]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (15, 14, 17, 18, 1, 0.0, '-6.0489343', 1, 0.0, '0.0', 1)

def test_zero2dihedral(dum_to):
    # The inverse of test_dihedral2zero where where the dihedral is in absence in the first top
    # should set to zero
    glu, glh, glu2glh, mapping = dum_to
    dihedral = glu2glh.content_dict['dihedrals'][73]
    assert (dihedral.i, dihedral.j, dihedral.k, dihedral.l, dihedral.func, dihedral.phase, dihedral.kd, dihedral.pn,
            dihedral.phaseB, dihedral.kdB, dihedral.pnB) == (15, 14, 17, 18, 1, 0.0, '0.0', 1, 0.0, '-6.0489343', 1)

def test_cmap():
    state_A = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'cmap/state_A.top'))
    state_B = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'cmap/state_B.top'))
    state, mapping = state_A.add_stateB(state_B, [None,], [219,])
    print(len(state.content_dict['cmaptypes'].content[0].data))
    assert len(state.content_dict['cmaptypes'].content[0].data) == len(state_A.content_dict['cmaptypes'].content[0].data)
    assert len(state.content_dict['cmaptypes'].content[0].data) == len(
        state_B.content_dict['cmaptypes'].content[0].data)

    #A number in state_C cmap is changed so that it is different
    state_A = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'cmap/state_A.top'))
    state_C = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'cmap/state_C.top'))
    with pytest.raises(NotImplementedError):
        state, mapping = state_A.add_stateB(state_C, [None, ], [219, ])

def test_alchem():
    # Test the case charlie's case where the state B has a atom
    # Which is later in the index
    state_A = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'example/expa.itp'))
    state_B = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'example/expb.itp'))
    new, mapping = state_A.add_stateB(state_B, [22, 23, 46, 32, 34, 36, 38, 39, 40],
                             [None, None, 25, None, None, None, None, None,
                              None, ])
    assert new.content_dict['atoms'][45].type == 'nh'
    assert new.content_dict['atoms'][45].typeB == 'f'

    mapping1, mapping2 = mapping
    for key in mapping1:
        assert mapping1[key] == key
    assert mapping2[22] == 24
    assert mapping2[26] == 27
    assert mapping2[31] == 33
    assert mapping2[25] == 46
