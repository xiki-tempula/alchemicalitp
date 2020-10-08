import os
import pytest
from pkg_resources import resource_filename

import MDAnalysis as mda

import alchemicalitp
from alchemicalitp.crd import merge_crd, extract_from_merged
from alchemicalitp.tutorial_files import (glh_top, glh_crd,
    glu_top, glu_crd,)

@pytest.fixture
def mapping():
    state_A = alchemicalitp.top.Topology(filename=glu_top)
    state_B = alchemicalitp.top.Topology(filename=glh_top)
    new, mapping = state_A.add_stateB(state_B,
                                      [19, None, ],
                                      [19, 20, ],)
    return mapping

def test_merge_crd(mapping):
    merged_crd = merge_crd(glu_crd, glh_crd, mapping, filename='merged.gro')
    os.remove('merged.gro')
    assert len(merged_crd.atoms) == len(mda.Universe(glh_crd).atoms)

def test_extract_from_merged(mapping):
    merged_crd = merge_crd(glu_crd, glh_crd, mapping, filename='merged.gro')
    extracted = extract_from_merged('merged.gro', mapping[0], filename='test.gro')
    os.remove('test.gro')
    os.remove('merged.gro')
    assert len(extracted.atoms) == len(mda.Universe(glu_crd).atoms)
