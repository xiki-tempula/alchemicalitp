'''
To test how does the charge conservation work for different cases
To make the initial pdb files
tleap
> source leaprc.protein.ff14SB
> lysine = sequence {ACE LYS NME}
> savepdb lysine, AKN.pdb
> valine = sequence {ACE VAL NME}
> savepdb valine, AVN.pdb
> glutamate = sequence {ACE GLU NME}
> savepdb glutamate, AEN.pdb

# fix the name to the amber14sb naming
for aa in AKN AVN AEN
do
gmx pdb2gmx -f $aa.pdb -o $aa.gro -water tip3p -ff amber14sb -ignh
done

# make the pmx files
export GMXLIB=/Users/XXX/GitHub/pmx/src/pmx/data/mutff
Do the K2Z, K2A, K2D, V2Y, E2D
pmx mutate -f AKN.gro -o AKZN.mutant.pdb -ff amber14sbmut
gmx pdb2gmx -f AKZN.mutant.pdb -o AKZN.mutant.gro
pmx gentop -p topol.top -o AKZN.top -ff amber14sbmut

'''
import alchemicalitp
import numpy as np
from pkg_resources import resource_filename

def test_glu2asp():
    topology = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'charge_conservation/AEDN.top'))
    top_A, top_B = topology.split_coul(charge_conservation='nearest')
    assert np.isclose(sum([float(atom.charge) for atom in top_B.content_dict['atoms'] if atom.charge]), 0, atol=0.001)

def test_val2tyr():
    topology = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'charge_conservation/AVYN.top'))
    top_A, top_B = topology.split_coul(charge_conservation='nearest')
    assert np.isclose(sum([float(atom.charge) for atom in top_B.content_dict['atoms'] if atom.charge]), 0, atol=0.001)

def test_lys2ala():
    topology = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'charge_conservation/AKAN.top'))
    top_A, top_B = topology.split_coul(charge_conservation='nearest')
    assert np.isclose(sum([float(atom.charge) for atom in top_B.content_dict['atoms'] if atom.charge]), 0, atol=0.001)

def test_lys2asp():
    topology = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'charge_conservation/AKDN.top'))
    top_A, top_B = topology.split_coul(charge_conservation='nearest')
    assert np.isclose(sum([float(atom.charge) for atom in top_B.content_dict['atoms'] if atom.charge]), 0, atol=0.001)

def test_lys2hip():
    topology = alchemicalitp.top.Topology(
        filename=resource_filename(__name__, 'charge_conservation/AKZN.top'))
    top_A, top_B = topology.split_coul(charge_conservation='nearest')
    assert np.isclose(sum([float(atom.charge) for atom in top_B.content_dict['atoms'] if atom.charge]), 0, atol=0.001)