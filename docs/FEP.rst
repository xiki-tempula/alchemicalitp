Generation of topology file for FEP calculations
================================================

For an alchemical transformation from state A, which, for example, is protonated
glutamate `GLH.top` to state B, deprotonated glutamate `GLU.top`. The topology
file from protonated glutamate (state A) to deprotonated glutamate (state B)
can be generated accoridng to the :ref:`rules <rules>` with the code. ::

    >>> import alchemicalitp
    >>> from alchemicalitp.crd import merge_crd
    >>> from alchemicalitp.tutorial_files import (glh_top, glh_crd,
    >>>                                           glu_top, glu_crd,)

    >>> state_A = alchemicalitp.top.Topology(filename=glh_top)
    >>> state_B = alchemicalitp.top.Topology(filename=glu_top)
    >>> new, mapping = state_A.add_stateB(state_B,
    >>>                                   [19, 20, ],
    >>>                                   [19, None,])
    >>> new.write('glh2glu.top')
    >>> new.write('glh2glu.itp')
    >>> merged_crd = merge_crd(glh_crd, glu_crd, mapping,
    >>>                        filename='merged.gro')

Where the `glh2glu.top` containes the *defaults* and *Atomtypes*, which
includes **dummy** atomtype, while `glh2glu.itp` contains *moleculetype*,
*Atoms*, *Bonds*, *Pairs*, *Angles* and *Dihedrals* with state B written for
performing alchemical transformation.

The merged coordinate file corresponding to the merged topology file is saved
as `merged.gro`.

Input Coordinate and Topology
-----------------------------
The topology and coordinate files are read through their filenames. ::

    >>> from alchemicalitp.tutorial_files import (glh_top, glh_crd,
    >>>                                           glu_top, glu_crd,)
    >>> print(glu_top)
    alchemicalitp/data/fep_example/GLU.top
    >>> print(glh_top)
    alchemicalitp/data/fep_example/GLH.top
    >>> state_A = alchemicalitp.top.Topology(filename=glh_top)
    >>> state_B = alchemicalitp.top.Topology(filename=glu_top)
    >>> print(glu_crd)
    alchemicalitp/data/fep_example/GLU.gro
    >>> print(glh_crd)
    alchemicalitp/data/fep_example/GLH.gro

Note that the topology parser doesn't support recursive loading yet and all the
lines starting with `#` cannot be parsed. Position restraint such as ::

    #ifdef POSRES
    #include "posre.itp"
    #endif

Should be removed before passing into the topology parser.

It is recommended that the user should included the *atomtypes* directive
in the topology file as **AlchemicalITP** would be able to concatenate the
*atomtypes* from both topologies and remove the duplicates.

**AlchemicalITP** will also tried to add a new dummy atomtype to the
*atomtypes* directive. If the *atomtypes* is not provided by the user, the user
should add the dummy atomtype to the *atomtypes* directive by
themselves. ::

    [ Atomtypes ]
    ;  name at.num mass     charge ptype sigma epsilon
    DUM     0      1.008000 0.0    A     0.0   0.0

Mapping between state A and state B
-----------------------------------
To link the second topology to the first topology, a mapping scheme is
required. The atoms that has the same *residue name*, *residues number* and the
same *atom name* would be assumed to be directly matched and is not required
to be written out explicitly. Examples: ::

    # The 7th atom from topology of protonated glutamate
    7          N      2    GLH      N      7 -0.51112000  14.010000   ; qtot -0.511120
    # The 7th atom from topology of deprotonated glutamate
    7          N      2    GLU      N      7 -0.70584000  14.010000   ; qtot -0.705840

Would not require explicit mapping.

For the cases, where one atom in the first topology should be alchemically
changed to another atom in the second topology, the mapping should be specified
explicitly. The 19th atom from the first topology should be mapped to the 19th
atom from the second topology. However, since they have different *atom atoms*,
explicit mapping is required. ::

    # The 19th atom from topology of protonated glutamate
    19         OH      2    GLH     OH     19 -0.53240000  16.000000   ; qtot -0.455710
    # The 7th atom from topology of deprotonated glutamate
    19         O3      2    GLU    OE2     19 -0.84375000  16.000000   ; qtot -0.941410

Note the original oxygen atom in the protonated glutamate is named *OE2*, I
have changed it to *OH* to show that if the atom name is different, explicit
mapping is required.

For the cases, where the atom is present in one topology but should be a dummy
in the other topology. Explicit mapping is also required. ::

    # The 20th atom from topology of protonated glutamate
    20         HO      2    GLH    HE2     20 0.40610000   1.008000   ; qtot -0.049610

Thus, we would map the 19th oxygen from the first topology to the 19th oxygen
from the second topology. The 20th hydrogen would be mapped to a dummy atom in
the state B, which gives ::

    >>> new, mapping = state_A.add_stateB(state_B,
    >>>                                   [19, 20, ],
    >>>                                   [19, None,])

I'm interested in implementing the `maximum common substructure <http://rdkit.org/docs/source/rdkit.Chem.MCS.html>`_
to replace this explicit mapping scheme, which is currently working under
progress.

For Splitting the Transformation into Three Parts
-------------------------------------------------

If both the growth and the annihilation of the atoms are required a
:ref:`three-stages transformation <mdps>` might be useful.

Validation of the Generated Topology
------------------------------------
I have tested the project thoroughly to make sure that the conversion is correct.
However, the user is also recommended to test their own system to make sure
that the conversion is correct.

To test if the generated topology is sensible, we need to check the potential
of the state A and state B from the generated topology and compare them
with the original topology. Since the discrepancy in potential comes from the
bonded potential from the dummy atoms, we need to energy minimise the merged
structure. ::

    >>> # Save the merged coordinate with a large box
    >>> merged_crd.dimensions = [100, 100, 100, 90, 90, 90]
    >>> merged_crd.atoms.write('merged.gro')

    >>> # Prepare the mdp file for energy minimisation
    >>> from alchemicalitp.tutorial_files import (mdp_em0,     mdp_em1,
    >>>                                           mdp_energy0, mdp_energy1)
    >>> import shutil
    >>> shutil.copy(mdp_em0, './')
    './minim0.mdp'
    >>> shutil.copy(mdp_em1, './')
    './minim1.mdp'
    >>> shutil.copy(mdp_energy0, './')
    './test0.mdp'
    >>> shutil.copy(mdp_energy1, './')
    './test1.mdp'

Add the following lines to the end of 'glh2glu.top'. ::

    #include "glh2glu.itp"

    [ system ]
    alchemicalitp

    [ molecules ]
    system12system1 1

and run energy minimisation. ::

    gmx grompp -f minim0.mdp -c merged.gro -o em_0.tpr -p glh2glu.top -maxwarn -1
    gmx mdrun -deffnm em_0
    gmx grompp -f minim1.mdp -c merged.gro -o em_1.tpr -p glh2glu.top -maxwarn -1
    gmx mdrun -deffnm em_1

Extract the energy minimisated state A and state B from the merged file. ::

    >>> from alchemicalitp.crd import extract_from_merged
    >>> extract_from_merged('em_0.gro', mapping[0], filename='stateA.gro')
    >>> extract_from_merged('em_1.gro', mapping[1], filename='stateB.gro')

Extract the potential of the state A and state B with their original topology.

    >>> shutil.copy(glh_top, './')
    './GLH.top'
    >>> shutil.copy(glu_top, './')
    './GLU.top'

Used another mdp file to extract energy through rerun. ::

    gmx grompp -f test0.mdp -c stateA.gro -o original_A.tpr -p GLH.top -maxwarn -1
    gmx mdrun -deffnm original_A -rerun stateA.gro
    gmx grompp -f test0.mdp -c stateB.gro -o original_B.tpr -p GLU.top -maxwarn -1
    gmx mdrun -deffnm original_B -rerun stateB.gro

    gmx grompp -f test0.mdp -c em_0.gro -o energy_A.tpr -p glh2glu.top -maxwarn -1
    gmx mdrun -deffnm energy_A -rerun em_0.gro
    gmx grompp -f test1.mdp -c em_1.gro -o energy_B.tpr -p glh2glu.top -maxwarn -1
    gmx mdrun -deffnm energy_B -rerun em_1.gro

Compare the potential from the original topology (`original_A.log` and
`original_A.log`) to the new potential (`energy_A.log` and `energy_B.log`).

.. list-table:: Potential Energy (kj/mol)
   :header-rows: 1

   * - Name
     - Bond
     - Angle
     - Proper Dih.
     - Improper Dih.
     - LJ (SR)
     - Coulomb (SR)
     - Potential
   * - Original State A
     - 5.46155
     - 16.7446
     - -58.7648
     - 1.01055
     - -2.06010
     - -407.300
     - -248.931
   * - New State A
     - 5.46155
     - 16.7446
     - -58.7648
     - 1.01055
     - -2.06010
     - -407.300
     - -248.931
   * - Original State B
     - 5.07483
     - 9.77376
     - -20.7308
     - 96.0576
     - -3.01319
     - -717.158
     - -98.1315
   * - New State B
     - 5.11319
     - 9.77407
     - -41.0334
     - 96.1632
     - -3.01318
     - -717.158
     - -118.395

The non-bonded interactions (*LJ (SR)* and *Coulomb (SR)*) should be the same.
A very small deviation might be observed for *bond* (< 1 kj/mol), a medium
deviation can be observed for the *angle* (< 10 kj/mol). Due to the flexibility
of the *dihedral* potential, a large deviation can be observed (<100 kj/mol).
However, as long as the energy difference do not exceed 1000 kj/mol, things
should be fine.