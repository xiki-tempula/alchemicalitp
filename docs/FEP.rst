Generation of topology file for FEP calculations
================================================

For an alchemical transformation from state A, which, for example, is protonated
glutamate `GLH.top` to state B, deprotonated glutamate `GLU.top`. The topology file
from protonated glutamate (state A) to deprotonated glutamate (state B) can
be generated with the code. ::

    >>> import alchemicalitp
    >>> from alchemicalitp.tutorial_files import glu, glh

    >>> state_A = alchemicalitp.top.Topology(filename=glh)
    >>> state_B = alchemicalitp.top.Topology(filename=glu)
    >>> new, mapping = state_A.add_stateB(state_B,
    >>>                                   [19, 20, ],
    >>>                                   [19, None,])
    >>> new.write('glh2glu.top')
    >>> new.write('glh2glu.itp')

Where the `glh2glu.top` containes the *defaults* and *Atomtypes*, which
includes **dummy** atomtype, while `glh2glu.itp` contains *moleculetype*,
*Atoms*, *Bonds*, *Pairs*, *Angles* and *Dihedrals* with state B written for
performing alchemical transformation.

Input Topology
--------------
The topology files are read through their filenames. ::

    >>> from alchemicalitp.tutorial_files import glu, glh
    >>> print(glu)
    alchemicalitp/data/fep_example/GLU.top
    >>> print(glh)
    alchemicalitp/data/fep_example/GLU.top
    >>> state_A = alchemicalitp.top.Topology(filename=glh)
    >>> state_B = alchemicalitp.top.Topology(filename=glu)

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
Separate topologies can be generated for splitting the alchemical transformation
into three parts (qon, vdw, qoff). ::

    >>> top_A, top_B = new.split_coul()
    >>> top_A.write('glh2glu.qoff_vdw.itp')
    >>> top_B.write('glh2glu.vdw_qon.itp')

Where the `glh2glu.qoff_vdw.itp` annihilates the partial charge of the atoms
which will be dummy in state B and the `glh2glu.vdw_qon.itp` recharges the
partial charge of the atoms which are dummy in state A.

Validation of the Generated Topology
------------------------------------
I have tested the project thoroughly to make sure that the conversion is correct.
However, the user is also recommended to test their own system to make sure
that the conversion is correct.

