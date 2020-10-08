.. _mdps:

Making a three stage transformation
===================================

If the alchemical transformation only involves the growth or the
If alchemical transofrmatiuon
It is recommended to split the tranformation into three stages **AlchemicalITP** is to make the topology file for the
transformation from one molecule (state A) to another molecule (state B).
Certain rules are applied, when making this transformation.
Separate topologies can be generated for splitting the alchemical transformation
into three parts (qon, vdw, qoff). ::

    >>> top_A, top_B = new.split_coul()
    >>> top_A.write('glh2glu.qoff_vdw.itp')
    >>> top_B.write('glh2glu.vdw_qon.itp')

Where the `glh2glu.qoff_vdw.itp` annihilates the partial charge of the atoms
which will be dummy in state B and the `glh2glu.vdw_qon.itp` recharges the
partial charge of the atoms which are dummy in state A.