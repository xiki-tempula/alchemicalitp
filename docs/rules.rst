.. _rules:

Rules in Making the Transformation
==================================

One usage of the **AlchemicalITP** is to make the topology file for the
transformation from one molecule (state A) to another molecule (state B).
Certain rules are applied, when making this transformation.

Atomtypes
---------

When the atomtype is given, **AlchemicalITP** will try to merge the *Atomtypes*
from both topologies and remove the duplicate *Atomtypes*.

An additional atomtype *dum* will be added to represent dummy atom with no LJ. ::

    [ Atomtypes ]
    ;  name at.num mass     charge ptype sigma epsilon
    DUM     0      1.008000 0.0    A     0.0   0.0

Atoms
-----

The atoms section will be treated differently based on the state of state B.

 - If the atom is the same in state B, the state B will be left blank
 - If the mass, partial charge or atomtype is changed in state B, state B will
   be set up explicitly
 - If the atom doesn't exist in state A or state B, a dummy atomtype with no
   partial change and the same mass as the other state will be set up

Bonded interactions
-------------------

The *bonds*, *pairs*, *angles* and *dihedrals* are set up, such that:

 - if an interaction is the same in state A and state B, state B will be left
   blank
 - If an interaction is present in both state A and state B, state B will be
   filled to reflect this change
 - If any interactions involve dummy atoms, the parameter when all the atoms
   are present will be used

cmap
----

If cmap is present, it can be read and write as it is but alchemical
transformation involving the cmap is not supported.

