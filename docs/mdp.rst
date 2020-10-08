.. _mdps:

Sequential Transformation
=========================

Atoms being annihilated
-----------------------
If the alchemical transformation only involves the growth **or** the annihilation
of atoms (there are only dummy atoms on one side of the transformation). A
single merged topology can be used. ::

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

In this case, the deprotonation of the glutamate means that a hydrogen atom is
annihilated in state B. To avoid electrostatic interactions between the dummy
atom and the rest of the molecule, the mdp files need to be set up such that
the partial charge of the hydrogen atom is turned off before changing the
hydrogen atom to a dummy with no Van der Waals interactions. ::

    free-energy              = yes
    separate-dhdl-file       = yes
    sc-alpha                 = 0.5
    sc-power                 = 1
    sc-sigma		         = 0.3
    init-lambda-state = 0
    coul-lambdas = 0.0 0.1 0.2 ... 0.9 1.0 1.0  1.0 1.0  ... 1.0  1.0
    fep-lambdas =  0.0 0.0 0.0 ... 0.0 0.0 0.05 0.1 0.15 ... 0.95 1.0
    nstdhdl                  = 1000
    calc-lambda-neighbors    = -1

Where a 0.1 spacing was used for coul-lambdas to change the partial charge and
0.05 spacing was used for changing everything else.

Growth of Atoms
---------------
If the transformation is inverted and instead of annihilating an atom, an atom
is being grown in state B, such as the protonation of the glutamate. ::

    >>> import alchemicalitp
    >>> from alchemicalitp.crd import merge_crd
    >>> from alchemicalitp.tutorial_files import (glh_top, glh_crd,
    >>>                                           glu_top, glu_crd,)

    >>> state_A = alchemicalitp.top.Topology(filename=glu_top)
    >>> state_B = alchemicalitp.top.Topology(filename=glh_top)
    >>> new, mapping = state_A.add_stateB(state_B,
    >>>                                   [19, None,],
    >>>                                   [19, 20, ],)
    >>> new.write('glu2glh.top')
    >>> new.write('glu2glh.itp')


The mdp files should be set up so that the hydrogen atom is recharged after
the Van der Waals transformation. ::

    free-energy              = yes
    separate-dhdl-file       = yes
    sc-alpha                 = 0.5
    sc-power                 = 1
    sc-sigma		         = 0.3
    init-lambda-state = 0
    coul-lambdas = 0.0 0.0  0.0 ... 0.0  0.0 0.1 0.2 ... 0.9 1.0
    fep-lambdas =  0.0 0.05 0.1 ... 0.95 1.0 1.0 1.0 ... 1.0 1.0
    nstdhdl                  = 1000
    calc-lambda-neighbors    = -1

Three-steps Transformation
--------------------------
However, if the atoms are being grown and annihilated at the same time, a
three-steps transformation is required.

1. The partial charge of the atoms, that will become dummy in state B, are
   annihilated.
2. Non-charge related transformation are performed.
3. The partial charge of the atoms, that are dummy in state A, are recharged.

Thus, separate topology files are needed for this three-steps transformation.::

    >>> top_A, top_B = new.split_coul()
    >>> top_A.write('glh2glu.qoff_vdw.itp')
    >>> top_B.write('glh2glu.vdw_qon.itp')

Splitting the topology file means that an intermediate partial charge stage is
introduced. By default, for the atoms that are present in both state A and
state B, the partial charge of this intermediate state is the average of the
partial charge from state A and state B.

However, due to the absence of the charge from dummy atoms, the total charge
of the molecule can be a non-integer value. Setting the *charge_conservation*
to `'nearest'` will round the total charge to the nearest integer. ::

    >>> top_A, top_B = new.split_coul(charge_conservation='nearest')
    >>> top_A.write('glh2glu.qoff_vdw.itp')
    >>> top_B.write('glh2glu.vdw_qon.itp')

The rounding is done by first taking the average of the partial charge between
state A and state B. The partial charge is then shifted in proportional to the
absolute difference between state A and state B to make the total charge into
an integer.

Depending on which part needs more sampling, the step 2 can be bundled with
either step 1 or step 3. For example, if the crystal structure of state A is
known, the step 2 can be bundled to step 3 to achieve more sampling in the
state B.

The mdp for `glh2glu.qoff_vdw.itp` will cover the step 1. ::

    free-energy              = yes
    separate-dhdl-file       = yes
    sc-alpha                 = 0.5
    sc-power                 = 1
    sc-sigma		         = 0.3
    init-lambda-state = 0
    coul-lambdas = 0.0 0.1 0.2 ... 0.9 1.0
    fep-lambdas =  0.0 0.0 0.0 ... 0.0 0.0
    nstdhdl                  = 1000
    calc-lambda-neighbors    = -1

Thus, the mdp for `glh2glu.vdw_qon.itp` will cover the step 2 and 3. ::

    free-energy              = yes
    separate-dhdl-file       = yes
    sc-alpha                 = 0.5
    sc-power                 = 1
    sc-sigma		         = 0.3
    init-lambda-state = 0
    coul-lambdas = 0.0 0.0  0.0 ... 0.0  0.0 0.1 0.2 ... 0.9 1.0
    fep-lambdas =  0.0 0.05 0.1 ... 0.95 1.0 1.0 1.0 ... 1.0 1.0
    nstdhdl                  = 1000
    calc-lambda-neighbors    = -1

On the other hand, if the crystal structure of the state B is known, the step
2 can be bundled with step 1 to achieve more sampling in state A.

Thus, the mdp for `glh2glu.qoff_vdw.itp` will cover the step 1 and 2. ::

    free-energy              = yes
    separate-dhdl-file       = yes
    sc-alpha                 = 0.5
    sc-power                 = 1
    sc-sigma		         = 0.3
    init-lambda-state = 0
    coul-lambdas = 0.0 0.1 0.2 ... 0.9 1.0 1.0  1.0 1.0  ... 1.0  1.0
    fep-lambdas =  0.0 0.0 0.0 ... 0.0 0.0 0.05 0.1 0.15 ... 0.95 1.0
    nstdhdl                  = 1000
    calc-lambda-neighbors    = -1

Thus, the mdp for `glh2glu.vdw_qon.itp` will cover the step 3. ::

    free-energy              = yes
    separate-dhdl-file       = yes
    sc-alpha                 = 0.5
    sc-power                 = 1
    sc-sigma		         = 0.3
    init-lambda-state = 0
    coul-lambdas = 0.0 0.1 0.2 ... 0.9 1.0
    fep-lambdas =  1.0 1.0 1.0 ... 1.0 1.0
    nstdhdl                  = 1000
    calc-lambda-neighbors    = -1
