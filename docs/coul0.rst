Generation of topology file for ABFE calculations
=================================================

Problem with the ordinary ABFE calculations
-------------------------------------------
The usual ABFE scheme of decoupling the coulombic and van der waals interactions
can be problematic for flexible ligand with a strong dipole. After decoupling
the coulombic interactions between the ligand and the protein, the strong dipole
in the ligand could will force the head of the ligand to interact with the
opposite charged tail of the ligand, which can give rise to the LINC error.

The restraint and charge annihilation
-------------------------------------
To solve this problem, we change the decoupling of the coulombic interactions
to the annihilation of the coulombic interactions by changing the partial
charge of the atom to 0.  ::

    [ moleculetype ]
    ;name            nrexcl
     ZA               3

    [ atoms ]
    ;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
         1   N3     1    ZA     N    1    -0.295300     14.01000    N3     0.     14.01000 ; qtot -0.295
         2    H     1    ZA    H1    2     0.258800      1.00800    H      0.      1.00800 ; qtot -0.037
         3    H     1    ZA    H2    3     0.258800      1.00800    H      0.      1.00800 ; qtot 0.222
         4    H     1    ZA    H3    4     0.258800      1.00800    H      0.      1.00800 ; qtot 0.481
         5   CT     1    ZA    CA    5     0.114500     12.01000    CT     0.     12.01000 ; qtot 0.596
         6   H1     1    ZA    HA    6     0.056000      1.00800    H1     0.      1.00800 ; qtot 0.652
         7   CT     1    ZA    CB    7    -0.174200     12.01000    CT     0.     12.01000 ; qtot 0.477
         8   HC     1    ZA   HB1    8     0.062900      1.00800    HC     0.      1.00800 ; qtot 0.540
         9   HC     1    ZA   HB2    9     0.062900      1.00800    HC     0.      1.00800 ; qtot 0.603
        10   HC     1    ZA   HB3   10     0.062900      1.00800    HC     0.      1.00800 ; qtot 0.666
        11    C     1    ZA     C   11     0.714501     12.01000     C     0.     12.01000 ; qtot 1.381
        12   O2     1    ZA     O   12    -0.690301     16.00000    O2     0.     16.00000 ; qtot 0.690
        13   O2     1    ZA   OXT   13    -0.690301     16.00000    O2     0.     16.00000 ; qtot 0.000

The topology file could be generated by ::

    >>> from alchemicalitp.tutorial_files import glh_top
    >>> top = alchemicalitp.top.Topology(filename=glh_top)
    >>> new_top = top.add_coul0()
    >>> new_top.write('ligand20.itp')

The restraint and charge annihilation could then be executed with these mdp
options. ::

    ;----------------------------------------------------
    ; FREE ENERGY CALCULATIONS
    ;----------------------------------------------------
    free-energy              = yes
    separate-dhdl-file       = yes
    init-lambda-state = 0
    bonded-lambdas = 0.0 0.01 0.025 0.05 0.075 0.1 0.2 0.35 0.5 0.75 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    coul-lambdas   = 0.0 0.0  0.0   0.0  0.0   0.0 0.0 0.0  0.0 0.0  0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
    nstdhdl                  = 1000
    calc-lambda-neighbors    = -1

The decoupling of van der waals
-------------------------------
The next step is the decoupling of the molecule as the annihilation of van der
waals interactions are difficult to converge. Given that this happens after
the charge annihilation, we would need a new topology file where the partial
charge is annihilated. ::

    [ moleculetype ]
    ;name            nrexcl
     ZA               3

    [ atoms ]
    ;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
         1   N3     1    ZA     N    1     0.     14.01000 ; qtot -0.295
         2    H     1    ZA    H1    2     0.      1.00800 ; qtot -0.037
         3    H     1    ZA    H2    3     0.      1.00800 ; qtot 0.222
         4    H     1    ZA    H3    4     0.      1.00800 ; qtot 0.481
         5   CT     1    ZA    CA    5     0.     12.01000 ; qtot 0.596
         6   H1     1    ZA    HA    6     0.      1.00800 ; qtot 0.652
         7   CT     1    ZA    CB    7     0.     12.01000 ; qtot 0.477
         8   HC     1    ZA   HB1    8     0.      1.00800 ; qtot 0.540
         9   HC     1    ZA   HB2    9     0.      1.00800 ; qtot 0.603
        10   HC     1    ZA   HB3   10     0.      1.00800 ; qtot 0.666
        11    C     1    ZA     C   11     0.     12.01000 ; qtot 1.381
        12   O2     1    ZA     O   12     0.     16.00000 ; qtot 0.690
        13   O2     1    ZA   OXT   13     0.     16.00000 ; qtot 0.000

The topology file can be generated with ::

    >>> from alchemicalitp.tutorial_files import glh_top
    >>> top = alchemicalitp.top.Topology(filename=glh_top)
    >>> new_top = top.add_coul0()
    >>> stateB = new_top.to_stateB()
    >>> new_top.write('ligand_0.itp')

The corresponding mdp file for the vdw decoupling is ::

    ;----------------------------------------------------
    ; FREE ENERGY CALCULATIONS
    ;----------------------------------------------------
    free-energy              = yes
    couple-moltype           = ZA
    couple-lambda0           = vdw
    couple-lambda1           = none
    separate-dhdl-file       = yes
    sc-alpha                 = 0.5
    sc-power                 = 1
    sc-sigma		             = 0.3
    init-lambda-state = 0
    bonded-lambdas = 1.0 1.0  1.0 1.0 1.0 1.0 1.0 1.0 1.0  1.0 1.0  1.0 1.0  1.0 1.0  1.0
    vdw-lambdas =    0.0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0
    nstdhdl                  = 1000
    calc-lambda-neighbors    = -1
