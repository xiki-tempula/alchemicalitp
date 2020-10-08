from pkg_resources import resource_filename

glu_top = resource_filename(__name__, 'data/fep_example/GLU.top')
glh_top = resource_filename(__name__, 'data/fep_example/GLH.top')
glu_crd = resource_filename(__name__, 'data/fep_example/GLU.gro')
glh_crd = resource_filename(__name__, 'data/fep_example/GLH.gro')
mdp_em0 = resource_filename(__name__, 'data/fep_example/minim0.mdp')
mdp_em1 = resource_filename(__name__, 'data/fep_example/minim1.mdp')
mdp_energy0 = resource_filename(__name__, 'data/fep_example/test0.mdp')
mdp_energy1 = resource_filename(__name__, 'data/fep_example/test1.mdp')
