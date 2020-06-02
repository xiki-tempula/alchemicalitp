# from state_A_B import Alchemistry, MakeLambdas
from top import Topology
import copy

def alchemical(top_A, top_B, top_A_list, top_B_list, itp, top=None):
    top_A = Topology(filename=top_A)
    top_B = Topology(filename=top_B)
    new_top = top_A.add_stateB(top_B, top_A_list, top_B_list)
    new_top.content_dict['moleculetype'].name = 'Protein'
    top_A, top_B = new_top.split_coul()
    top_A.write()
    top_B.write()
    # top_writer = itpWriter(top_A)
    # with open(itp, 'w') as f:
    #     top_writer.write_itp(f)
    # if top:
    #     with open(top, 'w') as f:
    #         top_writer.write_top(f)
#
# def to_coul0(filename):
#     top = Topology()
#     top.read(filename)
#     top = top.add_coul0()
#     MakeLambdas(top, lambdas={'coul-lambdas':[0,0.25,0.5,0.75,1.0]})

if __name__ == '__main__':
    alchemical('/Volumes/GoogleDrive/My Drive/Simulations/KDEL/ABFE/COO_CONH/KDEL_CONH.itp',
               '/Volumes/GoogleDrive/My Drive/Simulations/KDEL/ABFE/COO_CONH/KDEL_COO.itp',
               [None, ], [20, ], 'protein.itp', top='top.top')
    # alchemical('../example/GLH.top', '../example/GLU.top', [20, ], [None, ], 'protein.itp', top='top.top')
    # to_coul0('/Users/xiki_tempula/Downloads/ligand.itp')


