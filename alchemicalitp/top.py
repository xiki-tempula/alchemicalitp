import copy
import os
from .field import *
from .ff_io import itpWriter, Reader
from .state_A_B import Alchemistry


class Topology():
    def __init__(self, content_dict = None, content_list = None, name = None, filename=None):
        if filename:
            self.read(filename)
        else:
            if content_dict is None:
                content_dict = {}
            self.content_dict = content_dict
            self.name = name

    def read(self, filename):
        Reader(self, filename)

    def __deepcopy__(self, memo):
        new = Topology()
        new.name = self.name
        new.content_dict = copy.deepcopy(self.content_dict)
        return new

    def __repr__(self): # pragma: no cover
        return self.name

    def add_coul0(self):
        new_top = copy.deepcopy(self)
        for atom in new_top.content_dict['atoms']:
            atom.add_coul0()
        return new_top

    def assign_stateB_charge(self, charge_list):
        new_top = copy.deepcopy(self)
        new_top.content_dict['atoms'].merge_comment()
        for atom, charge in zip(new_top.content_dict['atoms'], charge_list):
            atom.copy2stateB()
            atom.chargeB = charge

        return new_top

    def to_stateB(self):
        new_top = copy.deepcopy(self)
        new_top.content_dict['atoms'].to_stateB()
        new_top.content_dict['bonds'].to_stateB()
        new_top.content_dict['pairs'].to_stateB()
        new_top.content_dict['angles'].to_stateB()
        new_top.content_dict['dihedrals'].to_stateB()
        return new_top

    def write(self, filename=None, format = None):
        if format is None:
            format = os.path.splitext(filename)[-1][1:]
        if format == 'itp':
            if filename is None:
                filename = self.name + '.itp'
            writer = itpWriter(self)
            with open(filename, 'w') as f:
                writer.write_itp(f)
        elif format == 'top':
            if filename is None:
                filename = self.name + '.top'
            writer = itpWriter(self)
            with open(filename, 'w') as f:
                writer.write_top(f)
        else:
            raise NameError('File extension {} not recognised'.format(format))

    # def return_intermediate(self, coul, vdw, bonded, restraint, mass):
    #     new_top = copy.deepcopy(self)
    #     for atom in new_top.content_dict['atoms']:
    #          atom.intermediate_coul(coul)
    #
    #     return new_top

    def remove_zero(self):
        self.content_dict['dihedrals'] = self.content_dict['dihedrals'].remove_zero()

    def split_coul(self, charge_conservation=None):
        # From state A to intermediate
        new_top_1 = copy.deepcopy(self)
        new_top_1.content_dict['atoms'].stateA2intermediate(charge_conservation)
        new_top_1.name = new_top_1.name + '2intermediate'
        # From state intermediate to B
        new_top_2 = copy.deepcopy(self)
        new_top_2.content_dict['atoms'].intermediate2stateB(charge_conservation)
        new_top_2.name = 'intermediate2' + new_top_2.name
        return new_top_1, new_top_2

    def add_stateB(self, top_B, top_A_list, top_B_list):
        alchemistry = Alchemistry(self, top_B, top_A_list, top_B_list)
        top_A_id_map, top_B_id_map = alchemistry.top_A_id_map, alchemistry.top_B_id_map
        return Topology(content_dict=alchemistry.content_dict, name=alchemistry.name), (top_A_id_map, top_B_id_map)


    # def convert_ff(self):
    #     for atomtype in self.content_dict['atomtypes']:
    #         atomtype._comment = self.name
    #
    #     bondtypes = Field('bondtypes')
    #     for bond in self.content_dict['bonds']:
    #         if not isinstance(bond, Comment):
    #             new_bond = copy.copy(bond)
    #             new_bond._comment = self.name + bond.to_str()
    #             bond._b0 = ''
    #             bond._kb = ''
    #             new_bond._i = self.content_dict['atoms'].atom_idx2attr(bond._i, '_type')
    #             new_bond._j = self.content_dict['atoms'].atom_idx2attr(bond._j, '_type')
    #             bondtypes.append(new_bond)
    #     bondtypes.uniqle()
    #     self.content_dict['bondtypes'] = bondtypes
    #
    #     angletypes = Field('angletypes')
    #     for angle in self.content_dict['angles']:
    #         if not isinstance(angle, Comment):
    #             new_angle = copy.copy(angle)
    #             new_angle._comment = self.name + angle.to_str()
    #             angle._th0 = ''
    #             angle._cth = ''
    #             new_angle._i = self.content_dict['atoms'].atom_idx2attr(angle._i, '_type')
    #             new_angle._j = self.content_dict['atoms'].atom_idx2attr(angle._j, '_type')
    #             new_angle._k = self.content_dict['atoms'].atom_idx2attr(angle._k, '_type')
    #             angletypes.append(new_angle)
    #     angletypes.uniqle()
    #     self.content_dict['angletypes'] = angletypes
    #
    #     dihedraltypes = Field('dihedraltypes')
    #     for dihedral in self.content_dict['dihedrals']:
    #         if not isinstance(dihedral, Comment):
    #             new_dihedral = copy.copy(dihedral)
    #             new_dihedral._comment = self.name + dihedral.to_str()
    #             dihedral._phase = ''
    #             dihedral._kd = ''
    #             dihedral._pn = ''
    #             if dihedral._func == '1':
    #                 dihedral._func = '9'
    #                 new_dihedral._func = '9'
    #             new_dihedral._i = self.content_dict['atoms'].atom_idx2attr(dihedral._i, '_type')
    #             new_dihedral._j = self.content_dict['atoms'].atom_idx2attr(dihedral._j, '_type')
    #             new_dihedral._k = self.content_dict['atoms'].atom_idx2attr(dihedral._k, '_type')
    #             new_dihedral._l = self.content_dict['atoms'].atom_idx2attr(dihedral._l, '_type')
    #
    #             dihedraltypes.append(new_dihedral)
    #     dihedraltypes.sort_dihedral()
    #     dihedraltypes.uniqle()
    #     self.content_dict['dihedrals'].uniqle()
    #     self.content_dict['dihedraltypes'] = dihedraltypes
    # #
    # def convert_ff_rtp(self, start_residue_suffix, end_residue_suffix):
    #     residue_list = []
    #     previous_residue = None
    #     atoms = Field('atoms')
    #     for line in self.content_list:
    #         if isinstance(line, Comment):
    #             if len(line._content.split()) == 7:
    #                 try:
    #                     _, resid, resname1, _, resname2, _, charge = line._content.split()
    #                     assert resname1 == resname2
    #                     resname = resname1
    #
    #                     if (not (previous_residue is None)) and (previous_residue != resname):
    #                         top = Topology()
    #                         bonds = Field('bonds')
    #                         impropers = Field('impropers')
    #                         top.content_dict = {'atoms': atoms, 'bonds': bonds, 'impropers': impropers}
    #                         top.name = previous_residue
    #                         residue_list.append(top)
    #                         previous_residue = resname
    #                         atoms = Field('atoms')
    #                     elif previous_residue is None:
    #                         previous_residue = resname
    #
    #                 except:
    #                     pass
    #
    #         elif isinstance(line, Atom):
    #             assert line._residue == previous_residue
    #             new_atom = copy.copy(line)
    #             new_atom._type, new_atom._atom = line._atom, line._type
    #             new_atom._resnr = ''
    #             new_atom._residue = ''
    #             new_atom._mass = ''
    #             new_atom._cgnr = ''
    #             atoms.append(new_atom)
    #     else:
    #         top = Topology()
    #         bonds = Field('bonds')
    #         impropers = Field('impropers')
    #         top.content_dict = {'atoms': atoms, 'bonds': bonds, 'impropers': impropers}
    #         top.name = previous_residue
    #         residue_list.append(top)
    #
    #     residue_start_index = []
    #     residue_end_index = []
    #     for residue in residue_list:
    #         residue_start_index.append(residue.content_dict['atoms'][0]._atom)
    #         residue_end_index.append(residue.content_dict['atoms'][-1]._atom)
    #
    #     for line in self.content_list:
    #         if isinstance(line, Bond):
    #             new_bond = copy.copy(line)
    #             new_bond._comment = self.name + line.to_str()
    #             new_bond._b0 = ''
    #             new_bond._kb = ''
    #             new_bond._i = self.content_dict['atoms'].atom_idx2attr(line._i, '_atom')
    #             new_bond._j = self.content_dict['atoms'].atom_idx2attr(line._j, '_atom')
    #             new_bond._func = ''
    #
    #             for residue_id, residue in enumerate(residue_list):
    #                 start_id = residue.content_dict['atoms'][0]._nr
    #                 end_id = residue.content_dict['atoms'][-1]._nr
    #                 if any([(getattr(line, idx) >= start_id and getattr(line, idx) <= end_id) for idx in ['_i', '_j']]):
    #                     residue_bond = copy.copy(new_bond)
    #                     # check intra molecule
    #                     if any([(getattr(line, idx) < start_id or getattr(line, idx) > end_id) for idx in ['_i', '_j']]):
    #                         if residue_id == 1:
    #                             # only add if it is the middle residue
    #                             for idx in ['_i', '_j']:
    #                                 if getattr(line, idx) < start_id:
    #                                     setattr(residue_bond, idx, '-' + getattr(new_bond, idx))
    #                                 elif getattr(line, idx) > end_id:
    #                                     setattr(residue_bond, idx, '+' + getattr(new_bond, idx))
    #                             residue.content_dict['bonds'].append(residue_bond)
    #                             break
    #                     else:
    #                         residue.content_dict['bonds'].append(residue_bond)
    #                         break
    #
    #         elif isinstance(line, Dihedral):
    #             if line._func == '4':
    #                 new_dihedral = copy.copy(line)
    #                 new_dihedral._comment = self.name + line.to_str()
    #                 new_dihedral._i = self.content_dict['atoms'].atom_idx2attr(line._i, '_atom')
    #                 new_dihedral._j = self.content_dict['atoms'].atom_idx2attr(line._j, '_atom')
    #                 new_dihedral._k = self.content_dict['atoms'].atom_idx2attr(line._k, '_atom')
    #                 new_dihedral._l = self.content_dict['atoms'].atom_idx2attr(line._l, '_atom')
    #                 new_dihedral._func = ''
    #                 new_dihedral._phase = ''
    #                 new_dihedral._kd = ''
    #                 new_dihedral._pn = ''
    #                 for residue_id, residue in enumerate(residue_list):
    #                     start_id = residue.content_dict['atoms'][0]._nr
    #                     end_id = residue.content_dict['atoms'][-1]._nr
    #                     # Try to only do the intra residue bond for the middle one
    #                     if any([(getattr(line, idx) >= start_id and getattr(line, idx) <= end_id) for idx in ['_i', '_j', '_k', '_l']]):
    #                         residue_dihedral = copy.copy(new_dihedral)
    #                         # check intra dihedral
    #                         if any([(getattr(line, idx) < start_id or getattr(line, idx) > end_id) for idx in
    #                                 ['_i', '_j', '_k', '_l']]):
    #                             if residue_id == 1:
    #                                 for idx in ['_i', '_j', '_k', '_l']:
    #                                     if getattr(line, idx) < start_id:
    #                                         setattr(residue_dihedral, idx, '-' + getattr(new_dihedral, idx))
    #                                     elif getattr(line, idx) > end_id:
    #                                         setattr(residue_dihedral, idx, '+' + getattr(new_dihedral, idx))
    #                                 residue.content_dict['impropers'].append(residue_dihedral)
    #                                 break
    #                         else:
    #                             residue.content_dict['impropers'].append(residue_dihedral)
    #                             break
    #
    #     for index, residue in enumerate(residue_list):
    #         residue.content_dict['atoms'].atom_reindex()
    #         residue.content_dict['bonds'].uniqle()
    #         residue.content_dict['impropers'].uniqle()
    #     return residue_list

