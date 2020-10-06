import copy
import warnings
from .field import Moleculetype, Atoms, Bonds, Pairs, Angles, Dihedrals, Cmaps
from .entry import Comment, Dummy_Atomtype

class Alchemistry():
    def __init__(self, top_A, top_B, top_A_list, top_B_list):
        self.content_dict = {}
        self.name = top_A.name + '2' + top_B.name
        self.top_A = top_A
        self.top_A.remove_zero()
        self.top_B = top_B
        self.top_B.remove_zero()
        self.top_A_list = top_A_list
        self.top_B_list = top_B_list

        self._merge_defaults()
        self._merge_atomtypes()
        self._add_dummy_atomtypes()
        self._merge_cmaptypes()
        self._create_moleculetype()
        self._merge_atoms()
        self.A_list = self._update_top_list(self.top_A_list, self.top_A_id_map)
        self.B_list = self._update_top_list(self.top_B_list, self.top_B_id_map)
        self._merge_bonds()
        self._merge_pairs()
        self._merge_angles()
        self._merge_dihedrals()
        self._merge_cmaps()

    def _merge_defaults(self):
        if 'defaults' in self.top_A.content_dict and 'defaults' in self.top_B.content_dict:
            assert self.top_A.content_dict['defaults'] == self.top_B.content_dict['defaults']
            self.content_dict['defaults'] = self.top_A.content_dict['defaults']

    def _merge_atomtypes(self):
        if 'atomtypes' in self.top_A.content_dict and 'atomtypes' in self.top_B.content_dict:
            self.content_dict['atomtypes'] = copy.copy(self.top_A.content_dict['atomtypes'])
            self.content_dict['atomtypes'].union(self.top_B.content_dict['atomtypes'])

    def _merge_cmaptypes(self):
        if 'cmaptypes' in self.top_A.content_dict or 'cmaptypes' in self.top_B.content_dict:
            self.content_dict['cmaptypes'] = copy.copy(self.top_A.content_dict['cmaptypes'])
            self.content_dict['cmaptypes'].union(self.top_B.content_dict['cmaptypes'])

    def _add_dummy_atomtypes(self):
        if 'atomtypes' in self.content_dict:
            self.content_dict['atomtypes'].append(Dummy_Atomtype())
        else:
            warnings.warn("No atomtypes directive found in the input file.\n" 
                          "Dummy atomtype added: \n{}".format(Dummy_Atomtype().to_str()), Warning)

    def _create_moleculetype(self):
        assert self.top_A.content_dict['moleculetype'].nrexcl == self.top_B.content_dict['moleculetype'].nrexcl
        name = self.top_A.content_dict['moleculetype'].name + '2' + self.top_B.content_dict['moleculetype'].name
        self.content_dict['moleculetype'] = Moleculetype()
        self.content_dict['moleculetype'].name = name
        self.content_dict['moleculetype'].nrexcl = self.top_A.content_dict['moleculetype'].nrexcl
        self.content_dict['moleculetype'].comments = ['{} comments'.format(self.top_A.name), ]
        self.content_dict['moleculetype'].comments.extend(self.top_A.content_dict['moleculetype'].comments)
        self.content_dict['moleculetype'].comments.append('{} comments'.format(self.top_B.name))
        self.content_dict['moleculetype'].comments.extend(self.top_B.content_dict['moleculetype'].comments)

    def _merge_atoms(self):
        self.top_A_id_map = {}
        self.top_B_id_map = {}

        # The index to be skipped
        skip_state_B = []

        field = Atoms()
        self.content_dict['atoms'] = field
        # track by atom id and index
        current_atom = 0
        top_A_idx = 0
        top_B_idx = 0

        # Ensure that the index can be mapped to atoms instead of comments
        self.top_A.content_dict['atoms'].merge_comment()
        self.top_B.content_dict['atoms'].merge_comment()
        top_A_length = len(self.top_A.content_dict['atoms'])
        top_B_length = len(self.top_B.content_dict['atoms'])

        while top_A_idx < top_A_length or top_B_idx < top_B_length:
            atom_A = self.top_A.content_dict['atoms'][top_A_idx]
            atom_B = self.top_B.content_dict['atoms'][top_B_idx]
            if (atom_A and atom_A.nr in self.top_A_list) or (atom_B and atom_B.nr in self.top_B_list):
                if atom_A and atom_A.nr in self.top_A_list:
                    atom_B_nr = self.top_B_list[self.top_A_list.index(top_A_idx+1)]
                    current_atom += 1
                    top_A_idx += 1
                    new_atom = copy.copy(atom_A)
                    new_atom.cgnr = new_atom.nr = current_atom
                    if atom_B_nr is None:
                        new_atom.typeB = 'DUM'
                        new_atom.massB = atom_A.mass
                        new_atom.chargeB = 0.0
                        field.append(new_atom)
                        self.top_A_id_map[atom_A.nr] = current_atom
                    else:
                        # Only increment the number when the atom B is not skipped previously
                        if atom_B_nr in skip_state_B:
                            pass
                        else:
                            top_B_idx += 1
                        new_atom.typeB = self.top_B.content_dict['atoms'].atom_idx2attr(atom_B_nr, 'type')
                        new_atom.massB = self.top_B.content_dict['atoms'].atom_idx2attr(atom_B_nr, 'mass')
                        new_atom.chargeB = self.top_B.content_dict['atoms'].atom_idx2attr(atom_B_nr, 'charge')
                        new_atom.comment = ' to '.join([atom_A.comment, atom_B.comment])
                        field.append(new_atom)
                        self.top_A_id_map[atom_A.nr] = current_atom
                        self.top_B_id_map[atom_B_nr] = current_atom
                elif atom_B and atom_B.nr in self.top_B_list:
                    if self.top_A_list[self.top_B_list.index(top_B_idx+1)] is None:
                        current_atom += 1
                        top_B_idx += 1
                        new_atom = copy.copy(atom_B)
                        new_atom.cgnr = new_atom.nr = current_atom
                        new_atom.typeB = new_atom.type
                        new_atom.massB = new_atom.mass
                        new_atom.chargeB = new_atom.charge
                        new_atom.type = 'DUM'
                        new_atom.charge = 0.0
                        field.append(new_atom)
                        self.top_B_id_map[atom_B.nr] = current_atom
                    else:
                        top_B_idx += 1
                        skip_state_B.append(top_B_idx)
            else:
                current_atom += 1
                new_atom = copy.copy(atom_A)
                new_atom.cgnr = new_atom.nr = current_atom
                new_atom.comment = ' to '.join([atom_A.comment, atom_B.comment])
                if atom_A.atom == atom_B.atom and atom_A.charge == atom_B.charge and atom_A.mass == atom_B.mass and \
                        atom_A.residue == atom_B.residue and atom_A.resnr == atom_B.resnr and atom_A.type == atom_B.type:
                    field.append(new_atom)
                    top_A_idx += 1
                    top_B_idx += 1
                    self.top_A_id_map[atom_A.nr] = current_atom
                    self.top_B_id_map[atom_B.nr] = current_atom

                elif atom_A.atom == atom_B.atom and atom_A.resnr == atom_B.resnr:

                    new_atom.typeB = atom_B.type
                    new_atom.massB = atom_B.mass
                    new_atom.chargeB = atom_B.charge
                    field.append(new_atom)
                    top_A_idx += 1
                    top_B_idx += 1
                    self.top_A_id_map[atom_A.nr] = current_atom
                    self.top_B_id_map[atom_B.nr] = current_atom
                else:
                    raise NameError('No know how to link {} to {}.'.format(atom_A, atom_B))

    def _update_top_list(self, old_list, mapping):
        new_list = []
        for index in old_list:
            if index:
                new_list.append(mapping[index])
            else:
                new_list.append(None)
        return new_list

    def _merge_entry(self, entry_A, entry_B):
        if entry_A is None:
            new_entry = copy.deepcopy(entry_B)
            new_entry.comment += 'Entry from topology 2'
        elif entry_B is None:
            new_entry = copy.deepcopy(entry_A)
            new_entry.comment += 'Entry from topology 1'
        else:
            new_entry = entry_A.add_stateB(entry_B)
        return new_entry
    def _merge_entries(self, top_A, top_B, field):
        top_A.update_idx(self.top_A_id_map)
        entries_A = top_A
        top_B.update_idx(self.top_B_id_map)
        entries_B = top_B

        entries_A.sort()
        entries_B.sort()

        top_A_idx = 0
        top_B_idx = 0
        top_A_length = len(entries_A)
        top_B_length = len(entries_B)

        while top_A_idx < top_A_length or top_B_idx < top_B_length:
            entry_A = entries_A[top_A_idx]
            entry_B = entries_B[top_B_idx]
            new_entry = self._merge_entry(entry_A, entry_B)
            if new_entry:
                top_A_idx += 1
                top_B_idx += 1
                field.append(new_entry)
            elif entry_A.idx_less(entry_B):
                top_A_idx += 1
                if entry_A.idx_in(self.A_list):
                    if entry_A.check_other_list(self.A_list, self.B_list):
                        entry_A.comment = 'Entry from topology 1'
                        field.append(entry_A)
                    else:
                        new_entry = self.create_dummy(entry_A, None)
                        field.append(new_entry)
                else:
                    new_entry = self.create_dummy(entry_A, None)
                    field.append(new_entry)
            else:
                top_B_idx += 1
                if entry_B.idx_in(self.B_list):
                    if entry_B.check_other_list(self.B_list, self.A_list):
                        entry_B.comment = 'Entry from topology 2'
                        field.append(entry_B)
                    else:
                        new_entry = self.create_dummy(None, entry_B)
                        field.append(new_entry)
                else:
                    new_entry = self.create_dummy(None, entry_B)
                    field.append(new_entry)

    def create_dummy(self, entry_A, entry_B):
        if entry_A is None:
            entry_A = entry_B.create_dummy()
        elif entry_B is None:
            entry_B = entry_A.create_dummy()
        dummy = entry_A.add_stateB(entry_B)
        if dummy:
            return dummy
        else:
            raise IndexError('Dummy {}/{} can not be created.'.format(entry_A, entry_B))


    def _merge_bonds(self):
        top_A_bond = self.top_A.content_dict['bonds']
        top_B_bond = self.top_B.content_dict['bonds']
        new_bond = Bonds()
        self._merge_entries(top_A_bond, top_B_bond, new_bond)
        self.content_dict['bonds'] = new_bond

    def _merge_pairs(self):
        top_A_pair = self.top_A.content_dict['pairs']
        top_B_pair = self.top_B.content_dict['pairs']
        new_pair = Pairs()
        self._merge_entries(top_A_pair, top_B_pair, new_pair)
        self.content_dict['pairs'] = new_pair

    def _merge_angles(self):
        top_A_angle = self.top_A.content_dict['angles']
        top_B_angle = self.top_B.content_dict['angles']
        new_angle = Angles()
        self._merge_entries(top_A_angle, top_B_angle, new_angle)
        self.content_dict['angles'] = new_angle

    def _merge_dihedrals(self):
        top_A_dihedral = self.top_A.content_dict['dihedrals']
        top_B_dihedral = self.top_B.content_dict['dihedrals']
        new_dihedral = Dihedrals()
        self._merge_entries(top_A_dihedral, top_B_dihedral, new_dihedral)
        self.content_dict['dihedrals'] = new_dihedral

    def _merge_cmaps(self):
        if 'cmap' in self.top_A.content_dict and 'cmap' in self.top_B.content_dict:
            top_A_cmap = self.top_A.content_dict['cmap']
            top_B_cmap = self.top_B.content_dict['cmap']
            new_cmap = Cmaps()
            self._merge_entries(top_A_cmap, top_B_cmap, new_cmap)
            self.content_dict['cmap'] = new_cmap
        elif (not 'cmap' in self.top_A.content_dict) and (not 'cmap' in self.top_B.content_dict):
            pass
        else:
            assert 'cmap' in self.top_A.content_dict == 'cmap' in self.top_B.content_dict, 'cmap need to present in both topology files'

class MakeLambdas():
    def __init__(self, top, lambdas):
        self._check_lambda(lambdas)
        self.tops = []
        self.top = top

        self._generate_tops()
        self._write_tops()

    def _check_lambda(self, lambdas):
        length = None
        for key, value in lambdas.items():
            if length is None:
                length = len(value)
            else:
                assert length == len(value), 'Lambda should have the same value'
        self.lambdas = lambdas
        self.lambdas_idx = range(length)

    def _idx2lambda(self, index):
        coul, vdw, bonded, restraint, mass = [None,] * 5
        if 'coul-lambdas' in self.lambdas:
            coul = self.lambdas['coul-lambdas'][index]
        if 'vdw-lambdas' in self.lambdas:
            vdw = self.lambdas['vdw-lambdas'][index]
        if 'bonded-lambdas' in self.lambdas:
            bonded = self.lambdas['bonded-lambdas'][index]
        if 'restraint-lambdas' in self.lambdas:
            restraint = self.lambdas['restraint-lambdas'][index]
        if 'mass-lambdas' in self.lambdas:
            mass = self.lambdas['mass-lambdas'][index]

        if 'fep-lambdas' in self.lambdas:
            if coul is None:
                coul = self.lambdas['fep-lambdas'][index]
            if vdw is None:
                vdw = self.lambdas['fep-lambdas'][index]
            if bonded is None:
                bonded = self.lambdas['fep-lambdas'][index]
            if restraint is None:
                restraint = self.lambdas['fep-lambdas'][index]
            if mass is None:
                mass = self.lambdas['fep-lambdas'][index]

        if coul is None:
            coul = 0
        if vdw is None:
            vdw = 0
        if bonded is None:
            bonded = 0
        if restraint is None:
            restraint = 0
        if mass is None:
            mass = 0

        return coul, vdw, bonded, restraint, mass

    def _generate_tops(self):
        for i in self.lambdas_idx:
            coul, vdw, bonded, restraint, mass = self._idx2lambda(i)
            new_top = self.top.return_intermediate(coul, vdw, bonded, restraint, mass)
            new_top.name = str(i) + '.' + self.top.name
            self.tops.append(new_top)

    def _write_tops(self):
        for i in self.lambdas_idx:
            top = self.tops[i]
            top.write()

