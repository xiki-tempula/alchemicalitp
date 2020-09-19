from .entry import Comment, Atomtype, Atom, Bond, Pair, Angle, Dihedral, Cmap
import copy
import numpy as np

class Field():
    def __init__(self, name):
        self.name = name
        self.content = []

    def __eq__(self, other):
        current = [line for line in self.content if not isinstance(line, Comment)]
        other = [line for line in other if not isinstance(line, Comment)]
        return current == other

    def append(self, line):
        self.content.append(line)

    def uniqle(self):
        content = []
        for line in self.content:
            if not line in content:
                content.append(line)
        self.content = content

    def __iter__(self):
        return self.content.__iter__()

    def __getitem__(self, index):
        try:
            return self.content.__getitem__(index)
        except IndexError:
            return None

    def __len__(self):
        return len(self.content)

    def to_str(self):
        output = []
        output.append('[ {} ]'.format(self.name))
        for line in self.content:
            output.append(line.to_str())
        return '\n'.join(output)

    def union(self, other):
        for line in other:
            if not line in self.content:
                self.content.append(line)

    def __deepcopy__(self, memo):
        newone = type(self)()
        newone.name = self.name
        newone.content = copy.deepcopy(self.content)
        return newone

    def atom_reindex(self):
        '''Alter the atom index to make sure that it can be used in rtp file'''
        for index, line in enumerate(self.content):
            if not isinstance(line, Comment):
                line.mass = index
                line.nr = ''

    def add_comment(self, line):
        self.content.append(Comment(line[1:]))

    def add_entry(self, line):
        if ';' in line:
            comment = line.split(';')[-1]
            line = line.split(';')[0]
            return line, comment
        else:
            comment = ''
            return line, comment

    def update_idx(self, mapping):
        for entry in self.content:
            entry.update_idx(mapping)


    def merge_comment(self):
        comment = None
        new_content = []
        for entry in self.content:
            if comment:
                entry.comment = entry.comment + comment.comment
                comment = None
            if isinstance(entry, Comment):
                comment = entry
            else:
                new_content.append(entry)
        self.content = new_content

    def to_stateB(self):
        for entry in self.content:
            entry.to_stateB()

class Default(Field):
    def __init__(self,):
        self.comments = []

    def edit(self, line):
        nbfunc, comb_rule, gen_pairs, fudgeLJ, fudgeQQ = line.split()
        self.nbfunc = nbfunc
        self.comb_rule = comb_rule
        self.gen_pairs = gen_pairs
        self.fudgeLJ = fudgeLJ
        self.fudgeQQ = fudgeQQ

    def __deepcopy__(self, memo):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def to_str(self):
        output = []
        output.append('[ defaults ]')
        output.append(';' + '\n;'.join(self.comments))
        output.append('{: <15} {: <15} {: <15} {: <12} {: <12}'.format(self.nbfunc, self.comb_rule, self.gen_pairs,
                                                                        self.fudgeLJ, self.fudgeQQ))
        return '\n'.join(output)

    def __eq__(self, other):
        if (self.nbfunc == other.nbfunc) and (self.comb_rule == other.comb_rule) and \
                (self.gen_pairs == other.gen_pairs) and (self.fudgeLJ == other.fudgeLJ) and \
                (self.fudgeQQ == other.fudgeQQ):
            return True
        else:
            return False

    def __repr__(self):
        return str((self.nbfunc, self.comb_rule, self.gen_pairs, self.fudgeLJ, self.fudgeQQ))

    def add_comment(self, line):
        self.comments.append(line[1:])

class Moleculetype(Field):
    def __init__(self):
        self.comments = []

    def __deepcopy__(self, memo):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def edit(self, line):
        name, nrexcl = line.split()
        self.name = name
        self.nrexcl = nrexcl

    def add_comment(self, line):
        self.comments.append(line[1:])

    def to_str(self):
        output = []
        output.append('[ moleculetype ]')
        output.append(';' + '\n;'.join(self.comments))
        output.append(self.name + ' ' + self.nrexcl)
        return '\n'.join(output)

    def __repr__(self):
        return str((self.name, self.nrexcl))

class Atomtypes(Field):
    def __init__(self):
        super().__init__('Atomtypes')

    def add_entry(self, line):
        line, comment = super().add_entry(line)
        name, at_num, mass, charge, ptype, sigma, epsilon = line.strip().split()
        self.content.append(Atomtype(name, at_num, mass, charge, ptype, sigma, epsilon, comment=comment))

class Cmaptypes(Field):
    def __init__(self):
        super().__init__('Cmaptypes')

    def add_entry(self, line):
        line, comment = super().add_entry(line)
        elements = line.strip().split()
        names = elements[:5]
        func = elements[5]
        if len(elements) > 6:
            rows, cols = elements[6], elements[7]
            data = elements[8:]
            assert len(data) == int(rows) * int(cols)
            self.content.append(Cmap(names, func, rows, cols, data, comment=comment))
        else:
            self.content.append(Cmap(names, func, comment=comment))

class Cmaps(Cmaptypes):
    def __init__(self):
        super(Cmaptypes, self).__init__('Cmap')

    def sort(self):
        self.merge_comment()
        self.content.sort(key = lambda x: [int(name) for name in x.names])

class Atoms(Field):
    def __init__(self):
        super().__init__('Atoms')

    def add_entry(self, line):
        line, comment = super().add_entry(line)
        if len(line.strip().split()) == 8:
            nr, type, resnr, residue, atom, cgnr, charge, mass = line.strip().split()
            self.content.append(Atom(nr, type, resnr, residue, atom, cgnr, charge, mass, comment=comment))
        else:
            nr, type, resnr, residue, atom, cgnr, charge, mass, typeB, chargeB, massB = line.strip().split()
            self.content.append(Atom(nr, type, resnr, residue, atom, cgnr, charge, mass, typeB, chargeB, massB, comment=comment))

    def atom_idx2attr(self, idx, attr):
        for line in self.content:
            if not isinstance(line, Comment):
                if line.nr == idx:
                    return getattr(line, attr)

    def _find_nearest_int_charge(self, state_A, state_B,
                                 state_inter_A, state_inter_B):
        '''Find the nearest interger charge for the middle state'''
        assert np.isclose(np.round(state_A) - state_A,
                          np.round(state_B) - state_B,
                          atol=0.001)
        offset = np.round(state_A) - state_A
        aim_charge = np.round((state_inter_A + state_inter_B) / 2
                              + offset)
        return aim_charge - offset

    def stateA2intermediate(self, charge_conservation):
        atoms = []
        for atom in self.content:
            if atom.typeB:
                atoms.append(atom)
        state_A = sum([float(atom.charge) for atom in atoms])
        state_B = sum([float(atom.chargeB) for atom in atoms])
        state_inter_A = sum(
            [float(atom.charge) for atom in atoms if atom.type[:3] != 'DUM' and atom.typeB[:3] != 'DUM'])
        state_inter_B = sum(
            [float(atom.chargeB) for atom in atoms if atom.type[:3] != 'DUM' and atom.typeB[:3] != 'DUM'])

        abs_change = sum([abs(float(atom.chargeB) - float(atom.charge)) for atom in atoms if atom.type[:3] != 'DUM' and atom.typeB[:3] != 'DUM'])
        for atom in atoms:
            if atom.type[:3] == 'DUM':
                # Makes sure that if it start from a dummy, it ends with a dummy as no vdw is modified
                atom.chargeB = 0
            elif atom.typeB[:3] == 'DUM':
                atom.chargeB = float(atom.chargeB)
            else:
                # make sure the charge is balanced
                if charge_conservation is None:
                    atom.chargeB = (float(atom.charge) + float(
                        atom.chargeB)) / 2
                elif charge_conservation == 'nearest':
                    aim_charge = self._find_nearest_int_charge(state_A, state_B,
                                 state_inter_A, state_inter_B)
                    current_charge = (state_inter_A + state_inter_B) / 2
                    atom.chargeB = (float(atom.charge) + float(atom.chargeB)) / 2 + \
                               abs((float(atom.chargeB) - float(atom.charge)) / 2) / (abs_change/2) * (aim_charge-current_charge)
        sum_charge = sum([atom.chargeB for atom in atoms])
        for atom in atoms:
            atom.str_charge()

    def intermediate2stateB(self, charge_conservation):
        atoms = []
        for atom in self.content:
            if atom.typeB:
                atoms.append(atom)
        state_A = sum([float(atom.charge) for atom in atoms])
        state_B =state_A = sum([float(atom.charge) for atom in atoms])
        state_B = sum([float(atom.chargeB) for atom in atoms])
        state_inter_A = sum(
            [float(atom.charge) for atom in atoms if atom.type[:3] != 'DUM' and atom.typeB[:3] != 'DUM'])
        state_inter_B = sum(
            [float(atom.chargeB) for atom in atoms if atom.type[:3] != 'DUM' and atom.typeB[:3] != 'DUM'])
        abs_change = sum([abs(float(atom.chargeB) - float(atom.charge)) for atom in atoms if atom.type[:3] != 'DUM' and atom.typeB[:3] != 'DUM'])

        for atom in atoms:
            if atom.type[:3] == 'DUM':
                atom.charge = float(atom.charge)
            elif atom.typeB[:3] == 'DUM':
                # Makes sure that if it end with a dummy, it start with a dummy as no vdw is modified
                atom.charge = 0
            else:
                if charge_conservation is None:
                    atom.charge = (float(atom.charge) + float(
                        atom.chargeB)) / 2
                elif charge_conservation == 'nearest':
                    aim_charge = self._find_nearest_int_charge(state_A, state_B,
                                 state_inter_A, state_inter_B)
                    current_charge = (state_inter_A + state_inter_B) / 2
                    # make sure the charge is balanced
                    atom.charge = (float(atom.charge) + float(atom.chargeB)) / 2 + \
                                   abs((float(atom.chargeB) - float(atom.charge)) / 2) / (abs_change/2) * (aim_charge-current_charge)

        for atom in atoms:
            atom.str_charge()

class Pairs(Field):
    def __init__(self):
        super().__init__('Pairs')

    def sort(self):
        self.merge_comment()
        self.content.sort(key = lambda x: (x.i, x.j))

    def add_entry(self, line):
        line, comment = super().add_entry(line)
        i, j, func = line.strip().split()
        self.content.append(Pair(i, j, func, comment=comment))

class Bonds(Pairs):
    def __init__(self):
        super(Pairs, self).__init__('Bonds')

    def add_entry(self, line):
        line, comment = super(Pairs, self).add_entry(line)
        func = line.strip().split()[2]
        if func == '1':
            if len(line.strip().split()) == 3:
                i, j, func, = line.strip().split()
                self.content.append(Bond(i, j, func, b0='', kb='', comment=comment))
            elif len(line.strip().split()) == 5:
                i, j, func, b0, kb = line.strip().split()
                self.content.append(Bond(i, j, func, b0=b0, kb=kb, comment=comment))
            else:
                i, j, func, b0, kb, b0B, kbB = line.strip().split()
                self.content.append(Bond(i, j, func, b0=b0, kb=kb, b0B=b0B, kbB=kbB, comment=comment))
        else:
            raise NotImplementedError('Function type {} not implemented yet'.format(func))

class Angles(Field):
    def __init__(self):
        super().__init__('Angles')

    def sort(self):
        self.merge_comment()
        self.content.sort(key = lambda x: (x.i, x.j, x.k))

    def add_entry(self, line):
        line, comment = super().add_entry(line)
        func = line.strip().split()[3]
        if func == '1':
            if len(line.strip().split()) == 4:
                i, j, k, func, = line.strip().split()
                self.content.append(Angle(i, j, k, func, th0='', cth='', comment=comment))
            elif len(line.strip().split()) == 6:
                i, j, k, func, th0, cth = line.strip().split()
                self.content.append(Angle(i, j, k, func, th0=th0, cth=cth, comment=comment))
            else:
                i, j, k, func, th0, cth, th0B, cthB = line.strip().split()
                self.content.append(Angle(i, j, k, func, th0=th0, cth=cth, th0B=th0B, cthB=cthB, comment=comment))
        else:
            raise NotImplementedError('Function type {} not implemented yet'.format(func))

class Dihedrals(Field):
    def __init__(self):
        super().__init__('Dihedrals')

    def sort(self):
        self.merge_comment()
        # Make sure that it is sorted by func, phase and pn
        self.content.sort(key = lambda x: (x.i, x.j, x.k, x.l, x.func, x.phase, x.pn))

    def remove_zero(self):
        new_content = []
        for dihedral in self.content:
            if dihedral.func in [1, 4, 9]:
                if np.abs(float(dihedral.kd)) > 0.000001:
                    new_content.append(dihedral)
            else:
                new_content.append(dihedral)
        new_field = Dihedrals()
        new_field.content = new_content
        return new_field


    def add_entry(self, line):
        line, comment = super().add_entry(line)
        func = line.strip().split()[4]
        if func in ['1', '4', '9']:
            if len(line.strip().split()) == 5:
                i, j, k, l, func= line.strip().split()
                self.content.append(Dihedral(i, j, k, l, func, phase="", kd='', pn="", comment=comment))
            elif len(line.strip().split()) == 8:
                i, j, k, l, func, phase, kd, pn = line.strip().split()
                self.content.append(Dihedral(i, j, k, l, func, phase=phase, kd=kd, pn=pn, comment=comment))
            else:
                i, j, k, l, func, phase, kd, pn, phaseB, kdB, pnB = line.strip().split()
                self.content.append(Dihedral(i, j, k, l, func, phase=phase, kd=kd, pn=pn, phaseB=phaseB, kdB=kdB, pnB=pnB, comment=comment))
        elif func == '3':
            if len(line.strip().split()) == 11:
                i, j, k, l, func, C0, C1, C2, C3, C4, C5 = line.strip().split()
                self.content.append(Dihedral(i, j, k, l, func, C0=C0, C1=C1, C2=C2, C3=C3, C4=C4, C5=C5, comment=comment))
            else:
                i, j, k, l, func, C0, C1, C2, C3, C4, C5, C0B, C1B, C2B, C3B, C4B, C5B = line.strip().split()
                self.content.append(
                    Dihedral(i, j, k, l, func, C0=C0, C1=C1, C2=C2, C3=C3, C4=C4, C5=C5,
                             C0B=C0B, C1B=C1B, C2B=C2B, C3B=C3B, C4B=C4B, C5B=C5B,
                             comment=comment))

        else:
            raise NotImplementedError('Function type {} not implemented yet'.format(func))






