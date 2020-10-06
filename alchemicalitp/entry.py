"""
entry.py
Handles each individual entry.

"""
import copy
import numpy as np
class EntryBase():
    def __deepcopy__(self, memo):
        return copy.copy(self)
    def add_coul0(self): # pragma: no cover
        pass
    def intermediate_coul(self, lam): # pragma: no cover
        pass

    def to_stateB(self):
        pass

class Comment(EntryBase):
    def __init__(self, comment):
        self.comment = comment
    def __deepcopy__(self, memo):
        return Comment(self.comment)
    def to_str(self):
        return '; ' + self.comment
    def __getattr__(self, item):
        return None
    def __repr__(self): # pragma: no cover
        return self.comment

    def update_idx(self, mapping, sort=True):
        pass

def join_comment(comment_A, comment_B):
    if comment_A and comment_B:
        return comment_A + ' to ' + comment_B
    elif comment_A:
        return comment_A
    elif comment_B:
        return comment_B
    else:
        return ''


class Atomtype(EntryBase):
    def __init__(self, name, at_num, mass, charge, ptype, sigma, epsilon, comment=''):
        self.name = name
        self.at_num = at_num
        self.mass = mass
        self.charge = charge
        self.ptype = ptype
        self.sigma = sigma
        self.epsilon = epsilon
        self.comment = comment

    def to_str(self):
        return '{: <15} {: <3} {: <11} {: <12} {: <6} {: <16} {: <16} ; {}'.format(self.name, self.at_num, self.mass,
                                                                self.charge, self.ptype, self.sigma,
                                                                self.epsilon, self.comment)

    def __eq__(self, other):
        if (self.name == other.name) and (self.at_num == other.at_num) and \
                (self.mass == other.mass) and (self.charge == other.charge) and \
                (self.ptype == other.ptype) and (self.sigma == other.sigma) and \
                (self.epsilon == other.epsilon):
            return True
        else:
            return False

    def __repr__(self): # pragma: no cover
        return str((self.name, self.at_num, self.mass, self.charge, self.ptype, self.sigma, self.epsilon))

class Cmap(EntryBase):
    def __init__(self, names, func, rows='', cols='', data='', comment=''):
        self.names = names
        self.func = func
        self.rows = rows
        self.cols = cols
        self.data = data
        self.comment = comment

    def to_str(self):
        name = '{: <7} {: <7} {: <7} {: <7} {: <7} {: <6} {: <6} {: <6}'.format(*self.names, self.func, self.rows, self.cols)
        if self.data:
            data = np.array_split(self.data, len(self.data)//10)
            data = [' '.join(['{: <16}'.format(num) for num in line]) for line in data]
            data = '\\\n'.join(data)
            return '{} ; {}'.format(name + '\\\n' + data, self.comment)
        else: # pragma: no cover
            return '{} ; {}'.format(name + self.data, self.comment)
    def __repr__(self): # pragma: no cover
        name = '{: <7} {: <7} {: <7} {: <7} {: <7}'.format(*self.names)
        func = '{: <6} {: <6} {: <6}'.format(self.func, self.rows, self.cols)
        return str(name, func)

    def update_idx(self, mapping):
        new_names = []
        for name in self.names:
            new_names.append(mapping[int(name)])
        self.names = new_names

    def equal_idx(self, other):
        if self.names == other.names:
            return True
        else:
            return False

    def add_stateB(self, other):
        new = copy.copy(self)
        if new == other:
            new.comment = join_comment(new.comment, other.comment)
            return new
        elif new.equal_idx(other):
            raise NotImplementedError('Cmap alchemcial transformation not supported')
        else:
            return False

    def __eq__(self, other):
        if (self.names, self.func, self.rows, self.cols, self.data) == (other.names, other.func, other.rows, other.cols, other.data):
            return True
        else:
            return False

class Dummy_Atomtype(Atomtype):
    def __init__(self):
        super().__init__('DUM', 0, '1.008000', 0.0, 'A', 0.0, 0.0)

class Atom(EntryBase):
    def __init__(self, nr, type, resnr, residue, atom, cgnr, charge, mass, typeB = '', chargeB = '', massB = '', comment=''):
        self.nr = int(nr)
        self.type = type
        self.resnr = int(resnr)
        self.residue = residue
        self.atom = atom
        self.cgnr = int(cgnr)
        self.charge = charge
        self.mass = mass
        self.typeB = typeB
        self.chargeB = chargeB
        self.massB = massB
        self.comment = comment

    def to_str(self):
        # not printing the B states
        if self.typeB:
            return '{: <5} {: <11} {: <7} {: <7} {: <7} {: <7} {: <12} {: <11} {: <11} {: <12} {: <11}; {}'.format(
                self.nr, self.type, self.resnr, self.residue, self.atom, self.cgnr, self.charge, self.mass,
                self.typeB, self.chargeB, self.massB, self.comment)
        else:
            return '{: <5} {: <11} {: <7} {: <7} {: <7} {: <7} {: <12} {: <11}; {}'.format(self.nr, self.type,
                                                                                           self.resnr, self.residue,
                                                                                           self.atom, self.cgnr,
                                                                                           self.charge, self.mass,
                                                                                           self.comment)
    def __repr__(self): # pragma: no cover
        return str((self.nr, self.type, self.resnr, self.residue, self.atom, self.cgnr, self.charge, self.mass))

    def add_coul0(self):
        '''Add a new B state with charge of 0'''
        if self.typeB == '':
            self.typeB = self.type
        if self.massB == '':
            self.massB = self.mass
        self.chargeB = 0

    def to_stateB(self):
        if self.typeB != '':
            self.type = self.typeB
            self.typeB = ''
        if self.massB != '':
            self.mass = self.massB
            self.massB = ''
        if self.chargeB != '':
            self.charge = self.chargeB
            self.chargeB = ''

    def str_charge(self):
        if type(self.charge) == float:
            self.charge = '{:.6f}'.format(self.charge)
        if type(self.chargeB) == float:
            self.chargeB = '{:.6f}'.format(self.chargeB)

    def __eq__(self, other):
        if (self.nr, self.type, self.resnr, self.residue, self.atom, self.cgnr, self.charge, self.mass,
            self.typeB, self.chargeB, self.massB) == (other.nr, other.type, other.resnr, other.residue, other.atom,
                                                      other.cgnr, other.charge, other.mass, other.typeB,
                                                      other.chargeB, other.massB):
            return True
        else:
            return False

class Pair(EntryBase):
    def __init__(self, i, j, func, comment=''):
        # Make sure that the entry is sorted
        if i > j:
            i, j = j, i

        self.i = int(i)
        self.j = int(j)
        self.func = int(func)
        self.comment = comment

    def to_str(self):
        return '{: <7} {: <7} {: <6}; {}'.format(self.i, self.j, self.func, self.comment)

    def idx_in(self, other):
        if self.i in other or self.j in other:
            return True
        else:
            return False

    def __eq__(self, other):
        if (self.i == other.i) and (self.j == other.j) and \
                (self.func == other.func):
            return True
        else:
            return False

    def idx_less(self, other):
        if (self.i, self.j) < (other.i, other.j):
            return True
        else:
            return False

    def __contains__(self, idx):
        return idx in [self.i, self.j]

    def __repr__(self): # pragma: no cover
        return str((self.i, self.j, self.func,))

    def equal_idx(self, other):
        if (self.i == other.i) and (self.j == other.j) and (self.func == other.func):
            return True
        else:
            return False

    def check_other_list(self, mapping_self, mapping_other):
        '''This function checks if the corresponding entry in the mapping_other is None'''
        if self.i in mapping_self and mapping_other[mapping_self.index(self.i)] is None:
            return True
        elif self.j in mapping_self and mapping_other[mapping_self.index(self.j)] is None:
            return True
        else:
            return False

    def update_idx(self, mapping, sort=True):
        i = mapping[self.i]
        j = mapping[self.j]
        if sort:
            if i > j:
                i, j = j, i
        self.i, self.j = i, j

    def add_stateB(self, other):
        if self == other:
            new = copy.copy(self)
            new.comment = join_comment(new.comment, other.comment)
            return self
        else:
            return False

class Bond(Pair):
    def __init__(self, i, j, func, comment='', **kwargs):
        # Make sure that the entry is sorted
        # b0 = '', kb = '', b0B = '', kbB = '',
        if i > j:
            i, j = j, i
        self.i = int(i)
        self.j = int(j)
        self.func = int(func)
        if self.func == 1:
            self.b0 = kwargs['b0']
            self.kb = kwargs['kb']
            if 'b0B' in kwargs:
                self.b0B = kwargs['b0B']
                self.kbB = kwargs['kbB']
            else:
                self.b0B = ''
                self.kbB = ''
        self.comment = comment

    def to_str(self):
        if self.func == 1:
            return '{: <7} {: <7} {: <6} {: <10} {: <14} {: <10} {: <14} ; {}'.format(self.i, self.j, self.func,
                                                                                          self.b0, self.kb,
                                                                                          self.b0B, self.kbB,
                                                                                          self.comment)

    def __eq__(self, other):
        if (self.i == other.i) and (self.j == other.j) and \
                (self.func == other.func) and (self.b0 == other.b0) and \
                (self.kb == other.kb):
            return True
        else:
            return False

    def add_stateB(self, other):
        new = copy.copy(self)
        if new == other:
            new.comment = join_comment(new.comment, other.comment)
            return new
        elif new.equal_idx(other):
            new.comment = join_comment(new.comment, other.comment)
            if new.func == 1:
                new.b0B = other.b0
                new.kbB = other.kb
            return new
        else:
            return False

    def __repr__(self):
        return str((self.i, self.j, self.func, self.b0, self.kb))

    def create_dummy(self):
        dummy = copy.copy(self)
        dummy.comment = 'Dummy'
        if dummy.func == 1:
            dummy.kb = '0.0'
        return dummy

    def to_stateB(self):
        if self.b0B != '':
            self.b0 = self.b0B
            self.b0B = ''
        if self.kbB != '':
            self.kb = self.kbB
            self.kbB = ''


class Angle(EntryBase):
    def __init__(self, i, j, k, func, comment='', **kwargs):
        # Make sure that the entry is sorted
        if i > k:
            i, k = k, i
        self.i = int(i)
        self.j = int(j)
        self.k = int(k)
        self.func = int(func)
        if self.func == 1:
            self.th0 = kwargs['th0']
            self.cth = kwargs['cth']
            if 'th0B' in kwargs:
                self.th0B = kwargs['th0B']
                self.cthB = kwargs['cthB']
            else:
                self.th0B = ''
                self.cthB = ''
        else:
            raise NotImplementedError('Function type {} not implemented yet'.format(self.func))

        self.comment = comment
    def to_str(self):
        if self.func == 1:
            return '{: <7} {: <7} {: <7} {: <6} {: <10} {: <14} {: <10} {: <14} ; {}'.format(self.i, self.j,
                                                                                                 self.k, self.func,
                                                                                                 self.th0, self.cth,
                                                                                                 self.th0B, self.cthB,
                                                                                                 self.comment)
    def __eq__(self, other):
        if (self.i == other.i) and (self.j == other.j) and \
                (self.func == other.func) and (self.k == other.k) and \
                (self.th0 == other.th0) and (self.cth == other.cth):
            return True
        else:
            return False
    def __repr__(self):
        return str((self.i, self.j, self.k, self.func, self.th0, self.cth))

    def update_idx(self, mapping, sort=True):
        i = mapping[self.i]
        j = mapping[self.j]
        k = mapping[self.k]
        if sort:
            if i > k:
                i, k = k, i
        self.i, self.j, self.k = i, j, k

    def idx_in(self, other):
        if self.i in other or self.j in other or self.k in other:
            return True
        else:
            return False

    def idx_less(self, other):
        if (self.i, self.j, self.k) < (other.i, other.j, other.k):
            return True
        else:
            return False

    def equal_idx(self, other):
        if (self.i == other.i) and (self.j == other.j) and (self.k == other.k) and (self.func == other.func):
            return True
        else:
            return False

    def check_other_list(self, mapping_self, mapping_other):
        '''This function checks if the corresponding entry in the mapping_other is None'''
        if self.i in mapping_self and mapping_other[mapping_self.index(self.i)] is None:
            return True
        elif self.j in mapping_self and mapping_other[mapping_self.index(self.j)] is None:
            return True
        elif self.k in mapping_self and mapping_other[mapping_self.index(self.k)] is None:
            return True
        else:
            return False

    def add_stateB(self, other):
        new = copy.copy(self)
        if new == other:
            new.comment = join_comment(new.comment, other.comment)
            return new
        elif new.equal_idx(other):
            new.comment = join_comment(new.comment, other.comment)
            if new.func == 1:
                new.th0B = other.th0
                new.cthB = other.cth
            return new
        else:
            return False

    def create_dummy(self):
        dummy = copy.copy(self)
        dummy.comment = 'Dummy'
        if dummy.func == 1:
            dummy.cth = '0.0'
        return dummy

    def to_stateB(self):
        if self.th0B != '':
            self.th0 = self.th0B
            self.th0B = ''
        if self.cthB != '':
            self.cth = self.cthB
            self.cthB = ''

class Dihedral(EntryBase):
    def __init__(self, i, j, k, l, func, comment='', **kwargs):
        # Make sure that the entry is sorted
        if i > l:
            i, j, k, l = l, k, j, i
        self.i = int(i)
        self.j = int(j)
        self.k = int(k)
        self.l = int(l)
        self.func = int(func)
        if self.func in [1, 4, 9]:
            if kwargs['phase']:
                self.phase = float(kwargs['phase'])
            else:
                self.phase = ''
            self.kd = kwargs['kd']
            if kwargs['pn']:
                self.pn = int(kwargs['pn'])
            else:
                self.pn = kwargs['pn']
            if 'phaseB' in kwargs:
                self.phaseB = float(kwargs['phaseB'])
                self.kdB = kwargs['kdB']
                self.pnB = int(kwargs['pnB'])
            else:
                self.phaseB = ''
                self.kdB = ''
                self.pnB = ''
        elif self.func == 3:
            self.C0 = float(kwargs['C0'])
            self.C1 = float(kwargs['C1'])
            self.C2 = float(kwargs['C2'])
            self.C3 = float(kwargs['C3'])
            self.C4 = float(kwargs['C4'])
            self.C5 = float(kwargs['C5'])
            if 'C0B' in kwargs:
                self.C0B = float(kwargs['C0B'])
                self.C1B = float(kwargs['C1B'])
                self.C2B = float(kwargs['C2B'])
                self.C3B = float(kwargs['C3B'])
                self.C4B = float(kwargs['C4B'])
                self.C5B = float(kwargs['C5B'])
            else:
                self.C0B = ''
                self.C1B = ''
                self.C2B = ''
                self.C3B = ''
                self.C4B = ''
                self.C5B = ''
        else:
            raise NotImplementedError('Function type {} not implemented yet'.format(self.func))

        self.comment = comment
    def to_str(self):
        if self.func in [1, 4, 9]:
            return '{: <7} {: <7} {: <7} {: <7} {: <6} {: <10} {: <14} {: <6} {: <10} {: <14} {: <6}; {}'.format(
                self.i, self.j, self.k, self.l, self.func,
                self.phase, self.kd, self.pn,
                self.phaseB, self.kdB, self.pnB,
                self.comment)
        elif self.func == 3:
            return '{: <7} {: <7} {: <7} {: <7} {: <6} {: <14} {: <14} {: <14} {: <14} {: <14} {: <14} {: <14} {: <14} {: <14} {: <14} {: <14} {: <14}; {}'.format(
                self.i, self.j, self.k, self.l, self.func,
                self.C0, self.C1, self.C2, self.C3, self.C4, self.C5,
                self.C0B, self.C1B, self.C2B, self.C3B, self.C4B, self.C5B,
                self.comment)

    def __repr__(self):
        if self.func in [1, 4, 9]:
            return str((self.i, self.j, self.k, self.l, self.func, self.phase, self.kd, self.pn))
        elif self.func == 3:
            return str((self.i, self.j, self.k, self.l, self.func, self.C0, self.C1, self.C2, self.C3, self.C4, self.C5))
    def __eq__(self, other):
        if self.func in [1, 4, 9]:
            if (self.func, self.kd, self.pn) == (other.func, other.kd, other.pn) and (
                    self.phase == other.phase or np.isclose(self.phase, other.phase, atol=0.01)):
                if (self.i, self.j, self.k, self.l) == (other.i, other.j, other.k, other.l):
                    return True
                elif (self.func == 4) and (self.i, self.j, self.k, self.l) == (other.i, other.l, other.k, other.j):
                    return True
                else:
                    return False
            else:
                return False
        elif self.func == 3:
            self_values = (self.C0, self.C1, self.C2, self.C3, self.C4, self.C5,
                           self.C0B, self.C1B, self.C2B, self.C3B, self.C4B, self.C5B)
            other_values = (other.C0, other.C1, other.C2, other.C3, other.C4, other.C5,
                           other.C0B, other.C1B, other.C2B, other.C3B, other.C4B, other.C5B)
            if (self.func) == (other.func) and \
                    (self_values == other_values or np.isclose(self_values, other_values, atol=0.01)) and \
                    (self.i, self.j, self.k, self.l) == (other.i, other.j, other.k, other.l):
                return True
            else:
                return False
        else:
            return False


    def __contains__(self, idx):
        return idx in [self.i, self.j, self.k, self.l]

    def update_idx(self, mapping, sort=True):
        i = mapping[self.i]
        j = mapping[self.j]
        k = mapping[self.k]
        l = mapping[self.l]
        if sort:
            if i > l:
                i, j, k, l = l, k, j, i
        self.i, self.j, self.k, self.l = i, j, k, l

    def idx_in(self, other):
        if self.i in other or self.j in other or self.k in other or self.l in other:
            return True
        else:
            return False

    def idx_less(self, other):
        if (self.i, self.j, self.k, self.l) < (other.i, other.j, other.k, other.l):
            return True
        else:
            return False

    def equal_idx(self, other):
        if (self.i == other.i) and (self.j == other.j) and (self.k == other.k) and (self.l == other.l) and (self.func == other.func):
            return True
        # improper
        elif (self.func == 4) and (self.i, self.j, self.k, self.l, self.func) == (other.i, other.l, other.k, other.j, self.func):
            return True
        else:
            return False

    def check_other_list(self, mapping_self, mapping_other):
        '''This function checks if the corresponding entry in the mapping_other is None'''
        if self.i in mapping_self and mapping_other[mapping_self.index(self.i)] is None:
            return True
        elif self.j in mapping_self and mapping_other[mapping_self.index(self.j)] is None:
            return True
        elif self.k in mapping_self and mapping_other[mapping_self.index(self.k)] is None:
            return True
        elif self.l in mapping_self and mapping_other[mapping_self.index(self.l)] is None:
            return True
        else:
            return False

    def add_stateB(self, other):
        new = copy.copy(self)
        if new == other:
            new.comment = join_comment(new.comment, other.comment)
            return new
        elif new.equal_idx(other):
            if new.func in [1, 4, 9] and new.pn == other.pn and np.isclose(self.phase, other.phase, atol=0.01):
                new.comment = join_comment(new.comment, other.comment)
                new.kdB = other.kd
                new.phaseB = other.phase
                new.pnB = other.pn
            else:
                return False
            return new
        else:
            return False

    def create_dummy(self):
        dummy = copy.copy(self)
        dummy.comment = 'Dummy'
        if dummy.func in [1, 4, 9]:
            dummy.kd = '0.0'
        return dummy

    def to_stateB(self):
        if self.func in [1, 4, 9]:
            if self.phaseB != '':
                self.phase = self.phaseB
                self.phaseB = ''
            if self.pnB != '':
                self.pn = self.pnB
                self.pnB = ''
            if self.kdB != '':
                self.kd = self.kdB
                self.kdB = ''
        else:
            if self.C0B != '':
                self.C0, self.C1, self.C2, self.C3, self.C4, self.C5 = self.C0B, self.C1B, self.C2B, self.C3B, self.C4B, self.C5B
                self.C0B, self.C1B, self.C2B, self.C3B, self.C4B, self.C5B = [0, ] * 6


