import os
from .field import Default, Moleculetype, Atomtypes, Atoms, Bonds, Pairs, Angles, Dihedrals, Cmaptypes, Cmaps
# from top import Topology

class Reader():
    def __init__(self, topology, filename):
        content_dict = {}

        with open(filename, 'r') as f:
            txt = f.read()

        # To change \ into none
        txt = txt.replace('\\\n', ' ')

        sections = txt.split('[')

        # Ensure at least one field exists
        content_dict['moleculetype'] = Moleculetype()
        for section in sections:
            if ']' in section:
                name = section[:section.index(']')].strip().lower()
                field = self._section_check(content_dict, name)
                if field is None:
                    print('Section {} is ignored'.format(name))
                    continue
                else:
                    content_dict[name] = field
                    content = section[section.index(']')+1:].strip()
                    for line in content.split('\n'):
                        if line and line[0] == ';':
                            field.add_comment(line)
                        elif name == 'defaults':
                            field.edit(line)
                        elif name == 'moleculetype':
                            field.edit(line)
                        else:
                            field.add_entry(line)

            else:
                # Assume the whole section is comments
                comments = section.split('\n')
                # Add the comments to the Moleculetype
                field = content_dict['moleculetype']
                for line in comments:
                    field.add_comment(line)
        topology.name = os.path.splitext(os.path.split(filename)[-1])[0]
        topology.content_dict = content_dict

    def _section_check(self, content_dict, name):
        '''Check for the section name and return the relevant section reader'''
        if name in content_dict:
            return content_dict[name]
        elif name == 'defaults':
            return Default()
        elif name == 'atomtypes':
            return Atomtypes()
        elif name == 'atoms':
            return Atoms()
        elif name == 'bonds':
            return Bonds()
        elif name == 'pairs':
            return Pairs()
        elif name == 'angles':
            return Angles()
        elif name == 'dihedrals':
            return Dihedrals()
        elif name == 'cmap':
            return Cmaps()
        elif name == 'cmaptypes':
            return Cmaptypes()
        else:
            return None

class itpWriter():
    def __init__(self, topology):
        self.top = topology
    def write_itp(self, f):
        f.write(self.top.content_dict['moleculetype'].to_str())
        f.write('\n')
        f.write(self.top.content_dict['atoms'].to_str())
        f.write('\n')
        f.write(self.top.content_dict['bonds'].to_str())
        f.write('\n')
        f.write(self.top.content_dict['pairs'].to_str())
        f.write('\n')
        f.write(self.top.content_dict['angles'].to_str())
        f.write('\n')
        f.write(self.top.content_dict['dihedrals'].to_str())
        f.write('\n')
        if 'cmap' in self.top.content_dict:
            f.write(self.top.content_dict['cmap'].to_str())
            f.write('\n')

    def write_top(self, f):
        if 'defaults' in self.top.content_dict:
            f.write(self.top.content_dict['defaults'].to_str())
            f.write('\n')
        f.write(self.top.content_dict['atomtypes'].to_str())
        f.write('\n')
        if 'cmaptypes' in self.top.content_dict:
            f.write(self.top.content_dict['cmaptypes'].to_str())
            f.write('\n')

# class rtpWriter():
#     def __init__(self, topology, f):
#         for comment in topology.content_list:
#             if isinstance(comment, Comment):
#                 f.write(comment.to_str() + '\n')
#         f.write('\n')
#         f.write('[ {} ]\n'.format(topology.name))
#
#         f.write(topology.content_dict['atoms'].to_str())
#         f.write('\n')
#         f.write(topology.content_dict['bonds'].to_str())
#         f.write('\n')
#         f.write(topology.content_dict['impropers'].to_str())
#         f.write('\n')
#
# class atpWriter():
#     def __init__(self, topology, f):
#         atomtypes = topology.content_dict['atomtypes']
#         for atom in atomtypes:
#             if not isinstance(atom, Comment):
#                 f.write('{: <15} {: <16}\n'.format(atom._name, atom._mass))
