from ff_io import Reader, FFWriter, itpWriter
from top import Topology


def merge(top_list, start_residue_suffix, end_residue_suffix):
    content_dict = {}
    content_list = []
    residue_list = []


    for top in top_list:
        top.convert_ff()
        residue_list.extend(top.convert_ff_rtp(start_residue_suffix, end_residue_suffix))
        if 'defaults' in content_dict:
            assert top.content_dict['defaults'] == content_dict['defaults']
        else:
            content_dict['defaults'] = top.content_dict['defaults']

        for section_name in ['atomtypes', 'bondtypes', 'angletypes', 'dihedraltypes']:
            if section_name in content_dict:
                content_dict[section_name].union(top.content_dict[section_name])
            else:
                content_dict[section_name] = top.content_dict[section_name]

    top = Topology()
    top.content_list = content_list
    top.content_dict = content_dict
    return top, residue_list

def unique_rtp(top_list, start_residue_suffix, end_residue_suffix):
    residue_list = []
    for top in top_list:
        if not top.name in [residue.name for residue in residue_list]:
            residue_list.append(top)
        else:
            for residue in residue_list:
                if top.name == residue.name:
                    if (top.content_dict['atoms'] == top.content_dict['atoms']) and \
                            (top.content_dict['bonds'] == top.content_dict['bonds']) and \
                            (top.content_dict['impropers'] == top.content_dict['impropers']):
                        pass
                    else:
                        raise NameError
                    break
    residuetypes = {}
    for residue in residue_list:
        try:
            if not residue.original_name in residuetypes:
                residuetypes[residue.original_name] = {}

            if residue.name[-len(start_residue_suffix):] == start_residue_suffix:
                residuetypes[residue.original_name]['start'] = residue.name
            elif residue.name[-len(end_residue_suffix):] == end_residue_suffix:
                residuetypes[residue.original_name]['end'] = residue.name
            else:
                raise NameError
        except AttributeError:
            residuetypes[residue.name] = None
    return residue_list, residuetypes
