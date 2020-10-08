from MDAnalysis.analysis.align import alignto
from MDAnalysis.core.universe import Merge
import MDAnalysis as mda

def merge_crd(crd1, crd2, mappings, filename=None):
    crd1 = mda.Universe(crd1)
    crd2 = mda.Universe(crd2)
    # The mapping contains two mapping files:
    # From molecule 1 to the merged molecule
    # And from molecule 2 to the merged molecule
    mapping1 = mappings[0]
    mapping2 = mappings[1]
    # Find the atom id in the merged molecule, which is in both molecule 1
    # and molecule 2
    common = set(mapping1.values()).intersection(mapping2.values())

    substructure_1 = []
    # Invert the mapping file so that it is from
    # Common merged id to molecule 1 id
    invert_mapping_1 = {mapping1[key]: key for key in mapping1
                        if mapping1[key] in common}
    # Get the atoms in molecule 1 according to the order in the Common merged id
    for id in sorted(list(invert_mapping_1.keys())):
        substructure_1.append(crd1.select_atoms('bynum {}'.format(
            invert_mapping_1[id])))

    substructure_2 = []
    # Invert the mapping file so that it is from
    # Common merged id to molecule 2 id
    invert_mapping_2 = {mapping2[key]: key for key in mapping2
                        if mapping2[key] in common}
    # Get the atoms in molecule 2 according to the order in the Common merged id
    for id in sorted(list(invert_mapping_2.keys())):
        substructure_2.append(crd2.select_atoms('bynum {}'.format(
            invert_mapping_2[id])))

    alignto(sum(substructure_2), sum(substructure_1), match_atoms=False)

    # Reorder the atomgroups in substructure_2
    substructure_2 = sorted(substructure_2, key=lambda x: x.ids)

    # align the original sutructure 2 to the common structure
    alignto(crd2, sum(substructure_2), select='bynum {}'.format(
        ' '.join([str(atom.ids[0]) for atom in substructure_2])))

    # Compile the dict from merged structure to (structure 1, structure 2)
    merged_length = max(max(mapping1.values()), max(mapping2.values()))
    merged_dict = {id: [None, None] for id in range(1, merged_length+1)}
    for id in mapping1:
        merged_dict[mapping1[id]][0] = id
    for id in mapping2:
        merged_dict[mapping2[id]][1] = id

    merged_structure = []
    for id in sorted(merged_dict.keys()):
        if merged_dict[id][0]:
            merged_structure.append(crd1.select_atoms('bynum {}'.format(
            merged_dict[id][0])))
        else:
            merged_structure.append(crd2.select_atoms('bynum {}'.format(
            merged_dict[id][1])))
    merged_structure = Merge(*merged_structure)
    merged_structure.atoms.ids = merged_structure.atoms.ix + 1
    if filename:
        merged_structure.atoms.write(filename)
    return merged_structure

def extract_from_merged(merged_crd, mapping, filename=None):
    merged_crd = mda.Universe(merged_crd)
    substructure = []
    for id in sorted(list(mapping.keys())):
        substructure.append(merged_crd.select_atoms('bynum {}'.format(
            mapping[id])))
    substructure = sum(substructure)
    substructure.atoms.ids = substructure.atoms.ix + 1
    if filename:
        substructure.atoms.write(filename)
    return substructure