"""Thin wrappers over the Olex2 API.

This is the bottom layer of the plugin: it imports nothing from the plugin's own
modules, so `selection`, `autoshape` and `octahedral_distortion` can all depend on
it without creating an import cycle.

Coordinate convention: every function here returns orthogonal coordinates as a
`tuple[float, float, float]`. Atom identity is an integer ORM tag; labels are
strings that may carry a symmetry suffix (`N2_$1`).
"""
import olex
import olexex
import olx


## SMALL HELPER FUNCTIONS
def get_orm_atoms() -> list:
    """Returns the atom list of the loaded refinement model.
    """
    return olexex.OlexRefinementModel().atoms()


def get_selected_atoms() -> str:
    # Gets the selection from Olex2 -  Returns a string with atom labels.
    # If no atoms are selected, returns ''
    selection = olex.f('sel()')
    print(selection)
    return selection


def get_id_from_label(atom_label, orm=None):
    """Returns the ORM tag for `atom_label`, or None if the atom is not in the model.

    Will strip symmetry suffixes like _$1 if no tag is found."""
    orm = get_orm_atoms() if orm is None else orm
    tag = next((atom['tag'] for atom in orm if atom['label'] == atom_label), None)

    if tag is None:
        clean_label = atom_label.split('_$')[0]  # strip symmetry suffix e.g. 'N2_$1' -> 'N2'
        tag = next((atom['tag'] for atom in orm if atom['label'] == clean_label), None)

    return tag


def get_label_from_id(atom_tag, orm=None):
    """Returns the ORM label for `atom_tag`, or None if the tag is not in the model."""
    orm = get_orm_atoms() if orm is None else orm
    return next((atom['label'] for atom in orm if atom['tag'] == atom_tag), None)


def get_xyz(atom_tag) -> tuple:
    """Returns the orthogonal coordinates of `atom_tag` as a tuple of floats."""
    crd = olx.xf.au.GetAtomCrd(atom_tag)
    xyz_string = olx.xf.au.Orthogonalise(crd).split(' ')
    return tuple(float(x) for x in xyz_string)


def get_part(atom_tag) -> int:
    """Returns the disorder part of `atom_tag`. Note this takes a tag, not a label."""
    return int(olx.xf.au.GetAtomPart(atom_tag))


def get_neighbours(atom_labels, orm=None):
    """Returns (per_atom_neighbours, unique_neighbours) for the given labels.

    Both elements are always lists, empty when nothing could be resolved. A
    neighbour is either an int tag, or a tuple whose first element is the tag and
    whose second element is the coordinate of a symmetry-generated image.
    """
    orm = get_orm_atoms() if orm is None else orm

    if not atom_labels or atom_labels == [""]:
        print("Could not find neighbours. No atoms selected.")
        return [], []

    neighbours_tags_list = []
    unique_neighbours = []
    for atom_label in atom_labels:
        # next finds the first occurrence in orm in which the label matches
        # with the sel and returns the atom's neighbours as a tuple of tags:
        neighbour_tags = next((atom['neighbours'] for atom in orm if atom['label'] == atom_label), ())
        # neighbour_tags is an empty tuple if the label wasn't found in the orm (e.g. a Q-peak)
        if not neighbour_tags:
            print(f'No connected atoms to {atom_label}.')

        neighbours_tags_list.append(neighbour_tags)

        for neighbour in neighbour_tags:
            if neighbour not in unique_neighbours:
                unique_neighbours.append(neighbour)

    return neighbours_tags_list, unique_neighbours


## DEBUG HELPERS - wired to the Timaeus-debug GUI panel.
def get_xyz_sel():
    """Prints and returns the coordinates of a single selected atom."""
    selection = olex.f('sel()')
    if selection == '':
        print('No atoms selected.')
        return None

    labels = selection.split(' ')
    if len(labels) != 1:
        print(f'Invalid atom selection: expected 1 atom, found {len(labels)}.')
        return None

    tag = get_id_from_label(labels[0])
    if tag is None:
        print(f'Could not find {labels[0]} in the orm.')
        return None

    print(selection)
    return get_xyz(tag)


def get_neighbours_on_sel():
    sel = olex.f('sel()')
    return get_neighbours(sel.split(' '))


def print_console_bs():
    sel = olex.f('sel()')
    print(f'Selection: {sel}')
    print(olex.f('xf.au.GetAtomCrd()'))
    print(olx.xf.au.GetAtomCrd())
    print(olex.f('Env()'))
    print(olex.f('Envi()'))


def print_orm():
    with open('orm.txt', 'w') as f:
        for line in get_orm_atoms():
            f.write(str(line) + '\n')


def test_selection_class():
    # Imported here rather than at module scope: `selection` depends on this
    # module, so a top-level import would create a cycle.
    from selection import AtomSelection, split_by_parts

    selection = AtomSelection(olex.f('sel()'))
    print(selection.labels)
    print('Add neighbours')
    selection.add_neighbours()
    print(selection.labels)
    print('Splitting in parts.')
    for p in split_by_parts(selection):
        print(f'{p.coords}\n{p.labels}')


def as_bool(value) -> bool:
    """Normalises the assorted true/false spellings that reach us from the GUI.

    Olex2 hands checkbox state back as a string, while a phil bool comes back as a
    real bool, so both have to be accepted at every param read.
    """
    if isinstance(value, str):
        return value.strip().lower() in ('true', '1', 'yes', 'on')
    return bool(value)
