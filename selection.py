"""Turning an Olex2 selection into something the measurement backends can consume."""
from typing import List

from helper_functions import (get_id_from_label, get_label_from_id, get_neighbours,
                              get_orm_atoms, get_part, get_xyz)


class MolecularStructure:
    """A set of labelled atoms handed to SHAPE / OctaDist / cosmochlore.

    When the structure is centered the central atom is at index 0.
    `coords` is always a list of (x, y, z) float tuples in orthogonal
    coordinates, one per entry in `labels`.
    """

    def __init__(self, coords, labels):
        self.coords = [tuple(float(c) for c in xyz) for xyz in coords]
        self.labels = [str(label) for label in labels]

        if len(self.coords) != len(self.labels):
            raise ValueError(f'Malformed structure: {len(self.labels)} labels '
                             f'but {len(self.coords)} coordinates.')

    def __len__(self):
        return len(self.labels)


class AtomSelection:
    """The current Olex2 selection, and the operations that grow or reshape it.

    `labels` holds the labels as Olex2 reports them, so entries taken straight
    from the selection may carry a symmetry suffix (N2_$1) while entries added
    by add_neighbours() are plain ORM labels. `tags`, `coords` and `parts` are
    always parallel to `labels`; every method here keeps all four in step.
    """

    def __init__(self, selection_string):
        self.orm_atoms = get_orm_atoms()

        self.labels = selection_string.split(' ') if selection_string else []
        self.tags = [get_id_from_label(label, self.orm_atoms) for label in self.labels]
        self.coords = [get_xyz(tag) for tag in self.tags]
        self.parts = [get_part(tag) for tag in self.tags]

    def __len__(self):
        return len(self.labels)

    def add_neighbours(self):
        """Appends the atoms bonded to each currently selected atom."""
        if not self.labels:
            print("Could not find neighbours. Selection is empty.")
            return

        for sel_label in self.labels.copy():
            neighbour_tags = next((atom['neighbours']
                                   for atom in self.orm_atoms
                                   if atom['label'] == sel_label),
                                  None)

            if not neighbour_tags:
                print(f'Could not find neighbours for {sel_label}.')
                continue

            unique_neighbours = []
            for neighbour_tag in neighbour_tags:
                if neighbour_tag not in unique_neighbours:
                    unique_neighbours.append(neighbour_tag)

            for neighbour in unique_neighbours:
                # A neighbour outside the ASU arrives as a tuple that already carries
                # the coordinates of its symmetry-generated image; one inside the ASU
                # is a bare tag whose coordinates we look up.
                if isinstance(neighbour, tuple):
                    tag = neighbour[0]
                    coord = tuple(float(c) for c in neighbour[1])
                else:
                    tag = neighbour
                    coord = get_xyz(neighbour)

                self.labels.append(get_label_from_id(tag, self.orm_atoms))
                self.tags.append(tag)
                self.coords.append(coord)
                self.parts.append(get_part(tag))

    def remove_duplicates(self, tolerance: int = 4):
        """Drops atoms that sit on a coordinate already present in the selection.

        Keyed on position rather than label: two symmetry-generated images share an
        ORM label but are distinct atoms, so a label-based test would discard real
        vertices.
        """
        seen = set()
        keep = []
        for i, coord in enumerate(self.coords):
            key = tuple(round(c, tolerance) for c in coord)
            if key in seen:
                continue
            seen.add(key)
            keep.append(i)

        self._keep_indices(keep)

    def merge_ligands(self):
        """Collapses each set of mutually bonded ligands into a single centroid.

        Used for pi-bonded ligands, where the coordination polyhedron should see the
        centroid of the bonded fragment rather than each of its atoms. Returns the
        fragments that were merged.
        """

        ####################
        # NEIGHBOUR SEARCH #
        ####################

        bonded_pairs = []
        # Iterate over the ligands (atoms added with add_neighbours())
        for label, tag in zip(self.labels[1:], self.tags[1:]):  # For each ligand

            # Get the neighbours of each ligand atom
            _, nei_uniques = get_neighbours([label], self.orm_atoms)

            # If the neighbour of the ligand is also a ligand, add the bonded pair to the list.
            for nei_tag in nei_uniques:
                if nei_tag in self.tags[1:]:
                    bonded_pairs.append((nei_tag, tag))

        # If there are no bonded pairs of ligands, stop the function.
        if not bonded_pairs:
            print('Nothing to merge.')
            return []

        # Create an adjacency list from the bonded pairs.
        adj = {}
        for a, b in bonded_pairs:  # Keys are all atoms, values are sets of their neighbours.
            # adj = {atom1: {atom2, atom3}, atom2: {atom1}, atom3: {atom1}}
            adj.setdefault(a, set()).add(b)  # If a is absent, add it with an empty set. To this set, add b
            adj.setdefault(b, set()).add(a)

        #####################
        # FRAGMENT BUILDING #
        #####################

        visited = set()  # Keeps track of atoms assigned to a fragment
        fragments = []   # Keeps list of final groups of fragments

        # Iterate over all atoms in the adjacency list:
        for atom in adj:
            if atom in visited:  # Skip if the atom was visited already
                continue

            stack = [atom]      # Add the atom to the stack (to-process list)
            fragment = set()    # Create a new fragment as an empty set

            while stack:  # While there are items in the stack
                current = stack.pop()       # Pop one atom
                if current in fragment:     # If the atom is already in the fragment, skip
                    continue
                fragment.add(current)       # Add the current atom to the fragment

                # Extend the stack with the adjacent atoms of the current one,
                # minus the ones already in the fragment
                stack.extend(adj[current] - fragment)

            visited |= fragment             # Mark all atoms as visited (|= merges sets)
            fragments.append(fragment)      # Adds the fragment to the final fragments list

        ##################
        # LIGAND MERGING #
        ##################

        tag_to_idx = {tag: i for i, tag in enumerate(self.tags)}
        centroid_at = {}  # index -> (label, coords) to place there
        skip = set()      # indices to drop

        for i, frag in enumerate(fragments):
            idxs = sorted(tag_to_idx[tag] for tag in frag)

            cx = sum(self.coords[j][0] for j in idxs) / len(idxs)
            cy = sum(self.coords[j][1] for j in idxs) / len(idxs)
            cz = sum(self.coords[j][2] for j in idxs) / len(idxs)

            first = idxs[0]
            centroid_at[first] = (f'Z{i}', (cx, cy, cz))
            skip.update(idxs[1:])  # keep 'first', drop the rest

        keep = [j for j in range(len(self.tags)) if j not in skip]
        self._keep_indices(keep)

        # Overwrite each surviving fragment atom with its fragment's centroid. The
        # part is inherited from the atom that was kept, so `parts` stays meaningful.
        for new_idx, old_idx in enumerate(keep):
            if old_idx in centroid_at:
                label, centroid = centroid_at[old_idx]
                self.labels[new_idx] = label
                self.coords[new_idx] = centroid
                self.tags[new_idx] = -(old_idx + 1)  # placeholder: no longer a real ORM atom

        return fragments

    def _keep_indices(self, keep):
        """Reduces every parallel list to `keep`, so they cannot drift out of step."""
        self.labels = [self.labels[i] for i in keep]
        self.tags = [self.tags[i] for i in keep]
        self.coords = [self.coords[i] for i in keep]
        self.parts = [self.parts[i] for i in keep]


def split_by_parts(selection: AtomSelection) -> List[MolecularStructure]:
    """Splits a disordered selection into one structure per disorder component.

    Atoms in part 0 are shared by every component, so each returned structure is
    part 0 plus one of the other parts. A selection with fewer than two disordered
    components is returned unchanged, as a single structure.
    """
    unique_parts = sorted(set(selection.parts))

    if len(unique_parts) <= 2:
        return [MolecularStructure(selection.coords, selection.labels)]

    # Separate atoms by part:
    by_part = {}
    for part, label, coord in zip(selection.parts, selection.labels, selection.coords):
        labels, coords = by_part.setdefault(part, ([], []))
        labels.append(label)
        coords.append(coord)

    base_labels, base_coords = by_part[unique_parts[0]]

    structures = []
    for part in unique_parts[1:]:
        labels, coords = by_part[part]
        structures.append(MolecularStructure(base_coords + coords, base_labels + labels))

    return structures
