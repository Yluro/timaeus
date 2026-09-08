import os

import olex
import olx
from olexFunctions import OlexFunctions

OV = OlexFunctions()

from reload_all import reload_all

from autoshape import (ShapeCalculation, can_find_shape_msg, find_shape,
                       print_shape_table, run_shape)
from constants import octadist_citation, shape21_citation
from helper_functions import (get_neighbours, get_neighbours_on_sel, get_selected_atoms,
                              get_xyz_sel, print_console_bs, print_orm, test_selection_class)
from octahedral_distortion import CalcDistortion
from selection import AtomSelection, split_by_parts

debug = bool(OV.GetParam("olex2.debug", False))

instance_path = OV.DataDir()

try:
    from_outside = False
    p_path = os.path.dirname(os.path.abspath(__file__))
except:
    from_outside = True
    p_path = os.path.dirname(os.path.abspath("__file__"))

l = open(os.sep.join([p_path, 'def.txt'])).readlines()
d = {}
for line in l:
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    d[line.split("=")[0].strip()] = line.split("=")[1].strip()

p_name = d['p_name']
p_htm = d['p_htm']
p_img = eval(d['p_img'])
p_scope = d['p_scope']

OV.SetVar('SymmetryMeasurements_plugin_path', p_path)

from PluginTools import PluginTools as PT


def as_bool(value) -> bool:
    """Normalises the assorted true/false spellings that reach us from the GUI.

    Olex2 hands checkbox state back as a string, while a phil bool comes back as a
    real bool, so both have to be accepted at every param read.
    """
    if isinstance(value, str):
        return value.strip().lower() in ('true', '1', 'yes', 'on')
    return bool(value)


# MAIN LOGIC FUNCTIONS.
def autoSHAPE():
    print('\n' + '-' * 50)
    print('Simple continuous Shape Analysis Using autoSHAPE')
    if not can_find_shape_msg():
        print('SHAPE executable not found in PATH.')
        return False

    sel_string = olex.f('sel()')
    if sel_string == '':
        print('Invalid atom selection: no atoms selected.')
        return False

    selection = AtomSelection(sel_string)

    if len(selection) > 1:
        # Multiple atoms selected: treat them as the vertices of a non-centered shape.
        selection.remove_duplicates()
        centered = False
    else:
        # A single atom selected: grow it into a coordination polyhedron.
        selection.add_neighbours()
        if as_bool(OV.GetParam('symmetrymeasurements.merge_ligands', False)):
            selection.merge_ligands()
        centered = True

    # Returns one structure per disorder component, or a single structure when the
    # selection is not disordered.
    structures = split_by_parts(selection)

    for i, structure in enumerate(structures):
        try:
            shape_measurement = ShapeCalculation(structure, f'{olx.FileName()}_{i}',
                                                 centered, ['%fullout'])
        except ValueError as e:
            print(f'Skipping part {i}: {e}')
            continue

        folder = shape_measurement.write_tab(olx.FilePath())
        if folder is None:
            continue

        for f in run_shape(folder):
            print_shape_table(os.path.join(folder, f'{f}.tab'))

    print(shape21_citation)
    return True


def autoOCTADIST():
    # Get the selected atoms.
    sel_string = olex.f('sel()')
    if sel_string == '':
        print('Invalid atom selection: no atoms selected.')
        return False

    selection = AtomSelection(sel_string)

    # Exit if the selection is anything other than the single central atom.
    if len(selection) != 1:
        print(f'Invalid atom selection: expected 1 atom, found {len(selection)}.')
        return False

    # Add coordinated atoms to the current selection.
    selection.add_neighbours()

    for structure in split_by_parts(selection):
        try:
            calculation = CalcDistortion(structure)
        except ValueError as e:
            print(f'Invalid polyhedra for {structure.labels[0]}: {e}')
            continue

        calculation.print_results(olx.FileName())
        calculation.draw_octahedron()

    # Citation
    print('\nThis calculations were made using a reimplementation of the OctaDist '
          'algorithm by David J. Harding et al.')
    print(octadist_citation)
    return True


def shape_status_html():
    where = find_shape()
    found = where is not None
    color = OV.GetParam('gui.green') if found else OV.GetParam('gui.grey')
    text = f'SHAPE executable found at: {where}' if found else 'Unable to find shape.exe in the system path.'
    return f"<font color='{color}'>{text}</font>"


class SymmetryMeasurements(PT):
    def __init__(self):
        super(SymmetryMeasurements, self).__init__()
        self.p_name = p_name
        self.p_path = p_path
        self.p_scope = p_scope
        self.p_htm = p_htm
        self.p_img = p_img
        self.deal_with_phil(operation='read')
        self.print_version_date()
        if not from_outside:
            self.setup_gui()

        # Main entry points.
        OV.registerFunction(autoSHAPE, True, "SymmetryMeasurements")
        OV.registerFunction(autoOCTADIST, True, "SymmetryMeasurements")
        OV.registerFunction(can_find_shape_msg, True, "SymmetryMeasurements")
        OV.registerFunction(shape_status_html, False, 'SymmetryMeasurements')

        # Debug panel helpers.
        OV.registerFunction(get_selected_atoms, True, "SymmetryMeasurements")
        OV.registerFunction(get_neighbours, True, "SymmetryMeasurements")
        OV.registerFunction(get_neighbours_on_sel, True, "SymmetryMeasurements")
        OV.registerFunction(get_xyz_sel, True, "SymmetryMeasurements")
        OV.registerFunction(print_console_bs, False, 'SymmetryMeasurements')
        OV.registerFunction(print_orm, False, 'SymmetryMeasurements')
        OV.registerFunction(test_selection_class, False, 'SymmetryMeasurements')
    # END Generated =======================================


SymmetryMeasurements_instance = SymmetryMeasurements()
print("Loading Symess modules.")
reload_all()
print("Symmetry Measurements by JSG loaded.")
