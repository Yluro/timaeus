import os

import olex
import olx
from olexFunctions import OlexFunctions

OV = OlexFunctions()

from reload_all import reload_all

from autoshape import (ShapeCalculation, can_find_shape_msg, find_shape,
                       print_shape_table, run_shape)
from constants import octadist_citation, shape21_citation
from cosmochlore import (CosmochloreError, build_cshm_args, build_csom_args,
                         build_odis_args, can_find_cosmochlore_msg, check_cosmochlore,
                         run_cosmochlore, trim_banner, write_xyz)
from helper_functions import (get_neighbours, get_neighbours_on_sel, get_selected_atoms,
                              get_xyz_sel, print_console_bs, print_orm, test_selection_class, as_bool)
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


def _prepare_structures(sel_string, merge=None):
    """Builds one MolecularStructure per disorder component from `sel_string`.

    A single selected atom is grown into a centered coordination polyhedron
    (its neighbours are added, and pi-bonded ligands are merged into a
    centroid if `merge`, or failing to the merge_ligands phil param.

    Multiple selected atoms are taken as-is as a non-centered shape.

    Returns (structures, centered), or (None, None) if the selection is empty.
    """
    if sel_string == '':
        print('Invalid atom selection: no atoms selected.')
        return None, None

    selection = AtomSelection(sel_string)

    if len(selection) > 1:
        # Multiple atoms selected: treat them as the vertices of a non-centered shape.
        selection.remove_duplicates()
        centered = False
    else:
        # A single atom selected: grow it into a coordination polyhedron.
        selection.add_neighbours()
        if merge is None:
            merge = as_bool(OV.GetParam('symmetrymeasurements.merge_ligands', False))
        if merge:
            selection.merge_ligands()
        centered = True

    # Returns one structure per disorder component, or a single structure when the strucutre is not disordered
    return split_by_parts(selection), centered


def _cosmochlore_exe_path():
    return OV.GetParam('symmetrymeasurements.cosmochlore.exe_path', '') or None


def _cosmochlore_workdir():
    """The scratch folder cosmochlore .xyz inputs and outputs are written to."""
    path = os.path.join(olx.FilePath(), 'cosmochlore')
    os.makedirs(path, exist_ok=True)
    return path


# MAIN LOGIC FUNCTIONS.
def autoSHAPE():
    print('\n' + '-' * 50)
    print('Simple continuous Shape Analysis Using autoSHAPE')
    if not can_find_shape_msg():
        print('SHAPE executable not found in PATH.')
        return False

    structures, centered = _prepare_structures(olex.f('sel()'))
    if structures is None:
        return False

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


_COSMOCHLORE_CREDIT = 'Computed using cosmochlore (https://github.com/Yluro/cosmochlore).'


def autoCSHM(shapes=None, user_shapes=None, table=False, ideal=False):
    """Continuous Shape Measures via cosmochlore, on the current selection.

    `shapes`: built-in reference-shape indices to restrict to (None = all
    applicable for the detected vertex count). `user_shapes`: paths to
    user-defined shape .yaml files. `table`/`ideal`: also write the
    corresponding cosmochlore output files next to the .xyz.
    """
    print('\n' + '-' * 50)
    print('Continuous Shape Measures using cosmochlore')
    try:
        exe = check_cosmochlore(_cosmochlore_exe_path())
    except CosmochloreError as e:
        print(e)
        return False

    structures, centered = _prepare_structures(olex.f('sel()'))
    if structures is None:
        return False

    workdir = _cosmochlore_workdir()
    ran_any = False
    for i, structure in enumerate(structures):
        xyz_path = os.path.join(workdir, f'{olx.FileName()}_{structure.labels[0]}_{i}.xyz')
        write_xyz(structure, xyz_path, comment=f'{olx.FileName()} part {i}')

        try:
            args = build_cshm_args(xyz_path, centered, shapes, user_shapes, table, ideal)
            print(trim_banner(run_cosmochlore(exe, args, cwd=workdir)))
            ran_any = True
        except CosmochloreError as e:
            print(f'cshm failed for part {i}: {e}')

    if ran_any:
        print(_COSMOCHLORE_CREDIT)
    return ran_any


def autoCSOM(point_groups=None, mode=None, vector=None, full=False, table=False,
            operated=False, samples=None, iterations=None, ignore_labels=False):
    """Continuous Symmetry Operation Measures via cosmochlore, on the current selection.

    point_groups: a space-separated string or list of Schoenflies, fails back to the symmetrymeasurements.cosmochlore.csom.point_groups phil param when not given.

    mode: centering mode (auto/first/centroid/manual), defaults to the matching phil param.

    `vector`: required 3-value centering vector when mode is 'manual'.
    """
    print('\n' + '-' * 50)
    print('Continuous Symmetry Operation Measures using cosmochlore')
    try:
        exe = check_cosmochlore(_cosmochlore_exe_path())
    except CosmochloreError as e:
        print(e)
        return False

    if point_groups is None:
        point_groups = OV.GetParam('symmetrymeasurements.cosmochlore.csom.point_groups', '')
    if isinstance(point_groups, str):
        point_groups = point_groups.split()

    if not point_groups:
        print('No point groups given. Pass e.g. '
              "spy.SymmetryMeasurements.autoCSOM('Oh D4h D3d'), or set "
              "symmetrymeasurements.cosmochlore.csom.point_groups.")
        return False

    if mode is None:
        mode = OV.GetParam('symmetrymeasurements.cosmochlore.csom.mode', 'auto')

    structures, centered = _prepare_structures(olex.f('sel()'))
    if structures is None:
        return False

    print('This may take a while depending on the number of point groups and atoms...')

    workdir = _cosmochlore_workdir()
    ran_any = False
    for i, structure in enumerate(structures):
        xyz_path = os.path.join(workdir, f'{olx.FileName()}_{structure.labels[0]}_{i}.xyz')
        write_xyz(structure, xyz_path, comment=f'{olx.FileName()} part {i}')

        try:
            args = build_csom_args(xyz_path, centered, point_groups, mode, vector, full,
                                   table, operated, samples, iterations, ignore_labels)
            print(trim_banner(run_cosmochlore(exe, args, cwd=workdir)))
            ran_any = True
        except CosmochloreError as e:
            print(f'csom failed for part {i}: {e}')

    if ran_any:
        print(_COSMOCHLORE_CREDIT)
    return ran_any


def autoODIS(full=False, table=False):
    """Octahedral distortion analysis via cosmochlore, on the current selection.

    Requires a single selected atom with exactly six neighbours (a 7-atom
    polyhedron once its neighbours are added).
    """
    print('\n' + '-' * 50)
    print('Octahedral distortion analysis using cosmochlore')
    try:
        exe = check_cosmochlore(_cosmochlore_exe_path())
    except CosmochloreError as e:
        print(e)
        return False

    sel_string = olex.f('sel()')
    if sel_string == '':
        print('Invalid atom selection: no atoms selected.')
        return False

    selection = AtomSelection(sel_string)
    if len(selection) != 1:
        print(f'Invalid atom selection: expected 1 atom, found {len(selection)}.')
        return False

    selection.add_neighbours()

    workdir = _cosmochlore_workdir()
    ran_any = False
    for i, structure in enumerate(split_by_parts(selection)):
        if len(structure) != 7:
            print(f'Skipping part {i}: expected 7 atoms (centre + six donors), '
                  f'found {len(structure)}.')
            continue

        xyz_path = os.path.join(workdir, f'{olx.FileName()}_{structure.labels[0]}_{i}.xyz')
        write_xyz(structure, xyz_path, comment=f'{olx.FileName()} part {i}')

        try:
            args = build_odis_args(xyz_path, full, table)
            print(trim_banner(run_cosmochlore(exe, args, cwd=workdir)))
            ran_any = True
        except CosmochloreError as e:
            print(f'odis failed for part {i}: {e}')

    if ran_any:
        print(_COSMOCHLORE_CREDIT)
    return ran_any


def compare_odis():
    """
    Cross-checks the Python OctaDist reimplementation against cosmochlore's
    odis on the current selection, printing both sets of results side by side.
    Meant for validating the two engines agree, not for routine use.
    """
    sel_string = olex.f('sel()')
    if sel_string == '':
        print('Invalid atom selection: no atoms selected.')
        return False

    selection = AtomSelection(sel_string)
    if len(selection) != 1:
        print(f'Invalid atom selection: expected 1 atom, found {len(selection)}.')
        return False
    selection.add_neighbours()

    try:
        exe = check_cosmochlore(_cosmochlore_exe_path())
    except CosmochloreError as e:
        print(f'cosmochlore unavailable, showing Python results only: {e}')
        exe = None

    workdir = _cosmochlore_workdir()
    # `key` names the field cosmochlore's odis CSV uses (see cosmochlore.read_odis_csv);
    # CalcDistortion names its mean-bond-distance attribute differently, hence the alias.
    rows = [('d_mean', 'mean_bond_distance', 'Mean d(M-X)', 'Ang'),
            ('zeta', 'zeta', 'Zeta', 'Ang'), ('delta', 'delta', 'Delta', ''),
            ('sigma', 'sigma', 'Sigma', 'deg'), ('theta', 'theta', 'Theta', 'deg'),
            ('tau', 'tau', 'Tau', 'deg'), ('mu', 'mu', 'Mu', 'Ang')]

    for i, structure in enumerate(split_by_parts(selection)):
        if len(structure) != 7:
            print(f'Skipping part {i}: expected 7 atoms, found {len(structure)}.')
            continue

        try:
            python_result = CalcDistortion(structure)
        except ValueError as e:
            print(f'Python OctaDist failed for part {i}: {e}')
            python_result = None

        cosmo_values = {}
        if exe:
            xyz_path = os.path.join(workdir, f'{olx.FileName()}_{structure.labels[0]}_{i}_cmp.xyz')
            write_xyz(structure, xyz_path, comment=f'{olx.FileName()} part {i}')
            try:
                run_cosmochlore(exe, build_odis_args(xyz_path, full=False, table=True), cwd=workdir)
                cosmo_values = read_odis_csv(xyz_path)
            except (CosmochloreError, OSError) as e:
                print(f'cosmochlore odis failed for part {i}: {e}')

        print(f'\n=== {structure.labels[0]}, part {i} ===')
        print(f'{"":<14}{"Python":>12}{"cosmochlore":>14}   {"":<5}')
        for key, py_attr, label, unit in rows:
            py_val = getattr(python_result, py_attr) if python_result is not None else float('nan')
            co_val = cosmo_values.get(key, float('nan'))
            print(f'{label:<14}{py_val:>12.4f}{co_val:>14.4f}   {unit:<5}')

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

        # cosmochlore entry points (Phase 1: console-only, no GUI yet).
        OV.registerFunction(autoCSHM, True, "SymmetryMeasurements")
        OV.registerFunction(autoCSOM, True, "SymmetryMeasurements")
        OV.registerFunction(autoODIS, True, "SymmetryMeasurements")
        OV.registerFunction(compare_odis, True, "SymmetryMeasurements")
        OV.registerFunction(can_find_cosmochlore_msg, False, "SymmetryMeasurements")

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
