import os
import re

import olex
import olx
from olexFunctions import OlexFunctions

OV = OlexFunctions()

from reload_all import reload_all

from autoshape import (ShapeCalculation, can_find_shape_msg, find_shape,
                       print_shape_table, run_shape)
from constants import octadist_citation, shape21_citation
import cosmochlore
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

OV.SetVar('Timaeus_plugin_path', p_path)

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
            merge = as_bool(OV.GetParam('timaeus.merge_ligands', False))
        if merge:
            selection.merge_ligands()
        centered = True

    # Returns one structure per disorder component, or a single structure when the strucutre is not disordered
    return split_by_parts(selection), centered


def _cosmochlore_exe_path():
    return OV.GetParam('timaeus.cosmochlore.exe_path', '') or None


def _cosmochlore_workdir():
    """The scratch folder cosmochlore .xyz inputs and outputs are written to."""
    path = os.path.join(olx.FilePath(), 'cosmochlore')
    os.makedirs(path, exist_ok=True)
    return path


def _user_shapes_dir():
    """Where autoCSHM()'s -r/--ref candidates are read from: a `user_shapes`
    folder living next to the plugin's own code, alongside SHAPE's `custom_shapes`."""
    return os.path.join(p_path, 'user_shapes')


def _safe_key(name: str) -> str:
    """Sanitises a shape's display name into a safe HTML control name / phil
    param suffix - the display name itself (a .yaml filename the user chose)
    isn't guaranteed to be free of spaces, quotes, or other characters that
    would break the generated markup or a dotted phil path."""
    return re.sub(r'[^A-Za-z0-9_-]', '_', name)


def _user_shape_param(name: str) -> str:
    return f'timaeus.cosmochlore.cshm.user_shape.{_safe_key(name)}'


def user_shapes_checkboxes_html():
    """Returns one checkbox per .yaml file found in the user_shapes folder, for
    picking which ones autoCSHM() passes to cshm's -r/--ref. Selection state is
    stored per file (see _user_shape_param) rather than as one phil-declared
    list, since the set of files - and so the set of possible param names - is
    only known at runtime; this means it is not phil-persisted across Olex2
    restarts, only for the current session.

    Mirrors the table/tr/td structure gui/snippets/input-checkbox-td expands
    to: the checkbox's own `label` attribute is not what renders visible text
    in Olex2's control - the label is a separate <td><b>...</b></td> cell next
    to it.
    """
    shapes = cosmochlore.list_user_shapes(_user_shapes_dir())
    if not shapes:
        return '<i>No .yaml files found.</i>'

    rows = []
    for name, _ in shapes:
        param = _user_shape_param(name)
        checked = 'true' if as_bool(OV.GetParam(param, False)) else 'false'
        label = name.replace('"', "'")
        rows.append(
            f'<tr><td><input type="checkbox" name="UserShape_{_safe_key(name)}" '
            f'checked="{checked}" '
            f'onclick="spy.SetParam(\'{param}\', html.GetState(\'~name~\'))">'
            f'</td><td align="left"><b>{label}</b></td></tr>'
        )
    return '<table cellpadding="0" cellspacing="0">' + ''.join(rows) + '</table>'


def _selected_user_shapes():
    """Full paths of the user_shapes .yaml files currently checked in the GUI."""
    return [path for name, path in cosmochlore.list_user_shapes(_user_shapes_dir())
           if as_bool(OV.GetParam(_user_shape_param(name), False))]


def open_user_shapes_folder():
    """Opens the user_shapes folder in the OS file browser. Creates it first
    if it doesn't exist yet, so the button always has somewhere to open."""
    path = _user_shapes_dir()
    os.makedirs(path, exist_ok=True)
    olx.Shell(path)


def reload_plugin():
    """Re-imports every plugin module from wherever the plugin is installed,
    so the Extras button doesn't depend on a hardcoded install path."""
    reload_all()


def open_plugin_folder():
    """Opens the plugin's own folder in the OS file browser."""
    olx.Shell(p_path)


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


def autoCSHM(shapes=None, user_shapes=None, table=None, ideal=None):
    """Continuous Shape Measures via cosmochlore, on the current selection.

    `shapes`: built-in reference-shape indices to restrict to (None = all
    applicable for the detected vertex count). `user_shapes`: paths to
    user-defined shape .yaml files; when not given, falls back to whichever
    files are checked in the user_shapes_checkboxes_html() GUI list. Passed to
    cshm's -r/--ref exactly as given - no vertex-count pre-filtering, cosmochlore
    reports a mismatch itself (aborting that part's run) rather than this
    having any fallback logic of its own to hide the error.
    `table`/`ideal`: also write the corresponding cosmochlore output files
    next to the .xyz; default to the matching phil params when not given.
    """
    print('\n' + '-' * 50)
    print('Continuous Shape Measures using cosmochlore')
    try:
        exe = cosmochlore.check_cosmochlore(_cosmochlore_exe_path())
    except cosmochlore.CosmochloreError as e:
        print(e)
        return False

    if user_shapes is None:
        user_shapes = _selected_user_shapes() or None
    if table is None:
        table = as_bool(OV.GetParam('timaeus.cosmochlore.cshm.table', False))
    if ideal is None:
        ideal = as_bool(OV.GetParam('timaeus.cosmochlore.cshm.ideal', False))

    structures, centered = _prepare_structures(olex.f('sel()'))
    if structures is None:
        return False

    workdir = _cosmochlore_workdir()
    ran_any = False
    for i, structure in enumerate(structures):
        xyz_path = os.path.join(workdir, f'{olx.FileName()}_{structure.labels[0]}_{i}.xyz')
        cosmochlore.write_xyz(structure, xyz_path, comment=f'{olx.FileName()} part {i}')

        try:
            args = cosmochlore.build_cshm_args(xyz_path, centered, shapes, user_shapes, table, ideal)
            print(cosmochlore.trim_banner(cosmochlore.run_cosmochlore(exe, args, cwd=workdir)))
            ran_any = True
        except cosmochlore.CosmochloreError as e:
            print(f'cshm failed for part {i}: {e}')

    if ran_any:
        print(_COSMOCHLORE_CREDIT)
    return ran_any


def autoCSOM(point_groups=None, mode=None, vector=None, full=None, table=None,
            operated=None, samples=None, iterations=None, ignore_labels=None):
    """Continuous Symmetry Operation Measures via cosmochlore, on the current selection.

    point_groups: a space-separated string or list of Schoenflies, falls back to the
    timaeus.cosmochlore.csom.point_groups phil param when not given.

    mode: centering mode (auto/first/centroid/manual), defaults to the matching phil param.

    `vector`: required 3-value centering vector when mode is 'manual'.

    `full`/`table`/`operated`/`ignore_labels`: default to the matching phil params.
    """
    print('\n' + '-' * 50)
    print('Continuous Symmetry Operation Measures using cosmochlore')
    try:
        exe = cosmochlore.check_cosmochlore(_cosmochlore_exe_path())
    except cosmochlore.CosmochloreError as e:
        print(e)
        return False

    if point_groups is None:
        point_groups = OV.GetParam('timaeus.cosmochlore.csom.point_groups', '')
    if isinstance(point_groups, str):
        point_groups = point_groups.split()

    if not point_groups:
        print('No point groups given. Pass e.g. '
              "spy.Timaeus.autoCSOM('Oh D4h D3d'), or set "
              "timaeus.cosmochlore.csom.point_groups.")
        return False

    if mode is None:
        mode = OV.GetParam('timaeus.cosmochlore.csom.mode', 'auto')
    # Olex2's combo control capitalises the value it hands back (e.g. 'auto'
    # -> 'Auto') regardless of the case used in the combo's own item list, so
    # normalise here rather than trust whatever case arrives from the GUI or a
    # console caller.
    mode = str(mode).strip().lower()

    if vector is None:
        vector_str = OV.GetParam('timaeus.cosmochlore.csom.vector', '')
        if vector_str.strip():
            try:
                vector = [float(v) for v in vector_str.split()]
            except ValueError:
                print(f'Could not parse the manual centering vector "{vector_str}": '
                      f'expected three numbers separated by spaces, e.g. "0.0 0.0 0.0".')
                return False
            if len(vector) != 3:
                print(f'The manual centering vector needs exactly 3 numbers, '
                      f'found {len(vector)} in "{vector_str}".')
                return False

    if mode == 'manual' and not vector:
        print('Centering mode is "manual" but no centering vector was given. Set it in '
              'the Cosmochlore section (x y z, space-separated), or pass '
              "spy.Timaeus.autoCSOM(vector=[x, y, z]).")
        return False

    if full is None:
        full = as_bool(OV.GetParam('timaeus.cosmochlore.csom.full', False))
    if table is None:
        table = as_bool(OV.GetParam('timaeus.cosmochlore.csom.table', False))
    if operated is None:
        operated = as_bool(OV.GetParam('timaeus.cosmochlore.csom.operated', False))
    if ignore_labels is None:
        ignore_labels = as_bool(OV.GetParam('timaeus.cosmochlore.csom.ignore_labels', False))

    structures, centered = _prepare_structures(olex.f('sel()'))
    if structures is None:
        return False

    print('This may take a while depending on the number of point groups and atoms...')

    workdir = _cosmochlore_workdir()
    ran_any = False
    for i, structure in enumerate(structures):
        xyz_path = os.path.join(workdir, f'{olx.FileName()}_{structure.labels[0]}_{i}.xyz')
        cosmochlore.write_xyz(structure, xyz_path, comment=f'{olx.FileName()} part {i}')

        try:
            args = cosmochlore.build_csom_args(xyz_path, centered, point_groups, mode, vector, full,
                                   table, operated, samples, iterations, ignore_labels)
            print(cosmochlore.trim_banner(cosmochlore.run_cosmochlore(exe, args, cwd=workdir)))
            ran_any = True
        except cosmochlore.CosmochloreError as e:
            print(f'csom failed for part {i}: {e}')

    if ran_any:
        print(_COSMOCHLORE_CREDIT)
    return ran_any


def autoODIS(full=None, table=None):
    """Octahedral distortion analysis via cosmochlore, on the current selection.

    Requires a single selected atom with exactly six neighbours (a 7-atom
    polyhedron once its neighbours are added). `full`/`table` default to the
    matching phil params when not given.
    """
    print('\n' + '-' * 50)
    print('Octahedral distortion analysis using cosmochlore')
    try:
        exe = cosmochlore.check_cosmochlore(_cosmochlore_exe_path())
    except cosmochlore.CosmochloreError as e:
        print(e)
        return False

    if full is None:
        full = as_bool(OV.GetParam('timaeus.cosmochlore.odis.full', False))
    if table is None:
        table = as_bool(OV.GetParam('timaeus.cosmochlore.odis.table', False))

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
        cosmochlore.write_xyz(structure, xyz_path, comment=f'{olx.FileName()} part {i}')

        try:
            args = cosmochlore.build_odis_args(xyz_path, full, table)
            print(cosmochlore.trim_banner(cosmochlore.run_cosmochlore(exe, args, cwd=workdir)))
            ran_any = True
        except cosmochlore.CosmochloreError as e:
            print(f'odis failed for part {i}: {e}')

    if ran_any:
        print(_COSMOCHLORE_CREDIT)
    return ran_any


def shape_status_html():
    where = find_shape()
    found = where is not None
    color = OV.GetParam('gui.green') if found else OV.GetParam('gui.grey')
    text = f'SHAPE executable found at: {where}' if found else 'Unable to find shape.exe in the system path.'
    return f"<font color='{color}'>{text}</font>"


def cosmochlore_status_html():
    try:
        exe = cosmochlore.check_cosmochlore(_cosmochlore_exe_path())
    except cosmochlore.CosmochloreError as e:
        return f"<font color='{OV.GetParam('gui.grey')}'>{e}</font>"

    version = '.'.join(str(v) for v in cosmochlore.get_version(exe))
    text = f'cosmochlore {version} found at: {exe}'
    return f"<font color='{OV.GetParam('gui.green')}'>{text}</font>"


class Timaeus(PT):
    def __init__(self):
        super(Timaeus, self).__init__()
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
        OV.registerFunction(autoSHAPE, True, "Timaeus")
        OV.registerFunction(autoOCTADIST, True, "Timaeus")
        OV.registerFunction(can_find_shape_msg, True, "Timaeus")
        OV.registerFunction(shape_status_html, False, 'Timaeus')

        # cosmochlore entry points.
        OV.registerFunction(autoCSHM, True, "Timaeus")
        OV.registerFunction(autoCSOM, True, "Timaeus")
        OV.registerFunction(autoODIS, True, "Timaeus")
        OV.registerFunction(cosmochlore.can_find_cosmochlore_msg, False, "Timaeus")
        OV.registerFunction(cosmochlore_status_html, False, "Timaeus")
        OV.registerFunction(user_shapes_checkboxes_html, False, "Timaeus")
        OV.registerFunction(open_user_shapes_folder, True, "Timaeus")

        # Extras panel.
        OV.registerFunction(reload_plugin, True, "Timaeus")
        OV.registerFunction(open_plugin_folder, True, "Timaeus")

        # Debug panel helpers.
        OV.registerFunction(get_selected_atoms, True, "Timaeus")
        OV.registerFunction(get_neighbours, True, "Timaeus")
        OV.registerFunction(get_neighbours_on_sel, True, "Timaeus")
        OV.registerFunction(get_xyz_sel, True, "Timaeus")
        OV.registerFunction(print_console_bs, False, 'Timaeus')
        OV.registerFunction(print_orm, False, 'Timaeus')
        OV.registerFunction(test_selection_class, False, 'Timaeus')
    # END Generated =======================================


Timaeus_instance = Timaeus()
print("Loading Timaeus modules.")
reload_all()
print("Timaeus by JSG loaded.")
