#Timaeus Plugin
Shape and symmetry analysis: a SHAPE 2.1 wrapper, an octahedral distortion calculator, and a cosmochlore integration (CShM, CSoM, ODis). URL[https://github.com/Yluro/cosmochlore]

#SHAPE 2.1
SHAPE 2.1 wrapper: Continuous Shape Measures against the built-in reference polyhedra.

#SHAPE Status
Whether **shape.exe**/**shape_2.1.bat** (SHAPE 2.1) was found on **PATH**.

#SHAPE 1
##SHAPE analysis
Runs SHAPE 2.1 on the current selection. One atom selected = centered shape (its neighbours are added); several atoms = non-centered shape. Results print to the console and save under **autoSHAPE/**.

#OctaDist
This plugin's own reimplementation of the OctaDist distortion parameters.

#OctaDist 1
##Octahedral distortion
Runs this plugin's own OctaDist reimplementation. Select the central atom of a 6-coordinate complex; its neighbours are added automatically. Saves a graph under **Oh_distortion/**.

#Options
General plugin-wide options.

#Options 1
##Merge pi-bonded ligands
When checked, pi-bonded ligand fragments are merged into one centroid before measuring (Cirera *et al.*, *Organometallics* 2005, 24, 1556).

#Timaeus Extras
Development and maintenance tools for this plugin.

#Extras 1
##Reload ALL
Reloads the plugin's Python code.

##Settings
Opens this plugin's phil parameter editor.

##Open Folder
Opens the plugin's install folder in the file browser.

#Timaeus Debug
Low-level diagnostics for the selection/refinement-model machinery.

#Debug 1
##Sel
Prints the current selection string.
##Neighbours
Prints the bonded neighbours of the selected atoms.
##Selection Class
Builds and splits an **AtomSelection**, printing each step.

#Debug 2
##XYZ-SEL
Prints the coordinates of a selected atom.
##ORM
Writes the refinement model's atom list to **orm.txt**.
##SHAPE
Prints whether **shape** was found on **PATH**.

#Cosmochlore
Shape and symmetry measures via the cosmochlore engine: CShM, CSoM, ODis. URL[https://github.com/Yluro/cosmochlore]

#Cosmochlore Status
Whether a **cosmochlore** executable (1.0.2+) was found. Every button below needs this to be green.

#CShM
Continuous Shape Measures: compares the selection to built-in and user-defined reference polyhedra. Same selection rule as SHAPE analysis.

#CShM 1
##Run CShM
Runs the measure and prints the results table.

#CShM 2
Options for this cshm run.

#CShM 3
##Results CSV
Also saves a **_cshm_table.csv** file.
##Ideal Polyhedra XYZ
Also saves a **_ideal.xyz** file with the aligned reference shapes.

#CShM 4
Custom reference shapes beyond the built-in 90.

#CShM 5
One checkbox per **.yaml** file in **user_shapes/**. Only checked ones are used; a mismatched vertex count is reported by cosmochlore itself.

#CShM 6
##Open user defined shapes folder
Opens the **user_shapes** folder, creating it if needed.

#CSoM
Continuous Symmetry Operation Measures: how well the selection matches given point groups. Slower than CShM/ODis, more so with more point groups
or atoms.

#CSoM 1
##Run CSoM
Runs the measure for every listed point group.

#CSoM 2
Options for this csom run.

#CSoM 3
##Point groups
Space-separated, case-sensitive Schoenflies symbols (e.g. **Oh D4h D3d**).
At least one is required.

#CSoM 4
##Centering mode
**auto** centers on the first atom (centered selection) or centroid (non-centered); 

**first**/**centroid** always use one or the other; 

**manual** uses the vector on the right.
##Manual vector (x y z)
Used only when the centering mode is set to manual, e.g. **5.1 3.3 -4.9**. In orthogonalised coordinates in Angstroms. 

#CSoM 5
##Summary CSV
Also saves a **_csom_table.csv** summary file.
##Per-operation details CSV
Also saves a per-operation **_details.csv** file, per point group.
##Operated .xyz/.mol2
Also saves the structure with each symmetry operation applied.

#CSoM 6
##Ignore atom labels
Match atoms by position only, ignoring element/label.

#ODis
Octahedral Distortion Parameters: a second, independent octahedral distortion implementation, via cosmochlore. Select the central atom of a 6-coordinate complex.

#ODis 1
##Run ODis
Runs the analysis and prints the results table.

#ODis 2
Options for this odis run.

#ODis 3
##Results CSV
Also saves a **_odis_table.csv** file.
##Full analysis
Also computes CShM (vs. octahedron/trigonal prism) and CSoM (vs. common distortions).
