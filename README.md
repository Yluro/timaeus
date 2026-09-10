# Timaeus

Timaeus is an [Olex2](https://www.olexsys.org/olex2/)$^1$ plugin that integrates several symmetry and shape analysis tools into a single graphical interface within Olex2.
### Features
- All analysis methods in Timaeus will read the atomic coordinates directly from the `OlexRefinementModel` — no external input files needed.
- Olex2 GUI panel.
- SHAPE 2.1 wrapper and output parser.
- Octahedral distortion parameters calculation.
- [cosmochlore](https://github.com/Yluro/cosmochlore) integration: Continuous Shape Measures (CShM), Continuous Symmetry Operation Measures (CSoM) and octahedral distortion (ODis), with user-defined reference shapes.
- Centroid merging for pi-bonded ligands
- Automatic disorder handling. If more than two parts are found in the selected strucutre, calculations will be run on each part separately. 


## Requirements
- Olex2 1.5.
- SHAPE 2.1 executable available on your system `PATH` (optional, needed for the SHAPE panel).
  - Download SHAPE 2.1 from the [Electronic Strucutre Group's webpage](https://www.ee.ub.edu/downloads/)
- [cosmochlore](https://github.com/Yluro/cosmochlore) 1.0.2 or newer on your `PATH` (optional, needed for the Cosmochlore panel).
  - The path can also be set explicitly in the plugin settings (`Extras` > `Settings`).

The plugin was developed/tested using a Windows 10/11 machine. The plugin should be system agnostic but please report any bugs found in any other operating systems. 


_**Note:** It is known that SHAPE 2.1 gives trouble in Mac machines with operating systems newer than 2022 and some Linux systems. Unfortunately, I can't do anything about that since the ESG hasn't published a precompiled SHAPE version since 2013. In the future I might introduce support for Cosymlib by ESG which is the updated version of their shape and symmetry measures program._

## Instalation
1. Download the source code from the lastest GitHub release.
2.  Go to the base directory of your Olex2 instalation. You can open the directory by typing `shell BaseDir()` in the Olex2 console.
3. Create a file called `plugins.xld` and write the following contents inside `plugins.xld`:
```xml
<Plugin
 <Timaeus>
>
```
4. Go to the `<BaseDir>\util\pyUtil\PluginLib\` folder and create a folder called `plugin-Timaeus`.
5. Extract the downloaded zip into that folder.
6. On restarting Olex2, a Timaeus panel should appear under the Tools tab.

_**Note:** I push development changes constantly to the master branch. You could git clone this repository to automatically keep the plugin updated. Usually, if changes are pushed it means that the plugin is in a usable state. But it does not guarantee that things won't break._

## AutoSHAPE
[SHAPE 2.1](https://www.ee.ub.edu/continuous-shape-and-symmetry-measures/)$^2$ is a software published by ESG used to calculate Continuous Shape Measures (CShM's). autoSHAPE are a collection of personal Python scripts I developed and used to run and parse SHAPE i/o files. Timaeus contains an implementation of autoSHAPE to:
- Generate the necessary `.dat` input files for SHAPE automatically.
- Run SHAPE and parse the resulting `.tab` output.
- Output a summary table from the `.out` and `.tab` files.
SM's autoSHAPE does not overwrite previous runs as it stores each run in a dedicated folder: `<FilePath>\autoSHAPE\<FileName>_<part>_<atoms>\<run>`.

### Usage

1. Open a structure in Olex2.
2. Select the atoms you want to include in the measurement.
 - If one atom is selected, the neighbouring atoms will be taken into account to form a centered shape (i.e. a coordination structure).
 - If multiple atoms are selected, the will be interpreted as a non centered shape (i.e. a borane cluster).
4. Run `spy.Timaeus.autoSHAPE()` from the Olex2 console or from the Tools/Timaeus panel.
5. Results are printed to the console and saved in `<FilePath>/autoSHAPE/`.

#### _New in version 0.2_
There is an option to merge pi-bonded ligands into a centroid. If checked, autoshape will interpret pi-bonded ligands as the average of the fragments as per Cirera _et al_$^2$ paper. 

## Octahedral Distortion Parameters.
Timaeus includes a reimplementation of the [OctaDist](https://octadist.github.io/)$^3$ algorithm. Unlike the original implementation, this version identifies opposite faces and vertices of an octahedron using topological criteria. It relies on constructing a convex hull that is topologically equivalent to an octahedron. As a result, the algorithm may fail when more than three vertices are coplanar, causing the convex hull to degenerate into a different polyhedral shape.

Several features of the original OctaDist program have been omitted to simplify integration with Olex2 and to remove redundant functionality.

The Octahedral Distortion Parameters module computes all distortion parameters reported by OctaDist, together with an additional $\tau$ and $\mu$ parameters that quantify the deviation of trans angles from the ideal 180 $^\circ$ and the deviation of the metal centre to the centroid of the octahedron, respectively. The calculated values are generally consistent with those produced by OctaDist, although discrepancies may occur for highly distorted octahedra or for structures approaching an ideal trigonal prism.


- Bond length distortion:
```math
\zeta = \sum_{i=1}^6 |d_i - d_{mean}|
```
where $d_i$ are the M-X bond distances and $d_{mean}$ is the mean M-X bond distance.
- Octahedral tilting parameter:
```math
\Delta = \frac{1}{6}\sum_{i=1}^6 \left(\frac{d_i - d_{mean}}{d_{mean}}\right)^2
```
where $d_i$ are the M-X bond distances and $d_{mean}$ is the mean M-X bond distance.
- Cis angle distortion:
```math
\Sigma = \sum_{i=1}^{12} |90^\circ - \phi_i|
```
where $\phi_i$ are the cis angles.
- Octahedral twisting distortion:
```math
\Theta = \sum_{i=1}^{24}  |60^\circ - \theta_i|
```
where $\theta_i$ are twisting angles between vectors of two opposite faces.

#### Extra parameters:
- Trans angle distortion:
```math
\tau = \sum_{i=1}^3 |180^\circ - \psi_i|
```
where $\psi_i$ are the trans angles.
- Deviation of the metal from the centroid:
```math
\vec{c} = \frac{1}{6}\sum^6_{i=1} \vec{r}_i
```
```math
\mu = |\vec{r}_M - \vec{c}|
```
where $\vec{r}_M$ is the position of the metal and $\vec{r}_i$ are the position of the donor atoms.



### Usage
1. Open a structure in Olex2.
2. Select the central atom of a 6-coordinate complex.
3. Run `spy.Timaeus.autoOCTADIST()` from the Olex2 console or from the Tools/Timaeus panel.
4. Results are printed in the console. A graph will saved in `<FilePath>/OH_distortion` showing the extracted octahedron. 


## Cosmochlore
[cosmochlore](https://github.com/Yluro/cosmochlore) is a separate Rust program that calculates shape and symmetry measures. Timaeus wraps its three subcommands; each writes its output next to the structure and prints the results table to the console. Every option below has a phil parameter, editable in `Extras` > `Settings` or from the panel itself.

| Function | Description |
| --- | --- |
| `spy.Timaeus.autoCSHM()` | Continuous Shape Measures against the built-in and user-defined reference polyhedra. Same selection rule as autoSHAPE. |
| `spy.Timaeus.autoCSOM()` | Continuous Symmetry Operation Measures against a list of point groups (space-separated Schoenflies symbols, e.g. `Oh D4h D3d`). Slower than CShM/ODis. |
| `spy.Timaeus.autoODIS()` | Octahedral distortion parameters, independent of the OctaDist reimplementation above. Select the central atom of a 6-coordinate complex. |

### User-defined shapes
Reference shapes beyond the built-in 90 can be added as `.yaml` files in the `user_shapes/` folder of the plugin (`Open user defined shapes folder` in the CShM panel). Each file gets a checkbox; only the checked ones are used. A mismatched vertex count is reported by cosmochlore itself.

## Known limitations/upcoming features.

| Status  | Features                                                                                                                 |
| ------------- |--------------------------------------------------------------------------------------------------------------------------|
| ✔️ | Suppport for disordered structures.                                                                                      |
| ✔️ | Centorid search for pi-bonded ligands. *Doesn't work with ligands outside ASU.                                           |
| ✔️ | Custom reference shapes. *Via cosmochlore, see `user_shapes/`.                                                           | 
| ✔️ | Smarter program logic (automatic coordination site detection, multiple selections, etc.).                                | 
| ✔️ | Support for other measurement programs (cosmochlore: CShM, CSoM, ODis).                                                  |
| ✔️ | Reimplementation of octahedral distortion parameters (Zeta, Sigma, Theta) — relevant for spin-crossover (SCO) complexes. |
| ✔️ | Non centered shapes.                                                                                                     |
| ✔️ | HTML UI.                                                                                                                 |

## License
- Not yet.

## Citations

1. _J. Appl. Cryst._ (2009). 42, 339–341. DOI: https://doi.org/10.1107/S0021889808042726
2. _Organometallics_ (2005), 24, **7**, 1556–1562. DOI: https://doi.org/10.1021/om049150z
3. _Dalton Trans._ (2021) 50, **3**, 1086–1096. DOI: https://doi.org/10.1039/d0dt03988h
