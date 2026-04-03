GTK4 app workspace

This is the new GTK4 rebuild target for GDIS.

Current state:
- A real GTK4 application scaffold exists here.
- The scaffold mirrors the legacy Linux single-window layout:
  - menu bar
  - toolbar
  - left model/property pane
  - right viewer pane
  - bottom status/log pane
- Legacy reference files remain available in `../legacy_snapshot/`.
- The app now builds locally on this Mac with GTK4 installed.
- `Open` uses a real GTK4 file dialog and adds the chosen file to the model list.
- `Save` now writes the active model back out using the loader/writer bridge.
- `Save As` now opens a real GTK4 save dialog even for already-opened files.
- `Close` now removes the active model from the session and selects a sensible neighbor.
- `Edit > New Model` now creates an empty model and opens the model editor immediately.
- Startup file arguments are accepted, so shell-expanded paths like `models/*` are no longer ignored.
- Legacy-derived model loading now works for:
  - `XYZ`
  - `PDB`
  - `ARC/CAR`
  - `CIF`
- Legacy-derived `ARC/CAR` handling now also preserves:
  - old BIOSYM-style atom records where field 5 is not literally `CORE`
  - `PBC=2D` headers
  - short 2D `PBC a b gamma` records
  - title lines before `!DATE`
- The viewer supports:
  - rotate by dragging
  - zoom by scroll wheel / trackpad
  - atom picking by click
  - atom / bond / cell / label toggles
- Most of the practical legacy selection modes in the sidebar are now live:
  - `Atoms`
  - `Atom Label`
  - `Atom FF Type`
  - `Elements`
  - `Elements in Molecule`
  - `Molecules`
  - `Molecule Fragments` via a 2-pick bonded path selection
- `Regions` is still shown as legacy reference, but remains disabled until region-labelled model data is ported into the GTK4 bridge.
- The viewer now honours the model image-range state when drawing periodic crystal images.
- Tool actions with real GTK4-native output now include:
  - `Measure`
  - `Surface`
  - `Diffraction`
- `Edit` now opens a model editor window with:
  - apply changes to the selected atom
  - add atom from the current fields
  - delete selected atom
  - delete the current selected set
  - delete picked atoms
  - add bond between the last 2 picked atoms
  - remove bond between the last 2 picked atoms
  - pick-driven bond add/remove modes from the viewer
  - clear pick history
- The `Images` page now has working crystal controls for:
  - periodic image ranges
  - `Confine To Cell`
  - `Force P1`
  - `Bake Supercell`
- Empty models now show an in-view hint telling you to add atoms from the editor and then save.
- `Iso-surfaces` now opens a concrete status window that explains the remaining engine work instead of doing nothing.

Local prerequisites:
- `pkg-config`
- GTK4 development files
- a C compiler

On macOS with Homebrew, the expected install is:
- `brew install gtk4 pkgconf`

Build:
```sh
make
```

Run:
```sh
make run
```

Run with startup files:
```sh
make run ARGS="../../models/deoxy.pdb"
make run ARGS="../../examples/water.xyz ../../models/deoxy.pdb"
./build/gdis-gtk4 ../../models/water.car
./build/gdis-gtk4 ../../models/deoxy.pdb
./build/gdis-gtk4 ../../models/*
```

How to use the current tools:
- `Edit > New Model` creates an empty working model and opens the model editor.
- Use `File > Save` to save the active model.
- Use `File > Save As` to export the current model to a new path/format.
- Use `File > Close` to remove the active model from the current session.
- Supported write targets currently follow the filename extension:
  - `.xyz`
  - `.pdb`
  - `.arc`
  - `.car`
  - `.cif`
- Click atoms in the viewer to build a pick history.
- `Tools > Measure` now opens a measurement tool window with explicit modes:
  - `Auto`
  - `Distance`
  - `Angle`
  - `Torsion`
- The main menu now restores more of the legacy interaction surface:
  - `File > New`
  - `Edit > Delete Selected`
  - `Edit > Select All`
  - `Edit > Invert Selection`
  - `View > Display Properties`
  - `View > Reset Model Images`
  - `Help > Manual`
- `Tools > Edit` now opens the model editor window.
- In `Edit`, the fields follow the currently selected atom in the viewer.
- The viewer selection mode buttons now change what a click selects.
- `Molecule Fragments` works as a 2-click path selection:
  - first click sets the fragment anchor
  - second click selects the bonded path to the new atom
- You can use the editor to:
  - change label
  - change element
  - change coordinates
  - add a new atom
  - delete a whole selected set
  - create/remove a bond from the last 2 picks
- `Pick Add Bond` and `Pick Remove Bond` restore a more legacy-like workflow:
  - click the mode button
  - pick 2 atoms in the viewer
  - the bond edit runs automatically
- Measurement results use the order you picked the atoms.
- Periodic measurement output now uses minimum-image geometry for periodic cells.
- The `Images` page now lets you:
  - change the stored image range for periodic models
  - confine atoms back into the active cell
  - force the crystal metadata to `P 1`
  - bake a true supercell into the loaded model
- `Tools > Surface` shows low-index planes and d-spacings for periodic models.
- `Tools > Diffraction` shows a lightweight powder X-ray preview for 3D periodic models.
- `Tools > Render` now opens a GTK4 display-properties window with viewer toggles, view presets, and image reset.
- A detailed legacy-vs-GTK4 parity matrix now lives in `../RESTORATION_AUDIT.md`.

Helpful checks:
```sh
make doctor
```

Notes:
- The current rebuild has moved beyond a shell:
  - file/model bridge is active
  - the viewer is active
  - several tool actions now produce real output
- The next deeper ports are the old full-feature dialogs and engines, especially:
  - region-labelled selection / region editing
  - slab/surface construction
  - full diffraction plotting
  - volumetric-grid loading for iso-surfaces
- Legacy file/model port targets are mapped in `PORTING.md`.
