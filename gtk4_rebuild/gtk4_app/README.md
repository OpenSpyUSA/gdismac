# GDIS GTK4 Rebuild

This directory contains the GTK4-based rebuild of GDIS for modern macOS and Linux desktops.

The goal of this rebuild is practical usability first:

- keep the legacy GDIS model-loading and chemistry workflows available
- restore the familiar single main window layout
- provide a native GTK4 application that builds cleanly on current systems
- make the macOS experience workable without requiring the old Linux-only stack

Legacy reference sources remain available in [../legacy_snapshot](../legacy_snapshot), and the current restoration status is tracked in [../RESTORATION_AUDIT.md](../RESTORATION_AUDIT.md).

## Current State

The GTK4 rebuild is no longer just a scaffold. It now includes:

- a real application window with menu bar, toolbar, sidebar, viewer, and status area
- native file open/save/save-as flows
- startup-file handling from the CLI
- restored structure loading for common legacy formats
- a working 3D viewer with rotate, zoom, picking, and display toggles
- native GTK4 tools for editing, measurement, diffraction, surface analysis, periodic table, animation, recording, executable paths, task management, and Qbox deck setup

On macOS, the rebuild also includes a native Qbox menu bridge so Qbox actions are reachable from the system menu bar.

## What Works

### File and Session

- `File > New` creates an empty model and opens the editor
- `File > Open` uses the native GTK4 file dialog
- `File > Save` and `File > Save As` write through the GTK4 model bridge
- `File > Close` removes the active model from the session
- shell startup paths are accepted, including multiple files

### Supported Model Formats

Current loader coverage includes:

- `XYZ`
- `PDB`
- `ARC`
- `CAR`
- `CIF`

Current write targets follow the filename extension:

- `.xyz`
- `.pdb`
- `.arc`
- `.car`
- `.cif`

### Viewer

The viewer currently supports:

- rotate by drag
- zoom by scroll or trackpad gesture
- atom picking by click
- selection-mode switching from the sidebar
- atom, bond, label, and cell display toggles
- periodic image-range drawing

### Editing and Analysis

The rebuild already provides practical native workflows for:

- atom label, element, and coordinate editing
- atom add/delete
- delete selected atoms
- bond add/remove
- pick-driven bond editing
- undo for model edits
- distance, angle, and torsion measurement
- periodic image controls
- low-index surface inspection
- lightweight powder diffraction setup and plot flow

### Computation Tools

`Tools > Computation > Qbox...` opens a native GTK4 Qbox window with:

- species mapping
- starter-deck generation
- practical setup controls for `ecutprec`, `dt`, spin, empty states, and charge mixing
- working-directory export
- pseudo and restart asset staging
- optional launcher prefix
- command helpers for `plot`, `spectrum`, `partial_charge`, `constraint`, `extforce`, and `set occ` / `print occ`
- direct `qbox` / `qb` launch support
- result/report session tracking

Other computation entries currently use the shared executable-path setup window:

- `GULP`
- `GAMESS`
- `Monty`
- `SIESTA`
- `VASP`
- `USPEX`

## What Is Still Incomplete

This rebuild is usable, but it is not yet full legacy parity.

Important remaining gaps include:

- some legacy interaction details such as full box-selection behavior and deeper selection semantics
- broader legacy notebook-style editing pages
- full slab/surface construction parity
- deeper plotting and analysis workflows
- complete external-code dialogs for non-Qbox backends
- full iso-surface engine parity
- some legacy feature families such as import/export branches, hidden-atom workflows, and specialized niche dialogs

For a more detailed feature-by-feature comparison, see [../RESTORATION_AUDIT.md](../RESTORATION_AUDIT.md).

## Prerequisites

You need:

- a C compiler
- `pkg-config`
- GTK4 development files

On macOS with Homebrew:

```sh
brew install gtk4 pkgconf
```

## Build

From this directory:

```sh
make
```

Helpful environment check:

```sh
make doctor
```

## Run

Run without a startup file:

```sh
make run
./run-gdis-gtk4
```

Run with one or more startup files:

```sh
make run ARGS="../../examples/water.xyz"
make run ARGS="../../models/deoxy.pdb ../../models/gibb.car"
./run-gdis-gtk4 ../../examples/water.xyz
./run-gdis-gtk4 ../../models/deoxy.pdb
./run-gdis-gtk4 ../../models/*
./build/gdis-gtk4 ../../examples/water.xyz
./build/gdis-gtk4 ../../models/deoxy.pdb
./build/gdis-gtk4 ../../models/*
```

Export a deterministic PNG directly from the GTK4 renderer:

```sh
./capture-gdis-gtk4 --model ../../models/water.car --output /tmp/water-view.png
./capture-gdis-gtk4 --model ../../models/water.car --output /tmp/water-iso.png --iso
./capture-gdis-gtk4 --model ../../models/water.car --output /tmp/water-shape.png --iso --color shape-index
make export-png MODEL=../../models/water.car OUT=/tmp/water-view.png
make export-png MODEL=../../models/water.car OUT=/tmp/water-iso.png ISO=1 COLOR=shape-index
```

Create and open a macOS app bundle (no manual binary path needed):

```sh
make app
make open-app
```

This creates:

```text
gtk4_rebuild/gtk4_app/build/GDIS.app
```

Create the app directly in the project `dist/` folder:

```sh
make dist-app
make open-dist-app
```

This creates:

```text
dist/GDIS.app
```

Create the full portable macOS release set:

```sh
make portable-release
```

This refreshes:

```text
dist/GDIS.app
dist/GDIS-macos-arm64.zip
dist/GDIS-macos-arm64.dmg
dist/GDIS Portable/
dist/GDIS-Portable-macos-arm64.zip
dist/GDIS-Portable-macos-arm64.dmg
```

## Typical Usage

1. Launch the app with or without a startup structure.
2. Use `File > Open` if you want to load another model into the current session.
3. Click in the viewer to select atoms and build pick history.
4. Use `Tools > Building > Editing...` for atom and bond edits.
5. Use `Tools > Analysis > Measurements...` for distance, angle, and torsion tools.
6. Use `Tools > Computation > Qbox...` if you want to prepare or execute a Qbox job.
7. Use `View > Executable paths...` to configure backend executables for the current session.

## macOS Notes

- The app builds as a normal GTK4 desktop binary on macOS.
- Qbox actions are also surfaced through a native macOS menu entry so they remain reachable from the system menu bar.
- The current layout is tuned for practical Mac usage, but the GTK4 rebuild still differs from the exact old Linux interface in a number of places.
- If desktop screenshots are unreliable, use `./capture-gdis-gtk4` or `make export-png ...` to export the viewer to PNG from inside the app itself. This avoids focus, scaling, and interrupted-window capture problems.

## Reference Files

- Legacy snapshot: [../legacy_snapshot](../legacy_snapshot)
- Restoration audit: [../RESTORATION_AUDIT.md](../RESTORATION_AUDIT.md)
- Porting notes: [../gtk4_app/PORTING.md](./PORTING.md)
