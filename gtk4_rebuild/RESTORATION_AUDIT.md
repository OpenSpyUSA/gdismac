# GTK4 Restoration Audit

Last checked: 2026-04-03

This audit compares the legacy Linux-oriented GDIS menu surface in:

- `legacy_snapshot/src/gui_main.c`

against the current GTK4 rebuild surface and bridge in:

- `gtk4_app/src/gdis_gtk4_window.c`
- `gtk4_app/src/gdis_model.c`
- `gtk4_app/src/gdis_reports.c`

Status legend:

- `Restored`: the GTK4 rebuild has a real working equivalent now
- `Partial`: the GTK4 rebuild has a meaningful equivalent, but not the full legacy workflow
- `Missing`: the legacy function is not exposed in the GTK4 rebuild
- `Blocked`: this needs deeper backend or renderer porting, not just menu wiring

## File

| Legacy feature | Legacy callback | GTK4 status | Notes |
| --- | --- | --- | --- |
| New | `edit_model_create` | `Restored` | `app.new-model` creates an empty model and opens the GTK4 editor. |
| Open | `file_load_dialog` | `Restored` | Native GTK4 file dialog loads through `gdis_model_load()`. |
| Save | `file_save_dialog` | `Restored` | Native GTK4 save path writes via `gdis_model_save()`. |
| Close | `tree_select_delete` | `Restored` | Current model can be removed from the GTK4 session. |
| Import / Geomview | `gui_import_geomview` | `Missing` | No GTK4 import branch yet. |
| Import / Project | `gui_import_project` | `Missing` | No GTK4 import branch yet. |
| Import / Graph | `gui_import_graph` | `Missing` | No GTK4 import branch yet. |
| Export / Canvas snapshot | `image_export_dialog` | `Missing` | No GTK4 canvas snapshot export dialog yet. |
| Export / Graph data | `analysis_export_dialog` | `Missing` | No GTK4 graph export path yet. |
| Quit | `gdis_exit_test` | `Restored` | Native GTK4 app quit action. |

## Edit

| Legacy feature | Legacy callback | GTK4 status | Notes |
| --- | --- | --- | --- |
| Undo | `undo_active` | `Missing` | No undo stack exists in the GTK4 model bridge yet. |
| Copy | `select_copy` | `Missing` | No GTK4 clipboard serialization for atom selections yet. |
| Paste | `select_paste` | `Missing` | No GTK4 clipboard import path yet. |
| Colour | `select_colour` | `Missing` | The GTK4 bridge has no per-atom display override state yet. |
| Delete selected | `select_delete` | `Restored` | GTK4 now exposes this directly and routes through `gdis_model_delete_atoms()`. |
| Select all | `select_all` | `Restored` | GTK4 now exposes full-model selection from the menu. |
| Invert selection | `select_invert` | `Restored` | GTK4 now exposes inverted atom-set selection from the menu. |
| Hide selected | `select_hide` | `Missing` | No hidden/display-mask state is stored on `GdisAtom` yet. |
| Hide unselected | `unselect_hide` | `Missing` | Same limitation as above. |
| Unhide all | `unhide_atoms` | `Missing` | Same limitation as above. |

## Tools / Visualization

| Legacy feature | Legacy callback | GTK4 status | Notes |
| --- | --- | --- | --- |
| Animation | `gui_animate_dialog` | `Missing` | No trajectory or animation dialog in the GTK4 app yet. |
| Iso-surfaces | `gui_isosurf_dialog` | `Blocked` | GTK4 only shows a status report. Real parity needs volumetric grid loaders, scalar-field storage, surface extraction, and mesh rendering. |
| Periodic table | `gui_gperiodic_dialog` | `Missing` | No GTK4 periodic table dialog yet. |

## Tools / Building

| Legacy feature | Legacy callback | GTK4 status | Notes |
| --- | --- | --- | --- |
| Editing | `gui_edit_dialog` | `Partial` | GTK4 editor supports label/element/coordinate edits, atom add/delete, bond add/remove, and pick-driven bond editing. The full legacy editing surface is still larger. |
| Dislocations | `gui_defect_dialog` | `Missing` | No GTK4 defect/dislocation workflow yet. |
| Docking | `gui_dock_dialog` | `Missing` | No GTK4 docking workflow yet. |
| Dynamics | `gui_mdi_dialog` | `Missing` | No GTK4 structure-building dynamics dialog yet. |
| Surfaces | `surface_dialog` | `Partial` | GTK4 computes low-index planes and d-spacings, but not the full slab-construction workflow from the legacy dialog. |
| Zmatrix | `gui_zmat_dialog` | `Missing` | No GTK4 Z-matrix editor yet. |

## Tools / Computation

| Legacy feature | Legacy callback | GTK4 status | Notes |
| --- | --- | --- | --- |
| Diffraction | `gui_diffract_dialog` | `Partial` | GTK4 now has a real setup dialog, calculation backend, export path, and plot window. It is still a narrower workflow than the legacy tool. |
| GULP | `gulp_dialog` | `Missing` | No GTK4 external-code runner dialog yet. |
| GAMESS | `gamess_dialog` | `Missing` | No GTK4 external-code runner dialog yet. |
| Monty | `monty_dialog` | `Missing` | No GTK4 external-code runner dialog yet. |
| SIESTA | `gui_siesta_dialog` | `Missing` | No GTK4 external-code runner dialog yet. |
| VASP | `gui_vasp_dialog` | `Missing` | No GTK4 external-code runner dialog yet. |
| USPEX | `gui_uspex_dialog` | `Missing` | No GTK4 external-code runner dialog yet. |

## Tools / Analysis

| Legacy feature | Legacy callback | GTK4 status | Notes |
| --- | --- | --- | --- |
| Dynamics analysis | `gui_analysis_dialog` | `Missing` | No GTK4 analysis dialog yet. |
| Measurements | `gui_measure_dialog` | `Partial` | GTK4 has explicit `Auto`, `Distance`, `Angle`, and `Torsion` modes, plus periodic minimum-image geometry support. It still lacks the legacy persistent measurement objects, search workflow, list controls, and editable measurement values. |
| Plots | `gui_plots_dialog` | `Missing` | No GTK4 plotting tool yet. |

## View

| Legacy feature | Legacy callback | GTK4 status | Notes |
| --- | --- | --- | --- |
| Display properties | `gui_render_dialog` | `Partial` | GTK4 now opens a native display-properties window with viewer toggles, view presets, and image reset. The deeper legacy render pages are still missing. |
| Reset model images | `space_image_widget_reset` | `Restored` | GTK4 now exposes the image-range reset directly from the menu. |
| Normal mode | `gui_mode_default` | `Missing` | No GTK4 mode-switch path yet. |
| Recording mode | `gui_mode_record` | `Missing` | No GTK4 recording mode path yet. |
| Task manager | `task_dialog` | `Missing` | No GTK4 task manager yet. |
| Grid manager | `gui_grid_dialog` | `Missing` | Optional legacy feature not ported. |
| Executable paths | `gui_setup_dialog` | `Missing` | No GTK4 setup dialog for external code paths yet. |

## Help

| Legacy feature | Legacy callback | GTK4 status | Notes |
| --- | --- | --- | --- |
| About | `gui_about_dialog` | `Restored` | GTK4 now opens a real summary window instead of only logging text. |
| Manual | `gui_help_dialog` | `Partial` | GTK4 now opens a native usage guide, but not the original full help/manual dialog. |

## Current Practical Parity

The GTK4 rebuild is already usable for this workflow:

- load and save structure models
- rotate, zoom, and pick atoms in the native GTK4 viewer
- switch models inside the same session
- select atoms using the major legacy selection modes
- edit atom labels, elements, coordinates, and bonds
- delete selected atoms or selected groups
- inspect distances, angles, and torsions
- work with periodic image ranges and crystal cell cleanup
- inspect low-index surface planes
- generate a lightweight powder diffraction plot for 3D periodic models

## Behavioral Mismatches Still Visible To Users

These are the places where a legacy user will still feel that the GTK4 rebuild behaves differently, even when a rough equivalent exists in the menus.

1. Mouse interaction model: legacy GDIS used primary-button drag for box-selection and supported `Shift`-click add/remove selection. GTK4 still uses primary-button drag for rotation and ordinary click for one-shot pick selection.
2. Selection and pick history are still coupled: ordinary selection clicks also mutate the shared pick history used by measurement and bond-edit workflows.
3. Measurements are still a reduced workflow: GTK4 computes real values, but it does not yet restore the legacy measurement list, search mode, select/delete/dump actions, or editable measurement objects.
4. Editing is still a subset: the GTK4 editor restores direct atom/bond editing, but not the wider legacy notebook surface such as transformations, regions, labelling/library pages, click-to-place atom editing, or midpoint bond deletion.
5. Surface is still a report, not the full slab-construction environment from legacy GDIS.
6. Diffraction is the closest advanced-tool match, but it still lacks the broader legacy graph workflow and all-frames style behavior.
7. Region selection is still unavailable because the GTK4 bridge does not yet load region-labelled model data.

## Highest-Impact Remaining Gaps

The next sensible restoration targets are:

1. Restore the legacy interaction model better: additive/toggle selection, cleaner separation of selection vs measurement picks, and box selection.
2. Hidden/unhidden atom display state, so `Hide selected`, `Hide unselected`, and `Unhide all` become real again.
3. A fuller GTK4 `Display properties` port that restores the deeper legacy render pages beyond the current everyday subset window.
4. A full surface-construction dialog instead of the current plane-ranking report.
5. Plot windows and the remaining analysis tools.
6. External-code dialogs and path setup for GULP, GAMESS, Monty, SIESTA, VASP, and USPEX.
7. The deeper iso-surface engine port.
