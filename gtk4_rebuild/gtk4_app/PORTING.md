GTK4 porting notes

Goal:
- rebuild the mac app as a native GTK4 single-window interface that matches the legacy Linux structure as closely as practical

Legacy reference:
- `../legacy_snapshot/src/gui_main.c`

Important legacy layout mapping:
- legacy menu bar -> GTK4 `GtkPopoverMenuBar`
- legacy toolbar icons -> GTK4 toolbar row of action buttons
- legacy `hpaned` main split -> GTK4 `GtkPaned` horizontal split
- legacy `mpane` left model pane -> GTK4 sidebar with model list + stack pages
- legacy `vpaned` viewer/status split -> GTK4 vertical `GtkPaned`
- legacy `display_box` -> GTK4 viewer placeholder for the future renderer
- legacy `tpane` -> GTK4 status/log text view

Planned next steps:
1. replace the placeholder viewer with a GTK4 rendering surface
2. port file loading into the GTK4 shell
3. connect model selection to real data
4. port editing/display/viewing panels incrementally
5. restore Linux-equivalent tools and dialogs in GTK4

Current scaffold boundary:
- UI shell only
- no chemistry/model engine integration yet
- no legacy OpenGL path reused yet
