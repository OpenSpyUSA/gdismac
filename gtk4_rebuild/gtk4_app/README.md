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

Helpful checks:
```sh
make doctor
```

Notes:
- This machine did not have GTK4 installed when the scaffold was created, so the scaffold was not compiled locally yet.
- The current shell is intentionally focused on the application structure and Linux-style layout first.
- Rendering, file parsing, and the model/data bridge are the next porting steps.
