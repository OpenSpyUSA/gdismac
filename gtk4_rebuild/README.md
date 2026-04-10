# GTK4 Rebuild Workspace

This directory keeps the modern GTK4 restoration work separate from the historical GDIS tree.

The goal is to make forward progress on modern macOS and Linux support without disturbing the original legacy source layout in the repository root.

## Structure

- `legacy_snapshot/`
  A frozen copy of the project tree from the point where the GTK4 rebuild began.
  Use this as reference material when checking legacy behavior and old workflows.
- `gtk4_app/`
  The active GTK4 application workspace, including the source code, build system, packaging scripts, and rebuild documentation.

The original project files at the repository root remain intact and are not moved or rewritten by this workspace.
