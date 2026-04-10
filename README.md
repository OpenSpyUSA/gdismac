# GDIS

GDIS is a long-running crystal and molecular structure viewer and builder originally developed by Sean Fleming and Andrew Rohl.

This repository now carries both:

- the original legacy GDIS source tree
- an actively restored GTK4 rebuild for modern macOS and Linux in [`gtk4_rebuild/gtk4_app`](./gtk4_rebuild/gtk4_app)

GDIS comes with ABSOLUTELY NO WARRANTY. It remains free software under the GPLv2 terms distributed with this repository.

## Project Status

The historical GTK2 application remains in the repository for reference and legacy workflows.

The current macOS-focused work happens in the GTK4 rebuild:

- native GTK4 application window and menus
- restored model loading, viewing, editing, and analysis workflows
- Apple Silicon portable app packaging in `dist/`
- bundled GTK runtime for easier macOS distribution

For current rebuild details, see:

- [`gtk4_rebuild/README.md`](./gtk4_rebuild/README.md)
- [`gtk4_rebuild/gtk4_app/README.md`](./gtk4_rebuild/gtk4_app/README.md)
- [`gtk4_rebuild/RESTORATION_AUDIT.md`](./gtk4_rebuild/RESTORATION_AUDIT.md)

## macOS Quick Start

If you want the modern macOS build, the easiest path is the GTK4 rebuild and the packaged release assets.

Prebuilt release assets are published on the [Releases page](https://github.com/OpenSpyUSA/gdismac/releases).

To build the portable macOS release locally on Apple Silicon:

```sh
brew install gtk4 pkgconf
make -C gtk4_rebuild/gtk4_app portable-release
```

This produces:

- `dist/GDIS.app`
- `dist/GDIS-macos-arm64.zip`
- `dist/GDIS-macos-arm64.dmg`
- `dist/GDIS Portable/`
- `dist/GDIS-Portable-macos-arm64.zip`
- `dist/GDIS-Portable-macos-arm64.dmg`

To launch the rebuilt app directly from the repo:

```sh
./gtk4_rebuild/gtk4_app/build/gdis-gtk4 ./models/deoxy.pdb
./gtk4_rebuild/gtk4_app/build/gdis-gtk4 ./models/gibb.car
```

## Repository Layout

- `src/`, `bin/`, and `install`
  The original legacy codebase and build flow.
- `gtk4_rebuild/legacy_snapshot/`
  Frozen reference copy of the legacy tree at the start of the GTK4 rebuild.
- `gtk4_rebuild/gtk4_app/`
  The active GTK4 application, packaging scripts, and rebuild documentation.
- `models/` and `examples/`
  Sample structures and test inputs used by both the legacy and GTK4 apps.
- `dist/`
  Generated macOS release artifacts.

## Legacy Source Build

If you want the historical GTK2 build flow instead of the GTK4 rebuild, use the original installer:

```sh
./install
```

This path expects a working C compiler plus GTK2-era dependencies such as `gtk+2` and `gtkglext`.

[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/arohl/gdis)

## Credits

Some GDIS features depend on the following projects:

1. [CDD](http://www.inf.ethz.ch/personal/fukudak/cdd_home/) for the computation of halfspace intersections.
   Copyright Komei Fukuda
2. [GPeriodic](http://www.frantz.fi/software/gperiodic.php) for periodic table display and editing.
   Copyright 1999 Kyle R. Burton
3. [SgInfo](http://cci.lbl.gov/sginfo/) for space group lookup support.
   Copyright 1994-96 Ralf W. Grosse-Kunstleve
4. [Brute force symmetry analyzer](http://www.cobalt.chem.ucalgary.ca/ps/symmetry/).
   Copyright 1996 S. Pachkovsky

Optional external tools that can enhance the broader legacy GDIS workflow include:

- [POVRay](http://www.povray.org)
- [ImageMagick](http://imagemagick.org)
- [GULP](http://nanochemistry.curtin.edu.au/gulp/)
- [GAMESS](http://www.msg.chem.iastate.edu/GAMESS/)
- [SIESTA](http://departments.icmab.es/leem/siesta/)
- [Monty](http://www.vsc.science.ru.nl/deij/monty.html)
- [VASP](http://www.vasp.at/)
- [USPEX](http://www.uspex-team.org/en/uspex/overview)
