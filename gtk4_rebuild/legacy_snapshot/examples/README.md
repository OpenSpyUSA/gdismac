# GDIS Quick Start

This folder contains a few small files you can open immediately in GDIS.

## Fastest Way To Try It

1. Launch the app:
   `/Users/erbai/Desktop/gdismac/dist/GDIS.app`
2. In GDIS, choose `File > Open...`
3. Open one of the files in this folder.

You can also open a file directly from Terminal:

```bash
open -n /Users/erbai/Desktop/gdismac/dist/GDIS.app --args /Users/erbai/Desktop/gdismac/examples/water.xyz
```

## Starter Files

- `water.xyz`
  Small molecule, easiest first test.
- `benzene.xyz`
  Slightly larger molecule, useful for rotate/zoom/select practice.
- `water_motion.ani`
  Multi-frame XYZ animation. Load this to try animation support.
- `rocksalt_demo.cif`
  Small periodic crystal example.

## Good Built-In Repo Examples

The project already ships with a larger `models/` folder. These are good next steps:

- `/Users/erbai/Desktop/gdismac/models/methane.gin`
- `/Users/erbai/Desktop/gdismac/models/burk1.cif`
- `/Users/erbai/Desktop/gdismac/models/deoxy.pdb`
- `/Users/erbai/Desktop/gdismac/models/gibb_opt_bulk.arc`

## Basic Controls

- Left click: select an atom
- Left drag: box-select atoms
- Right drag: rotate the model
- Shift + middle-drag: zoom
- Click a model name in the left pane to make it active

## Nice First Things To Try

- Open `water.xyz`, then right-drag to rotate it.
- Open `rocksalt_demo.cif`, then look at the left pane and active-model info.
- Open `water_motion.ani`, then use the animation tools/dialogs to inspect frames.
- Open one of the `models/*.gin` files if you want to explore the computational side of GDIS.
