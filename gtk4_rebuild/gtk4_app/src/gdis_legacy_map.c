#include <stddef.h>

#include "gdis_legacy_map.h"

const char *const gdis_legacy_sidebar_pages[] = {
  "Content",
  "Editing",
  "Display",
  "Images",
  "Symmetry",
  "Viewing",
  NULL
};

const char *const gdis_legacy_selection_modes[] = {
  "Select : Atoms",
  "Select : Atom Label",
  "Select : Atom FF Type",
  "Select : Elements",
  "Select : Elements in Molecule",
  "Select : Molecules",
  "Select : Molecule Fragments",
  "Select : Regions",
  NULL
};

const char *const gdis_legacy_toolbar_actions[][2] = {
  {"Open", "open"},
  {"Save", "save"},
  {"New", "new-model"},
  {"Edit", "edit"},
  {"Render", "render"},
  {"Measure", "measure"},
  {"Iso", "isosurface"},
  {"Surface", "surface"},
  {"Diffraction", "diffraction"},
  {"Qbox", "qbox"},
  {"Reset View", "reset-view"},
  {NULL, NULL}
};

const char *const gdis_legacy_model_samples[] = {
  "../../examples/water.xyz",
  "../../examples/benzene.xyz",
  "../../examples/rocksalt_demo.cif",
  "../../models/deoxy.pdb",
  "../../models/1_C2H4_HOOH_Tifer.arc",
  NULL
};
