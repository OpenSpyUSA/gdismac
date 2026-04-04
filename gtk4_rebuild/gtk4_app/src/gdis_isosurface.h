#pragma once

#include <glib.h>

#include "gdis_model.h"

G_BEGIN_DECLS

typedef enum
{
  GDIS_ISOSURFACE_MODE_MOLECULAR = 0,
  GDIS_ISOSURFACE_MODE_PROMOLECULE,
  GDIS_ISOSURFACE_MODE_HIRSHFELD,
  GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY
} GdisIsosurfaceMode;

typedef struct
{
  gdouble position[3];
  gdouble normal[3];
} GdisIsoVertex;

typedef struct
{
  GdisIsoVertex vertices[3];
} GdisIsoTriangle;

typedef struct
{
  GdisIsosurfaceMode mode;
  gdouble grid_size;
  gdouble blur;
  gdouble isovalue;
  gdouble color_rgba[4];
  GArray *triangles;
  gchar *summary;
} GdisIsoSurface;

typedef struct
{
  GdisIsosurfaceMode mode;
  gdouble grid_size;
  gdouble blur;
  gdouble isovalue;
  const GArray *selected_atoms;
} GdisIsosurfaceSettings;

void gdis_isosurface_settings_init(GdisIsosurfaceSettings *settings);
const char *gdis_isosurface_mode_label(GdisIsosurfaceMode mode);
gboolean gdis_isosurface_mode_supported(GdisIsosurfaceMode mode);
gboolean gdis_isosurface_generate(const GdisModel *model,
                                  const GdisIsosurfaceSettings *settings,
                                  GdisIsoSurface **surface_out,
                                  GError **error);
void gdis_isosurface_free(GdisIsoSurface *surface);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GdisIsoSurface, gdis_isosurface_free)

G_END_DECLS
