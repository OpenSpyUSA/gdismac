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

typedef enum
{
  GDIS_ISOSURFACE_COLOR_DEFAULT = 0,
  GDIS_ISOSURFACE_COLOR_AFM,
  GDIS_ISOSURFACE_COLOR_ELECTROSTATIC,
  GDIS_ISOSURFACE_COLOR_CURVEDNESS,
  GDIS_ISOSURFACE_COLOR_SHAPE_INDEX,
  GDIS_ISOSURFACE_COLOR_DE
} GdisIsosurfaceColorMode;

typedef struct
{
  gdouble position[3];
  gdouble normal[3];
  gdouble color_rgba[4];
  gdouble property_value;
} GdisIsoVertex;

typedef struct
{
  GdisIsoVertex vertices[3];
} GdisIsoTriangle;

typedef struct
{
  GdisIsosurfaceMode mode;
  GdisIsosurfaceColorMode color_mode;
  gdouble grid_size;
  gdouble blur;
  gdouble isovalue;
  gdouble color_rgba[4];
  GArray *triangles;
  gboolean used_color_fallback;
  gdouble property_min;
  gdouble property_max;
  gchar *summary;
} GdisIsoSurface;

typedef struct
{
  GdisIsosurfaceMode mode;
  GdisIsosurfaceColorMode color_mode;
  gdouble grid_size;
  gdouble blur;
  gdouble isovalue;
  gboolean electrostatic_autoscale;
  gdouble electrostatic_min;
  gdouble electrostatic_max;
  const GArray *selected_atoms;
} GdisIsosurfaceSettings;

void gdis_isosurface_settings_init(GdisIsosurfaceSettings *settings);
const char *gdis_isosurface_mode_label(GdisIsosurfaceMode mode);
const char *gdis_isosurface_color_mode_label(GdisIsosurfaceColorMode mode);
GdisIsosurfaceColorMode gdis_isosurface_default_color_mode(GdisIsosurfaceMode mode);
gboolean gdis_isosurface_mode_supported(GdisIsosurfaceMode mode);
gboolean gdis_isosurface_color_mode_supported(GdisIsosurfaceMode mode,
                                              GdisIsosurfaceColorMode color_mode);
gboolean gdis_isosurface_generate(const GdisModel *model,
                                  const GdisIsosurfaceSettings *settings,
                                  GdisIsoSurface **surface_out,
                                  GError **error);
void gdis_isosurface_free(GdisIsoSurface *surface);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GdisIsoSurface, gdis_isosurface_free)

G_END_DECLS
