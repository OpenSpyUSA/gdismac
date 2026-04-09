#include "gdis_isosurface.h"
#include "gdis_elements.h"

#include <math.h>
#include <string.h>

typedef struct
{
  gdouble point[3];
  gdouble value;
} GdisIsoSample;

typedef struct
{
  gdouble *pseudocharges;
  gdouble property_min;
  gdouble property_max;
  gboolean property_valid;
  gboolean used_color_fallback;
  gboolean used_heuristic_electrostatic;
  GdisIsosurfaceColorMode effective_color_mode;
} GdisIsoColoringInfo;


static const guint gdis_iso_tetrahedra[6][4] = {
  {0, 5, 1, 6},
  {0, 1, 2, 6},
  {0, 2, 3, 6},
  {0, 3, 7, 6},
  {0, 7, 4, 6},
  {0, 4, 5, 6}
};

static const guint gdis_iso_corner_offsets[8][3] = {
  {0, 0, 0},
  {1, 0, 0},
  {1, 1, 0},
  {0, 1, 0},
  {0, 0, 1},
  {1, 0, 1},
  {1, 1, 1},
  {0, 1, 1}
};

static void gdis_iso_vec3_copy(gdouble dest[3], const gdouble src[3]);
static void gdis_iso_vec3_set(gdouble dest[3], gdouble x, gdouble y, gdouble z);
static void gdis_iso_vec3_subtract(const gdouble a[3], const gdouble b[3], gdouble out[3]);
static gdouble gdis_iso_vec3_dot(const gdouble a[3], const gdouble b[3]);
static void gdis_iso_vec3_cross(const gdouble a[3], const gdouble b[3], gdouble out[3]);
static gdouble gdis_iso_vec3_length(const gdouble vector[3]);
static gboolean gdis_iso_vec3_normalize(gdouble vector[3]);
static gboolean gdis_iso_atom_array_contains(const GArray *array, guint atom_index);
static gboolean gdis_iso_mode_uses_selection(GdisIsosurfaceMode mode);
static gboolean gdis_iso_mode_supports_afm(GdisIsosurfaceMode mode);
static gdouble gdis_iso_atomic_density(const GdisElementInfo *element,
                                       gdouble blur,
                                       gdouble dist2);
static gdouble gdis_iso_density_sum(const GdisModel *model,
                                    const GArray *selected_atoms,
                                    gboolean selected_only,
                                    gdouble blur,
                                    const gdouble point[3]);
static gdouble gdis_iso_field_value(const GdisModel *model,
                                    const GdisIsosurfaceSettings *settings,
                                    const gdouble point[3]);
static void gdis_iso_field_gradient(const GdisModel *model,
                                    const GdisIsosurfaceSettings *settings,
                                    gdouble delta,
                                    const gdouble point[3],
                                    gdouble gradient_out[3]);
static gboolean gdis_iso_choose_grid(const GdisModel *model,
                                     const GdisIsosurfaceSettings *settings,
                                     gdouble min_bound[3],
                                     gdouble max_bound[3],
                                     guint dims[3],
                                     gdouble *step_out,
                                     GError **error);
static gboolean gdis_iso_append_triangle(GArray *triangles,
                                         const GdisModel *model,
                                         const GdisIsosurfaceSettings *settings,
                                         gdouble gradient_delta,
                                         const gdouble a[3],
                                         const gdouble b[3],
                                         const gdouble c[3]);
static void gdis_iso_hsv_to_rgb(gdouble hue,
                                gdouble saturation,
                                gdouble value,
                                gdouble rgb_out[3]);
static void gdis_iso_color_copy(gdouble dest[4], const gdouble src[4]);
static void gdis_iso_set_rgba(gdouble rgba[4],
                              gdouble r,
                              gdouble g,
                              gdouble b,
                              gdouble a);
static void gdis_iso_default_mode_color(GdisIsosurfaceMode mode,
                                        gdouble rgba_out[4]);
static void gdis_iso_touch_color(const GdisModel *model,
                                 const GdisIsosurfaceSettings *settings,
                                 const gdouble point[3],
                                 gdouble rgba_out[4]);
static gdouble gdis_iso_shape_property(const GdisModel *model,
                                       const GdisIsosurfaceSettings *settings,
                                       const gdouble point[3],
                                       gboolean curvedness,
                                       gboolean *valid_out);
static gdouble gdis_iso_de_property(const GdisModel *model,
                                    const GdisIsosurfaceSettings *settings,
                                    const gdouble point[3],
                                    gboolean *valid_out);
static gdouble gdis_iso_electronegativity(const GdisElementInfo *element);
static gdouble *gdis_iso_build_pseudocharges(const GdisModel *model);
static gdouble gdis_iso_electrostatic_property(const GdisModel *model,
                                               const GdisIsosurfaceSettings *settings,
                                               const gdouble *pseudocharges,
                                               const gdouble point[3],
                                               gboolean *valid_out);
static gdouble gdis_iso_vertex_property(const GdisModel *model,
                                        const GdisIsosurfaceSettings *settings,
                                        GdisIsosurfaceColorMode color_mode,
                                        const gdouble *pseudocharges,
                                        const gdouble point[3],
                                        gboolean *valid_out);
static void gdis_iso_apply_property_color(GdisIsosurfaceColorMode color_mode,
                                          gdouble value,
                                          gdouble min_value,
                                          gdouble max_value,
                                          gdouble rgba_out[4]);
static GdisIsosurfaceColorMode gdis_iso_resolve_color_mode(GdisIsosurfaceMode mode,
                                                           GdisIsosurfaceColorMode requested_mode,
                                                           gboolean *used_fallback_out);
static void gdis_iso_apply_coloring(const GdisModel *model,
                                    const GdisIsosurfaceSettings *settings,
                                    GdisIsoSurface *surface,
                                    GdisIsoColoringInfo *info);
static void gdis_iso_interpolate_vertex(const GdisIsoSample *inside,
                                        const GdisIsoSample *outside,
                                        gdouble out[3]);
static gboolean gdis_iso_polygonize_tetra(const GdisIsoSample samples[4],
                                          GArray *triangles,
                                          const GdisModel *model,
                                          const GdisIsosurfaceSettings *settings,
                                          gdouble gradient_delta);
static void gdis_iso_compute_area_volume(const GArray *triangles,
                                         gdouble *area_out,
                                         gdouble *volume_out);

void
gdis_isosurface_settings_init(GdisIsosurfaceSettings *settings)
{
  g_return_if_fail(settings != NULL);

  settings->mode = GDIS_ISOSURFACE_MODE_MOLECULAR;
  settings->color_mode = GDIS_ISOSURFACE_COLOR_DEFAULT;
  settings->grid_size = 0.30;
  settings->blur = 0.30;
  settings->isovalue = 0.08;
  settings->electrostatic_autoscale = TRUE;
  settings->electrostatic_min = -0.50;
  settings->electrostatic_max = 0.50;
  settings->selected_atoms = NULL;
}

const char *
gdis_isosurface_mode_label(GdisIsosurfaceMode mode)
{
  switch (mode)
    {
    case GDIS_ISOSURFACE_MODE_MOLECULAR:
      return "Molecular surface";
    case GDIS_ISOSURFACE_MODE_PROMOLECULE:
      return "Promolecule isosurface";
    case GDIS_ISOSURFACE_MODE_HIRSHFELD:
      return "Hirshfeld surface";
    case GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY:
      return "Electron density";
    default:
      return "Iso-surface";
    }
}

const char *
gdis_isosurface_color_mode_label(GdisIsosurfaceColorMode mode)
{
  switch (mode)
    {
    case GDIS_ISOSURFACE_COLOR_DEFAULT:
      return "Default";
    case GDIS_ISOSURFACE_COLOR_AFM:
      return "AFM";
    case GDIS_ISOSURFACE_COLOR_ELECTROSTATIC:
      return "Electrostatic";
    case GDIS_ISOSURFACE_COLOR_CURVEDNESS:
      return "Curvedness";
    case GDIS_ISOSURFACE_COLOR_SHAPE_INDEX:
      return "Shape Index";
    case GDIS_ISOSURFACE_COLOR_DE:
      return "De";
    default:
      return "Default";
    }
}

GdisIsosurfaceColorMode
gdis_isosurface_default_color_mode(GdisIsosurfaceMode mode)
{
  if (mode == GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY)
    return GDIS_ISOSURFACE_COLOR_DEFAULT;
  return GDIS_ISOSURFACE_COLOR_DEFAULT;
}

gboolean
gdis_isosurface_mode_supported(GdisIsosurfaceMode mode)
{
  return mode == GDIS_ISOSURFACE_MODE_MOLECULAR ||
         mode == GDIS_ISOSURFACE_MODE_PROMOLECULE ||
         mode == GDIS_ISOSURFACE_MODE_HIRSHFELD ||
         mode == GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY;
}

gboolean
gdis_isosurface_color_mode_supported(GdisIsosurfaceMode mode,
                                     GdisIsosurfaceColorMode color_mode)
{
  switch (color_mode)
    {
    case GDIS_ISOSURFACE_COLOR_DEFAULT:
    case GDIS_ISOSURFACE_COLOR_ELECTROSTATIC:
    case GDIS_ISOSURFACE_COLOR_CURVEDNESS:
    case GDIS_ISOSURFACE_COLOR_SHAPE_INDEX:
      return TRUE;
    case GDIS_ISOSURFACE_COLOR_AFM:
      return gdis_iso_mode_supports_afm(mode);
    case GDIS_ISOSURFACE_COLOR_DE:
      return mode == GDIS_ISOSURFACE_MODE_HIRSHFELD;
    default:
      return FALSE;
    }
}

void
gdis_isosurface_free(GdisIsoSurface *surface)
{
  if (!surface)
    return;

  g_clear_pointer(&surface->triangles, g_array_unref);
  g_clear_pointer(&surface->summary, g_free);
  g_free(surface);
}

gboolean
gdis_isosurface_generate(const GdisModel *model,
                         const GdisIsosurfaceSettings *settings,
                         GdisIsoSurface **surface_out,
                         GError **error)
{
  GdisIsoSurface *surface;
  GdisIsoColoringInfo coloring;
  gdouble min_bound[3];
  gdouble max_bound[3];
  guint dims[3];
  gdouble step;
  gdouble gradient_delta;
  GArray *triangles;
  gdouble *grid_values;
  gsize grid_point_count;
  guint nx;
  guint ny;
  guint nz;
  guint px;
  guint py;
  guint pz;
  guint focus_atom_count;
  gdouble area;
  gdouble volume;
  g_autofree gchar *selection_summary = NULL;
  g_autofree gchar *notes = NULL;

  g_return_val_if_fail(surface_out != NULL, FALSE);
  *surface_out = NULL;
  memset(&coloring, 0, sizeof(coloring));

  if (!model)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Iso-surface generation needs an active model.");
      return FALSE;
    }

  if (!model->atoms || model->atoms->len == 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Iso-surface generation needs at least one atom.");
      return FALSE;
    }

  if (!settings)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Iso-surface settings were missing.");
      return FALSE;
    }

  if (!gdis_isosurface_mode_supported(settings->mode))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "%s is not available in this app.",
                  gdis_isosurface_mode_label(settings->mode));
      return FALSE;
    }

  focus_atom_count = settings->selected_atoms ? settings->selected_atoms->len : 0u;
  if (gdis_iso_mode_uses_selection(settings->mode) && focus_atom_count == 0u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "%s needs a current atom selection first.",
                  gdis_isosurface_mode_label(settings->mode));
      return FALSE;
    }

  if (!gdis_iso_choose_grid(model, settings, min_bound, max_bound, dims, &step, error))
    return FALSE;

  nx = dims[0];
  ny = dims[1];
  nz = dims[2];
  px = nx + 1u;
  py = ny + 1u;
  pz = nz + 1u;
  grid_point_count = (gsize) px * (gsize) py * (gsize) pz;
  gradient_delta = MAX(step * 0.35, 0.05);

  grid_values = g_new(gdouble, grid_point_count);
  for (guint iz = 0; iz < pz; iz++)
    {
      for (guint iy = 0; iy < py; iy++)
        {
          for (guint ix = 0; ix < px; ix++)
            {
              gdouble point[3];
              gsize index;

              point[0] = min_bound[0] + (gdouble) ix * step;
              point[1] = min_bound[1] + (gdouble) iy * step;
              point[2] = min_bound[2] + (gdouble) iz * step;
              index = ((gsize) iz * py + iy) * px + ix;
              grid_values[index] = gdis_iso_field_value(model, settings, point);
            }
        }
    }

  triangles = g_array_new(FALSE, FALSE, sizeof(GdisIsoTriangle));
  for (guint iz = 0; iz < nz; iz++)
    {
      for (guint iy = 0; iy < ny; iy++)
        {
          for (guint ix = 0; ix < nx; ix++)
            {
              GdisIsoSample cube_samples[8];

              for (guint corner = 0; corner < 8; corner++)
                {
                  guint sx;
                  guint sy;
                  guint sz;
                  gsize point_index;

                  sx = ix + gdis_iso_corner_offsets[corner][0];
                  sy = iy + gdis_iso_corner_offsets[corner][1];
                  sz = iz + gdis_iso_corner_offsets[corner][2];
                  cube_samples[corner].point[0] = min_bound[0] + (gdouble) sx * step;
                  cube_samples[corner].point[1] = min_bound[1] + (gdouble) sy * step;
                  cube_samples[corner].point[2] = min_bound[2] + (gdouble) sz * step;
                  point_index = ((gsize) sz * py + sy) * px + sx;
                  cube_samples[corner].value = grid_values[point_index];
                }

              for (guint tetra = 0; tetra < G_N_ELEMENTS(gdis_iso_tetrahedra); tetra++)
                {
                  GdisIsoSample tetra_samples[4];

                  for (guint i = 0; i < 4; i++)
                    tetra_samples[i] = cube_samples[gdis_iso_tetrahedra[tetra][i]];

                  if (!gdis_iso_polygonize_tetra(tetra_samples,
                                                 triangles,
                                                 model,
                                                 settings,
                                                 gradient_delta))
                    {
                      g_array_unref(triangles);
                      g_free(grid_values);
                      g_set_error(error,
                                  GDIS_MODEL_ERROR,
                                  GDIS_MODEL_ERROR_FAILED,
                                  "Iso-surface generation exceeded the current triangle budget.");
                      return FALSE;
                    }
                }
            }
        }
    }
  g_free(grid_values);

  if (triangles->len == 0)
    {
      g_array_unref(triangles);
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "No iso-surface triangles were generated. Try a smaller grid size or different isovalue.");
      return FALSE;
    }

  surface = g_new0(GdisIsoSurface, 1);
  surface->mode = settings->mode;
  surface->color_mode = settings->color_mode;
  surface->grid_size = step;
  surface->blur = settings->blur;
  surface->isovalue = settings->isovalue;
  surface->triangles = triangles;
  gdis_iso_default_mode_color(settings->mode, surface->color_rgba);
  gdis_iso_apply_coloring(model, settings, surface, &coloring);
  surface->color_mode = coloring.effective_color_mode;
  surface->used_color_fallback = coloring.used_color_fallback;
  surface->property_min = coloring.property_min;
  surface->property_max = coloring.property_max;
  gdis_iso_compute_area_volume(triangles, &area, &volume);
  if (gdis_iso_mode_uses_selection(settings->mode))
    selection_summary = g_strdup_printf("%u selected atom%s",
                                        focus_atom_count,
                                        focus_atom_count == 1u ? "" : "s");
  else
    selection_summary = g_strdup("Whole model");
  surface->summary = g_strdup_printf(
    "%s\n"
    "Model: %s\n"
    "Selection: %s\n"
    "Colour method: %s\n"
    "Triangles: %u\n"
    "Estimated area: %.3f A^2\n"
    "Estimated volume: %.3f A^3\n"
    "Grid step: %.3f A\n"
    "Blur: %.3f\n"
    "Isovalue: %.4f\n",
    gdis_isosurface_mode_label(settings->mode),
    model->basename ? model->basename : "(model)",
    selection_summary,
    gdis_isosurface_color_mode_label(surface->color_mode),
    triangles->len,
    area,
    volume,
    step,
    settings->blur,
    settings->isovalue);
  if (coloring.property_valid &&
      (surface->color_mode == GDIS_ISOSURFACE_COLOR_AFM ||
       surface->color_mode == GDIS_ISOSURFACE_COLOR_ELECTROSTATIC ||
       surface->color_mode == GDIS_ISOSURFACE_COLOR_CURVEDNESS ||
       surface->color_mode == GDIS_ISOSURFACE_COLOR_SHAPE_INDEX ||
       surface->color_mode == GDIS_ISOSURFACE_COLOR_DE))
    {
      gchar *with_range;

      with_range = g_strdup_printf("%sProperty range: %.4f -> %.4f\n",
                                   surface->summary,
                                   surface->property_min,
                                   surface->property_max);
      g_free(surface->summary);
      surface->summary = with_range;
    }
  if (surface->used_color_fallback || coloring.used_heuristic_electrostatic)
    {
      GString *note_builder;

      note_builder = g_string_new("");
      if (surface->used_color_fallback)
        g_string_append_printf(note_builder,
                               "Colour note: %s is not available for %s, so GTK4 used %s.\n",
                               gdis_isosurface_color_mode_label(settings->color_mode),
                               gdis_isosurface_mode_label(settings->mode),
                               gdis_isosurface_color_mode_label(surface->color_mode));
      if (coloring.used_heuristic_electrostatic)
        g_string_append(note_builder,
                        "Colour note: electrostatic colouring uses the GTK4 bonded electronegativity field.\n");
      notes = g_string_free(note_builder, FALSE);
    }
  if (settings->mode == GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY)
    {
      gchar *with_note;

      with_note = g_strconcat(surface->summary,
                              "Mode note: electron density uses the GTK4 analytic whole-model density field.\n",
                              NULL);
      g_free(surface->summary);
      surface->summary = with_note;
    }
  if (notes)
    {
      gchar *with_notes;

      with_notes = g_strconcat(surface->summary, notes, NULL);
      g_free(surface->summary);
      surface->summary = with_notes;
    }
  g_free(coloring.pseudocharges);

  *surface_out = surface;
  return TRUE;
}

static gboolean
gdis_iso_mode_supports_afm(GdisIsosurfaceMode mode)
{
  return mode == GDIS_ISOSURFACE_MODE_MOLECULAR ||
         mode == GDIS_ISOSURFACE_MODE_PROMOLECULE;
}

static void
gdis_iso_hsv_to_rgb(gdouble hue,
                    gdouble saturation,
                    gdouble value,
                    gdouble rgb_out[3])
{
  gdouble chroma;
  gdouble hue_prime;
  gdouble x;
  gdouble m;
  gdouble r1;
  gdouble g1;
  gdouble b1;

  g_return_if_fail(rgb_out != NULL);

  hue = fmod(hue, 360.0);
  if (hue < 0.0)
    hue += 360.0;

  chroma = value * saturation;
  hue_prime = hue / 60.0;
  x = chroma * (1.0 - fabs(fmod(hue_prime, 2.0) - 1.0));
  r1 = g1 = b1 = 0.0;

  if (hue_prime < 1.0)
    {
      r1 = chroma;
      g1 = x;
    }
  else if (hue_prime < 2.0)
    {
      r1 = x;
      g1 = chroma;
    }
  else if (hue_prime < 3.0)
    {
      g1 = chroma;
      b1 = x;
    }
  else if (hue_prime < 4.0)
    {
      g1 = x;
      b1 = chroma;
    }
  else if (hue_prime < 5.0)
    {
      r1 = x;
      b1 = chroma;
    }
  else
    {
      r1 = chroma;
      b1 = x;
    }

  m = value - chroma;
  rgb_out[0] = r1 + m;
  rgb_out[1] = g1 + m;
  rgb_out[2] = b1 + m;
}

static void
gdis_iso_color_copy(gdouble dest[4], const gdouble src[4])
{
  dest[0] = src[0];
  dest[1] = src[1];
  dest[2] = src[2];
  dest[3] = src[3];
}

static void
gdis_iso_set_rgba(gdouble rgba[4],
                  gdouble r,
                  gdouble g,
                  gdouble b,
                  gdouble a)
{
  rgba[0] = CLAMP(r, 0.0, 1.0);
  rgba[1] = CLAMP(g, 0.0, 1.0);
  rgba[2] = CLAMP(b, 0.0, 1.0);
  rgba[3] = CLAMP(a, 0.0, 1.0);
}

static void
gdis_iso_default_mode_color(GdisIsosurfaceMode mode,
                            gdouble rgba_out[4])
{
  switch (mode)
    {
    case GDIS_ISOSURFACE_MODE_PROMOLECULE:
      gdis_iso_set_rgba(rgba_out, 0.87, 0.72, 0.26, 0.48);
      break;
    case GDIS_ISOSURFACE_MODE_HIRSHFELD:
      gdis_iso_set_rgba(rgba_out, 0.88, 0.46, 0.66, 0.44);
      break;
    case GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY:
      gdis_iso_set_rgba(rgba_out, 0.20, 0.20, 0.65, 0.42);
      break;
    case GDIS_ISOSURFACE_MODE_MOLECULAR:
    default:
      gdis_iso_set_rgba(rgba_out, 0.32, 0.78, 0.92, 0.44);
      break;
    }
}

static void
gdis_iso_touch_color(const GdisModel *model,
                     const GdisIsosurfaceSettings *settings,
                     const gdouble point[3],
                     gdouble rgba_out[4])
{
  gdouble fallback[4];
  const GdisAtom *touch_atom;
  guint touch_count;
  gdouble touch_factor;

  gdis_iso_default_mode_color(settings->mode, fallback);
  if (settings->mode == GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY)
    {
      gdis_iso_color_copy(rgba_out, fallback);
      return;
    }

  touch_atom = NULL;
  touch_count = 0u;
  touch_factor = exp(0.6666667 * MAX(settings->blur, 0.05));
  for (guint i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      const GdisElementInfo *element;
      gdouble dx;
      gdouble dy;
      gdouble dz;
      gdouble distance;
      gdouble radius;

      atom = g_ptr_array_index(model->atoms, i);
      element = gdis_element_lookup(atom->element);
      dx = point[0] - atom->position[0];
      dy = point[1] - atom->position[1];
      dz = point[2] - atom->position[2];
      distance = sqrt(dx * dx + dy * dy + dz * dz);
      radius = (element ? element->vdw_radius : 1.70) * touch_factor;
      if (distance < radius)
        {
          touch_atom = atom;
          touch_count++;
          if (touch_count > 1u)
            break;
        }
    }

  if (touch_count == 1u && touch_atom)
    {
      const GdisElementInfo *element;

      element = gdis_element_lookup(touch_atom->element);
      if (element)
        {
          gdis_iso_set_rgba(rgba_out,
                            element->color_rgb[0],
                            element->color_rgb[1],
                            element->color_rgb[2],
                            fallback[3]);
          return;
        }
    }

  if (model->atoms->len > 0u)
    {
      const GdisAtom *nearest_atom;
      gdouble nearest_distance2;

      nearest_atom = NULL;
      nearest_distance2 = G_MAXDOUBLE;
      for (guint i = 0; i < model->atoms->len; i++)
        {
          const GdisAtom *atom;
          gdouble dx;
          gdouble dy;
          gdouble dz;
          gdouble distance2;

          atom = g_ptr_array_index(model->atoms, i);
          dx = point[0] - atom->position[0];
          dy = point[1] - atom->position[1];
          dz = point[2] - atom->position[2];
          distance2 = dx * dx + dy * dy + dz * dz;
          if (distance2 < nearest_distance2)
            {
              nearest_distance2 = distance2;
              nearest_atom = atom;
            }
        }
      if (nearest_atom)
        {
          const GdisElementInfo *element;

          element = gdis_element_lookup(nearest_atom->element);
          if (element)
            {
              gdis_iso_set_rgba(rgba_out,
                                element->color_rgb[0],
                                element->color_rgb[1],
                                element->color_rgb[2],
                                fallback[3]);
              return;
            }
        }
    }

  gdis_iso_color_copy(rgba_out, fallback);
}

static gdouble
gdis_iso_shape_property(const GdisModel *model,
                        const GdisIsosurfaceSettings *settings,
                        const gdouble point[3],
                        gboolean curvedness,
                        gboolean *valid_out)
{
  const gdouble delta = 0.5;
  const gdouble one_over_2_delta = 0.5 / delta;
  const gdouble d2 = 1.0 / (delta * delta);
  const gdouble d4 = d2 / 4.0;
  gdouble f;
  gdouble fx;
  gdouble fy;
  gdouble fz;
  gdouble f_x;
  gdouble f_y;
  gdouble f_z;
  gdouble fxy;
  gdouble fx_y;
  gdouble f_xy;
  gdouble f_x_y;
  gdouble fxz;
  gdouble fx_z;
  gdouble f_xz;
  gdouble f_x_z;
  gdouble fyz;
  gdouble fy_z;
  gdouble f_yz;
  gdouble f_y_z;
  gdouble h[3][3];
  gdouble grad[3];
  gdouble grad_mag;
  gdouble normal[3];
  gdouble tangent1[3];
  gdouble tangent2[3];
  gdouble helper[3];
  gdouble ht1[3];
  gdouble ht2[3];
  gdouble m00;
  gdouble m01;
  gdouble m11;
  gdouble theta;
  gdouble sin_theta;
  gdouble cos_theta;
  gdouble eig0;
  gdouble eig1;
  gdouble h1;
  gdouble h2;
  gdouble sample[3];

  if (valid_out)
    *valid_out = FALSE;

  gdis_iso_vec3_copy(sample, point);
  f = gdis_iso_field_value(model, settings, sample);
  sample[0] = point[0] + delta;
  sample[1] = point[1];
  sample[2] = point[2];
  fx = gdis_iso_field_value(model, settings, sample);
  sample[0] = point[0] - delta;
  f_x = gdis_iso_field_value(model, settings, sample);
  sample[0] = point[0];
  sample[1] = point[1] + delta;
  fy = gdis_iso_field_value(model, settings, sample);
  sample[1] = point[1] - delta;
  f_y = gdis_iso_field_value(model, settings, sample);
  sample[1] = point[1];
  sample[2] = point[2] + delta;
  fz = gdis_iso_field_value(model, settings, sample);
  sample[2] = point[2] - delta;
  f_z = gdis_iso_field_value(model, settings, sample);

  sample[0] = point[0] + delta;
  sample[1] = point[1] + delta;
  sample[2] = point[2];
  fxy = gdis_iso_field_value(model, settings, sample);
  sample[1] = point[1] - delta;
  fx_y = gdis_iso_field_value(model, settings, sample);
  sample[0] = point[0] - delta;
  sample[1] = point[1] + delta;
  f_xy = gdis_iso_field_value(model, settings, sample);
  sample[1] = point[1] - delta;
  f_x_y = gdis_iso_field_value(model, settings, sample);

  sample[0] = point[0] + delta;
  sample[1] = point[1];
  sample[2] = point[2] + delta;
  fxz = gdis_iso_field_value(model, settings, sample);
  sample[2] = point[2] - delta;
  fx_z = gdis_iso_field_value(model, settings, sample);
  sample[0] = point[0] - delta;
  sample[2] = point[2] + delta;
  f_xz = gdis_iso_field_value(model, settings, sample);
  sample[2] = point[2] - delta;
  f_x_z = gdis_iso_field_value(model, settings, sample);

  sample[0] = point[0];
  sample[1] = point[1] + delta;
  sample[2] = point[2] + delta;
  fyz = gdis_iso_field_value(model, settings, sample);
  sample[2] = point[2] - delta;
  fy_z = gdis_iso_field_value(model, settings, sample);
  sample[1] = point[1] - delta;
  sample[2] = point[2] + delta;
  f_yz = gdis_iso_field_value(model, settings, sample);
  sample[2] = point[2] - delta;
  f_y_z = gdis_iso_field_value(model, settings, sample);

  h[0][0] = (fx + f_x - 2.0 * f) * d2;
  h[1][1] = (fy + f_y - 2.0 * f) * d2;
  h[2][2] = (fz + f_z - 2.0 * f) * d2;
  h[0][1] = h[1][0] = (fxy + f_x_y - f_xy - fx_y) * d4;
  h[0][2] = h[2][0] = (fxz + f_x_z - f_xz - fx_z) * d4;
  h[1][2] = h[2][1] = (fyz + f_y_z - f_yz - fy_z) * d4;

  grad[0] = (fx - f_x) * one_over_2_delta;
  grad[1] = (fy - f_y) * one_over_2_delta;
  grad[2] = (fz - f_z) * one_over_2_delta;
  grad_mag = gdis_iso_vec3_length(grad);
  if (grad_mag <= 1.0e-12)
    return 0.0;

  gdis_iso_vec3_copy(normal, grad);
  gdis_iso_vec3_normalize(normal);
  gdis_iso_vec3_set(helper, 0.0, 0.0, 1.0);
  if (fabs(gdis_iso_vec3_dot(normal, helper)) > 0.90)
    gdis_iso_vec3_set(helper, 0.0, 1.0, 0.0);
  gdis_iso_vec3_cross(normal, helper, tangent1);
  if (!gdis_iso_vec3_normalize(tangent1))
    return 0.0;
  gdis_iso_vec3_cross(normal, tangent1, tangent2);
  if (!gdis_iso_vec3_normalize(tangent2))
    return 0.0;

  ht1[0] = h[0][0] * tangent1[0] + h[0][1] * tangent1[1] + h[0][2] * tangent1[2];
  ht1[1] = h[1][0] * tangent1[0] + h[1][1] * tangent1[1] + h[1][2] * tangent1[2];
  ht1[2] = h[2][0] * tangent1[0] + h[2][1] * tangent1[1] + h[2][2] * tangent1[2];
  ht2[0] = h[0][0] * tangent2[0] + h[0][1] * tangent2[1] + h[0][2] * tangent2[2];
  ht2[1] = h[1][0] * tangent2[0] + h[1][1] * tangent2[1] + h[1][2] * tangent2[2];
  ht2[2] = h[2][0] * tangent2[0] + h[2][1] * tangent2[1] + h[2][2] * tangent2[2];
  m00 = gdis_iso_vec3_dot(ht1, tangent1);
  m01 = gdis_iso_vec3_dot(ht1, tangent2);
  m11 = gdis_iso_vec3_dot(ht2, tangent2);

  theta = 0.5 * atan2(2.0 * m01, m11 - m00);
  sin_theta = sin(theta);
  cos_theta = cos(theta);
  eig0 = cos_theta * cos_theta * m00 + sin_theta * sin_theta * m11 - 2.0 * sin_theta * cos_theta * m01;
  eig1 = sin_theta * sin_theta * m00 + cos_theta * cos_theta * m11 + 2.0 * sin_theta * cos_theta * m01;
  h1 = -eig0 / grad_mag;
  h2 = -eig1 / grad_mag;

  if (valid_out)
    *valid_out = TRUE;
  if (curvedness)
    {
      gdouble mean_square;

      mean_square = (h1 * h1 + h2 * h2) / 2.0;
      if (mean_square <= 1.0e-18)
        mean_square = 1.0e-18;
      return log(mean_square) / G_PI;
    }

  if (fabs(h1 - h2) <= 1.0e-12)
    return 0.0;
  return 2.0 * atan2(h1 + h2, MAX(h1, h2) - MIN(h1, h2)) / G_PI;
}

static gdouble
gdis_iso_de_property(const GdisModel *model,
                     const GdisIsosurfaceSettings *settings,
                     const gdouble point[3],
                     gboolean *valid_out)
{
  gdouble min_distance;
  gboolean found;

  if (valid_out)
    *valid_out = FALSE;
  if (!settings->selected_atoms || settings->selected_atoms->len == 0u)
    return 0.0;

  min_distance = G_MAXDOUBLE;
  found = FALSE;
  for (guint i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      gdouble dx;
      gdouble dy;
      gdouble dz;
      gdouble distance;

      if (gdis_iso_atom_array_contains(settings->selected_atoms, i))
        continue;

      atom = g_ptr_array_index(model->atoms, i);
      dx = point[0] - atom->position[0];
      dy = point[1] - atom->position[1];
      dz = point[2] - atom->position[2];
      distance = sqrt(dx * dx + dy * dy + dz * dz);
      if (distance < min_distance)
        min_distance = distance;
      found = TRUE;
    }

  if (valid_out)
    *valid_out = found;
  return found ? min_distance : 0.0;
}

static gdouble
gdis_iso_electronegativity(const GdisElementInfo *element)
{
  if (!element)
    return 2.5;

  switch (element->atomic_number)
    {
    case 1u: return 2.20;
    case 5u: return 2.04;
    case 6u: return 2.55;
    case 7u: return 3.04;
    case 8u: return 3.44;
    case 9u: return 3.98;
    case 14u: return 1.90;
    case 15u: return 2.19;
    case 16u: return 2.58;
    case 17u: return 3.16;
    case 35u: return 2.96;
    case 53u: return 2.66;
    default:
      break;
    }

  switch (element->family)
    {
    case GDIS_ELEMENT_FAMILY_HALOGEN:
      return 3.1;
    case GDIS_ELEMENT_FAMILY_NONMETAL:
      return 2.7;
    case GDIS_ELEMENT_FAMILY_METALLOID:
      return 2.1;
    case GDIS_ELEMENT_FAMILY_POST_TRANSITION_METAL:
      return 1.7;
    case GDIS_ELEMENT_FAMILY_TRANSITION_METAL:
      return 1.8;
    case GDIS_ELEMENT_FAMILY_ALKALINE_EARTH_METAL:
      return 1.1;
    case GDIS_ELEMENT_FAMILY_ALKALI_METAL:
      return 0.9;
    default:
      return 2.0;
    }
}

static gdouble *
gdis_iso_build_pseudocharges(const GdisModel *model)
{
  gdouble *charges;

  g_return_val_if_fail(model != NULL, NULL);

  charges = g_new0(gdouble, model->atoms->len);
  for (guint i = 0; i < model->bonds->len; i++)
    {
      const GdisBond *bond;
      const GdisAtom *atom_a;
      const GdisAtom *atom_b;
      const GdisElementInfo *element_a;
      const GdisElementInfo *element_b;
      gdouble chi_a;
      gdouble chi_b;
      gdouble transfer;

      bond = &g_array_index(model->bonds, GdisBond, i);
      if (bond->atom_index_a >= model->atoms->len || bond->atom_index_b >= model->atoms->len)
        continue;

      atom_a = g_ptr_array_index(model->atoms, bond->atom_index_a);
      atom_b = g_ptr_array_index(model->atoms, bond->atom_index_b);
      element_a = gdis_element_lookup(atom_a->element);
      element_b = gdis_element_lookup(atom_b->element);
      chi_a = gdis_iso_electronegativity(element_a);
      chi_b = gdis_iso_electronegativity(element_b);
      transfer = CLAMP(0.18 * (chi_b - chi_a), -0.75, 0.75);
      charges[bond->atom_index_a] += transfer;
      charges[bond->atom_index_b] -= transfer;
    }

  return charges;
}

static gdouble
gdis_iso_electrostatic_property(const GdisModel *model,
                                const GdisIsosurfaceSettings *settings,
                                const gdouble *pseudocharges,
                                const gdouble point[3],
                                gboolean *valid_out)
{
  gdouble potential;
  gdouble softness;

  if (valid_out)
    *valid_out = FALSE;

  if (!model || !pseudocharges)
    return 0.0;

  potential = 0.0;
  softness = MAX(0.30, settings->blur + 0.20);
  for (guint i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      gdouble dx;
      gdouble dy;
      gdouble dz;
      gdouble dist2;

      atom = g_ptr_array_index(model->atoms, i);
      dx = point[0] - atom->position[0];
      dy = point[1] - atom->position[1];
      dz = point[2] - atom->position[2];
      dist2 = dx * dx + dy * dy + dz * dz;
      potential += pseudocharges[i] / sqrt(dist2 + softness * softness);
    }

  if (valid_out)
    *valid_out = TRUE;
  return potential;
}

static gdouble
gdis_iso_vertex_property(const GdisModel *model,
                         const GdisIsosurfaceSettings *settings,
                         GdisIsosurfaceColorMode color_mode,
                         const gdouble *pseudocharges,
                         const gdouble point[3],
                         gboolean *valid_out)
{
  if (valid_out)
    *valid_out = FALSE;

  switch (color_mode)
    {
    case GDIS_ISOSURFACE_COLOR_AFM:
      if (valid_out)
        *valid_out = TRUE;
      return point[2];
    case GDIS_ISOSURFACE_COLOR_CURVEDNESS:
      return gdis_iso_shape_property(model, settings, point, TRUE, valid_out);
    case GDIS_ISOSURFACE_COLOR_SHAPE_INDEX:
      return gdis_iso_shape_property(model, settings, point, FALSE, valid_out);
    case GDIS_ISOSURFACE_COLOR_DE:
      return gdis_iso_de_property(model, settings, point, valid_out);
    case GDIS_ISOSURFACE_COLOR_ELECTROSTATIC:
      return gdis_iso_electrostatic_property(model, settings, pseudocharges, point, valid_out);
    case GDIS_ISOSURFACE_COLOR_DEFAULT:
    default:
      if (valid_out)
        *valid_out = FALSE;
      return 0.0;
    }
}

static void
gdis_iso_apply_property_color(GdisIsosurfaceColorMode color_mode,
                              gdouble value,
                              gdouble min_value,
                              gdouble max_value,
                              gdouble rgba_out[4])
{
  gdouble normalized;
  gdouble rgb[3];

  if (fabs(max_value - min_value) <= 1.0e-12)
    normalized = 0.5;
  else
    normalized = CLAMP((value - min_value) / (max_value - min_value), 0.0, 1.0);

  switch (color_mode)
    {
    case GDIS_ISOSURFACE_COLOR_AFM:
      gdis_iso_set_rgba(rgba_out,
                        0.7 * normalized + 0.2,
                        0.8 * normalized,
                        0.5 * pow(normalized, 4.0),
                        0.46);
      return;
    case GDIS_ISOSURFACE_COLOR_ELECTROSTATIC:
      if (value < 0.0)
        {
          gdouble scale;

          scale = (fabs(min_value) <= 1.0e-12) ? 0.0 : CLAMP(1.0 - value / min_value, 0.0, 1.0);
          gdis_iso_set_rgba(rgba_out, 1.0, scale, scale, 0.46);
        }
      else
        {
          gdouble scale;

          scale = (fabs(max_value) <= 1.0e-12) ? 0.0 : CLAMP(1.0 - value / max_value, 0.0, 1.0);
          gdis_iso_set_rgba(rgba_out, scale, scale, 1.0, 0.46);
        }
      return;
    case GDIS_ISOSURFACE_COLOR_CURVEDNESS:
    case GDIS_ISOSURFACE_COLOR_SHAPE_INDEX:
    case GDIS_ISOSURFACE_COLOR_DE:
      gdis_iso_hsv_to_rgb(normalized * 240.0, 1.0, 1.0, rgb);
      gdis_iso_set_rgba(rgba_out, rgb[0], rgb[1], rgb[2], 0.46);
      return;
    case GDIS_ISOSURFACE_COLOR_DEFAULT:
    default:
      gdis_iso_set_rgba(rgba_out, 0.7, 0.7, 0.7, 0.46);
      return;
    }
}

static GdisIsosurfaceColorMode
gdis_iso_resolve_color_mode(GdisIsosurfaceMode mode,
                            GdisIsosurfaceColorMode requested_mode,
                            gboolean *used_fallback_out)
{
  GdisIsosurfaceColorMode resolved;

  resolved = requested_mode;
  if (!gdis_isosurface_color_mode_supported(mode, resolved))
    {
      resolved = GDIS_ISOSURFACE_COLOR_DEFAULT;
      if (used_fallback_out)
        *used_fallback_out = TRUE;
    }
  else if (used_fallback_out)
    {
      *used_fallback_out = FALSE;
    }

  return resolved;
}

static void
gdis_iso_apply_coloring(const GdisModel *model,
                        const GdisIsosurfaceSettings *settings,
                        GdisIsoSurface *surface,
                        GdisIsoColoringInfo *info)
{
  GdisIsosurfaceColorMode color_mode;
  gboolean needs_range_scan;

  g_return_if_fail(model != NULL);
  g_return_if_fail(settings != NULL);
  g_return_if_fail(surface != NULL);
  g_return_if_fail(info != NULL);

  color_mode = gdis_iso_resolve_color_mode(settings->mode,
                                           settings->color_mode,
                                           &info->used_color_fallback);
  info->effective_color_mode = color_mode;
  info->property_valid = FALSE;
  info->property_min = 0.0;
  info->property_max = 0.0;

  if (color_mode == GDIS_ISOSURFACE_COLOR_ELECTROSTATIC)
    {
      info->pseudocharges = gdis_iso_build_pseudocharges(model);
      info->used_heuristic_electrostatic = TRUE;
    }

  needs_range_scan = (color_mode == GDIS_ISOSURFACE_COLOR_AFM ||
                      color_mode == GDIS_ISOSURFACE_COLOR_DE ||
                      (color_mode == GDIS_ISOSURFACE_COLOR_ELECTROSTATIC &&
                       settings->electrostatic_autoscale));
  if (color_mode == GDIS_ISOSURFACE_COLOR_CURVEDNESS)
    {
      info->property_min = -4.0;
      info->property_max = 0.4;
      info->property_valid = TRUE;
    }
  else if (color_mode == GDIS_ISOSURFACE_COLOR_SHAPE_INDEX)
    {
      info->property_min = -1.0;
      info->property_max = 1.0;
      info->property_valid = TRUE;
    }
  else if (color_mode == GDIS_ISOSURFACE_COLOR_ELECTROSTATIC &&
           !settings->electrostatic_autoscale)
    {
      info->property_min = settings->electrostatic_min;
      info->property_max = settings->electrostatic_max;
      info->property_valid = TRUE;
    }

  if (needs_range_scan)
    {
      for (guint i = 0; i < surface->triangles->len; i++)
        {
          GdisIsoTriangle *triangle;

          triangle = &g_array_index(surface->triangles, GdisIsoTriangle, i);
          for (guint v = 0; v < 3; v++)
            {
              gboolean valid;
              gdouble value;

              value = gdis_iso_vertex_property(model,
                                               settings,
                                               color_mode,
                                               info->pseudocharges,
                                               triangle->vertices[v].position,
                                               &valid);
              if (!valid)
                continue;

              if (!info->property_valid)
                {
                  info->property_min = value;
                  info->property_max = value;
                  info->property_valid = TRUE;
                }
              else
                {
                  info->property_min = MIN(info->property_min, value);
                  info->property_max = MAX(info->property_max, value);
                }
            }
        }
    }

  for (guint i = 0; i < surface->triangles->len; i++)
    {
      GdisIsoTriangle *triangle;

      triangle = &g_array_index(surface->triangles, GdisIsoTriangle, i);
      for (guint v = 0; v < 3; v++)
        {
          gboolean valid;
          gdouble value;

          triangle->vertices[v].property_value = 0.0;
          if (color_mode == GDIS_ISOSURFACE_COLOR_DEFAULT)
            {
              gdis_iso_touch_color(model,
                                   settings,
                                   triangle->vertices[v].position,
                                   triangle->vertices[v].color_rgba);
              triangle->vertices[v].color_rgba[3] = surface->color_rgba[3];
              continue;
            }

          value = gdis_iso_vertex_property(model,
                                           settings,
                                           color_mode,
                                           info->pseudocharges,
                                           triangle->vertices[v].position,
                                           &valid);
          triangle->vertices[v].property_value = value;
          if (!valid || !info->property_valid)
            {
              gdis_iso_touch_color(model,
                                   settings,
                                   triangle->vertices[v].position,
                                   triangle->vertices[v].color_rgba);
              triangle->vertices[v].color_rgba[3] = surface->color_rgba[3];
              continue;
            }

          gdis_iso_apply_property_color(color_mode,
                                        value,
                                        info->property_min,
                                        info->property_max,
                                        triangle->vertices[v].color_rgba);
          triangle->vertices[v].color_rgba[3] = surface->color_rgba[3];
        }
    }
}

static void
gdis_iso_vec3_copy(gdouble dest[3], const gdouble src[3])
{
  dest[0] = src[0];
  dest[1] = src[1];
  dest[2] = src[2];
}

static void
gdis_iso_vec3_set(gdouble dest[3], gdouble x, gdouble y, gdouble z)
{
  dest[0] = x;
  dest[1] = y;
  dest[2] = z;
}

static void
gdis_iso_vec3_subtract(const gdouble a[3], const gdouble b[3], gdouble out[3])
{
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}

static gdouble
gdis_iso_vec3_dot(const gdouble a[3], const gdouble b[3])
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void
gdis_iso_vec3_cross(const gdouble a[3], const gdouble b[3], gdouble out[3])
{
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

static gdouble
gdis_iso_vec3_length(const gdouble vector[3])
{
  return sqrt(gdis_iso_vec3_dot(vector, vector));
}

static gboolean
gdis_iso_vec3_normalize(gdouble vector[3])
{
  gdouble length;

  length = gdis_iso_vec3_length(vector);
  if (length <= 1.0e-12)
    return FALSE;

  vector[0] /= length;
  vector[1] /= length;
  vector[2] /= length;
  return TRUE;
}

static gboolean
gdis_iso_atom_array_contains(const GArray *array, guint atom_index)
{
  if (!array)
    return FALSE;

  for (guint i = 0; i < array->len; i++)
    {
      if (g_array_index(array, guint, i) == atom_index)
        return TRUE;
    }

  return FALSE;
}

static gboolean
gdis_iso_mode_uses_selection(GdisIsosurfaceMode mode)
{
  return mode == GDIS_ISOSURFACE_MODE_PROMOLECULE ||
         mode == GDIS_ISOSURFACE_MODE_HIRSHFELD;
}

static gdouble
gdis_iso_atomic_density(const GdisElementInfo *element,
                        gdouble blur,
                        gdouble dist2)
{
  gdouble sigma;
  gdouble weight;

  sigma = MAX(0.18, ((element ? element->vdw_radius : 1.70) * 0.38) + blur * 0.22);
  weight = (gdouble) (element ? element->atomic_number : 6u);
  return weight * exp(-(dist2 / (sigma * sigma)));
}

static gdouble
gdis_iso_density_sum(const GdisModel *model,
                     const GArray *selected_atoms,
                     gboolean selected_only,
                     gdouble blur,
                     const gdouble point[3])
{
  gdouble density;

  density = 0.0;
  for (guint i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      const GdisElementInfo *element;
      gdouble dx;
      gdouble dy;
      gdouble dz;
      gdouble dist2;

      if (selected_only && !gdis_iso_atom_array_contains(selected_atoms, i))
        continue;

      atom = g_ptr_array_index(model->atoms, i);
      element = gdis_element_lookup(atom->element);
      dx = point[0] - atom->position[0];
      dy = point[1] - atom->position[1];
      dz = point[2] - atom->position[2];
      dist2 = dx * dx + dy * dy + dz * dz;
      density += gdis_iso_atomic_density(element, blur, dist2);
    }

  return density;
}

static gdouble
gdis_iso_field_value(const GdisModel *model,
                     const GdisIsosurfaceSettings *settings,
                     const gdouble point[3])
{
  GdisIsosurfaceMode mode;
  gdouble value;

  g_return_val_if_fail(settings != NULL, -G_MAXDOUBLE);

  mode = settings->mode;
  value = (mode == GDIS_ISOSURFACE_MODE_MOLECULAR) ? -G_MAXDOUBLE : 0.0;

  for (guint i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      const GdisElementInfo *element;
      gdouble dx;
      gdouble dy;
      gdouble dz;
      gdouble dist2;

      atom = g_ptr_array_index(model->atoms, i);
      element = gdis_element_lookup(atom->element);
      dx = point[0] - atom->position[0];
      dy = point[1] - atom->position[1];
      dz = point[2] - atom->position[2];
      dist2 = dx * dx + dy * dy + dz * dz;

      if (mode == GDIS_ISOSURFACE_MODE_MOLECULAR)
        {
          gdouble radius;
          gdouble distance;
          gdouble signed_distance;

          radius = (element ? element->vdw_radius : 1.70) + settings->blur;
          distance = sqrt(dist2);
          signed_distance = radius - distance;
          if (signed_distance > value)
            value = signed_distance;
        }
    }

  if (mode == GDIS_ISOSURFACE_MODE_MOLECULAR)
    return value;

  if (mode == GDIS_ISOSURFACE_MODE_PROMOLECULE)
    {
      value = gdis_iso_density_sum(model,
                                   settings->selected_atoms,
                                   TRUE,
                                   settings->blur,
                                   point);
      return value - settings->isovalue;
    }

  if (mode == GDIS_ISOSURFACE_MODE_HIRSHFELD)
    {
      gdouble focus_density;
      gdouble total_density;

      focus_density = gdis_iso_density_sum(model,
                                           settings->selected_atoms,
                                           TRUE,
                                           settings->blur,
                                           point);
      total_density = gdis_iso_density_sum(model,
                                           settings->selected_atoms,
                                           FALSE,
                                           settings->blur,
                                           point);
      if (total_density <= 1.0e-12)
        return -settings->isovalue;
      return (focus_density / total_density) - settings->isovalue;
    }

  if (mode == GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY)
    {
      value = gdis_iso_density_sum(model,
                                   settings->selected_atoms,
                                   FALSE,
                                   settings->blur,
                                   point);
      return value - settings->isovalue;
    }

  return -settings->isovalue;
}

static void
gdis_iso_field_gradient(const GdisModel *model,
                        const GdisIsosurfaceSettings *settings,
                        gdouble delta,
                        const gdouble point[3],
                        gdouble gradient_out[3])
{
  gdouble sample[3];

  gdis_iso_vec3_copy(sample, point);
  sample[0] = point[0] + delta;
  gradient_out[0] = gdis_iso_field_value(model, settings, sample);
  sample[0] = point[0] - delta;
  gradient_out[0] -= gdis_iso_field_value(model, settings, sample);

  gdis_iso_vec3_copy(sample, point);
  sample[1] = point[1] + delta;
  gradient_out[1] = gdis_iso_field_value(model, settings, sample);
  sample[1] = point[1] - delta;
  gradient_out[1] -= gdis_iso_field_value(model, settings, sample);

  gdis_iso_vec3_copy(sample, point);
  sample[2] = point[2] + delta;
  gradient_out[2] = gdis_iso_field_value(model, settings, sample);
  sample[2] = point[2] - delta;
  gradient_out[2] -= gdis_iso_field_value(model, settings, sample);
}

static gboolean
gdis_iso_choose_grid(const GdisModel *model,
                     const GdisIsosurfaceSettings *settings,
                     gdouble min_bound[3],
                     gdouble max_bound[3],
                     guint dims[3],
                     gdouble *step_out,
                     GError **error)
{
  gdouble padding;
  gdouble step;
  guint64 cube_budget;
  gboolean use_selected_bounds;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(settings != NULL, FALSE);
  g_return_val_if_fail(min_bound != NULL, FALSE);
  g_return_val_if_fail(max_bound != NULL, FALSE);
  g_return_val_if_fail(dims != NULL, FALSE);
  g_return_val_if_fail(step_out != NULL, FALSE);

  padding = 2.5 + settings->blur * 2.0;
  if (settings->mode == GDIS_ISOSURFACE_MODE_PROMOLECULE)
    padding += 1.5;
  else if (settings->mode == GDIS_ISOSURFACE_MODE_HIRSHFELD)
    padding += 3.0;
  else if (settings->mode == GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY)
    padding += 2.0;
  use_selected_bounds = gdis_iso_mode_uses_selection(settings->mode) &&
                        settings->selected_atoms &&
                        settings->selected_atoms->len > 0u;

  for (guint axis = 0; axis < 3; axis++)
    {
      min_bound[axis] = G_MAXDOUBLE;
      max_bound[axis] = -G_MAXDOUBLE;
    }

  for (guint i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      const GdisElementInfo *element;
      gdouble radius;

      if (use_selected_bounds &&
          !gdis_iso_atom_array_contains(settings->selected_atoms, i))
        continue;

      atom = g_ptr_array_index(model->atoms, i);
      element = gdis_element_lookup(atom->element);
      radius = (element ? element->vdw_radius : 1.70) + padding;

      for (guint axis = 0; axis < 3; axis++)
        {
          min_bound[axis] = MIN(min_bound[axis], atom->position[axis] - radius);
          max_bound[axis] = MAX(max_bound[axis], atom->position[axis] + radius);
        }
    }

  step = MAX(settings->grid_size, 0.12);
  cube_budget = 900000;

  for (;;)
    {
      guint64 total_cubes;

      total_cubes = 1u;
      for (guint axis = 0; axis < 3; axis++)
        {
          dims[axis] = MAX(2u, (guint) ceil((max_bound[axis] - min_bound[axis]) / step));
          total_cubes *= dims[axis];
        }

      if (total_cubes <= cube_budget)
        break;

      step *= 1.15;
      if (step > 2.5)
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_FAILED,
                      "Iso-surface grid would be too large for the current model.");
          return FALSE;
        }
    }

  *step_out = step;
  return TRUE;
}

static gboolean
gdis_iso_append_triangle(GArray *triangles,
                         const GdisModel *model,
                         const GdisIsosurfaceSettings *settings,
                         gdouble gradient_delta,
                         const gdouble a[3],
                         const gdouble b[3],
                         const gdouble c[3])
{
  GdisIsoTriangle triangle;
  gdouble ab[3];
  gdouble ac[3];
  gdouble face_normal[3];
  gdouble average_normal[3];
  gdouble area;

  if (triangles->len >= 120000)
    return FALSE;

  gdis_iso_vec3_copy(triangle.vertices[0].position, a);
  gdis_iso_vec3_copy(triangle.vertices[1].position, b);
  gdis_iso_vec3_copy(triangle.vertices[2].position, c);

  for (guint i = 0; i < 3; i++)
    {
      gdis_iso_field_gradient(model,
                              settings,
                              gradient_delta,
                              triangle.vertices[i].position,
                              triangle.vertices[i].normal);
      if (!gdis_iso_vec3_normalize(triangle.vertices[i].normal))
        gdis_iso_vec3_set(triangle.vertices[i].normal, 0.0, 0.0, 1.0);

      triangle.vertices[i].normal[0] = -triangle.vertices[i].normal[0];
      triangle.vertices[i].normal[1] = -triangle.vertices[i].normal[1];
      triangle.vertices[i].normal[2] = -triangle.vertices[i].normal[2];
    }

  gdis_iso_vec3_subtract(b, a, ab);
  gdis_iso_vec3_subtract(c, a, ac);
  gdis_iso_vec3_cross(ab, ac, face_normal);
  area = gdis_iso_vec3_length(face_normal);
  if (area <= 1.0e-8)
    return TRUE;

  average_normal[0] = triangle.vertices[0].normal[0] +
                      triangle.vertices[1].normal[0] +
                      triangle.vertices[2].normal[0];
  average_normal[1] = triangle.vertices[0].normal[1] +
                      triangle.vertices[1].normal[1] +
                      triangle.vertices[2].normal[1];
  average_normal[2] = triangle.vertices[0].normal[2] +
                      triangle.vertices[1].normal[2] +
                      triangle.vertices[2].normal[2];

  if (gdis_iso_vec3_dot(face_normal, average_normal) < 0.0)
    {
      GdisIsoVertex swap;

      swap = triangle.vertices[1];
      triangle.vertices[1] = triangle.vertices[2];
      triangle.vertices[2] = swap;
    }

  g_array_append_val(triangles, triangle);
  return TRUE;
}

static void
gdis_iso_interpolate_vertex(const GdisIsoSample *inside,
                            const GdisIsoSample *outside,
                            gdouble out[3])
{
  gdouble fraction;
  gdouble denominator;

  denominator = inside->value - outside->value;
  if (fabs(denominator) < 1.0e-12)
    fraction = 0.5;
  else
    fraction = inside->value / denominator;
  fraction = CLAMP(fraction, 0.0, 1.0);

  out[0] = inside->point[0] + (outside->point[0] - inside->point[0]) * fraction;
  out[1] = inside->point[1] + (outside->point[1] - inside->point[1]) * fraction;
  out[2] = inside->point[2] + (outside->point[2] - inside->point[2]) * fraction;
}

static gboolean
gdis_iso_polygonize_tetra(const GdisIsoSample samples[4],
                          GArray *triangles,
                          const GdisModel *model,
                          const GdisIsosurfaceSettings *settings,
                          gdouble gradient_delta)
{
  guint inside_indices[4];
  guint outside_indices[4];
  guint inside_count;
  guint outside_count;

  inside_count = 0;
  outside_count = 0;
  for (guint i = 0; i < 4; i++)
    {
      if (samples[i].value >= 0.0)
        inside_indices[inside_count++] = i;
      else
        outside_indices[outside_count++] = i;
    }

  if (inside_count == 0 || inside_count == 4)
    return TRUE;

  if (inside_count == 1 || inside_count == 3)
    {
      gdouble v0[3];
      gdouble v1[3];
      gdouble v2[3];

      if (inside_count == 1)
        {
          gdis_iso_interpolate_vertex(&samples[inside_indices[0]], &samples[outside_indices[0]], v0);
          gdis_iso_interpolate_vertex(&samples[inside_indices[0]], &samples[outside_indices[1]], v1);
          gdis_iso_interpolate_vertex(&samples[inside_indices[0]], &samples[outside_indices[2]], v2);
        }
      else
        {
          gdis_iso_interpolate_vertex(&samples[inside_indices[0]], &samples[outside_indices[0]], v0);
          gdis_iso_interpolate_vertex(&samples[inside_indices[1]], &samples[outside_indices[0]], v1);
          gdis_iso_interpolate_vertex(&samples[inside_indices[2]], &samples[outside_indices[0]], v2);
        }

      return gdis_iso_append_triangle(triangles,
                                      model,
                                      settings,
                                      gradient_delta,
                                      v0,
                                      v1,
                                      v2);
    }

  if (inside_count == 2)
    {
      gdouble v0[3];
      gdouble v1[3];
      gdouble v2[3];
      gdouble v3[3];

      gdis_iso_interpolate_vertex(&samples[inside_indices[0]], &samples[outside_indices[0]], v0);
      gdis_iso_interpolate_vertex(&samples[inside_indices[0]], &samples[outside_indices[1]], v1);
      gdis_iso_interpolate_vertex(&samples[inside_indices[1]], &samples[outside_indices[0]], v2);
      gdis_iso_interpolate_vertex(&samples[inside_indices[1]], &samples[outside_indices[1]], v3);

      if (!gdis_iso_append_triangle(triangles,
                                    model,
                                    settings,
                                    gradient_delta,
                                    v0,
                                    v1,
                                    v2))
        return FALSE;

      return gdis_iso_append_triangle(triangles,
                                      model,
                                      settings,
                                      gradient_delta,
                                      v1,
                                      v3,
                                      v2);
    }

  return TRUE;
}

static void
gdis_iso_compute_area_volume(const GArray *triangles,
                             gdouble *area_out,
                             gdouble *volume_out)
{
  gdouble area;
  gdouble volume;

  area = 0.0;
  volume = 0.0;
  if (triangles)
    {
      for (guint i = 0; i < triangles->len; i++)
        {
          const GdisIsoTriangle *triangle;
          gdouble ab[3];
          gdouble ac[3];
          gdouble cross[3];

          triangle = &g_array_index(triangles, GdisIsoTriangle, i);
          gdis_iso_vec3_subtract(triangle->vertices[1].position,
                                 triangle->vertices[0].position,
                                 ab);
          gdis_iso_vec3_subtract(triangle->vertices[2].position,
                                 triangle->vertices[0].position,
                                 ac);
          gdis_iso_vec3_cross(ab, ac, cross);
          area += 0.5 * gdis_iso_vec3_length(cross);
          volume += gdis_iso_vec3_dot(triangle->vertices[0].position, cross) / 6.0;
        }
    }

  if (area_out)
    *area_out = area;
  if (volume_out)
    *volume_out = fabs(volume);
}
