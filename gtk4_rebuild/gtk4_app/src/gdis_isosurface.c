#include "gdis_isosurface.h"
#include "gdis_elements.h"

#include <math.h>
#include <string.h>

typedef struct
{
  gdouble point[3];
  gdouble value;
} GdisIsoSample;


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
  settings->grid_size = 0.40;
  settings->blur = 0.35;
  settings->isovalue = 0.08;
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

gboolean
gdis_isosurface_mode_supported(GdisIsosurfaceMode mode)
{
  return mode == GDIS_ISOSURFACE_MODE_MOLECULAR ||
         mode == GDIS_ISOSURFACE_MODE_PROMOLECULE ||
         mode == GDIS_ISOSURFACE_MODE_HIRSHFELD ||
         mode == GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY;
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

  g_return_val_if_fail(surface_out != NULL, FALSE);
  *surface_out = NULL;

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
                  "%s is not restored yet in the GTK4 rebuild.",
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
  surface->grid_size = step;
  surface->blur = settings->blur;
  surface->isovalue = settings->isovalue;
  surface->triangles = triangles;
  if (settings->mode == GDIS_ISOSURFACE_MODE_PROMOLECULE)
    {
      surface->color_rgba[0] = 0.87;
      surface->color_rgba[1] = 0.72;
      surface->color_rgba[2] = 0.26;
      surface->color_rgba[3] = 0.48;
    }
  else if (settings->mode == GDIS_ISOSURFACE_MODE_HIRSHFELD)
    {
      surface->color_rgba[0] = 0.88;
      surface->color_rgba[1] = 0.46;
      surface->color_rgba[2] = 0.66;
      surface->color_rgba[3] = 0.42;
    }
  else if (settings->mode == GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY)
    {
      surface->color_rgba[0] = 0.40;
      surface->color_rgba[1] = 0.86;
      surface->color_rgba[2] = 0.54;
      surface->color_rgba[3] = 0.40;
    }
  else
    {
      surface->color_rgba[0] = 0.32;
      surface->color_rgba[1] = 0.78;
      surface->color_rgba[2] = 0.92;
      surface->color_rgba[3] = 0.42;
    }
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
    "Triangles: %u\n"
    "Estimated area: %.3f A^2\n"
    "Estimated volume: %.3f A^3\n"
    "Grid step: %.3f A\n"
    "Blur: %.3f\n"
    "Isovalue: %.4f\n",
    gdis_isosurface_mode_label(settings->mode),
    model->basename ? model->basename : "(model)",
    selection_summary,
    triangles->len,
    area,
    volume,
    step,
    settings->blur,
    settings->isovalue);
  if (settings->mode == GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY)
    {
      gchar *with_note;

      with_note = g_strconcat(surface->summary,
                              "Mode note: GTK4 currently restores this as an analytic whole-model density approximation, not a volumetric grid import.\n",
                              NULL);
      g_free(surface->summary);
      surface->summary = with_note;
    }

  *surface_out = surface;
  return TRUE;
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
