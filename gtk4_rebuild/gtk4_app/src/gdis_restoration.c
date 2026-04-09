#include "gdis_restoration.h"

#include <glib/gstdio.h>
#include <math.h>
#include <string.h>

#include "gdis_reports.h"

static GArray *gdis_build_scope_array(const GdisModel *model,
                                      const GArray *selected_atoms,
                                      gboolean use_selection,
                                      GError **error);
static gboolean gdis_scope_contains(const GArray *scope, guint atom_index);
static gint gdis_scope_find_index(const GArray *scope, guint atom_index);
static gboolean gdis_atoms_are_bonded(const GdisModel *model,
                                      guint atom_index_a,
                                      guint atom_index_b);
static gdouble gdis_distance_sq(const gdouble a[3], const gdouble b[3]);
static void gdis_vec3_copy(gdouble dest[3], const gdouble src[3]);
static void gdis_vec3_set(gdouble vector[3], gdouble x, gdouble y, gdouble z);
static void gdis_vec3_subtract(const gdouble a[3], const gdouble b[3], gdouble out[3]);
static void gdis_vec3_add_scaled(gdouble vector[3], const gdouble addend[3], gdouble scale);
static gdouble gdis_vec3_dot(const gdouble a[3], const gdouble b[3]);
static void gdis_vec3_cross(const gdouble a[3], const gdouble b[3], gdouble out[3]);
static gdouble gdis_vec3_length(const gdouble vector[3]);
static gboolean gdis_vec3_normalize(gdouble vector[3]);
static void gdis_rotate_vector_axis_angle(const gdouble vector[3],
                                          const gdouble axis[3],
                                          gdouble angle,
                                          gdouble out[3]);
static void gdis_rotate_euler_xyz(const gdouble vector[3],
                                  gdouble angle_x,
                                  gdouble angle_y,
                                  gdouble angle_z,
                                  gdouble out[3]);
static gboolean gdis_measure_value(const GdisModel *model,
                                   GdisMeasureMode mode,
                                   const guint *atom_indices,
                                   guint count,
                                   gdouble *value_out);
static guint gdis_pick_distance_reference(const GdisModel *model,
                                          const GArray *scope,
                                          guint scope_position);
static guint gdis_pick_angle_reference(const GdisModel *model,
                                       const GArray *scope,
                                       guint scope_position,
                                       guint distance_ref);
static guint gdis_pick_torsion_reference(const GdisModel *model,
                                         const GArray *scope,
                                         guint scope_position,
                                         guint distance_ref,
                                         guint angle_ref);
static gboolean gdis_write_xyz_subset(FILE *stream,
                                      const GdisModel *model,
                                      const GArray *scope,
                                      const gchar *title);
static gboolean gdis_write_xyz_pose(FILE *stream,
                                    const GdisModel *model,
                                    const GArray *scope,
                                    const gdouble *positions,
                                    const gchar *title);
static void gdis_sort_scope_by_connectivity(const GdisModel *model, GArray *scope);
static void gdis_compute_scope_centroid(const GdisModel *model,
                                        const GArray *scope,
                                        gdouble centroid[3]);
static gdouble gdis_measure_scope_radius(const GdisModel *model,
                                         const GArray *scope,
                                         const gdouble centroid[3]);
static void gdis_pick_orthogonal_axis(const gdouble axis[3], gdouble orthogonal[3]);
static gboolean gdis_zmatrix_place_with_two_refs(const gdouble ref_b[3],
                                                 const gdouble ref_c[3],
                                                 gdouble distance,
                                                 gdouble angle_degrees,
                                                 gdouble position_out[3]);
static gboolean gdis_zmatrix_place_with_three_refs(const gdouble ref_b[3],
                                                   const gdouble ref_c[3],
                                                   const gdouble ref_d[3],
                                                   gdouble distance,
                                                   gdouble angle_degrees,
                                                   gdouble torsion_degrees,
                                                   gdouble position_out[3]);
static void gdis_zmatrix_build_world_basis(const gdouble *original_positions,
                                           guint scope_len,
                                           gdouble origin[3],
                                           gdouble axis_x[3],
                                           gdouble axis_y[3],
                                           gdouble axis_z[3]);

void
gdis_dislocation_settings_init(GdisDislocationSettings *settings)
{
  g_return_if_fail(settings != NULL);

  memset(settings, 0, sizeof(*settings));
  settings->type = GDIS_DISLOCATION_TYPE_SCREW;
  settings->line_direction[2] = 1.0;
  settings->burgers_vector[2] = 1.0;
  settings->poisson_ratio = 0.33;
}

const char *
gdis_dislocation_type_label(GdisDislocationType type)
{
  switch (type)
    {
    case GDIS_DISLOCATION_TYPE_EDGE:
      return "Edge (experimental)";
    case GDIS_DISLOCATION_TYPE_SCREW:
    default:
      return "Screw";
    }
}

gboolean
gdis_model_apply_dislocation(GdisModel *model,
                             const GdisDislocationSettings *settings,
                             gchar **summary_out,
                             GError **error)
{
  GArray *scope;
  gdouble axis[3];
  gdouble burgers[3];
  gdouble origin[3];
  gdouble basis_u[3];
  gdouble basis_v[3];
  gdouble basis_ref[3];
  guint affected;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(settings != NULL, FALSE);

  scope = gdis_build_scope_array(model, settings->selected_atoms, TRUE, error);
  if (!scope)
    return FALSE;

  gdis_vec3_copy(axis, settings->line_direction);
  if (!gdis_vec3_normalize(axis))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Dislocation line direction must not be the zero vector.");
      g_array_free(scope, TRUE);
      return FALSE;
    }

  gdis_vec3_copy(burgers, settings->burgers_vector);
  gdis_vec3_copy(origin, settings->origin);

  if (fabs(axis[2]) < 0.9)
    gdis_vec3_set(basis_ref, 0.0, 0.0, 1.0);
  else
    gdis_vec3_set(basis_ref, 0.0, 1.0, 0.0);
  gdis_vec3_cross(axis, basis_ref, basis_u);
  if (!gdis_vec3_normalize(basis_u))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not build a stable dislocation basis from the line direction.");
      g_array_free(scope, TRUE);
      return FALSE;
    }
  gdis_vec3_cross(axis, basis_u, basis_v);
  gdis_vec3_normalize(basis_v);

  affected = 0u;
  gdis_model_discard_frames(model);
  for (guint i = 0; i < scope->len; i++)
    {
      guint atom_index;
      GdisAtom *atom;
      gdouble rel[3];
      gdouble parallel_component;
      gdouble radial_vector[3];
      gdouble x;
      gdouble y;
      gdouble r2;

      atom_index = g_array_index(scope, guint, i);
      if (atom_index >= model->atoms->len)
        continue;

      atom = g_ptr_array_index(model->atoms, atom_index);
      gdis_vec3_subtract(atom->position, origin, rel);
      parallel_component = gdis_vec3_dot(rel, axis);
      gdis_vec3_copy(radial_vector, rel);
      gdis_vec3_add_scaled(radial_vector, axis, -parallel_component);
      x = gdis_vec3_dot(radial_vector, basis_u);
      y = gdis_vec3_dot(radial_vector, basis_v);
      r2 = x * x + y * y;

      if (settings->radius > 0.0 && r2 > settings->radius * settings->radius)
        continue;
      if (r2 < 1.0e-12)
        continue;

      if (settings->type == GDIS_DISLOCATION_TYPE_SCREW)
        {
          gdouble burgers_parallel;
          gdouble phi;
          gdouble displacement;

          burgers_parallel = gdis_vec3_dot(burgers, axis);
          if (fabs(burgers_parallel) < 1.0e-10)
            burgers_parallel = gdis_vec3_length(burgers);
          phi = atan2(y, x);
          displacement = (burgers_parallel / (2.0 * G_PI)) * phi;
          gdis_vec3_add_scaled(atom->position, axis, displacement);
        }
      else
        {
          gdouble burgers_in_plane[3];
          gdouble burgers_mag;
          gdouble nu;
          gdouble ux;
          gdouble uy;

          gdis_vec3_copy(burgers_in_plane, burgers);
          gdis_vec3_add_scaled(burgers_in_plane, axis, -gdis_vec3_dot(burgers, axis));
          burgers_mag = gdis_vec3_length(burgers_in_plane);
          if (burgers_mag < 1.0e-10)
            burgers_mag = 1.0;
          if (!gdis_vec3_normalize(burgers_in_plane))
            gdis_vec3_copy(burgers_in_plane, basis_u);

          gdis_vec3_copy(basis_u, burgers_in_plane);
          gdis_vec3_cross(axis, basis_u, basis_v);
          gdis_vec3_normalize(basis_v);
          x = gdis_vec3_dot(radial_vector, basis_u);
          y = gdis_vec3_dot(radial_vector, basis_v);
          r2 = MAX(x * x + y * y, 1.0e-10);
          nu = CLAMP(settings->poisson_ratio, 0.01, 0.49);

          ux = (burgers_mag / (2.0 * G_PI)) *
               (atan2(y, x) + (x * y) / (2.0 * (1.0 - nu) * r2));
          uy = -(burgers_mag / (2.0 * G_PI)) *
               ((((1.0 - 2.0 * nu) / (4.0 * (1.0 - nu))) * log(r2)) +
                ((x * x - y * y) / (4.0 * (1.0 - nu) * r2)));
          gdis_vec3_add_scaled(atom->position, basis_u, ux);
          gdis_vec3_add_scaled(atom->position, basis_v, uy);
        }

      affected++;
    }

  if (model->explicit_bond_count == 0u)
    gdis_model_reset_inferred_bonds(model);

  if (summary_out)
    {
      *summary_out = g_strdup_printf(
        "%s dislocation applied.\n"
        "Affected atoms: %u of %u\n"
        "Line direction: %.4f %.4f %.4f\n"
        "Burgers vector: %.4f %.4f %.4f\n"
        "Origin: %.4f %.4f %.4f\n"
        "Radius filter: %s\n\n"
        "%s",
        gdis_dislocation_type_label(settings->type),
        affected,
        scope->len,
        axis[0], axis[1], axis[2],
        burgers[0], burgers[1], burgers[2],
        origin[0], origin[1], origin[2],
        settings->radius > 0.0 ? "enabled" : "disabled",
        settings->type == GDIS_DISLOCATION_TYPE_EDGE ?
        "The edge mode is a GTK4-native isotropic approximation." :
        "The screw mode is a GTK4-native displacement-field restoration.");
    }

  g_array_free(scope, TRUE);
  return TRUE;
}

void
gdis_docking_settings_init(GdisDockingSettings *settings)
{
  g_return_if_fail(settings != NULL);

  memset(settings, 0, sizeof(*settings));
  settings->grid_counts[0] = 3u;
  settings->grid_counts[1] = 3u;
  settings->grid_counts[2] = 1u;
  settings->rotation_counts[0] = 1u;
  settings->rotation_counts[1] = 1u;
  settings->rotation_counts[2] = 4u;
  settings->translation_span[0] = 6.0;
  settings->translation_span[1] = 6.0;
  settings->translation_span[2] = 2.0;
  settings->preview_limit = 24u;
}

gboolean
gdis_generate_docking_project(const GdisModel *model,
                              const GdisDockingSettings *settings,
                              gchar **summary_out,
                              GError **error)
{
  GArray *scope;
  gchar *project_name;
  gchar *project_dir;
  gchar *fragment_path;
  gchar *poses_path;
  gchar *config_path;
  FILE *fragment_stream;
  FILE *poses_stream;
  FILE *config_stream;
  gdouble centroid[3];
  guint64 total_poses;
  guint64 preview_written;
  gboolean success;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(settings != NULL, FALSE);

  scope = gdis_build_scope_array(model, settings->selected_atoms, TRUE, error);
  if (!scope)
    return FALSE;

  if (!settings->output_dir || settings->output_dir[0] == '\0')
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Docking project generation needs an output directory.");
      g_array_free(scope, TRUE);
      return FALSE;
    }

  project_name = g_strdup((settings->project_name && settings->project_name[0] != '\0')
                          ? settings->project_name
                          : "dock_project");
  project_dir = g_build_filename(settings->output_dir, project_name, NULL);
  fragment_path = g_build_filename(project_dir, "fragment.xyz", NULL);
  poses_path = g_build_filename(project_dir, "poses.xyz", NULL);
  config_path = g_build_filename(project_dir, "project.txt", NULL);
  fragment_stream = NULL;
  poses_stream = NULL;
  config_stream = NULL;
  success = FALSE;

  if (g_mkdir_with_parents(project_dir, 0755) != 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not create docking project directory '%s'.",
                  project_dir);
      goto docking_done;
    }

  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;
  for (guint i = 0; i < scope->len; i++)
    {
      const GdisAtom *atom;

      atom = g_ptr_array_index(model->atoms, g_array_index(scope, guint, i));
      centroid[0] += atom->position[0];
      centroid[1] += atom->position[1];
      centroid[2] += atom->position[2];
    }
  centroid[0] /= (gdouble) scope->len;
  centroid[1] /= (gdouble) scope->len;
  centroid[2] /= (gdouble) scope->len;

  total_poses = (guint64) MAX(settings->grid_counts[0], 1u) *
                (guint64) MAX(settings->grid_counts[1], 1u) *
                (guint64) MAX(settings->grid_counts[2], 1u) *
                (guint64) MAX(settings->rotation_counts[0], 1u) *
                (guint64) MAX(settings->rotation_counts[1], 1u) *
                (guint64) MAX(settings->rotation_counts[2], 1u);

  fragment_stream = fopen(fragment_path, "w");
  if (!fragment_stream || !gdis_write_xyz_subset(fragment_stream, model, scope, "Docking fragment"))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not write the docking fragment file.");
      goto docking_done;
    }
  fclose(fragment_stream);
  fragment_stream = NULL;

  poses_stream = fopen(poses_path, "w");
  if (!poses_stream)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not open the docking pose file for writing.");
      goto docking_done;
    }

  preview_written = 0u;
  for (guint gx = 0; gx < MAX(settings->grid_counts[0], 1u); gx++)
    {
      for (guint gy = 0; gy < MAX(settings->grid_counts[1], 1u); gy++)
        {
          for (guint gz = 0; gz < MAX(settings->grid_counts[2], 1u); gz++)
            {
              gdouble translate[3];

              translate[0] = (MAX(settings->grid_counts[0], 1u) == 1u)
                             ? 0.0
                             : (-0.5 * settings->translation_span[0]) +
                               settings->translation_span[0] *
                               ((gdouble) gx / (gdouble) (MAX(settings->grid_counts[0], 1u) - 1u));
              translate[1] = (MAX(settings->grid_counts[1], 1u) == 1u)
                             ? 0.0
                             : (-0.5 * settings->translation_span[1]) +
                               settings->translation_span[1] *
                               ((gdouble) gy / (gdouble) (MAX(settings->grid_counts[1], 1u) - 1u));
              translate[2] = (MAX(settings->grid_counts[2], 1u) == 1u)
                             ? 0.0
                             : (-0.5 * settings->translation_span[2]) +
                               settings->translation_span[2] *
                               ((gdouble) gz / (gdouble) (MAX(settings->grid_counts[2], 1u) - 1u));

              for (guint rx = 0; rx < MAX(settings->rotation_counts[0], 1u); rx++)
                {
                  for (guint ry = 0; ry < MAX(settings->rotation_counts[1], 1u); ry++)
                    {
                      for (guint rz = 0; rz < MAX(settings->rotation_counts[2], 1u); rz++)
                        {
                          gdouble angle_x;
                          gdouble angle_y;
                          gdouble angle_z;
                          g_autofree gchar *title = NULL;
                          gdouble *positions;

                          if (preview_written >= MAX(settings->preview_limit, 1u))
                            continue;

                          angle_x = 2.0 * G_PI *
                                    ((gdouble) rx / (gdouble) MAX(settings->rotation_counts[0], 1u));
                          angle_y = 2.0 * G_PI *
                                    ((gdouble) ry / (gdouble) MAX(settings->rotation_counts[1], 1u));
                          angle_z = 2.0 * G_PI *
                                    ((gdouble) rz / (gdouble) MAX(settings->rotation_counts[2], 1u));
                          positions = g_new(gdouble, scope->len * 3u);

                          for (guint atom_slot = 0; atom_slot < scope->len; atom_slot++)
                            {
                              const GdisAtom *atom;
                              gdouble rel[3];
                              gdouble rotated[3];

                              atom = g_ptr_array_index(model->atoms,
                                                       g_array_index(scope, guint, atom_slot));
                              gdis_vec3_subtract(atom->position, centroid, rel);
                              gdis_rotate_euler_xyz(rel, angle_x, angle_y, angle_z, rotated);
                              positions[3u * atom_slot] = centroid[0] + translate[0] + rotated[0];
                              positions[3u * atom_slot + 1u] = centroid[1] + translate[1] + rotated[1];
                              positions[3u * atom_slot + 2u] = centroid[2] + translate[2] + rotated[2];
                            }

                          title = g_strdup_printf("pose %06u  shift=(%.3f %.3f %.3f)  rot=(%.1f %.1f %.1f)",
                                                  (guint) preview_written + 1u,
                                                  translate[0], translate[1], translate[2],
                                                  angle_x * (180.0 / G_PI),
                                                  angle_y * (180.0 / G_PI),
                                                  angle_z * (180.0 / G_PI));
                          if (!gdis_write_xyz_pose(poses_stream, model, scope, positions, title))
                            {
                              g_free(positions);
                              g_set_error(error,
                                          GDIS_MODEL_ERROR,
                                          GDIS_MODEL_ERROR_IO,
                                          "Could not write the docking pose ensemble.");
                              goto docking_done;
                            }

                          g_free(positions);
                          preview_written++;
                        }
                    }
                }
            }
        }
    }
  fclose(poses_stream);
  poses_stream = NULL;

  config_stream = fopen(config_path, "w");
  if (!config_stream)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not write the docking project summary.");
      goto docking_done;
    }
  fprintf(config_stream,
          "Docking project: %s\n"
          "Source model: %s\n"
          "Selected atoms: %u\n"
          "Grid counts: %u %u %u\n"
          "Rotation counts: %u %u %u\n"
          "Translation span: %.4f %.4f %.4f\n"
          "Total sampled poses: %llu\n"
          "Preview poses written: %llu\n",
          project_name,
          model->basename ? model->basename : "model",
          scope->len,
          MAX(settings->grid_counts[0], 1u),
          MAX(settings->grid_counts[1], 1u),
          MAX(settings->grid_counts[2], 1u),
          MAX(settings->rotation_counts[0], 1u),
          MAX(settings->rotation_counts[1], 1u),
          MAX(settings->rotation_counts[2], 1u),
          settings->translation_span[0],
          settings->translation_span[1],
          settings->translation_span[2],
          (unsigned long long) total_poses,
          (unsigned long long) preview_written);
  fclose(config_stream);
  config_stream = NULL;

  if (summary_out)
    {
      *summary_out = g_strdup_printf(
        "Docking project generated.\n"
        "Directory: %s\n"
        "Selected atoms: %u\n"
        "Total sampled poses: %llu\n"
        "Preview poses written: %llu\n\n"
        "Files:\n"
        "  fragment.xyz\n"
        "  poses.xyz\n"
        "  project.txt\n\n"
        "GTK4 restores the legacy rigid-fragment pose sampling workflow here without depending on the old GULP batch runner.",
        project_dir,
        scope->len,
        (unsigned long long) total_poses,
        (unsigned long long) preview_written);
    }

  success = TRUE;

docking_done:
  if (fragment_stream)
    fclose(fragment_stream);
  if (poses_stream)
    fclose(poses_stream);
  if (config_stream)
    fclose(config_stream);
  g_free(project_name);
  g_free(project_dir);
  g_free(fragment_path);
  g_free(poses_path);
  g_free(config_path);
  g_array_free(scope, TRUE);
  return success;
}

void
gdis_mdi_settings_init(GdisMdiSettings *settings)
{
  g_return_if_fail(settings != NULL);

  memset(settings, 0, sizeof(*settings));
  settings->box_dim = 4u;
  settings->random_rotate = TRUE;
}

gboolean
gdis_generate_mdi_model(const GPtrArray *source_models,
                        const GdisMdiSettings *settings,
                        GdisModel **model_out,
                        gchar **summary_out,
                        GError **error)
{
  GdisModel *model;
  guint component_count;
  guint solvent_index;
  guint box_dim;
  guint64 total_sites;
  guint64 solute_total;
  guint64 solvent_total;
  gdouble lattice_sep;
  gdouble box_len;
  guint *assignment;
  GArray *occupied_sites;
  guint copied_atoms;
  guint copied_molecules;

  g_return_val_if_fail(source_models != NULL, FALSE);
  g_return_val_if_fail(settings != NULL, FALSE);
  g_return_val_if_fail(model_out != NULL, FALSE);

  *model_out = NULL;
  component_count = source_models->len;
  if (component_count < 2u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Load solvent and solute molecules first.");
      return FALSE;
    }

  if (!settings->component_counts || settings->component_count != component_count)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "MD initializer component counts are out of sync with the loaded model list.");
      return FALSE;
    }

  box_dim = MAX(settings->box_dim, 3u);
  solvent_index = MIN(settings->solvent_index, component_count - 1u);
  total_sites = (guint64) box_dim * (guint64) box_dim * (guint64) box_dim;
  solute_total = 0u;
  lattice_sep = 0.0;

  for (guint i = 0; i < component_count; i++)
    {
      const GdisModel *source_model;
      guint requested;
      GArray *scope;
      gdouble centroid[3];
      gdouble radius;

      source_model = g_ptr_array_index((GPtrArray *) source_models, i);
      if (!source_model || !source_model->atoms || source_model->atoms->len == 0u)
        continue;

      requested = (i == solvent_index)
        ? 1u
        : settings->component_counts[i];
      if (i != solvent_index)
        solute_total += requested;
      if (requested == 0u)
        continue;

      scope = g_array_sized_new(FALSE, FALSE, sizeof(guint), source_model->atoms->len);
      for (guint atom_index = 0; atom_index < source_model->atoms->len; atom_index++)
        g_array_append_val(scope, atom_index);
      gdis_compute_scope_centroid(source_model, scope, centroid);
      radius = gdis_measure_scope_radius(source_model, scope, centroid);
      lattice_sep = MAX(lattice_sep, 2.0 * radius + 1.0);
      g_array_free(scope, TRUE);
    }

  if (solute_total >= total_sites)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Too many solute molecules were requested for a %u x %u x %u box.",
                  box_dim,
                  box_dim,
                  box_dim);
      return FALSE;
    }

  solvent_total = total_sites - solute_total;
  if (solvent_total < 9u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Too many solute molecules were requested. The solvent needs at least 9 remaining lattice sites.");
      return FALSE;
    }

  if (!(lattice_sep > 0.0))
    lattice_sep = 4.0;
  box_len = lattice_sep * (gdouble) box_dim;
  assignment = g_new0(guint, (gsize) total_sites);
  occupied_sites = g_array_new(FALSE, FALSE, sizeof(guint));
  for (guint64 pos = 0; pos < total_sites; pos++)
    assignment[pos] = solvent_index;

  for (guint component_index = 0; component_index < component_count; component_index++)
    {
      guint requested;

      if (component_index == solvent_index)
        continue;

      requested = settings->component_counts[component_index];
      for (guint placed = 0; placed < requested; placed++)
        {
          guint chosen_pos = 0u;
          gint best_min_distance = G_MININT;
          GArray *candidates;

          candidates = g_array_new(FALSE, FALSE, sizeof(guint));
          for (guint iz = 0; iz < box_dim; iz++)
            {
              for (guint iy = 0; iy < box_dim; iy++)
                {
                  for (guint ix = 0; ix < box_dim; ix++)
                    {
                      guint pos;
                      gint min_distance_sq;

                      pos = ix + box_dim * (iy + box_dim * iz);
                      min_distance_sq = G_MAXINT;
                      if (occupied_sites->len == 0u)
                        min_distance_sq = (gint) total_sites;
                      else
                        {
                          for (guint s = 0; s < occupied_sites->len; s++)
                            {
                              guint occupied_pos;
                              guint occupied_ix;
                              guint occupied_iy;
                              guint occupied_iz;
                              gint dx;
                              gint dy;
                              gint dz;
                              gint distance_sq;

                              occupied_pos = g_array_index(occupied_sites, guint, s);
                              occupied_ix = occupied_pos % box_dim;
                              occupied_iy = (occupied_pos / box_dim) % box_dim;
                              occupied_iz = occupied_pos / (box_dim * box_dim);
                              dx = abs((gint) occupied_ix - (gint) ix);
                              dy = abs((gint) occupied_iy - (gint) iy);
                              dz = abs((gint) occupied_iz - (gint) iz);
                              if (dx > (gint) box_dim / 2)
                                dx = (gint) box_dim - dx;
                              if (dy > (gint) box_dim / 2)
                                dy = (gint) box_dim - dy;
                              if (dz > (gint) box_dim / 2)
                                dz = (gint) box_dim - dz;
                              distance_sq = dx * dx + dy * dy + dz * dz;
                              min_distance_sq = MIN(min_distance_sq, distance_sq);
                            }
                        }

                      if (min_distance_sq > best_min_distance)
                        {
                          best_min_distance = min_distance_sq;
                          g_array_set_size(candidates, 0u);
                          g_array_append_val(candidates, pos);
                        }
                      else if (min_distance_sq == best_min_distance)
                        {
                          g_array_append_val(candidates, pos);
                        }
                    }
                }
            }

          if (candidates->len == 0u)
            {
              g_array_free(candidates, TRUE);
              g_array_free(occupied_sites, TRUE);
              g_free(assignment);
              g_set_error(error,
                          GDIS_MODEL_ERROR,
                          GDIS_MODEL_ERROR_FAILED,
                          "No candidate lattice sites were available for component placement.");
              return FALSE;
            }

          chosen_pos = g_array_index(candidates,
                                     guint,
                                     g_random_int_range(0, (gint) candidates->len));
          assignment[chosen_pos] = component_index;
          g_array_append_val(occupied_sites, chosen_pos);
          g_array_free(candidates, TRUE);
        }
    }

  model = gdis_model_create_empty("MDI model.cif", GDIS_MODEL_FORMAT_CIF);
  if (!model)
    {
      g_array_free(occupied_sites, TRUE);
      g_free(assignment);
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not allocate the MDI model.");
      return FALSE;
    }

  model->periodic = TRUE;
  model->periodicity = 3u;
  model->cell_lengths[0] = box_len;
  model->cell_lengths[1] = box_len;
  model->cell_lengths[2] = box_len;
  model->cell_angles[0] = 90.0;
  model->cell_angles[1] = 90.0;
  model->cell_angles[2] = 90.0;
  gdis_model_reset_image_limits(model);

  copied_atoms = 0u;
  copied_molecules = 0u;
  for (guint iz = 0; iz < box_dim; iz++)
    {
      for (guint iy = 0; iy < box_dim; iy++)
        {
          for (guint ix = 0; ix < box_dim; ix++)
            {
              guint pos;
              guint component_index;
              GdisModel *source_model;
              GArray *scope;
              gdouble centroid[3];
              gdouble translate[3];
              gdouble angle_x;
              gdouble angle_y;
              gdouble angle_z;

              pos = ix + box_dim * (iy + box_dim * iz);
              component_index = assignment[pos];
              source_model = g_ptr_array_index((GPtrArray *) source_models, component_index);
              if (!source_model || !source_model->atoms || source_model->atoms->len == 0u)
                continue;

              scope = g_array_sized_new(FALSE, FALSE, sizeof(guint), source_model->atoms->len);
              for (guint atom_index = 0; atom_index < source_model->atoms->len; atom_index++)
                g_array_append_val(scope, atom_index);
              gdis_compute_scope_centroid(source_model, scope, centroid);
              g_array_free(scope, TRUE);

              translate[0] = (lattice_sep * 0.5) + lattice_sep * (gdouble) ix;
              translate[1] = (lattice_sep * 0.5) + lattice_sep * (gdouble) iy;
              translate[2] = (lattice_sep * 0.5) + lattice_sep * (gdouble) iz;
              if (settings->random_rotate)
                {
                  angle_x = g_random_double_range(0.0, 2.0 * G_PI);
                  angle_y = g_random_double_range(0.0, 2.0 * G_PI);
                  angle_z = g_random_double_range(0.0, 2.0 * G_PI);
                }
              else
                {
                  angle_x = 0.0;
                  angle_y = 0.0;
                  angle_z = 0.0;
                }

              for (guint atom_index = 0; atom_index < source_model->atoms->len; atom_index++)
                {
                  const GdisAtom *source_atom;
                  GdisAtom *dest_atom;
                  gdouble relative[3];
                  gdouble rotated[3];

                  source_atom = g_ptr_array_index(source_model->atoms, atom_index);
                  relative[0] = source_atom->position[0] - centroid[0];
                  relative[1] = source_atom->position[1] - centroid[1];
                  relative[2] = source_atom->position[2] - centroid[2];
                  gdis_rotate_euler_xyz(relative, angle_x, angle_y, angle_z, rotated);

                  if (!gdis_model_add_atom(model,
                                           source_atom->label,
                                           source_atom->element,
                                           source_atom->ff_type,
                                           source_atom->region,
                                           translate[0] + rotated[0],
                                           translate[1] + rotated[1],
                                           translate[2] + rotated[2],
                                           error))
                    {
                      gdis_model_free(model);
                      g_array_free(occupied_sites, TRUE);
                      g_free(assignment);
                      return FALSE;
                    }

                  dest_atom = g_ptr_array_index(model->atoms, model->atoms->len - 1u);
                  dest_atom->occupancy = source_atom->occupancy;
                  dest_atom->charge = source_atom->charge;
                  dest_atom->has_charge = source_atom->has_charge;
                  copied_atoms++;
                }

              copied_molecules++;
            }
        }
    }

  g_array_free(occupied_sites, TRUE);
  g_free(assignment);

  g_free(model->path);
  g_free(model->basename);
  model->path = g_strdup_printf("MDI-model-%ux%ux%u.cif", box_dim, box_dim, box_dim);
  model->basename = g_strdup_printf("MDI model %ux%ux%u", box_dim, box_dim, box_dim);
  g_clear_pointer(&model->title, g_free);
  model->title = g_strdup_printf("MDI model from %s",
                                 ((GdisModel *) g_ptr_array_index((GPtrArray *) source_models, solvent_index))->basename);

  if (summary_out)
    {
      GString *summary;

      summary = g_string_new(NULL);
      g_string_append_printf(summary,
                             "MD initializer complete.\n"
                             "Solvent: %s\n"
                             "Box: %u x %u x %u lattice sites\n"
                             "Lattice spacing: %.3f A\n"
                             "Cell length: %.3f A\n"
                             "Placed molecules: %u\n"
                             "Placed atoms: %u\n\n"
                             "Component counts:\n"
                             "  %s (solvent): %llu\n",
                             ((GdisModel *) g_ptr_array_index((GPtrArray *) source_models, solvent_index))->basename,
                             box_dim,
                             box_dim,
                             box_dim,
                             lattice_sep,
                             box_len,
                             copied_molecules,
                             copied_atoms,
                             ((GdisModel *) g_ptr_array_index((GPtrArray *) source_models, solvent_index))->basename,
                             (unsigned long long) solvent_total);
      for (guint i = 0; i < component_count; i++)
        {
          GdisModel *source_model;

          if (i == solvent_index || settings->component_counts[i] == 0u)
            continue;
          source_model = g_ptr_array_index((GPtrArray *) source_models, i);
          g_string_append_printf(summary,
                                 "  %s: %u\n",
                                 source_model ? source_model->basename : "model",
                                 settings->component_counts[i]);
        }
      g_string_append(summary,
                      "\nGTK4 now restores the legacy lattice-fill MD initializer workflow for structure packing.");
      *summary_out = g_string_free(summary, FALSE);
    }

  *model_out = model;
  return TRUE;
}

gboolean
gdis_build_zmatrix_rows(const GdisModel *model,
                        const GArray *selected_atoms,
                        gboolean use_selection,
                        GArray **scope_out,
                        GArray **rows_out,
                        GError **error)
{
  GArray *scope;
  GArray *rows;

  g_return_val_if_fail(model != NULL, FALSE);

  scope = gdis_build_scope_array(model, selected_atoms, use_selection, error);
  if (!scope)
    return FALSE;

  rows = g_array_sized_new(FALSE, FALSE, sizeof(GdisZmatrixRow), scope->len);
  for (guint i = 0; i < scope->len; i++)
    {
      GdisZmatrixRow row;
      guint pair[2];
      guint trio[3];
      guint quartet[4];

      memset(&row, 0, sizeof(row));
      row.atom_index = g_array_index(scope, guint, i);
      row.distance_ref = gdis_pick_distance_reference(model, scope, i);
      row.angle_ref = (row.distance_ref != G_MAXUINT)
                      ? gdis_pick_angle_reference(model, scope, i, row.distance_ref)
                      : G_MAXUINT;
      row.torsion_ref = (row.angle_ref != G_MAXUINT)
                        ? gdis_pick_torsion_reference(model, scope, i, row.distance_ref, row.angle_ref)
                        : G_MAXUINT;

      if (row.distance_ref != G_MAXUINT)
        {
          pair[0] = row.atom_index;
          pair[1] = row.distance_ref;
          gdis_measure_value(model, GDIS_MEASURE_MODE_DISTANCE, pair, 2u, &row.distance);
        }
      if (row.angle_ref != G_MAXUINT)
        {
          trio[0] = row.atom_index;
          trio[1] = row.distance_ref;
          trio[2] = row.angle_ref;
          gdis_measure_value(model, GDIS_MEASURE_MODE_ANGLE, trio, 3u, &row.angle);
        }
      if (row.torsion_ref != G_MAXUINT)
        {
          quartet[0] = row.atom_index;
          quartet[1] = row.distance_ref;
          quartet[2] = row.angle_ref;
          quartet[3] = row.torsion_ref;
          gdis_measure_value(model, GDIS_MEASURE_MODE_TORSION, quartet, 4u, &row.torsion);
        }

      g_array_append_val(rows, row);
    }

  if (scope_out)
    *scope_out = scope;
  else
    g_array_free(scope, TRUE);

  if (rows_out)
    *rows_out = rows;
  else
    g_array_free(rows, TRUE);

  return TRUE;
}

char *
gdis_format_zmatrix_rows(const GdisModel *model,
                         const GArray *scope,
                         const GArray *rows,
                         gboolean use_selection)
{
  GString *report;

  g_return_val_if_fail(model != NULL, NULL);
  g_return_val_if_fail(scope != NULL, NULL);
  g_return_val_if_fail(rows != NULL, NULL);

  report = g_string_new("");
  g_string_append_printf(report,
                         "ZMATRIX editor\n"
                         "Model: %s\n"
                         "Scope: %s\n"
                         "Atoms: %u\n"
                         "Distance units: Angstrom\n"
                         "Angle units: degrees\n\n",
                         model->basename ? model->basename : "model",
                         use_selection ? "current selection" : "full model",
                         scope->len);

  for (guint i = 0; i < rows->len; i++)
    {
      const GdisZmatrixRow *row;
      const GdisAtom *atom;
      guint distance_ref;
      guint angle_ref;
      guint torsion_ref;
      gdouble distance_value;
      gdouble angle_value;
      gdouble torsion_value;

      row = &g_array_index(rows, GdisZmatrixRow, i);
      atom = g_ptr_array_index(model->atoms, row->atom_index);

      distance_ref = (row->distance_ref != G_MAXUINT)
                     ? (guint) (gdis_scope_find_index(scope, row->distance_ref) + 1)
                     : 0u;
      angle_ref = (row->angle_ref != G_MAXUINT)
                  ? (guint) (gdis_scope_find_index(scope, row->angle_ref) + 1)
                  : 0u;
      torsion_ref = (row->torsion_ref != G_MAXUINT)
                    ? (guint) (gdis_scope_find_index(scope, row->torsion_ref) + 1)
                    : 0u;
      distance_value = (row->distance_ref != G_MAXUINT) ? row->distance : 0.0;
      angle_value = (row->angle_ref != G_MAXUINT) ? row->angle : 0.0;
      torsion_value = (row->torsion_ref != G_MAXUINT) ? row->torsion : 0.0;

      g_string_append_printf(
        report,
        "[%u]  %-4s  %d %d %d  %-9.4f  %-9.4f  %-9.4f\n",
        i + 1u,
        atom->label && atom->label[0] != '\0' ? atom->label : atom->element,
        (gint) distance_ref,
        (gint) angle_ref,
        (gint) torsion_ref,
        distance_value,
        angle_value,
        torsion_value);
    }

  g_string_append(report,
                  "\nUse the row editor below to change numeric values, then choose Recompute geometry to apply them.\n"
                  "Build zmatrix from selection follows the current selected atom set after GTK4 reorders it into a connectivity-friendly path.\n");
  return g_string_free(report, FALSE);
}

gboolean
gdis_apply_zmatrix_rows(GdisModel *model,
                        const GArray *scope,
                        const GArray *rows,
                        gchar **summary_out,
                        GError **error)
{
  gdouble *original_positions;
  gdouble *local_positions;
  gdouble origin[3];
  gdouble axis_x[3];
  gdouble axis_y[3];
  gdouble axis_z[3];

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(scope != NULL, FALSE);
  g_return_val_if_fail(rows != NULL, FALSE);

  if (scope->len == 0u || rows->len == 0u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "The Z-matrix scope is empty.");
      return FALSE;
    }
  if (scope->len != rows->len)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "The Z-matrix scope and row list are out of sync.");
      return FALSE;
    }

  original_positions = g_new0(gdouble, 3u * scope->len);
  local_positions = g_new0(gdouble, 3u * scope->len);

  for (guint i = 0; i < scope->len; i++)
    {
      const GdisAtom *atom;

      atom = g_ptr_array_index(model->atoms, g_array_index(scope, guint, i));
      original_positions[3u * i] = atom->position[0];
      original_positions[3u * i + 1u] = atom->position[1];
      original_positions[3u * i + 2u] = atom->position[2];
    }

  gdis_zmatrix_build_world_basis(original_positions, scope->len, origin, axis_x, axis_y, axis_z);

  for (guint i = 0; i < rows->len; i++)
    {
      const GdisZmatrixRow *row;
      gdouble *position;

      row = &g_array_index(rows, GdisZmatrixRow, i);
      position = local_positions + 3u * i;

      if (row->atom_index != g_array_index(scope, guint, i))
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_FAILED,
                      "The Z-matrix row order no longer matches the current scope.");
          g_free(local_positions);
          g_free(original_positions);
          return FALSE;
        }

      if (i == 0u)
        continue;

      if (row->distance_ref == G_MAXUINT)
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_FAILED,
                      "Row %u is missing its bond reference.",
                      i + 1u);
          g_free(local_positions);
          g_free(original_positions);
          return FALSE;
        }

      {
        gint distance_index;

        distance_index = gdis_scope_find_index(scope, row->distance_ref);
        if (distance_index < 0 || distance_index >= (gint) i)
          {
            g_set_error(error,
                        GDIS_MODEL_ERROR,
                        GDIS_MODEL_ERROR_FAILED,
                        "Row %u has an invalid bond reference.",
                        i + 1u);
            g_free(local_positions);
            g_free(original_positions);
            return FALSE;
          }

        if (row->angle_ref == G_MAXUINT)
          {
            const gdouble *ref_b = local_positions + 3u * (guint) distance_index;

            position[0] = ref_b[0] + MAX(row->distance, 1.0e-6);
            position[1] = ref_b[1];
            position[2] = ref_b[2];
            continue;
          }

        {
          gint angle_index;

          angle_index = gdis_scope_find_index(scope, row->angle_ref);
          if (angle_index < 0 || angle_index >= (gint) i)
            {
              g_set_error(error,
                          GDIS_MODEL_ERROR,
                          GDIS_MODEL_ERROR_FAILED,
                          "Row %u has an invalid angle reference.",
                          i + 1u);
              g_free(local_positions);
              g_free(original_positions);
              return FALSE;
            }

          if (row->torsion_ref == G_MAXUINT)
            {
              if (!gdis_zmatrix_place_with_two_refs(local_positions + 3u * (guint) distance_index,
                                                    local_positions + 3u * (guint) angle_index,
                                                    MAX(row->distance, 1.0e-6),
                                                    row->angle,
                                                    position))
                {
                  g_set_error(error,
                              GDIS_MODEL_ERROR,
                              GDIS_MODEL_ERROR_FAILED,
                              "Row %u could not be positioned from its bond and angle references.",
                              i + 1u);
                  g_free(local_positions);
                  g_free(original_positions);
                  return FALSE;
                }
              continue;
            }

          {
            gint torsion_index;

            torsion_index = gdis_scope_find_index(scope, row->torsion_ref);
            if (torsion_index < 0 || torsion_index >= (gint) i)
              {
                g_set_error(error,
                            GDIS_MODEL_ERROR,
                            GDIS_MODEL_ERROR_FAILED,
                            "Row %u has an invalid torsion reference.",
                            i + 1u);
                g_free(local_positions);
                g_free(original_positions);
                return FALSE;
              }

            if (!gdis_zmatrix_place_with_three_refs(local_positions + 3u * (guint) distance_index,
                                                    local_positions + 3u * (guint) angle_index,
                                                    local_positions + 3u * (guint) torsion_index,
                                                    MAX(row->distance, 1.0e-6),
                                                    row->angle,
                                                    row->torsion,
                                                    position))
              {
                g_set_error(error,
                            GDIS_MODEL_ERROR,
                            GDIS_MODEL_ERROR_FAILED,
                            "Row %u could not be positioned from its bond, angle, and torsion references.",
                            i + 1u);
                g_free(local_positions);
                g_free(original_positions);
                return FALSE;
              }
          }
        }
      }
    }

  for (guint i = 0; i < scope->len; i++)
    {
      GdisAtom *atom;
      const gdouble *local_position;

      atom = g_ptr_array_index(model->atoms, g_array_index(scope, guint, i));
      local_position = local_positions + 3u * i;
      atom->position[0] = origin[0];
      atom->position[1] = origin[1];
      atom->position[2] = origin[2];
      gdis_vec3_add_scaled(atom->position, axis_x, local_position[0]);
      gdis_vec3_add_scaled(atom->position, axis_y, local_position[1]);
      gdis_vec3_add_scaled(atom->position, axis_z, local_position[2]);
    }

  gdis_model_reset_inferred_bonds(model);

  if (summary_out)
    {
      *summary_out = g_strdup_printf(
        "Recomputed geometry from %u Z-matrix row%s.\n"
        "Scope: %s\n"
        "GTK4 now applies the edited internal coordinates directly into the live model.",
        scope->len,
        scope->len == 1u ? "" : "s",
        scope->len == model->atoms->len ? "whole model" : "current selection");
    }

  g_free(local_positions);
  g_free(original_positions);
  return TRUE;
}

char *
gdis_build_zmatrix_report(const GdisModel *model,
                          const GArray *selected_atoms,
                          gboolean use_selection,
                          GError **error)
{
  GArray *scope;
  GArray *rows;
  gchar *report;

  g_return_val_if_fail(model != NULL, NULL);

  if (!gdis_build_zmatrix_rows(model, selected_atoms, use_selection, &scope, &rows, error))
    return NULL;

  report = gdis_format_zmatrix_rows(model, scope, rows, use_selection);
  g_array_free(rows, TRUE);
  g_array_free(scope, TRUE);
  return report;
}

static GArray *
gdis_build_scope_array(const GdisModel *model,
                       const GArray *selected_atoms,
                       gboolean use_selection,
                       GError **error)
{
  GArray *scope;

  g_return_val_if_fail(model != NULL, NULL);
  g_return_val_if_fail(model->atoms != NULL, NULL);

  scope = g_array_new(FALSE, FALSE, sizeof(guint));
  if (use_selection && selected_atoms && selected_atoms->len > 0u)
    {
      for (guint i = 0; i < selected_atoms->len; i++)
        {
          guint atom_index;

          atom_index = g_array_index(selected_atoms, guint, i);
          if (atom_index >= model->atoms->len || gdis_scope_contains(scope, atom_index))
            continue;
          g_array_append_val(scope, atom_index);
        }
    }
  else
    {
      for (guint i = 0; i < model->atoms->len; i++)
        g_array_append_val(scope, i);
    }

  if (scope->len == 0u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  use_selection
                    ? "The current selection is empty."
                    : "The model does not contain any atoms.");
      g_array_free(scope, TRUE);
      return NULL;
    }

  gdis_sort_scope_by_connectivity(model, scope);
  return scope;
}

static void
gdis_sort_scope_by_connectivity(const GdisModel *model, GArray *scope)
{
  GArray *sorted;
  GArray *remaining;
  guint current_atom;

  g_return_if_fail(model != NULL);
  g_return_if_fail(scope != NULL);

  if (scope->len < 3u)
    return;

  sorted = g_array_sized_new(FALSE, FALSE, sizeof(guint), scope->len);
  remaining = g_array_sized_new(FALSE, FALSE, sizeof(guint), scope->len);
  g_array_append_vals(remaining, scope->data, scope->len);

  current_atom = g_array_index(remaining, guint, 0u);
  g_array_append_val(sorted, current_atom);
  g_array_remove_index(remaining, 0u);

  while (remaining->len > 0u)
    {
      gint bonded_to_current = -1;
      gint bonded_to_scope = -1;
      gint nearest_any = 0;
      gdouble bonded_to_current_dist = G_MAXDOUBLE;
      gdouble bonded_to_scope_dist = G_MAXDOUBLE;
      gdouble nearest_any_dist = G_MAXDOUBLE;
      const GdisAtom *current_atom_record;

      current_atom_record = g_ptr_array_index(model->atoms, current_atom);
      for (guint i = 0; i < remaining->len; i++)
        {
          guint candidate_atom;
          const GdisAtom *candidate_record;
          gdouble distance_sq;
          gboolean bonded_to_previous_scope;

          candidate_atom = g_array_index(remaining, guint, i);
          candidate_record = g_ptr_array_index(model->atoms, candidate_atom);
          distance_sq = gdis_distance_sq(current_atom_record->position, candidate_record->position);
          if (distance_sq < nearest_any_dist)
            {
              nearest_any_dist = distance_sq;
              nearest_any = (gint) i;
            }

          bonded_to_previous_scope = FALSE;
          if (gdis_atoms_are_bonded(model, current_atom, candidate_atom) &&
              distance_sq < bonded_to_current_dist)
            {
              bonded_to_current_dist = distance_sq;
              bonded_to_current = (gint) i;
            }

          for (guint j = 0; j < sorted->len; j++)
            {
              guint previous_atom;

              previous_atom = g_array_index(sorted, guint, j);
              if (gdis_atoms_are_bonded(model, previous_atom, candidate_atom))
                {
                  bonded_to_previous_scope = TRUE;
                  break;
                }
            }

          if (bonded_to_previous_scope && distance_sq < bonded_to_scope_dist)
            {
              bonded_to_scope_dist = distance_sq;
              bonded_to_scope = (gint) i;
            }
        }

      if (bonded_to_current >= 0)
        current_atom = g_array_index(remaining, guint, (guint) bonded_to_current);
      else if (bonded_to_scope >= 0)
        current_atom = g_array_index(remaining, guint, (guint) bonded_to_scope);
      else
        current_atom = g_array_index(remaining, guint, (guint) nearest_any);

      g_array_append_val(sorted, current_atom);
      g_array_remove_index(remaining,
                           bonded_to_current >= 0 ? (guint) bonded_to_current :
                           bonded_to_scope >= 0 ? (guint) bonded_to_scope :
                           (guint) nearest_any);
    }

  g_array_set_size(scope, 0u);
  g_array_append_vals(scope, sorted->data, sorted->len);
  g_array_free(sorted, TRUE);
  g_array_free(remaining, TRUE);
}

static void
gdis_compute_scope_centroid(const GdisModel *model,
                            const GArray *scope,
                            gdouble centroid[3])
{
  g_return_if_fail(model != NULL);
  g_return_if_fail(scope != NULL);
  g_return_if_fail(centroid != NULL);

  gdis_vec3_set(centroid, 0.0, 0.0, 0.0);
  if (scope->len == 0u)
    return;

  for (guint i = 0; i < scope->len; i++)
    {
      const GdisAtom *atom;

      atom = g_ptr_array_index(model->atoms, g_array_index(scope, guint, i));
      centroid[0] += atom->position[0];
      centroid[1] += atom->position[1];
      centroid[2] += atom->position[2];
    }

  centroid[0] /= (gdouble) scope->len;
  centroid[1] /= (gdouble) scope->len;
  centroid[2] /= (gdouble) scope->len;
}

static gdouble
gdis_measure_scope_radius(const GdisModel *model,
                          const GArray *scope,
                          const gdouble centroid[3])
{
  gdouble radius;

  g_return_val_if_fail(model != NULL, 0.0);
  g_return_val_if_fail(scope != NULL, 0.0);
  g_return_val_if_fail(centroid != NULL, 0.0);

  radius = 0.0;
  for (guint i = 0; i < scope->len; i++)
    {
      const GdisAtom *atom;
      gdouble rel[3];

      atom = g_ptr_array_index(model->atoms, g_array_index(scope, guint, i));
      gdis_vec3_subtract(atom->position, centroid, rel);
      radius = MAX(radius, gdis_vec3_length(rel));
    }

  return radius;
}

static gboolean
gdis_scope_contains(const GArray *scope, guint atom_index)
{
  return gdis_scope_find_index(scope, atom_index) >= 0;
}

static gint
gdis_scope_find_index(const GArray *scope, guint atom_index)
{
  g_return_val_if_fail(scope != NULL, -1);

  for (guint i = 0; i < scope->len; i++)
    {
      if (g_array_index(scope, guint, i) == atom_index)
        return (gint) i;
    }

  return -1;
}

static gboolean
gdis_atoms_are_bonded(const GdisModel *model, guint atom_index_a, guint atom_index_b)
{
  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(model->bonds != NULL, FALSE);

  for (guint i = 0; i < model->bonds->len; i++)
    {
      const GdisBond *bond;

      bond = &g_array_index(model->bonds, GdisBond, i);
      if ((bond->atom_index_a == atom_index_a && bond->atom_index_b == atom_index_b) ||
          (bond->atom_index_a == atom_index_b && bond->atom_index_b == atom_index_a))
        return TRUE;
    }

  return FALSE;
}

static gdouble
gdis_distance_sq(const gdouble a[3], const gdouble b[3])
{
  const gdouble dx = a[0] - b[0];
  const gdouble dy = a[1] - b[1];
  const gdouble dz = a[2] - b[2];

  return dx * dx + dy * dy + dz * dz;
}

static void
gdis_vec3_copy(gdouble dest[3], const gdouble src[3])
{
  memcpy(dest, src, sizeof(gdouble) * 3u);
}

static void
gdis_vec3_set(gdouble vector[3], gdouble x, gdouble y, gdouble z)
{
  vector[0] = x;
  vector[1] = y;
  vector[2] = z;
}

static void
gdis_vec3_subtract(const gdouble a[3], const gdouble b[3], gdouble out[3])
{
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}

static void
gdis_vec3_add_scaled(gdouble vector[3], const gdouble addend[3], gdouble scale)
{
  vector[0] += addend[0] * scale;
  vector[1] += addend[1] * scale;
  vector[2] += addend[2] * scale;
}

static gdouble
gdis_vec3_dot(const gdouble a[3], const gdouble b[3])
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void
gdis_vec3_cross(const gdouble a[3], const gdouble b[3], gdouble out[3])
{
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

static gdouble
gdis_vec3_length(const gdouble vector[3])
{
  return sqrt(gdis_vec3_dot(vector, vector));
}

static gboolean
gdis_vec3_normalize(gdouble vector[3])
{
  gdouble length;

  length = gdis_vec3_length(vector);
  if (length < 1.0e-12)
    return FALSE;

  vector[0] /= length;
  vector[1] /= length;
  vector[2] /= length;
  return TRUE;
}

static void
gdis_rotate_vector_axis_angle(const gdouble vector[3],
                              const gdouble axis[3],
                              gdouble angle,
                              gdouble out[3])
{
  gdouble unit_axis[3];
  gdouble cross[3];
  gdouble axis_projection;
  gdouble cosine;
  gdouble sine;

  gdis_vec3_copy(unit_axis, axis);
  if (!gdis_vec3_normalize(unit_axis))
    {
      gdis_vec3_copy(out, vector);
      return;
    }

  cosine = cos(angle);
  sine = sin(angle);
  axis_projection = gdis_vec3_dot(unit_axis, vector);
  gdis_vec3_cross(unit_axis, vector, cross);

  out[0] = vector[0] * cosine + cross[0] * sine + unit_axis[0] * axis_projection * (1.0 - cosine);
  out[1] = vector[1] * cosine + cross[1] * sine + unit_axis[1] * axis_projection * (1.0 - cosine);
  out[2] = vector[2] * cosine + cross[2] * sine + unit_axis[2] * axis_projection * (1.0 - cosine);
}

static void
gdis_rotate_euler_xyz(const gdouble vector[3],
                      gdouble angle_x,
                      gdouble angle_y,
                      gdouble angle_z,
                      gdouble out[3])
{
  gdouble temp[3];
  gdouble axis_x[3] = {1.0, 0.0, 0.0};
  gdouble axis_y[3] = {0.0, 1.0, 0.0};
  gdouble axis_z[3] = {0.0, 0.0, 1.0};

  gdis_rotate_vector_axis_angle(vector, axis_x, angle_x, temp);
  gdis_rotate_vector_axis_angle(temp, axis_y, angle_y, out);
  gdis_rotate_vector_axis_angle(out, axis_z, angle_z, temp);
  gdis_vec3_copy(out, temp);
}

static gboolean
gdis_measure_value(const GdisModel *model,
                   GdisMeasureMode mode,
                   const guint *atom_indices,
                   guint count,
                   gdouble *value_out)
{
  guint used_indices[4] = {0u, 0u, 0u, 0u};
  guint used_count = 0u;
  GdisMeasureMode resolved_mode = mode;
  GError *error = NULL;
  gboolean ok;

  ok = gdis_measure_calculate(model,
                              mode,
                              atom_indices,
                              count,
                              used_indices,
                              &used_count,
                              &resolved_mode,
                              value_out,
                              &error);
  g_clear_error(&error);
  return ok;
}

static guint
gdis_pick_distance_reference(const GdisModel *model,
                             const GArray *scope,
                             guint scope_position)
{
  const GdisAtom *atom;
  guint atom_index;
  guint best_bonded;
  guint best_any;
  gdouble best_bonded_dist;
  gdouble best_any_dist;

  g_return_val_if_fail(model != NULL, G_MAXUINT);
  g_return_val_if_fail(scope != NULL, G_MAXUINT);

  if (scope_position == 0u)
    return G_MAXUINT;

  atom_index = g_array_index(scope, guint, scope_position);
  atom = g_ptr_array_index(model->atoms, atom_index);
  best_bonded = G_MAXUINT;
  best_any = G_MAXUINT;
  best_bonded_dist = G_MAXDOUBLE;
  best_any_dist = G_MAXDOUBLE;

  for (guint i = 0; i < scope_position; i++)
    {
      guint candidate_index;
      const GdisAtom *candidate;
      gdouble dist_sq;

      candidate_index = g_array_index(scope, guint, i);
      candidate = g_ptr_array_index(model->atoms, candidate_index);
      dist_sq = gdis_distance_sq(atom->position, candidate->position);

      if (dist_sq < best_any_dist)
        {
          best_any_dist = dist_sq;
          best_any = candidate_index;
        }

      if (gdis_atoms_are_bonded(model, atom_index, candidate_index) &&
          dist_sq < best_bonded_dist)
        {
          best_bonded_dist = dist_sq;
          best_bonded = candidate_index;
        }
    }

  return best_bonded != G_MAXUINT ? best_bonded : best_any;
}

static guint
gdis_pick_angle_reference(const GdisModel *model,
                          const GArray *scope,
                          guint scope_position,
                          guint distance_ref)
{
  const GdisAtom *distance_atom;
  guint best_bonded;
  guint best_any;
  gdouble best_bonded_dist;
  gdouble best_any_dist;

  g_return_val_if_fail(model != NULL, G_MAXUINT);
  g_return_val_if_fail(scope != NULL, G_MAXUINT);

  if (scope_position < 2u)
    return G_MAXUINT;

  distance_atom = g_ptr_array_index(model->atoms, distance_ref);
  best_bonded = G_MAXUINT;
  best_any = G_MAXUINT;
  best_bonded_dist = G_MAXDOUBLE;
  best_any_dist = G_MAXDOUBLE;

  for (guint i = 0; i < scope_position; i++)
    {
      guint candidate_index;
      const GdisAtom *candidate;
      gdouble dist_sq;

      candidate_index = g_array_index(scope, guint, i);
      if (candidate_index == distance_ref)
        continue;

      candidate = g_ptr_array_index(model->atoms, candidate_index);
      dist_sq = gdis_distance_sq(distance_atom->position, candidate->position);

      if (dist_sq < best_any_dist)
        {
          best_any_dist = dist_sq;
          best_any = candidate_index;
        }

      if (gdis_atoms_are_bonded(model, distance_ref, candidate_index) &&
          dist_sq < best_bonded_dist)
        {
          best_bonded_dist = dist_sq;
          best_bonded = candidate_index;
        }
    }

  return best_bonded != G_MAXUINT ? best_bonded : best_any;
}

static guint
gdis_pick_torsion_reference(const GdisModel *model,
                            const GArray *scope,
                            guint scope_position,
                            guint distance_ref,
                            guint angle_ref)
{
  const GdisAtom *angle_atom;
  guint best_bonded;
  guint best_any;
  gdouble best_bonded_dist;
  gdouble best_any_dist;

  g_return_val_if_fail(model != NULL, G_MAXUINT);
  g_return_val_if_fail(scope != NULL, G_MAXUINT);

  if (scope_position < 3u)
    return G_MAXUINT;

  angle_atom = g_ptr_array_index(model->atoms, angle_ref);
  best_bonded = G_MAXUINT;
  best_any = G_MAXUINT;
  best_bonded_dist = G_MAXDOUBLE;
  best_any_dist = G_MAXDOUBLE;

  for (guint i = 0; i < scope_position; i++)
    {
      guint candidate_index;
      const GdisAtom *candidate;
      gdouble dist_sq;

      candidate_index = g_array_index(scope, guint, i);
      if (candidate_index == distance_ref || candidate_index == angle_ref)
        continue;

      candidate = g_ptr_array_index(model->atoms, candidate_index);
      dist_sq = gdis_distance_sq(angle_atom->position, candidate->position);

      if (dist_sq < best_any_dist)
        {
          best_any_dist = dist_sq;
          best_any = candidate_index;
        }

      if (gdis_atoms_are_bonded(model, angle_ref, candidate_index) &&
          dist_sq < best_bonded_dist)
        {
          best_bonded_dist = dist_sq;
          best_bonded = candidate_index;
        }
    }

  return best_bonded != G_MAXUINT ? best_bonded : best_any;
}

static void
gdis_pick_orthogonal_axis(const gdouble axis[3], gdouble orthogonal[3])
{
  gdouble reference[3] = {1.0, 0.0, 0.0};

  if (fabs(axis[0]) > 0.85)
    {
      reference[0] = 0.0;
      reference[1] = 1.0;
    }
  if (fabs(axis[0]) > 0.85 && fabs(axis[1]) > 0.85)
    {
      reference[1] = 0.0;
      reference[2] = 1.0;
    }

  gdis_vec3_cross(axis, reference, orthogonal);
  if (!gdis_vec3_normalize(orthogonal))
    gdis_vec3_set(orthogonal, 0.0, 1.0, 0.0);
}

static gboolean
gdis_zmatrix_place_with_two_refs(const gdouble ref_b[3],
                                 const gdouble ref_c[3],
                                 gdouble distance,
                                 gdouble angle_degrees,
                                 gdouble position_out[3])
{
  gdouble axis[3];
  gdouble normal[3] = {0.0, 0.0, 1.0};
  gdouble tangent[3];
  const gdouble angle = angle_degrees * (G_PI / 180.0);

  gdis_vec3_subtract(ref_c, ref_b, axis);
  if (!gdis_vec3_normalize(axis))
    return FALSE;

  if (fabs(gdis_vec3_dot(axis, normal)) > 0.95)
    gdis_vec3_set(normal, 0.0, 1.0, 0.0);

  gdis_vec3_cross(normal, axis, tangent);
  if (!gdis_vec3_normalize(tangent))
    return FALSE;

  gdis_vec3_copy(position_out, ref_b);
  gdis_vec3_add_scaled(position_out, axis, distance * cos(angle));
  gdis_vec3_add_scaled(position_out, tangent, distance * sin(angle));
  return TRUE;
}

static gboolean
gdis_zmatrix_place_with_three_refs(const gdouble ref_b[3],
                                   const gdouble ref_c[3],
                                   const gdouble ref_d[3],
                                   gdouble distance,
                                   gdouble angle_degrees,
                                   gdouble torsion_degrees,
                                   gdouble position_out[3])
{
  gdouble axis_bc[3];
  gdouble axis_cd[3];
  gdouble normal[3];
  gdouble binormal[3];
  const gdouble angle = angle_degrees * (G_PI / 180.0);
  const gdouble torsion = -torsion_degrees * (G_PI / 180.0);

  gdis_vec3_subtract(ref_c, ref_b, axis_bc);
  if (!gdis_vec3_normalize(axis_bc))
    return FALSE;

  gdis_vec3_subtract(ref_d, ref_c, axis_cd);
  gdis_vec3_cross(axis_cd, axis_bc, normal);
  if (!gdis_vec3_normalize(normal))
    gdis_pick_orthogonal_axis(axis_bc, normal);

  gdis_vec3_cross(axis_bc, normal, binormal);
  if (!gdis_vec3_normalize(binormal))
    return FALSE;

  gdis_vec3_copy(position_out, ref_b);
  gdis_vec3_add_scaled(position_out, axis_bc, distance * cos(angle));
  gdis_vec3_add_scaled(position_out, binormal, distance * sin(angle) * cos(torsion));
  gdis_vec3_add_scaled(position_out, normal, distance * sin(angle) * sin(torsion));
  return TRUE;
}

static void
gdis_zmatrix_build_world_basis(const gdouble *original_positions,
                               guint scope_len,
                               gdouble origin[3],
                               gdouble axis_x[3],
                               gdouble axis_y[3],
                               gdouble axis_z[3])
{
  gdouble candidate[3];

  g_return_if_fail(original_positions != NULL);
  g_return_if_fail(scope_len > 0u);

  gdis_vec3_set(origin,
                original_positions[0],
                original_positions[1],
                original_positions[2]);

  if (scope_len >= 2u)
    {
      gdis_vec3_set(axis_x,
                    original_positions[3u] - original_positions[0],
                    original_positions[4u] - original_positions[1],
                    original_positions[5u] - original_positions[2]);
    }
  else
    {
      gdis_vec3_set(axis_x, 1.0, 0.0, 0.0);
    }

  if (!gdis_vec3_normalize(axis_x))
    gdis_vec3_set(axis_x, 1.0, 0.0, 0.0);

  if (scope_len >= 3u)
    {
      gdis_vec3_set(candidate,
                    original_positions[6u] - original_positions[0],
                    original_positions[7u] - original_positions[1],
                    original_positions[8u] - original_positions[2]);
      gdis_vec3_cross(axis_x, candidate, axis_z);
    }
  else
    {
      gdis_vec3_set(axis_z, 0.0, 0.0, 0.0);
    }

  if (!gdis_vec3_normalize(axis_z))
    {
      gdis_pick_orthogonal_axis(axis_x, axis_y);
      gdis_vec3_cross(axis_x, axis_y, axis_z);
      gdis_vec3_normalize(axis_z);
      gdis_vec3_cross(axis_z, axis_x, axis_y);
      gdis_vec3_normalize(axis_y);
      return;
    }

  gdis_vec3_cross(axis_z, axis_x, axis_y);
  if (!gdis_vec3_normalize(axis_y))
    {
      gdis_pick_orthogonal_axis(axis_x, axis_y);
      gdis_vec3_cross(axis_x, axis_y, axis_z);
      gdis_vec3_normalize(axis_z);
      gdis_vec3_cross(axis_z, axis_x, axis_y);
      gdis_vec3_normalize(axis_y);
    }
}

static gboolean
gdis_write_xyz_subset(FILE *stream,
                      const GdisModel *model,
                      const GArray *scope,
                      const gchar *title)
{
  g_return_val_if_fail(stream != NULL, FALSE);
  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(scope != NULL, FALSE);

  fprintf(stream, "%u\n%s\n", scope->len, title ? title : "subset");
  for (guint i = 0; i < scope->len; i++)
    {
      const GdisAtom *atom;

      atom = g_ptr_array_index(model->atoms, g_array_index(scope, guint, i));
      fprintf(stream,
              "%-2s % .8f % .8f % .8f\n",
              atom->element ? atom->element : "X",
              atom->position[0],
              atom->position[1],
              atom->position[2]);
    }

  return TRUE;
}

static gboolean
gdis_write_xyz_pose(FILE *stream,
                    const GdisModel *model,
                    const GArray *scope,
                    const gdouble *positions,
                    const gchar *title)
{
  g_return_val_if_fail(stream != NULL, FALSE);
  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(scope != NULL, FALSE);
  g_return_val_if_fail(positions != NULL, FALSE);

  fprintf(stream, "%u\n%s\n", scope->len, title ? title : "pose");
  for (guint i = 0; i < scope->len; i++)
    {
      const GdisAtom *atom;

      atom = g_ptr_array_index(model->atoms, g_array_index(scope, guint, i));
      fprintf(stream,
              "%-2s % .8f % .8f % .8f\n",
              atom->element ? atom->element : "X",
              positions[3u * i],
              positions[3u * i + 1u],
              positions[3u * i + 2u]);
    }

  return TRUE;
}
