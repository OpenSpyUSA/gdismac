#include "gdis_model.h"

#include <ctype.h>
#include <math.h>
#include <string.h>

typedef struct
{
  gchar *label;
  gchar *element;
  gdouble coords[3];
  gdouble occupancy;
  gboolean fractional;
} GdisCifAtomRecord;

typedef struct
{
  gchar *title;
  gchar *space_group;
  gboolean periodic;
  guint periodicity;
  gdouble cell_lengths[3];
  gdouble cell_angles[3];
  GPtrArray *atoms;
} GdisModelFrame;

typedef struct
{
  GdisModel *model;
  GHashTable *species_symbols;
  gchar *current_species_name;
  gchar *current_atom_name;
  gchar *current_atom_species;
  GString *text;
  gdouble current_atom_position[3];
  gdouble a_vec[3];
  gdouble b_vec[3];
  gdouble c_vec[3];
  gboolean inside_species_symbol;
  gboolean inside_atom_position;
  gboolean have_current_atom_position;
  gboolean have_cell;
  gboolean looks_like_qbox;
} GdisQboxXmlParseState;

typedef struct
{
  const char *symbol;
  gdouble radius;
} GdisRadiusEntry;

static GdisModel *gdis_model_new(const char *path, GdisModelFormat format);
static void gdis_model_clear_contents(GdisModel *model);
static GdisModelFrame *gdis_model_frame_new(void);
static void gdis_model_frame_free(gpointer data);
static GdisModelFrame *gdis_model_frame_clone(const GdisModelFrame *frame);
static void gdis_model_frame_apply(GdisModel *model, const GdisModelFrame *frame);
static gboolean gdis_model_frame_bonds_need_reset(const GdisModel *model);
static gboolean gdis_model_load_xyz(GdisModel *model, const gchar *contents, GError **error);
static gboolean gdis_model_load_pdb(GdisModel *model, const gchar *contents, GError **error);
static gboolean gdis_model_load_arc_like(GdisModel *model,
                                         const gchar *contents,
                                         GdisModelFormat format,
                                         GError **error);
static gboolean gdis_model_load_cif(GdisModel *model, const gchar *contents, GError **error);
static gboolean gdis_model_load_qbox_xml(GdisModel *model, const gchar *contents, GError **error);
static gboolean gdis_model_load_gulp_input_like(GdisModel *model,
                                                const gchar *contents,
                                                GdisModelFormat format,
                                                GError **error);
static gboolean gdis_model_load_gulp_output(GdisModel *model, const gchar *contents, GError **error);
static gboolean gdis_model_load_xtl(GdisModel *model, const gchar *contents, GError **error);
static gboolean gdis_model_load_qe_input(GdisModel *model, const gchar *contents, GError **error);
static gboolean gdis_model_load_qe_output(GdisModel *model, const gchar *contents, GError **error);
static gboolean gdis_model_load_gamess_input(GdisModel *model, const gchar *contents, GError **error);
static gchar *gdis_model_write_xyz(const GdisModel *model);
static gchar *gdis_model_write_pdb(const GdisModel *model);
static gchar *gdis_model_write_arc_like(const GdisModel *model, GdisModelFormat format);
static gchar *gdis_model_write_cif(const GdisModel *model, GError **error);
static gboolean gdis_model_infer_bonds(GdisModel *model);
static gboolean gdis_model_add_bond(GdisModel *model,
                                    guint atom_index_a,
                                    guint atom_index_b,
                                    guint8 order,
                                    gboolean inferred);
static void gdis_model_make_all_bonds_explicit(GdisModel *model);
static gint gdis_model_find_bond_index(const GdisModel *model,
                                       guint atom_index_a,
                                       guint atom_index_b);
static void gdis_model_update_identity(GdisModel *model,
                                       const char *path,
                                       GdisModelFormat format);
static void gdis_model_refresh_counts(GdisModel *model);
static void gdis_model_finalize_metadata(GdisModel *model);
static GdisAtom *gdis_atom_new(const gchar *label,
                               const gchar *element,
                               const gchar *ff_type,
                               gdouble x,
                               gdouble y,
                               gdouble z,
                               gdouble occupancy,
                               gint region,
                               guint serial);
static void gdis_atom_free(gpointer data);
static void gdis_cif_atom_record_free(gpointer data);
static gboolean gdis_model_has_valid_cell(const GdisModel *model);
static gboolean gdis_build_cell_matrix_from_parameters(const gdouble lengths[3],
                                                       const gdouble angles[3],
                                                       gdouble matrix[9],
                                                       gdouble inverse[9]);
static gboolean gdis_model_build_cell_matrix(const GdisModel *model, gdouble matrix[9], gdouble inverse[9]);
static gboolean gdis_model_invert_matrix3(const gdouble matrix[9], gdouble inverse[9]);
static void gdis_frac_to_cart(const gdouble matrix[9], const gdouble frac[3], gdouble cart[3]);
static void gdis_cart_to_frac(const gdouble inverse[9], const gdouble cart[3], gdouble frac[3]);
static void gdis_model_init_image_limits(GdisModel *model);
static void gdis_model_wrap_atoms_to_cell(GPtrArray *atoms,
                                          guint periodicity,
                                          const gdouble matrix[9],
                                          const gdouble inverse[9]);
static gboolean gdis_model_repack_periodic_components(GdisModel *model,
                                                      const gdouble matrix[9],
                                                      const gdouble inverse[9]);
static void gdis_model_normalize_loaded_periodic_geometry(GdisModel *model);
static guint *gdis_model_build_component_ids(const GdisModel *model, guint *component_count_out);
static void gdis_model_minimum_image_delta(const GdisModel *model,
                                           const gdouble matrix[9],
                                           const gdouble inverse[9],
                                           const gdouble from[3],
                                           const gdouble to[3],
                                           gdouble delta[3]);
static gdouble gdis_model_distance2(const GdisModel *model,
                                    const gdouble a[3],
                                    const gdouble b[3],
                                    const gdouble matrix[9],
                                    const gdouble inverse[9],
                                    gboolean use_periodic_image);
static gchar *gdis_strdup_strip(const gchar *text);
static gchar *gdis_copy_field(const gchar *line, gsize start, gsize width);
static gchar *gdis_normalize_element_symbol(const gchar *text);
static gchar *gdis_pdb_guess_element(const gchar *atom_name);
static gboolean gdis_arc_like_type_is_marvin_label(const gchar *text);
static gint gdis_arc_like_region_from_type(const gchar *text);
static gboolean gdis_try_parse_double(const gchar *text, gdouble *value);
static gboolean gdis_try_parse_double_relaxed(const gchar *text, gdouble *value);
static gboolean gdis_try_parse_uint(const gchar *text, guint *value);
static guint gdis_collect_doubles_from_line(const gchar *line,
                                            gdouble *values,
                                            guint max_values);
static gboolean gdis_collect_three_doubles_from_tokens(gchar **tokens,
                                                       gint token_count,
                                                       gint start_index,
                                                       gdouble out[3]);
static gchar **gdis_split_simple(const gchar *line, gint *count);
static GPtrArray *gdis_tokenize_cif(const gchar *line);
static gboolean gdis_cif_parse_value_line(const gchar *line, gchar **tag_out, gchar **value_out);
static gdouble gdis_parse_cif_number(const gchar *text);
static gboolean gdis_cif_is_control_line(const gchar *line);
static gdouble gdis_lookup_covalent_radius(const gchar *element);
static const gchar *gdis_xml_local_name(const gchar *name);
static gboolean gdis_parse_vector3_text(const gchar *text, gdouble vector[3]);
static void gdis_vec3_set(gdouble vector[3], gdouble x, gdouble y, gdouble z);
static void gdis_vec3_copy(gdouble dest[3], const gdouble src[3]);
static void gdis_vec3_add_scaled(gdouble vector[3], const gdouble addend[3], gdouble scale);
static gdouble gdis_vec3_dot(const gdouble a[3], const gdouble b[3]);
static void gdis_vec3_cross(const gdouble a[3], const gdouble b[3], gdouble out[3]);
static gdouble gdis_vec3_length(const gdouble vector[3]);
static gboolean gdis_vec3_normalize(gdouble vector[3]);
static gint gdis_int_gcd(gint left, gint right);
static void gdis_int_vector_reduce(gint vector[3]);
static gboolean gdis_surface_build_inplane_indices(gint h,
                                                   gint k,
                                                   gint l,
                                                   gint first[3],
                                                   gint second[3]);
static gboolean gdis_surface_build_cell_from_vectors(const gdouble a_vec[3],
                                                     const gdouble b_vec[3],
                                                     const gdouble c_vec[3],
                                                     gdouble lengths[3],
                                                     gdouble angles[3]);
static void gdis_qbox_xml_start_element(GMarkupParseContext *context,
                                        const gchar         *element_name,
                                        const gchar        **attribute_names,
                                        const gchar        **attribute_values,
                                        gpointer             user_data,
                                        GError             **parse_error);
static void gdis_qbox_xml_end_element(GMarkupParseContext *context,
                                      const gchar         *element_name,
                                      gpointer             user_data,
                                      GError             **parse_error);
static void gdis_qbox_xml_text(GMarkupParseContext *context,
                               const gchar         *text,
                               gsize                text_len,
                               gpointer             user_data,
                               GError             **parse_error);
static gchar *gdis_model_surface_default_path(const GdisModel *source,
                                              gint h,
                                              gint k,
                                              gint l);

GQuark
gdis_model_error_quark(void)
{
  return g_quark_from_static_string("gdis-model-error-quark");
}

GdisModel *
gdis_model_create_empty(const char *path, GdisModelFormat format)
{
  GdisModel *model;

  g_return_val_if_fail(path != NULL, NULL);

  model = gdis_model_new(path, format);
  model->title = g_strdup(model->basename ? model->basename : "Untitled");
  gdis_model_finalize_metadata(model);
  return model;
}

GdisModel *
gdis_model_load(const char *path, GError **error)
{
  GdisModel *model;
  GdisModelFormat format;
  gchar *contents;
  gsize length;
  GError *local_error;
  gboolean ok;

  g_return_val_if_fail(path != NULL, NULL);

  format = gdis_model_format_from_path(path);
  if (format == GDIS_MODEL_FORMAT_UNKNOWN)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_UNSUPPORTED_FORMAT,
                  "Could not determine a supported model format from '%s'.",
                  path);
      return NULL;
    }

  contents = NULL;
  length = 0;
  local_error = NULL;
  if (!g_file_get_contents(path, &contents, &length, &local_error))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not read '%s': %s",
                  path,
                  local_error ? local_error->message : "unknown error");
      g_clear_error(&local_error);
      return NULL;
    }

  model = gdis_model_new(path, format);
  ok = FALSE;

  switch (format)
    {
    case GDIS_MODEL_FORMAT_XYZ:
      ok = gdis_model_load_xyz(model, contents, error);
      break;
    case GDIS_MODEL_FORMAT_PDB:
      ok = gdis_model_load_pdb(model, contents, error);
      break;
    case GDIS_MODEL_FORMAT_ARC:
    case GDIS_MODEL_FORMAT_CAR:
      ok = gdis_model_load_arc_like(model, contents, format, error);
      break;
    case GDIS_MODEL_FORMAT_CIF:
      ok = gdis_model_load_cif(model, contents, error);
      break;
    case GDIS_MODEL_FORMAT_QBOX_XML:
      ok = gdis_model_load_qbox_xml(model, contents, error);
      break;
    case GDIS_MODEL_FORMAT_GULP_INPUT:
      ok = gdis_model_load_gulp_input_like(model, contents, format, error);
      break;
    case GDIS_MODEL_FORMAT_GULP_OUTPUT:
      ok = gdis_model_load_gulp_output(model, contents, error);
      break;
    case GDIS_MODEL_FORMAT_XTL:
      ok = gdis_model_load_xtl(model, contents, error);
      break;
    case GDIS_MODEL_FORMAT_QE_INPUT:
      ok = gdis_model_load_qe_input(model, contents, error);
      break;
    case GDIS_MODEL_FORMAT_QE_OUTPUT:
      ok = gdis_model_load_qe_output(model, contents, error);
      break;
    case GDIS_MODEL_FORMAT_GAMESS_INPUT:
      ok = gdis_model_load_gamess_input(model, contents, error);
      break;
    case GDIS_MODEL_FORMAT_UNKNOWN:
    default:
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_UNSUPPORTED_FORMAT,
                  "Unsupported model format for '%s'.",
                  path);
      break;
    }

  g_free(contents);

  if (!ok)
    {
      gdis_model_free(model);
      return NULL;
    }

  model->explicit_bond_count = model->bonds->len;
  if (model->explicit_bond_count == 0)
    gdis_model_infer_bonds(model);
  gdis_model_normalize_loaded_periodic_geometry(model);

  gdis_model_finalize_metadata(model);

  if (model->atom_count == 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No atoms were parsed from '%s'.",
                  path);
      gdis_model_free(model);
      return NULL;
    }

  return model;
}

GdisModel *
gdis_model_clone(const GdisModel *model)
{
  GdisModel *copy;
  guint i;

  g_return_val_if_fail(model != NULL, NULL);

  copy = g_new0(GdisModel, 1);
  copy->path = g_strdup(model->path);
  copy->basename = g_strdup(model->basename);
  copy->title = g_strdup(model->title);
  copy->format_label = g_strdup(model->format_label);
  copy->space_group = g_strdup(model->space_group);
  copy->format = model->format;
  copy->periodic = model->periodic;
  copy->periodicity = model->periodicity;
  memcpy(copy->cell_lengths, model->cell_lengths, sizeof(copy->cell_lengths));
  memcpy(copy->cell_angles, model->cell_angles, sizeof(copy->cell_angles));
  memcpy(copy->image_limits, model->image_limits, sizeof(copy->image_limits));
  copy->atom_count = model->atom_count;
  copy->bond_count = model->bond_count;
  copy->explicit_bond_count = model->explicit_bond_count;
  copy->atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  copy->bonds = g_array_sized_new(FALSE,
                                  FALSE,
                                  sizeof(GdisBond),
                                  model->bonds ? model->bonds->len : 0u);

  if (model->atoms)
    {
      for (i = 0; i < model->atoms->len; i++)
        {
          const GdisAtom *atom;

          atom = g_ptr_array_index(model->atoms, i);
          g_ptr_array_add(copy->atoms,
                          gdis_atom_new(atom->label,
                                        atom->element,
                                        atom->ff_type,
                                        atom->position[0],
                                        atom->position[1],
                                        atom->position[2],
                                        atom->occupancy,
                                        atom->region,
                                        atom->serial));
        }
    }

  if (model->bonds && model->bonds->len > 0)
    g_array_append_vals(copy->bonds, model->bonds->data, model->bonds->len);

  if (model->frames && model->frames->len > 0)
    {
      copy->frames = g_ptr_array_new_with_free_func(gdis_model_frame_free);
      for (i = 0; i < model->frames->len; i++)
        g_ptr_array_add(copy->frames,
                        gdis_model_frame_clone(g_ptr_array_index(model->frames, i)));
    }
  copy->current_frame_index = model->current_frame_index;

  return copy;
}

gboolean
gdis_model_copy_from(GdisModel *dest, const GdisModel *src, GError **error)
{
  GdisModel *copy;

  g_return_val_if_fail(dest != NULL, FALSE);
  g_return_val_if_fail(src != NULL, FALSE);

  if (dest == src)
    return TRUE;

  copy = gdis_model_clone(src);
  if (!copy)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not clone the model for undo restoration.");
      return FALSE;
    }

  gdis_model_clear_contents(dest);
  *dest = *copy;
  g_free(copy);
  return TRUE;
}

void
gdis_model_free(GdisModel *model)
{
  if (!model)
    return;

  gdis_model_clear_contents(model);
  g_free(model);
}

static void
gdis_model_clear_contents(GdisModel *model)
{
  if (!model)
    return;

  g_clear_pointer(&model->path, g_free);
  g_clear_pointer(&model->basename, g_free);
  g_clear_pointer(&model->title, g_free);
  g_clear_pointer(&model->format_label, g_free);
  g_clear_pointer(&model->space_group, g_free);

  if (model->atoms)
    g_ptr_array_free(model->atoms, TRUE);
  if (model->bonds)
    g_array_free(model->bonds, TRUE);
  if (model->frames)
    g_ptr_array_free(model->frames, TRUE);
  model->atoms = NULL;
  model->bonds = NULL;
  model->frames = NULL;
}

gboolean
gdis_model_save(GdisModel *model, const char *path, GError **error)
{
  GdisModelFormat format;
  g_autofree gchar *contents = NULL;
  GError *local_error;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(path != NULL, FALSE);

  format = gdis_model_format_from_path(path);
  if (format == GDIS_MODEL_FORMAT_UNKNOWN)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_UNSUPPORTED_FORMAT,
                  "Could not determine a writable format from '%s'.",
                  path);
      return FALSE;
    }

  switch (format)
    {
    case GDIS_MODEL_FORMAT_XYZ:
      contents = gdis_model_write_xyz(model);
      break;
    case GDIS_MODEL_FORMAT_PDB:
      contents = gdis_model_write_pdb(model);
      break;
    case GDIS_MODEL_FORMAT_ARC:
    case GDIS_MODEL_FORMAT_CAR:
      contents = gdis_model_write_arc_like(model, format);
      break;
    case GDIS_MODEL_FORMAT_CIF:
      contents = gdis_model_write_cif(model, error);
      break;
    case GDIS_MODEL_FORMAT_QBOX_XML:
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_UNSUPPORTED_FORMAT,
                  "Saving Qbox XML is not supported from the GTK4 rebuild yet.");
      break;
    case GDIS_MODEL_FORMAT_GULP_INPUT:
    case GDIS_MODEL_FORMAT_GULP_OUTPUT:
    case GDIS_MODEL_FORMAT_XTL:
    case GDIS_MODEL_FORMAT_QE_INPUT:
    case GDIS_MODEL_FORMAT_QE_OUTPUT:
    case GDIS_MODEL_FORMAT_GAMESS_INPUT:
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_UNSUPPORTED_FORMAT,
                  "Saving %s is not supported from the GTK4 rebuild yet.",
                  gdis_model_format_label(format));
      break;
    case GDIS_MODEL_FORMAT_UNKNOWN:
    default:
      break;
    }

  if (!contents)
    return FALSE;

  local_error = NULL;
  if (!g_file_set_contents(path, contents, -1, &local_error))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not write '%s': %s",
                  path,
                  local_error ? local_error->message : "unknown error");
      g_clear_error(&local_error);
      return FALSE;
    }

  gdis_model_update_identity(model, path, format);
  gdis_model_finalize_metadata(model);
  return TRUE;
}

gboolean
gdis_model_delete_atoms(GdisModel *model,
                        const guint *atom_indices,
                        guint count,
                        GError **error)
{
  gboolean *remove;
  guint *mapping;
  GPtrArray *new_atoms;
  GArray *new_bonds;
  guint removed_count;
  guint old_len;
  guint i;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(model->atoms != NULL, FALSE);
  g_return_val_if_fail(model->bonds != NULL, FALSE);

  if (!atom_indices || count == 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "No atoms were selected for deletion.");
      return FALSE;
    }

  old_len = model->atoms->len;
  if (old_len == 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "This model does not contain any atoms.");
      return FALSE;
    }

  remove = g_new0(gboolean, old_len);
  mapping = g_new0(guint, old_len);
  removed_count = 0;

  for (i = 0; i < count; i++)
    {
      guint index;

      index = atom_indices[i];
      if (index >= old_len || remove[index])
        continue;

      remove[index] = TRUE;
      removed_count++;
    }

  if (removed_count == 0)
    {
      g_free(remove);
      g_free(mapping);
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "The selected atoms were not valid for this model.");
      return FALSE;
    }

  gdis_model_discard_frames(model);

  new_atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  for (i = 0; i < old_len; i++)
    {
      GdisAtom *atom;

      atom = g_ptr_array_index(model->atoms, i);
      if (remove[i])
        continue;

      mapping[i] = new_atoms->len;
      atom->serial = new_atoms->len + 1;
      g_ptr_array_add(new_atoms, atom);
      model->atoms->pdata[i] = NULL;
    }

  g_ptr_array_free(model->atoms, TRUE);
  model->atoms = new_atoms;

  new_bonds = g_array_new(FALSE, FALSE, sizeof(GdisBond));
  if (model->explicit_bond_count > 0)
    {
      for (i = 0; i < model->bonds->len; i++)
        {
          GdisBond bond;

          bond = g_array_index(model->bonds, GdisBond, i);
          if (bond.atom_index_a >= old_len || bond.atom_index_b >= old_len)
            continue;
          if (remove[bond.atom_index_a] || remove[bond.atom_index_b])
            continue;

          bond.atom_index_a = mapping[bond.atom_index_a];
          bond.atom_index_b = mapping[bond.atom_index_b];
          g_array_append_val(new_bonds, bond);
        }
      model->explicit_bond_count = new_bonds->len;
    }
  else
    {
      model->explicit_bond_count = 0;
    }

  g_array_free(model->bonds, TRUE);
  model->bonds = new_bonds;

  if (model->explicit_bond_count == 0)
    gdis_model_infer_bonds(model);

  gdis_model_finalize_metadata(model);

  g_free(remove);
  g_free(mapping);
  return TRUE;
}

gboolean
gdis_model_add_atom(GdisModel *model,
                    const char *label,
                    const char *element,
                    const char *ff_type,
                    gint region,
                    gdouble x,
                    gdouble y,
                    gdouble z,
                    GError **error)
{
  g_autofree gchar *normalized_element = NULL;
  g_autofree gchar *normalized_ff_type = NULL;
  g_autofree gchar *final_label = NULL;
  GdisAtom *atom;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(model->atoms != NULL, FALSE);
  g_return_val_if_fail(model->bonds != NULL, FALSE);

  normalized_element = gdis_normalize_element_symbol(element);
  if (!normalized_element || normalized_element[0] == '\0')
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "A valid element symbol is required when adding an atom.");
      return FALSE;
    }

  final_label = g_strdup((label && *label) ? label : normalized_element);
  normalized_ff_type = g_strdup((ff_type && *ff_type) ? ff_type : normalized_element);
  atom = gdis_atom_new(final_label,
                       normalized_element,
                       normalized_ff_type,
                       x,
                       y,
                       z,
                       1.0,
                       region,
                       model->atoms->len + 1);
  gdis_model_discard_frames(model);
  g_ptr_array_add(model->atoms, atom);

  gdis_model_refresh_counts(model);
  if (model->explicit_bond_count == 0)
    gdis_model_reset_inferred_bonds(model);
  gdis_model_finalize_metadata(model);
  return TRUE;
}

gboolean
gdis_model_update_atom(GdisModel *model,
                       guint atom_index,
                       const char *label,
                       const char *element,
                       const char *ff_type,
                       gint region,
                       gdouble x,
                       gdouble y,
                       gdouble z,
                       GError **error)
{
  GdisAtom *atom;
  g_autofree gchar *normalized_element = NULL;
  g_autofree gchar *normalized_ff_type = NULL;
  gchar *new_label;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(model->atoms != NULL, FALSE);
  g_return_val_if_fail(model->bonds != NULL, FALSE);

  if (atom_index >= model->atoms->len)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "The selected atom index is out of range.");
      return FALSE;
    }

  normalized_element = gdis_normalize_element_symbol(element);
  if (!normalized_element || normalized_element[0] == '\0')
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "A valid element symbol is required when editing an atom.");
      return FALSE;
    }

  atom = g_ptr_array_index(model->atoms, atom_index);
  gdis_model_discard_frames(model);
  new_label = g_strdup((label && *label) ? label : normalized_element);
  normalized_ff_type = g_strdup((ff_type && *ff_type) ? ff_type : normalized_element);

  g_free(atom->label);
  g_free(atom->element);
  atom->label = new_label;
  atom->element = g_strdup(normalized_element);
  g_free(atom->ff_type);
  atom->ff_type = g_strdup(normalized_ff_type);
  atom->position[0] = x;
  atom->position[1] = y;
  atom->position[2] = z;
  atom->region = region;

  if (model->explicit_bond_count == 0)
    gdis_model_reset_inferred_bonds(model);
  gdis_model_finalize_metadata(model);
  return TRUE;
}

gboolean
gdis_model_add_explicit_bond(GdisModel *model,
                             guint atom_index_a,
                             guint atom_index_b,
                             guint8 order,
                             GError **error)
{
  gint bond_index;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(model->atoms != NULL, FALSE);
  g_return_val_if_fail(model->bonds != NULL, FALSE);

  if (atom_index_a >= model->atoms->len || atom_index_b >= model->atoms->len)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Both atoms must be valid before creating a bond.");
      return FALSE;
    }

  if (atom_index_a == atom_index_b)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "A bond cannot connect an atom to itself.");
      return FALSE;
    }

  gdis_model_make_all_bonds_explicit(model);
  gdis_model_discard_frames(model);
  bond_index = gdis_model_find_bond_index(model, atom_index_a, atom_index_b);
  if (bond_index >= 0)
    {
      GdisBond *bond;

      bond = &g_array_index(model->bonds, GdisBond, bond_index);
      bond->order = MAX(order, 1u);
      bond->inferred = FALSE;
      return TRUE;
    }

  if (!gdis_model_add_bond(model, atom_index_a, atom_index_b, MAX(order, 1u), FALSE))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Failed to add the requested bond.");
      return FALSE;
    }

  model->explicit_bond_count = model->bonds->len;
  gdis_model_finalize_metadata(model);
  return TRUE;
}

gboolean
gdis_model_remove_bond(GdisModel *model,
                       guint atom_index_a,
                       guint atom_index_b,
                       GError **error)
{
  gint bond_index;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(model->atoms != NULL, FALSE);
  g_return_val_if_fail(model->bonds != NULL, FALSE);

  if (atom_index_a >= model->atoms->len || atom_index_b >= model->atoms->len)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Both atoms must be valid before removing a bond.");
      return FALSE;
    }

  gdis_model_make_all_bonds_explicit(model);
  gdis_model_discard_frames(model);
  bond_index = gdis_model_find_bond_index(model, atom_index_a, atom_index_b);
  if (bond_index < 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "No bond exists between the selected atoms.");
      return FALSE;
    }

  g_array_remove_index(model->bonds, bond_index);
  model->explicit_bond_count = model->bonds->len;
  gdis_model_finalize_metadata(model);
  return TRUE;
}

void
gdis_model_get_image_limits(const GdisModel *model, gint limits_out[6])
{
  g_return_if_fail(model != NULL);
  g_return_if_fail(limits_out != NULL);

  memcpy(limits_out, model->image_limits, sizeof(model->image_limits));
}

gboolean
gdis_model_set_image_limits(GdisModel *model,
                            const gint limits[6],
                            GError **error)
{
  guint dims;
  guint axis;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(limits != NULL, FALSE);

  dims = model->periodic ? model->periodicity : 0u;
  if (dims == 0 && model->periodic)
    dims = 3;
  dims = MIN(dims, 3u);

  for (axis = 0; axis < 3; axis++)
    {
      gint neg_limit;
      gint pos_limit;

      neg_limit = limits[2 * axis];
      pos_limit = limits[2 * axis + 1];
      if (neg_limit < 0 || pos_limit < 1)
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_FAILED,
                      "Image limits must satisfy: negative >= 0 and positive >= 1.");
          return FALSE;
        }

      if (axis < dims)
        {
          model->image_limits[2 * axis] = neg_limit;
          model->image_limits[2 * axis + 1] = pos_limit;
        }
      else
        {
          model->image_limits[2 * axis] = 0;
          model->image_limits[2 * axis + 1] = 1;
        }
    }

  return TRUE;
}

void
gdis_model_reset_image_limits(GdisModel *model)
{
  g_return_if_fail(model != NULL);

  gdis_model_init_image_limits(model);
}

gboolean
gdis_model_confine_atoms_to_cell(GdisModel *model, GError **error)
{
  gdouble matrix[9];
  gdouble inverse[9];

  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->periodic)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Confine to cell requires a periodic model.");
      return FALSE;
    }

  if (!gdis_model_build_cell_matrix(model, matrix, inverse))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Confine to cell requires a valid unit cell.");
      return FALSE;
    }

  gdis_model_wrap_atoms_to_cell(model->atoms,
                                model->periodicity,
                                matrix,
                                inverse);

  gdis_model_discard_frames(model);
  if (model->explicit_bond_count == 0)
    gdis_model_reset_inferred_bonds(model);
  gdis_model_finalize_metadata(model);
  return TRUE;
}

gboolean
gdis_model_confine_molecules_to_cell(GdisModel *model, GError **error)
{
  gdouble matrix[9];
  gdouble inverse[9];

  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->periodic)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Confine molecules to cell requires a periodic model.");
      return FALSE;
    }

  if (!gdis_model_build_cell_matrix(model, matrix, inverse))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Confine molecules to cell requires a valid unit cell.");
      return FALSE;
    }

  if (!gdis_model_repack_periodic_components(model, matrix, inverse))
    gdis_model_wrap_atoms_to_cell(model->atoms,
                                  model->periodicity,
                                  matrix,
                                  inverse);

  gdis_model_discard_frames(model);
  if (model->explicit_bond_count == 0)
    gdis_model_reset_inferred_bonds(model);
  gdis_model_finalize_metadata(model);
  return TRUE;
}

gboolean
gdis_model_force_p1(GdisModel *model, GError **error)
{
  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->periodic)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Force to P1 requires a periodic model.");
      return FALSE;
    }

  g_clear_pointer(&model->space_group, g_free);
  model->space_group = g_strdup("P 1");
  if (model->periodicity == 0)
    model->periodicity = 3;
  gdis_model_discard_frames(model);
  gdis_model_finalize_metadata(model);
  return TRUE;
}

gboolean
gdis_model_make_supercell_from_image_limits(GdisModel *model, GError **error)
{
  gdouble matrix[9];
  gdouble inverse[9];
  gdouble origin_shift[3];
  guint dims;
  guint repeat_a;
  guint repeat_b;
  guint repeat_c;
  guint i;

  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->periodic)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Make supercell from images requires a periodic model.");
      return FALSE;
    }

  dims = model->periodicity;
  if (dims == 0)
    dims = 3;
  dims = MIN(dims, 3u);

  repeat_a = (guint) MAX(1, model->image_limits[0] + model->image_limits[1]);
  repeat_b = (guint) MAX(1, model->image_limits[2] + model->image_limits[3]);
  repeat_c = (guint) MAX(1, model->image_limits[4] + model->image_limits[5]);

  if (!gdis_model_build_cell_matrix(model, matrix, inverse))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Make supercell from images requires a valid unit cell.");
      return FALSE;
    }

  (void) inverse;

  origin_shift[0] = -(gdouble) model->image_limits[0] * matrix[0] -
                    -(gdouble) 0 * matrix[1];
  origin_shift[1] = 0.0;
  origin_shift[2] = 0.0;
  origin_shift[0] = -((gdouble) model->image_limits[0] * matrix[0] +
                      (gdouble) model->image_limits[2] * matrix[1] +
                      (gdouble) model->image_limits[4] * matrix[2]);
  origin_shift[1] = -((gdouble) model->image_limits[0] * matrix[3] +
                      (gdouble) model->image_limits[2] * matrix[4] +
                      (gdouble) model->image_limits[4] * matrix[5]);
  origin_shift[2] = -((gdouble) model->image_limits[0] * matrix[6] +
                      (gdouble) model->image_limits[2] * matrix[7] +
                      (gdouble) model->image_limits[4] * matrix[8]);

  for (i = 0; i < model->atoms->len; i++)
    {
      GdisAtom *atom;

      atom = g_ptr_array_index(model->atoms, i);
      atom->position[0] += origin_shift[0];
      atom->position[1] += origin_shift[1];
      atom->position[2] += origin_shift[2];
    }

  if (dims < 3)
    repeat_c = 1u;
  if (dims < 2)
    repeat_b = 1u;

  return gdis_model_make_supercell(model, repeat_a, repeat_b, repeat_c, error);
}

gboolean
gdis_model_make_supercell(GdisModel *model,
                          guint repeat_a,
                          guint repeat_b,
                          guint repeat_c,
                          GError **error)
{
  gdouble matrix[9];
  gdouble inverse[9];
  gdouble a_vec[3];
  gdouble b_vec[3];
  gdouble c_vec[3];
  guint dims;
  guint64 total_cells;
  guint64 total_atoms;
  GPtrArray *old_atoms;
  GArray *old_bonds;
  GPtrArray *new_atoms;
  GArray *new_bonds;
  guint old_atom_count;
  guint old_explicit_count;
  guint ia;
  guint ib;
  guint ic;

  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->periodic)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Make supercell requires a periodic model.");
      return FALSE;
    }

  if (repeat_a == 0 || repeat_b == 0 || repeat_c == 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Supercell repeats must all be at least 1.");
      return FALSE;
    }

  dims = model->periodicity;
  if (dims == 0)
    dims = 3;
  dims = MIN(dims, 3u);

  if ((dims < 1 && repeat_a != 1) ||
      (dims < 2 && repeat_b != 1) ||
      (dims < 3 && repeat_c != 1))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Cannot repeat along non-periodic directions.");
      return FALSE;
    }

  if (!gdis_model_build_cell_matrix(model, matrix, inverse))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Make supercell requires a valid unit cell.");
      return FALSE;
    }

  if (repeat_a == 1 && repeat_b == 1 && repeat_c == 1)
    return TRUE;

  gdis_model_discard_frames(model);
  old_atoms = model->atoms;
  old_bonds = model->bonds;
  old_atom_count = old_atoms ? old_atoms->len : 0u;
  old_explicit_count = MIN(model->explicit_bond_count,
                           old_bonds ? old_bonds->len : 0u);

  total_cells = (guint64) repeat_a * (guint64) repeat_b * (guint64) repeat_c;
  total_atoms = total_cells * (guint64) old_atom_count;
  if (total_atoms > (guint64) G_MAXUINT)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Requested supercell is too large.");
      return FALSE;
    }

  a_vec[0] = matrix[0];
  a_vec[1] = matrix[3];
  a_vec[2] = matrix[6];
  b_vec[0] = matrix[1];
  b_vec[1] = matrix[4];
  b_vec[2] = matrix[7];
  c_vec[0] = matrix[2];
  c_vec[1] = matrix[5];
  c_vec[2] = matrix[8];

  new_atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  for (ia = 0; ia < repeat_a; ia++)
    {
      for (ib = 0; ib < repeat_b; ib++)
        {
          for (ic = 0; ic < repeat_c; ic++)
            {
              gdouble shift[3];
              guint atom_index;

              shift[0] = (gdouble) ia * a_vec[0] + (gdouble) ib * b_vec[0] + (gdouble) ic * c_vec[0];
              shift[1] = (gdouble) ia * a_vec[1] + (gdouble) ib * b_vec[1] + (gdouble) ic * c_vec[1];
              shift[2] = (gdouble) ia * a_vec[2] + (gdouble) ib * b_vec[2] + (gdouble) ic * c_vec[2];

              for (atom_index = 0; atom_index < old_atom_count; atom_index++)
                {
                  const GdisAtom *src;
                  GdisAtom *copy;

                  src = g_ptr_array_index(old_atoms, atom_index);
                  copy = gdis_atom_new(src->label,
                                       src->element,
                                       src->ff_type,
                                       src->position[0] + shift[0],
                                       src->position[1] + shift[1],
                                       src->position[2] + shift[2],
                                       src->occupancy,
                                       src->region,
                                       new_atoms->len + 1);
                  g_ptr_array_add(new_atoms, copy);
                }
            }
        }
    }

  new_bonds = g_array_new(FALSE, FALSE, sizeof(GdisBond));
  if (old_explicit_count > 0 && old_atom_count > 0)
    {
      guint64 cell_index;

      for (cell_index = 0; cell_index < total_cells; cell_index++)
        {
          guint64 cell_atom_offset;
          guint bond_index;

          cell_atom_offset = cell_index * (guint64) old_atom_count;
          for (bond_index = 0; bond_index < old_explicit_count; bond_index++)
            {
              GdisBond src_bond;
              GdisBond mapped_bond;

              src_bond = g_array_index(old_bonds, GdisBond, bond_index);
              if (src_bond.atom_index_a >= old_atom_count ||
                  src_bond.atom_index_b >= old_atom_count)
                continue;

              mapped_bond = src_bond;
              mapped_bond.atom_index_a = (guint) (cell_atom_offset + src_bond.atom_index_a);
              mapped_bond.atom_index_b = (guint) (cell_atom_offset + src_bond.atom_index_b);
              mapped_bond.inferred = FALSE;
              g_array_append_val(new_bonds, mapped_bond);
            }
        }
    }

  model->atoms = new_atoms;
  model->bonds = new_bonds;
  model->explicit_bond_count = (old_explicit_count > 0) ? new_bonds->len : 0u;

  g_ptr_array_free(old_atoms, TRUE);
  g_array_free(old_bonds, TRUE);

  if (dims >= 1)
    model->cell_lengths[0] *= (gdouble) repeat_a;
  if (dims >= 2)
    model->cell_lengths[1] *= (gdouble) repeat_b;
  if (dims >= 3)
    model->cell_lengths[2] *= (gdouble) repeat_c;

  if (model->explicit_bond_count == 0)
    gdis_model_infer_bonds(model);
  else
    gdis_model_refresh_counts(model);
  gdis_model_init_image_limits(model);
  gdis_model_finalize_metadata(model);
  return TRUE;
}

gboolean
gdis_model_build_surface_slab(const GdisModel *source,
                              gint h,
                              gint k,
                              gint l,
                              gdouble shift,
                              guint region_a,
                              guint region_b,
                              guint repeat_a,
                              guint repeat_b,
                              gdouble vacuum,
                              GdisModel **surface_out,
                              gchar **summary_out,
                              GError **error)
{
  GdisModel *surface;
  GHashTable *seen;
  gchar *path;
  gdouble source_matrix[9];
  gdouble source_inverse[9];
  gdouble source_a[3];
  gdouble source_b[3];
  gdouble source_c[3];
  gdouble normal[3];
  gdouble first_vec[3];
  gdouble second_vec[3];
  gdouble slab_a[3];
  gdouble slab_b[3];
  gdouble slab_c[3];
  gdouble slab_actual_matrix[9];
  gdouble slab_actual_inverse[9];
  gdouble slab_canonical_matrix[9];
  gdouble lengths[3];
  gdouble angles[3];
  gdouble d_spacing;
  gdouble slab_thickness;
  gdouble cell_depth;
  gdouble base_projection;
  gdouble start_projection;
  gdouble split_projection;
  gdouble end_projection;
  gdouble tolerance;
  gint first_indices[3];
  gint second_indices[3];
  guint total_layers;
  gint range_a;
  gint range_b;
  gint range_c;
  guint source_atom_count;

  g_return_val_if_fail(surface_out != NULL, FALSE);
  *surface_out = NULL;

  if (summary_out)
    *summary_out = NULL;

  if (!source)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Surface construction needs an active model.");
      return FALSE;
    }

  if (h == 0 && k == 0 && l == 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Surface construction needs a non-zero Miller index.");
      return FALSE;
    }

  if (!source->periodic || source->periodicity < 3)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Surface construction currently targets fully 3D periodic crystal models.");
      return FALSE;
    }

  if (repeat_a == 0 || repeat_b == 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Surface repeats must both be at least 1.");
      return FALSE;
    }

  if (!gdis_model_build_cell_matrix(source, source_matrix, source_inverse))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Surface construction requires a valid source unit cell.");
      return FALSE;
    }

  source_a[0] = source_matrix[0];
  source_a[1] = source_matrix[3];
  source_a[2] = source_matrix[6];
  source_b[0] = source_matrix[1];
  source_b[1] = source_matrix[4];
  source_b[2] = source_matrix[7];
  source_c[0] = source_matrix[2];
  source_c[1] = source_matrix[5];
  source_c[2] = source_matrix[8];

  normal[0] = source_inverse[0] * h + source_inverse[3] * k + source_inverse[6] * l;
  normal[1] = source_inverse[1] * h + source_inverse[4] * k + source_inverse[7] * l;
  normal[2] = source_inverse[2] * h + source_inverse[5] * k + source_inverse[8] * l;
  d_spacing = gdis_vec3_length(normal);
  if (d_spacing <= 1.0e-10)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not resolve a valid surface normal for the selected Miller index.");
      return FALSE;
    }
  d_spacing = 1.0 / d_spacing;
  if (!gdis_vec3_normalize(normal))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not normalize the surface normal.");
      return FALSE;
    }

  if (!gdis_surface_build_inplane_indices(h, k, l, first_indices, second_indices))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not construct in-plane lattice vectors for this Miller index.");
      return FALSE;
    }

  gdis_vec3_set(first_vec, 0.0, 0.0, 0.0);
  gdis_vec3_add_scaled(first_vec, source_a, (gdouble) first_indices[0]);
  gdis_vec3_add_scaled(first_vec, source_b, (gdouble) first_indices[1]);
  gdis_vec3_add_scaled(first_vec, source_c, (gdouble) first_indices[2]);

  gdis_vec3_set(second_vec, 0.0, 0.0, 0.0);
  gdis_vec3_add_scaled(second_vec, source_a, (gdouble) second_indices[0]);
  gdis_vec3_add_scaled(second_vec, source_b, (gdouble) second_indices[1]);
  gdis_vec3_add_scaled(second_vec, source_c, (gdouble) second_indices[2]);

  if (gdis_vec3_length(first_vec) <= 1.0e-8 ||
      gdis_vec3_length(second_vec) <= 1.0e-8)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "The derived surface plane vectors were degenerate.");
      return FALSE;
    }

  {
    gdouble handedness[3];

    gdis_vec3_cross(first_vec, second_vec, handedness);
    if (gdis_vec3_dot(handedness, normal) < 0.0)
      {
        for (guint axis = 0; axis < 3; axis++)
          second_vec[axis] = -second_vec[axis];
      }
  }

  gdis_vec3_copy(slab_a, first_vec);
  for (guint axis = 0; axis < 3; axis++)
    slab_a[axis] *= (gdouble) repeat_a;

  gdis_vec3_copy(slab_b, second_vec);
  for (guint axis = 0; axis < 3; axis++)
    slab_b[axis] *= (gdouble) repeat_b;

  total_layers = MAX(region_a + region_b, 1u);
  slab_thickness = (gdouble) total_layers * d_spacing;
  cell_depth = slab_thickness + MAX(vacuum, 0.0);
  gdis_vec3_copy(slab_c, normal);
  for (guint axis = 0; axis < 3; axis++)
    slab_c[axis] *= cell_depth;

  if (!gdis_surface_build_cell_from_vectors(slab_a, slab_b, slab_c, lengths, angles))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not convert the surface lattice vectors into a valid cell.");
      return FALSE;
    }

  slab_actual_matrix[0] = slab_a[0];
  slab_actual_matrix[1] = slab_b[0];
  slab_actual_matrix[2] = slab_c[0];
  slab_actual_matrix[3] = slab_a[1];
  slab_actual_matrix[4] = slab_b[1];
  slab_actual_matrix[5] = slab_c[1];
  slab_actual_matrix[6] = slab_a[2];
  slab_actual_matrix[7] = slab_b[2];
  slab_actual_matrix[8] = slab_c[2];
  if (!gdis_model_invert_matrix3(slab_actual_matrix, slab_actual_inverse))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not invert the constructed surface cell.");
      return FALSE;
    }

  path = gdis_model_surface_default_path(source, h, k, l);
  surface = gdis_model_new(path, GDIS_MODEL_FORMAT_CIF);
  g_free(path);

  g_clear_pointer(&surface->title, g_free);
  surface->title = g_strdup_printf("%s surface (%d %d %d)",
                                   source->title ? source->title : source->basename,
                                   h,
                                   k,
                                   l);
  g_clear_pointer(&surface->space_group, g_free);
  surface->space_group = g_strdup("P 1");
  surface->periodic = TRUE;
  surface->periodicity = 2;
  memcpy(surface->cell_lengths, lengths, sizeof(lengths));
  memcpy(surface->cell_angles, angles, sizeof(angles));

  if (!gdis_model_build_cell_matrix(surface, slab_canonical_matrix, NULL))
    {
      gdis_model_free(surface);
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not build the canonical slab cell for display.");
      return FALSE;
    }

  base_projection = G_MAXDOUBLE;
  source_atom_count = source->atoms ? source->atoms->len : 0u;
  for (guint atom_index = 0; atom_index < source_atom_count; atom_index++)
    {
      const GdisAtom *atom;
      gdouble projection;

      atom = g_ptr_array_index(source->atoms, atom_index);
      projection = gdis_vec3_dot(normal, atom->position);
      if (projection < base_projection)
        base_projection = projection;
    }

  if (base_projection == G_MAXDOUBLE)
    {
      gdis_model_free(surface);
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Surface construction could not find any atoms in the source model.");
      return FALSE;
    }

  start_projection = base_projection + shift * d_spacing;
  split_projection = start_projection + ((gdouble) region_a * d_spacing);
  end_projection = start_projection + slab_thickness;
  tolerance = MAX(1.0e-4, 0.02 * d_spacing);

  range_a = (gint) (repeat_a * ABS(first_indices[0]) +
                    repeat_b * ABS(second_indices[0]) +
                    total_layers + 4u);
  range_b = (gint) (repeat_a * ABS(first_indices[1]) +
                    repeat_b * ABS(second_indices[1]) +
                    total_layers + 4u);
  range_c = (gint) (repeat_a * ABS(first_indices[2]) +
                    repeat_b * ABS(second_indices[2]) +
                    total_layers + 4u);

  seen = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
  for (gint ia = -range_a; ia <= range_a; ia++)
    {
      for (gint ib = -range_b; ib <= range_b; ib++)
        {
          for (gint ic = -range_c; ic <= range_c; ic++)
            {
              for (guint atom_index = 0; atom_index < source_atom_count; atom_index++)
                {
                  const GdisAtom *atom;
                  gdouble candidate[3];
                  gdouble frac[3];
                  gdouble projection;
                  gdouble cart[3];
                  gint region;
                  g_autofree gchar *key = NULL;

                  atom = g_ptr_array_index(source->atoms, atom_index);
                  gdis_vec3_copy(candidate, atom->position);
                  gdis_vec3_add_scaled(candidate, source_a, (gdouble) ia);
                  gdis_vec3_add_scaled(candidate, source_b, (gdouble) ib);
                  gdis_vec3_add_scaled(candidate, source_c, (gdouble) ic);

                  projection = gdis_vec3_dot(normal, candidate);
                  if (projection < start_projection - tolerance ||
                      projection > end_projection + tolerance)
                    continue;

                  gdis_cart_to_frac(slab_actual_inverse, candidate, frac);
                  frac[0] -= floor(frac[0]);
                  frac[1] -= floor(frac[1]);
                  frac[2] = (projection - start_projection + MAX(vacuum, 0.0) * 0.5) / cell_depth;

                  if (frac[0] >= 1.0 - 1.0e-6)
                    frac[0] = 0.0;
                  if (frac[1] >= 1.0 - 1.0e-6)
                    frac[1] = 0.0;
                  if (frac[2] < -1.0e-6 || frac[2] > 1.0 + 1.0e-6)
                    continue;
                  frac[2] = CLAMP(frac[2], 0.0, 1.0 - 1.0e-6);

                  key = g_strdup_printf("%.5f|%.5f|%.5f|%s",
                                        floor(frac[0] * 100000.0 + 0.5) / 100000.0,
                                        floor(frac[1] * 100000.0 + 0.5) / 100000.0,
                                        floor(frac[2] * 100000.0 + 0.5) / 100000.0,
                                        atom->element ? atom->element : "X");
                  if (!g_hash_table_add(seen, g_steal_pointer(&key)))
                    continue;

                  gdis_frac_to_cart(slab_canonical_matrix, frac, cart);
                  region = atom->region;
                  if (region < 0)
                    region = (projection < split_projection - tolerance || region_b == 0u) ? 1 : 2;
                  g_ptr_array_add(surface->atoms,
                                  gdis_atom_new(atom->label,
                                                atom->element,
                                                atom->ff_type,
                                                cart[0],
                                                cart[1],
                                                cart[2],
                                                atom->occupancy,
                                                region,
                                                surface->atoms->len + 1));
                }
            }
        }
    }
  g_hash_table_unref(seen);

  if (!surface->atoms || surface->atoms->len == 0)
    {
      gdis_model_free(surface);
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Surface construction did not capture any atoms. Try a smaller shift or larger slab depth.");
      return FALSE;
    }

  surface->explicit_bond_count = 0;
  gdis_model_infer_bonds(surface);
  gdis_model_init_image_limits(surface);
  gdis_model_finalize_metadata(surface);

  if (summary_out)
    {
      *summary_out = g_strdup_printf(
        "Built slab from %s\n"
        "Miller index: (%d %d %d)\n"
        "Dhkl: %.4f A\n"
        "Layers: %u + %u\n"
        "In-plane repeats: %u x %u\n"
        "Vacuum: %.2f A\n"
        "Result: %u atoms, cell %.4f %.4f %.4f A\n",
        source->basename,
        h,
        k,
        l,
        d_spacing,
        region_a,
        region_b,
        repeat_a,
        repeat_b,
        MAX(vacuum, 0.0),
        surface->atom_count,
        surface->cell_lengths[0],
        surface->cell_lengths[1],
        surface->cell_lengths[2]);
    }

  *surface_out = surface;
  return TRUE;
}

GdisModelFormat
gdis_model_format_from_path(const char *path)
{
  const gchar *dot;
  gchar *lower;
  GdisModelFormat format;

  if (!path)
    return GDIS_MODEL_FORMAT_UNKNOWN;

  dot = strrchr(path, '.');
  if (!dot || dot[1] == '\0')
    return GDIS_MODEL_FORMAT_UNKNOWN;

  lower = g_ascii_strdown(dot + 1, -1);
  format = GDIS_MODEL_FORMAT_UNKNOWN;

  if (g_str_equal(lower, "xyz"))
    format = GDIS_MODEL_FORMAT_XYZ;
  else if (g_str_equal(lower, "pdb") || g_str_equal(lower, "ent"))
    format = GDIS_MODEL_FORMAT_PDB;
  else if (g_str_equal(lower, "arc"))
    format = GDIS_MODEL_FORMAT_ARC;
  else if (g_str_equal(lower, "car"))
    format = GDIS_MODEL_FORMAT_CAR;
  else if (g_str_equal(lower, "cif"))
    format = GDIS_MODEL_FORMAT_CIF;
  else if (g_str_equal(lower, "xml"))
    {
      g_autofree gchar *basename = g_path_get_basename(path);
      g_autofree gchar *basename_lower = g_ascii_strdown(basename ? basename : "", -1);

      if (g_str_has_prefix(basename_lower, "vasprun"))
        format = GDIS_MODEL_FORMAT_UNKNOWN;
      else
        format = GDIS_MODEL_FORMAT_QBOX_XML;
    }
  else if (g_str_equal(lower, "gin") || g_str_equal(lower, "res"))
    format = GDIS_MODEL_FORMAT_GULP_INPUT;
  else if (g_str_equal(lower, "got") || g_str_equal(lower, "gout"))
    format = GDIS_MODEL_FORMAT_GULP_OUTPUT;
  else if (g_str_equal(lower, "xtl"))
    format = GDIS_MODEL_FORMAT_XTL;
  else if (g_str_equal(lower, "qein"))
    format = GDIS_MODEL_FORMAT_QE_INPUT;
  else if (g_str_equal(lower, "qeout"))
    format = GDIS_MODEL_FORMAT_QE_OUTPUT;
  else if (g_str_equal(lower, "inp"))
    format = GDIS_MODEL_FORMAT_GAMESS_INPUT;

  g_free(lower);
  return format;
}

const char *
gdis_model_format_label(GdisModelFormat format)
{
  switch (format)
    {
    case GDIS_MODEL_FORMAT_XYZ:
      return "XYZ";
    case GDIS_MODEL_FORMAT_PDB:
      return "PDB";
    case GDIS_MODEL_FORMAT_ARC:
      return "BIOSYM ARC";
    case GDIS_MODEL_FORMAT_CAR:
      return "BIOSYM CAR";
    case GDIS_MODEL_FORMAT_CIF:
      return "CIF";
    case GDIS_MODEL_FORMAT_QBOX_XML:
      return "Qbox XML";
    case GDIS_MODEL_FORMAT_GULP_INPUT:
      return "GULP Input";
    case GDIS_MODEL_FORMAT_GULP_OUTPUT:
      return "GULP Output";
    case GDIS_MODEL_FORMAT_XTL:
      return "XTL";
    case GDIS_MODEL_FORMAT_QE_INPUT:
      return "Quantum ESPRESSO Input";
    case GDIS_MODEL_FORMAT_QE_OUTPUT:
      return "Quantum ESPRESSO Output";
    case GDIS_MODEL_FORMAT_GAMESS_INPUT:
      return "GAMESS Input";
    case GDIS_MODEL_FORMAT_UNKNOWN:
    default:
      return "Unknown";
    }
}

gboolean
gdis_model_reset_inferred_bonds(GdisModel *model)
{
  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(model->bonds != NULL, FALSE);

  if (model->bonds->len > model->explicit_bond_count)
    g_array_set_size(model->bonds, model->explicit_bond_count);

  gdis_model_refresh_counts(model);

  if (model->explicit_bond_count == 0)
    return gdis_model_infer_bonds(model);

  return TRUE;
}

guint
gdis_model_get_frame_count(const GdisModel *model)
{
  g_return_val_if_fail(model != NULL, 0u);

  if (model->frames && model->frames->len > 0)
    return model->frames->len;

  return model->atoms && model->atoms->len > 0 ? 1u : 0u;
}

guint
gdis_model_get_current_frame_index(const GdisModel *model)
{
  g_return_val_if_fail(model != NULL, 0u);

  return model->current_frame_index;
}

const char *
gdis_model_get_frame_title(const GdisModel *model, guint frame_index)
{
  GdisModelFrame *frame;

  g_return_val_if_fail(model != NULL, NULL);

  if (!model->frames || frame_index >= model->frames->len)
    {
      if (frame_index == 0u)
        return model->title;
      return NULL;
    }

  frame = g_ptr_array_index(model->frames, frame_index);
  return frame ? frame->title : NULL;
}

gboolean
gdis_model_set_frame_index(GdisModel *model, guint frame_index, GError **error)
{
  GdisModelFrame *frame;

  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->frames || model->frames->len == 0u)
    {
      if (frame_index == 0u)
        {
          model->current_frame_index = 0u;
          return TRUE;
        }

      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "The current model does not contain multiple animation frames.");
      return FALSE;
    }

  if (frame_index >= model->frames->len)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Requested frame %u is outside the available range 1..%u.",
                  frame_index + 1u,
                  model->frames->len);
      return FALSE;
    }

  frame = g_ptr_array_index(model->frames, frame_index);
  if (!frame)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Requested frame %u could not be loaded.",
                  frame_index + 1u);
      return FALSE;
    }

  gdis_model_frame_apply(model, frame);
  if (model->explicit_bond_count == 0u)
    gdis_model_reset_inferred_bonds(model);
  else
    gdis_model_finalize_metadata(model);
  model->current_frame_index = frame_index;
  return TRUE;
}

void
gdis_model_discard_frames(GdisModel *model)
{
  g_return_if_fail(model != NULL);

  if (model->frames)
    g_ptr_array_set_size(model->frames, 0u);
  g_clear_pointer(&model->frames, g_ptr_array_unref);
  model->current_frame_index = 0u;
}

static GdisModel *
gdis_model_new(const char *path, GdisModelFormat format)
{
  GdisModel *model;

  model = g_new0(GdisModel, 1);
  model->path = g_strdup(path);
  model->basename = g_path_get_basename(path);
  model->format = format;
  model->format_label = g_strdup(gdis_model_format_label(format));
  model->atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  model->bonds = g_array_new(FALSE, FALSE, sizeof(GdisBond));
  model->frames = NULL;
  model->current_frame_index = 0u;
  model->cell_angles[0] = 90.0;
  model->cell_angles[1] = 90.0;
  model->cell_angles[2] = 90.0;
  gdis_model_init_image_limits(model);

  return model;
}

static GdisModelFrame *
gdis_model_frame_new(void)
{
  GdisModelFrame *frame;

  frame = g_new0(GdisModelFrame, 1);
  frame->atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  frame->cell_angles[0] = 90.0;
  frame->cell_angles[1] = 90.0;
  frame->cell_angles[2] = 90.0;
  return frame;
}

static void
gdis_model_frame_free(gpointer data)
{
  GdisModelFrame *frame;

  frame = data;
  if (!frame)
    return;

  g_clear_pointer(&frame->title, g_free);
  g_clear_pointer(&frame->space_group, g_free);
  g_clear_pointer(&frame->atoms, g_ptr_array_unref);
  g_free(frame);
}

static GdisModelFrame *
gdis_model_frame_clone(const GdisModelFrame *frame)
{
  GdisModelFrame *copy;

  g_return_val_if_fail(frame != NULL, NULL);

  copy = gdis_model_frame_new();
  copy->title = g_strdup(frame->title);
  copy->space_group = g_strdup(frame->space_group);
  copy->periodic = frame->periodic;
  copy->periodicity = frame->periodicity;
  memcpy(copy->cell_lengths, frame->cell_lengths, sizeof(copy->cell_lengths));
  memcpy(copy->cell_angles, frame->cell_angles, sizeof(copy->cell_angles));

  if (frame->atoms)
    {
      for (guint i = 0; i < frame->atoms->len; i++)
        {
          const GdisAtom *atom;

          atom = g_ptr_array_index(frame->atoms, i);
          g_ptr_array_add(copy->atoms,
                          gdis_atom_new(atom->label,
                                        atom->element,
                                        atom->ff_type,
                                        atom->position[0],
                                        atom->position[1],
                                        atom->position[2],
                                        atom->occupancy,
                                        atom->region,
                                        atom->serial));
        }
    }

  return copy;
}

static gboolean
gdis_model_frame_bonds_need_reset(const GdisModel *model)
{
  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->bonds || !model->atoms)
    return FALSE;

  for (guint i = 0; i < model->bonds->len; i++)
    {
      const GdisBond *bond;

      bond = &g_array_index(model->bonds, GdisBond, i);
      if (bond->atom_index_a >= model->atoms->len ||
          bond->atom_index_b >= model->atoms->len)
        return TRUE;
    }

  return FALSE;
}

static void
gdis_model_frame_apply(GdisModel *model, const GdisModelFrame *frame)
{
  GPtrArray *atoms;

  g_return_if_fail(model != NULL);
  g_return_if_fail(frame != NULL);

  atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  if (frame->atoms)
    {
      for (guint i = 0; i < frame->atoms->len; i++)
        {
          const GdisAtom *atom;

          atom = g_ptr_array_index(frame->atoms, i);
          g_ptr_array_add(atoms,
                          gdis_atom_new(atom->label,
                                        atom->element,
                                        atom->ff_type,
                                        atom->position[0],
                                        atom->position[1],
                                        atom->position[2],
                                        atom->occupancy,
                                        atom->region,
                                        atom->serial));
        }
    }

  if (model->atoms)
    g_ptr_array_free(model->atoms, TRUE);
  model->atoms = atoms;

  g_clear_pointer(&model->title, g_free);
  model->title = g_strdup(frame->title);
  g_clear_pointer(&model->space_group, g_free);
  model->space_group = g_strdup(frame->space_group);
  model->periodic = frame->periodic;
  model->periodicity = frame->periodicity;
  memcpy(model->cell_lengths, frame->cell_lengths, sizeof(model->cell_lengths));
  memcpy(model->cell_angles, frame->cell_angles, sizeof(model->cell_angles));

  if (gdis_model_frame_bonds_need_reset(model))
    {
      if (model->bonds)
        g_array_set_size(model->bonds, 0u);
      model->explicit_bond_count = 0u;
    }

  gdis_model_finalize_metadata(model);
}

static void
gdis_model_init_image_limits(GdisModel *model)
{
  g_return_if_fail(model != NULL);

  model->image_limits[0] = 0;
  model->image_limits[1] = 1;
  model->image_limits[2] = 0;
  model->image_limits[3] = 1;
  model->image_limits[4] = 0;
  model->image_limits[5] = 1;
}

static void
gdis_model_wrap_atoms_to_cell(GPtrArray *atoms,
                              guint periodicity,
                              const gdouble matrix[9],
                              const gdouble inverse[9])
{
  guint dims;

  g_return_if_fail(atoms != NULL);
  g_return_if_fail(matrix != NULL);
  g_return_if_fail(inverse != NULL);

  dims = periodicity;
  if (dims == 0u)
    dims = 3u;
  dims = MIN(dims, 3u);

  for (guint i = 0; i < atoms->len; i++)
    {
      GdisAtom *atom;
      gdouble frac[3];

      atom = g_ptr_array_index(atoms, i);
      gdis_cart_to_frac(inverse, atom->position, frac);

      for (guint axis = 0; axis < dims; axis++)
        frac[axis] -= floor(frac[axis]);

      gdis_frac_to_cart(matrix, frac, atom->position);
    }
}

static gboolean
gdis_model_repack_periodic_components(GdisModel *model,
                                      const gdouble matrix[9],
                                      const gdouble inverse[9])
{
  static const gdouble compact_span_tolerance = 1.05;
  gdouble *unwrapped_positions;
  guint *component_ids;
  guint component_count;
  guint dims;
  guint atom_count;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(matrix != NULL, FALSE);
  g_return_val_if_fail(inverse != NULL, FALSE);

  if (!model->atoms || model->atoms->len == 0u)
    return TRUE;

  atom_count = model->atoms->len;
  dims = model->periodicity;
  if (dims == 0u)
    dims = 3u;
  dims = MIN(dims, 3u);

  component_ids = gdis_model_build_component_ids(model, &component_count);
  if (!component_ids)
    return FALSE;

  unwrapped_positions = g_new0(gdouble, atom_count * 3u);
  for (guint i = 0; i < atom_count; i++)
    {
      const GdisAtom *atom;

      atom = g_ptr_array_index(model->atoms, i);
      unwrapped_positions[3u * i] = atom->position[0];
      unwrapped_positions[3u * i + 1u] = atom->position[1];
      unwrapped_positions[3u * i + 2u] = atom->position[2];
    }

  for (guint component = 0; component < component_count; component++)
    {
      GQueue queue = G_QUEUE_INIT;
      gboolean *visited;
      gboolean found_start;
      gboolean component_is_compact;
      guint start_index;
      guint count;
      gdouble centroid[3];
      gdouble centroid_frac[3];
      gdouble wrapped_frac[3];
      gdouble shift_frac[3];
      gdouble shift_cart[3];
      gdouble min_frac[3];
      gdouble max_frac[3];

      visited = g_new0(gboolean, atom_count);
      found_start = FALSE;
      start_index = 0u;
      count = 0u;
      centroid[0] = 0.0;
      centroid[1] = 0.0;
      centroid[2] = 0.0;
      min_frac[0] = G_MAXDOUBLE;
      min_frac[1] = G_MAXDOUBLE;
      min_frac[2] = G_MAXDOUBLE;
      max_frac[0] = -G_MAXDOUBLE;
      max_frac[1] = -G_MAXDOUBLE;
      max_frac[2] = -G_MAXDOUBLE;

      for (guint i = 0; i < atom_count; i++)
        {
          if (component_ids[i] == component)
            {
              start_index = i;
              found_start = TRUE;
              break;
            }
        }

      if (!found_start)
        {
          g_free(visited);
          continue;
        }

      visited[start_index] = TRUE;
      g_queue_push_tail(&queue, GUINT_TO_POINTER(start_index));

      while (!g_queue_is_empty(&queue))
        {
          guint atom_index;

          atom_index = GPOINTER_TO_UINT(g_queue_pop_head(&queue));
          for (guint bond_index = 0; bond_index < model->bonds->len; bond_index++)
            {
              const GdisBond *bond;
              guint neighbor;

              bond = &g_array_index(model->bonds, GdisBond, bond_index);
              if (bond->atom_index_a == atom_index)
                neighbor = bond->atom_index_b;
              else if (bond->atom_index_b == atom_index)
                neighbor = bond->atom_index_a;
              else
                continue;

              if (neighbor >= atom_count ||
                  component_ids[neighbor] != component ||
                  visited[neighbor])
                continue;

              {
                const GdisAtom *from_atom;
                const GdisAtom *to_atom;
                gdouble delta[3];

                from_atom = g_ptr_array_index(model->atoms, atom_index);
                to_atom = g_ptr_array_index(model->atoms, neighbor);
                gdis_model_minimum_image_delta(model,
                                               matrix,
                                               inverse,
                                               from_atom->position,
                                               to_atom->position,
                                               delta);
                unwrapped_positions[3u * neighbor] = unwrapped_positions[3u * atom_index] + delta[0];
                unwrapped_positions[3u * neighbor + 1u] = unwrapped_positions[3u * atom_index + 1u] + delta[1];
                unwrapped_positions[3u * neighbor + 2u] = unwrapped_positions[3u * atom_index + 2u] + delta[2];
              }

              visited[neighbor] = TRUE;
              g_queue_push_tail(&queue, GUINT_TO_POINTER(neighbor));
            }
        }

      for (guint i = 0; i < atom_count; i++)
        {
          gdouble frac[3];

          if (component_ids[i] != component)
            continue;

          centroid[0] += unwrapped_positions[3u * i];
          centroid[1] += unwrapped_positions[3u * i + 1u];
          centroid[2] += unwrapped_positions[3u * i + 2u];
          gdis_cart_to_frac(inverse, &unwrapped_positions[3u * i], frac);
          for (guint axis = 0; axis < dims; axis++)
            {
              min_frac[axis] = MIN(min_frac[axis], frac[axis]);
              max_frac[axis] = MAX(max_frac[axis], frac[axis]);
            }
          count++;
        }

      if (count == 0u)
        {
          g_free(visited);
          continue;
        }

      centroid[0] /= (gdouble) count;
      centroid[1] /= (gdouble) count;
      centroid[2] /= (gdouble) count;

      component_is_compact = TRUE;
      for (guint axis = 0; axis < dims; axis++)
        {
          if ((max_frac[axis] - min_frac[axis]) > compact_span_tolerance)
            {
              component_is_compact = FALSE;
              break;
            }
        }

      if (component_is_compact)
        {
          gdis_cart_to_frac(inverse, centroid, centroid_frac);
          wrapped_frac[0] = centroid_frac[0];
          wrapped_frac[1] = centroid_frac[1];
          wrapped_frac[2] = centroid_frac[2];
          for (guint axis = 0; axis < dims; axis++)
            wrapped_frac[axis] -= floor(wrapped_frac[axis]);

          shift_frac[0] = wrapped_frac[0] - centroid_frac[0];
          shift_frac[1] = wrapped_frac[1] - centroid_frac[1];
          shift_frac[2] = wrapped_frac[2] - centroid_frac[2];
          gdis_frac_to_cart(matrix, shift_frac, shift_cart);

          for (guint i = 0; i < atom_count; i++)
            {
              GdisAtom *atom;

              if (component_ids[i] != component)
                continue;

              atom = g_ptr_array_index(model->atoms, i);
              atom->position[0] = unwrapped_positions[3u * i] + shift_cart[0];
              atom->position[1] = unwrapped_positions[3u * i + 1u] + shift_cart[1];
              atom->position[2] = unwrapped_positions[3u * i + 2u] + shift_cart[2];
            }
        }
      else
        {
          for (guint i = 0; i < atom_count; i++)
            {
              GdisAtom *atom;
              gdouble frac[3];

              if (component_ids[i] != component)
                continue;

              gdis_cart_to_frac(inverse, &unwrapped_positions[3u * i], frac);
              for (guint axis = 0; axis < dims; axis++)
                frac[axis] -= floor(frac[axis]);

              atom = g_ptr_array_index(model->atoms, i);
              gdis_frac_to_cart(matrix, frac, atom->position);
            }
        }

      g_free(visited);
    }

  g_free(unwrapped_positions);
  g_free(component_ids);

  return TRUE;
}

static void
gdis_model_normalize_loaded_periodic_geometry(GdisModel *model)
{
  gdouble matrix[9];
  gdouble inverse[9];

  g_return_if_fail(model != NULL);

  if (!model->periodic ||
      !model->atoms ||
      model->atoms->len == 0u ||
      !gdis_model_build_cell_matrix(model, matrix, inverse))
    return;

  if (!gdis_model_repack_periodic_components(model, matrix, inverse))
    gdis_model_wrap_atoms_to_cell(model->atoms,
                                  model->periodicity,
                                  matrix,
                                  inverse);

  if (model->frames)
    {
      for (guint i = 0; i < model->frames->len; i++)
        {
          GdisModelFrame *frame;

          frame = g_ptr_array_index(model->frames, i);
          if (!frame ||
              !frame->periodic ||
              !frame->atoms ||
              frame->atoms->len == 0u ||
              !gdis_build_cell_matrix_from_parameters(frame->cell_lengths,
                                                      frame->cell_angles,
                                                      matrix,
                                                      inverse))
            continue;

          gdis_model_wrap_atoms_to_cell(frame->atoms,
                                        frame->periodicity,
                                        matrix,
                                        inverse);
        }
    }
}

static guint *
gdis_model_build_component_ids(const GdisModel *model, guint *component_count_out)
{
  guint atom_count;
  guint *component_ids;
  guint component;
  guint i;

  g_return_val_if_fail(model != NULL, NULL);

  atom_count = model->atoms ? model->atoms->len : 0u;
  if (component_count_out)
    *component_count_out = 0;
  if (atom_count == 0)
    return NULL;

  component_ids = g_new(guint, atom_count);
  for (i = 0; i < atom_count; i++)
    component_ids[i] = G_MAXUINT;

  component = 0;
  for (i = 0; i < atom_count; i++)
    {
      GQueue queue = G_QUEUE_INIT;

      if (component_ids[i] != G_MAXUINT)
        continue;

      component_ids[i] = component;
      g_queue_push_tail(&queue, GUINT_TO_POINTER(i));
      while (!g_queue_is_empty(&queue))
        {
          guint atom_index;
          guint bond_index;

          atom_index = GPOINTER_TO_UINT(g_queue_pop_head(&queue));
          for (bond_index = 0; bond_index < model->bonds->len; bond_index++)
            {
              const GdisBond *bond;
              guint neighbor;

              bond = &g_array_index(model->bonds, GdisBond, bond_index);
              if (bond->atom_index_a == atom_index)
                neighbor = bond->atom_index_b;
              else if (bond->atom_index_b == atom_index)
                neighbor = bond->atom_index_a;
              else
                continue;

              if (neighbor >= atom_count || component_ids[neighbor] != G_MAXUINT)
                continue;

              component_ids[neighbor] = component;
              g_queue_push_tail(&queue, GUINT_TO_POINTER(neighbor));
            }
        }
      component++;
    }

  if (component_count_out)
    *component_count_out = component;
  return component_ids;
}

static void
gdis_model_minimum_image_delta(const GdisModel *model,
                               const gdouble matrix[9],
                               const gdouble inverse[9],
                               const gdouble from[3],
                               const gdouble to[3],
                               gdouble delta[3])
{
  gdouble from_frac[3];
  gdouble to_frac[3];
  gdouble diff_frac[3];
  guint dims;
  guint axis;

  g_return_if_fail(model != NULL);
  g_return_if_fail(matrix != NULL);
  g_return_if_fail(inverse != NULL);
  g_return_if_fail(from != NULL);
  g_return_if_fail(to != NULL);
  g_return_if_fail(delta != NULL);

  dims = model->periodicity;
  if (dims == 0)
    dims = 3;
  dims = MIN(dims, 3u);

  gdis_cart_to_frac(inverse, from, from_frac);
  gdis_cart_to_frac(inverse, to, to_frac);

  for (axis = 0; axis < 3; axis++)
    diff_frac[axis] = to_frac[axis] - from_frac[axis];

  for (axis = 0; axis < dims; axis++)
    diff_frac[axis] -= floor(diff_frac[axis] + 0.5);

  gdis_frac_to_cart(matrix, diff_frac, delta);
}

static gboolean
gdis_model_load_xyz(GdisModel *model, const gchar *contents, GError **error)
{
  gchar **lines;
  guint line_index;
  GPtrArray *frames;
  gboolean ok;

  lines = g_strsplit(contents, "\n", -1);
  frames = g_ptr_array_new_with_free_func(gdis_model_frame_free);
  ok = FALSE;
  line_index = 0u;

  while (lines[line_index] != NULL)
    {
      GdisModelFrame *frame;
      gchar *endptr;
      glong expected_atoms;
      guint parsed_atoms;

      g_strchomp(lines[line_index]);
      if (g_strstrip(lines[line_index])[0] == '\0')
        {
          line_index++;
          continue;
        }

      expected_atoms = g_ascii_strtoll(lines[line_index], &endptr, 10);
      if (endptr == lines[line_index] || expected_atoms < 0)
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_PARSE,
                      "XYZ atom count header is invalid at line %u in '%s'.",
                      line_index + 1u,
                      model->path);
          goto xyz_done;
        }

      if (lines[line_index + 1u] == NULL)
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_PARSE,
                      "XYZ file '%s' ended before the frame title line.",
                      model->path);
          goto xyz_done;
        }

      frame = gdis_model_frame_new();
      g_strchomp(lines[line_index + 1u]);
      frame->title = gdis_strdup_strip(lines[line_index + 1u]);
      if (!frame->title || frame->title[0] == '\0')
        {
          g_clear_pointer(&frame->title, g_free);
          frame->title = g_strdup_printf("%s frame %u",
                                         model->basename ? model->basename : "trajectory",
                                         frames->len + 1u);
        }

      line_index += 2u;
      parsed_atoms = 0u;
      while (lines[line_index] != NULL && parsed_atoms < (guint) expected_atoms)
        {
          gint token_count;
          gchar **tokens;
          gchar *element;
          gdouble x;
          gdouble y;
          gdouble z;

          g_strchomp(lines[line_index]);
          if (g_strstrip(lines[line_index])[0] == '\0')
            {
              line_index++;
              continue;
            }

          tokens = gdis_split_simple(lines[line_index], &token_count);
          if (token_count < 4)
            {
              g_strfreev(tokens);
              gdis_model_frame_free(frame);
              g_set_error(error,
                          GDIS_MODEL_ERROR,
                          GDIS_MODEL_ERROR_PARSE,
                          "XYZ atom line %u is incomplete in '%s'.",
                          line_index + 1u,
                          model->path);
              goto xyz_done;
            }

          if (!gdis_try_parse_double(tokens[1], &x) ||
              !gdis_try_parse_double(tokens[2], &y) ||
              !gdis_try_parse_double(tokens[3], &z))
            {
              g_strfreev(tokens);
              gdis_model_frame_free(frame);
              g_set_error(error,
                          GDIS_MODEL_ERROR,
                          GDIS_MODEL_ERROR_PARSE,
                          "XYZ atom line %u has invalid coordinates in '%s'.",
                          line_index + 1u,
                          model->path);
              goto xyz_done;
            }

          element = gdis_normalize_element_symbol(tokens[0]);
          g_ptr_array_add(frame->atoms,
                          gdis_atom_new(tokens[0],
                                        element,
                                        element,
                                        x,
                                        y,
                                        z,
                                        1.0,
                                        -1,
                                        parsed_atoms + 1u));
          g_free(element);
          g_strfreev(tokens);
          parsed_atoms++;
          line_index++;
        }

      if (parsed_atoms != (guint) expected_atoms)
        {
          gdis_model_frame_free(frame);
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_PARSE,
                      "XYZ atom count mismatch in '%s': expected %ld, parsed %u.",
                      model->path,
                      expected_atoms,
                      parsed_atoms);
          goto xyz_done;
        }

      g_ptr_array_add(frames, frame);
    }

  if (frames->len == 0u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "XYZ file '%s' is empty.",
                  model->path);
      goto xyz_done;
    }

  gdis_model_frame_apply(model, g_ptr_array_index(frames, 0u));
  model->current_frame_index = 0u;
  if (frames->len > 1u)
    {
      model->frames = frames;
      frames = NULL;
    }
  ok = TRUE;

xyz_done:
  if (frames)
    g_ptr_array_free(frames, TRUE);
  g_strfreev(lines);
  return ok;
}

static gboolean
gdis_model_load_pdb(GdisModel *model, const gchar *contents, GError **error)
{
  gchar **lines;
  GHashTable *serial_to_index;
  GHashTable *seen_bonds;
  guint i;

  lines = g_strsplit(contents, "\n", -1);
  serial_to_index = g_hash_table_new(g_direct_hash, g_direct_equal);
  seen_bonds = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);

  for (i = 0; lines[i] != NULL; i++)
    {
      gchar *line;

      line = lines[i];
      g_strchomp(line);

      if (g_str_has_prefix(line, "ENDMDL") || g_str_has_prefix(line, "END"))
        break;

      if (g_str_has_prefix(line, "TITLE"))
        {
          gchar *chunk;
          gchar *joined;

          chunk = gdis_copy_field(line, 10, MAX(0, (gint) strlen(line) - 10));
          g_strstrip(chunk);
          if (chunk[0] != '\0')
            {
              if (!model->title || model->title[0] == '\0')
                {
                  g_clear_pointer(&model->title, g_free);
                  model->title = g_strdup(chunk);
                }
              else
                {
                  joined = g_strjoin(" ", model->title, chunk, NULL);
                  g_free(model->title);
                  model->title = joined;
                }
            }
          g_free(chunk);
          continue;
        }

      if (g_str_has_prefix(line, "CRYST1"))
        {
          gchar *field;

          field = gdis_copy_field(line, 6, 9);
          model->cell_lengths[0] = g_ascii_strtod(field, NULL);
          g_free(field);

          field = gdis_copy_field(line, 15, 9);
          model->cell_lengths[1] = g_ascii_strtod(field, NULL);
          g_free(field);

          field = gdis_copy_field(line, 24, 9);
          model->cell_lengths[2] = g_ascii_strtod(field, NULL);
          g_free(field);

          field = gdis_copy_field(line, 33, 7);
          model->cell_angles[0] = g_ascii_strtod(field, NULL);
          g_free(field);

          field = gdis_copy_field(line, 40, 7);
          model->cell_angles[1] = g_ascii_strtod(field, NULL);
          g_free(field);

          field = gdis_copy_field(line, 47, 7);
          model->cell_angles[2] = g_ascii_strtod(field, NULL);
          g_free(field);

          field = gdis_copy_field(line, 55, 11);
          g_strstrip(field);
          if (field[0] != '\0')
            {
              g_clear_pointer(&model->space_group, g_free);
              model->space_group = g_strdup(field);
            }
          g_free(field);

          if (model->cell_lengths[0] > 1.0 &&
              model->cell_lengths[1] > 1.0 &&
              model->cell_lengths[2] > 1.0)
            {
              model->periodic = TRUE;
              model->periodicity = 3;
            }
          continue;
        }

      if (g_str_has_prefix(line, "ATOM  ") || g_str_has_prefix(line, "HETATM"))
        {
          gchar *serial_text;
          gchar *atom_name;
          gchar *element_text;
          gchar *label;
          gchar *element;
          gchar *field;
          gdouble x;
          gdouble y;
          gdouble z;
          gdouble occupancy;
          guint serial;

          serial_text = gdis_copy_field(line, 6, 5);
          if (!gdis_try_parse_uint(serial_text, &serial))
            serial = model->atoms->len + 1;
          g_free(serial_text);

          atom_name = gdis_copy_field(line, 12, 4);
          label = gdis_strdup_strip(atom_name);

          element_text = NULL;
          if (strlen(line) >= 78)
            element_text = gdis_copy_field(line, 76, 2);

          if (element_text)
            g_strstrip(element_text);

          if (!element_text || element_text[0] == '\0')
            element = gdis_pdb_guess_element(atom_name);
          else
            element = gdis_normalize_element_symbol(element_text);

          field = gdis_copy_field(line, 30, 8);
          x = g_ascii_strtod(field, NULL);
          g_free(field);

          field = gdis_copy_field(line, 38, 8);
          y = g_ascii_strtod(field, NULL);
          g_free(field);

          field = gdis_copy_field(line, 46, 8);
          z = g_ascii_strtod(field, NULL);
          g_free(field);

          occupancy = 1.0;
          if (strlen(line) >= 60)
            {
              field = gdis_copy_field(line, 54, 6);
              if (g_strstrip(field)[0] != '\0')
                occupancy = g_ascii_strtod(field, NULL);
              g_free(field);
            }

          g_ptr_array_add(model->atoms,
                          gdis_atom_new(label[0] != '\0' ? label : atom_name,
                                        element,
                                        element,
                                        x,
                                        y,
                                        z,
                                        occupancy,
                                        -1,
                                        serial));
          g_hash_table_insert(serial_to_index,
                              GUINT_TO_POINTER(serial),
                              GUINT_TO_POINTER(model->atoms->len));

          g_free(atom_name);
          g_free(element_text);
          g_free(label);
          g_free(element);
          continue;
        }

      if (g_str_has_prefix(line, "CONECT"))
        {
          gint token_count;
          gchar **tokens;
          guint source_serial;
          gpointer source_index_ptr;
          gint j;

          tokens = gdis_split_simple(line, &token_count);
          if (token_count < 2 || !gdis_try_parse_uint(tokens[1], &source_serial))
            {
              g_strfreev(tokens);
              continue;
            }

          source_index_ptr = g_hash_table_lookup(serial_to_index, GUINT_TO_POINTER(source_serial));
          if (!source_index_ptr)
            {
              g_strfreev(tokens);
              continue;
            }

          for (j = 2; j < token_count; j++)
            {
              guint target_serial;
              gpointer target_index_ptr;
              guint source_index;
              guint target_index;
              guint min_index;
              guint max_index;
              gchar *bond_key;

              if (!gdis_try_parse_uint(tokens[j], &target_serial))
                continue;

              target_index_ptr = g_hash_table_lookup(serial_to_index, GUINT_TO_POINTER(target_serial));
              if (!target_index_ptr)
                continue;

              source_index = GPOINTER_TO_UINT(source_index_ptr) - 1;
              target_index = GPOINTER_TO_UINT(target_index_ptr) - 1;
              min_index = MIN(source_index, target_index);
              max_index = MAX(source_index, target_index);
              bond_key = g_strdup_printf("%u:%u", min_index, max_index);

              if (g_hash_table_add(seen_bonds, bond_key))
                gdis_model_add_bond(model, min_index, max_index, 1, FALSE);
            }

          g_strfreev(tokens);
        }
    }

  g_hash_table_destroy(seen_bonds);
  g_hash_table_destroy(serial_to_index);
  g_strfreev(lines);

  if (model->atoms->len == 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No ATOM/HETATM records were found in '%s'.",
                  model->path);
      return FALSE;
    }

  return TRUE;
}

static gboolean
gdis_model_load_arc_like(GdisModel *model,
                         const gchar *contents,
                         GdisModelFormat format,
                         GError **error)
{
  gchar **lines;
  GPtrArray *frames;
  GdisModelFrame *current_frame;
  gchar *header_title;
  gboolean header_pbc_2d;
  gboolean header_periodic;
  guint header_periodicity;
  guint i;
  gboolean ok;

  lines = g_strsplit(contents, "\n", -1);
  frames = g_ptr_array_new_with_free_func(gdis_model_frame_free);
  current_frame = NULL;
  header_title = NULL;
  header_pbc_2d = FALSE;
  header_periodic = FALSE;
  header_periodicity = 0u;
  ok = FALSE;

  for (i = 0; lines[i] != NULL; i++)
    {
      gchar *line;
      gint token_count;
      gchar **tokens;

      line = lines[i];
      g_strchomp(line);
      g_strstrip(line);

      if (line[0] == '\0')
        continue;

      if (!current_frame)
        {
          if (g_ascii_strncasecmp(line, "!BIOSYM", 7) == 0 ||
              g_ascii_strncasecmp(line, "!MSI", 4) == 0)
            continue;

          if (g_ascii_strncasecmp(line, "PBC=", 4) == 0)
            {
              const gchar *value;

              value = line + 4;
              while (*value && g_ascii_isspace(*value))
                value++;

              if (g_ascii_strncasecmp(value, "2D", 2) == 0)
                {
                  header_periodic = TRUE;
                  header_periodicity = 2u;
                  header_pbc_2d = TRUE;
                }
              else if (g_ascii_strncasecmp(value, "ON", 2) == 0)
                {
                  header_periodic = TRUE;
                  header_periodicity = 3u;
                  header_pbc_2d = FALSE;
                }
              else if (g_ascii_strncasecmp(value, "OFF", 3) == 0)
                {
                  header_periodic = FALSE;
                  header_periodicity = 0u;
                  header_pbc_2d = FALSE;
                }
              continue;
            }

          if (g_ascii_strncasecmp(line, "!DATE", 5) == 0)
            {
              current_frame = gdis_model_frame_new();
              current_frame->periodic = header_periodic;
              current_frame->periodicity = header_periodicity;
              current_frame->title = g_strdup((header_title && header_title[0] != '\0')
                                              ? header_title
                                              : (model->basename ? model->basename : "model"));
              current_frame->cell_angles[0] = 90.0;
              current_frame->cell_angles[1] = 90.0;
              current_frame->cell_angles[2] = 90.0;
              continue;
            }

          if (line[0] != '!')
            {
              g_clear_pointer(&header_title, g_free);
              header_title = gdis_strdup_strip(line);
            }
          continue;
        }

      if (g_ascii_strncasecmp(line, "end", 3) == 0)
        {
          if (current_frame->atoms->len > 0u)
            {
              if (!current_frame->title || current_frame->title[0] == '\0')
                {
                  g_clear_pointer(&current_frame->title, g_free);
                  current_frame->title = g_strdup_printf("%s frame %u",
                                                         model->basename ? model->basename : "model",
                                                         frames->len + 1u);
                }
              g_ptr_array_add(frames, current_frame);
            }
          else
            {
              gdis_model_frame_free(current_frame);
            }
          current_frame = NULL;
          continue;
        }

      if (g_ascii_strncasecmp(line, "PBC", 3) == 0)
        {
          const gchar *paren_start;
          const gchar *paren_end;

          tokens = gdis_split_simple(line, &token_count);
          paren_start = strchr(line, '(');
          if (paren_start && current_frame)
            {
              gchar *space_group;

              paren_end = strchr(paren_start + 1, ')');
              if (paren_end && paren_end > paren_start + 1)
                {
                  space_group = g_strndup(paren_start + 1, paren_end - paren_start - 1);
                  g_strstrip(space_group);
                  if (space_group[0] != '\0')
                    {
                      g_clear_pointer(&current_frame->space_group, g_free);
                      current_frame->space_group = space_group;
                    }
                  else
                    {
                      g_free(space_group);
                    }
                }
            }

          if (token_count >= 7)
            {
              current_frame->cell_lengths[0] = g_ascii_strtod(tokens[1], NULL);
              current_frame->cell_lengths[1] = g_ascii_strtod(tokens[2], NULL);
              current_frame->cell_lengths[2] = g_ascii_strtod(tokens[3], NULL);
              current_frame->cell_angles[0] = g_ascii_strtod(tokens[4], NULL);
              current_frame->cell_angles[1] = g_ascii_strtod(tokens[5], NULL);
              current_frame->cell_angles[2] = g_ascii_strtod(tokens[6], NULL);
              current_frame->periodic = TRUE;
              current_frame->periodicity =
                (header_pbc_2d || current_frame->cell_lengths[2] <= 0.0) ? 2u : 3u;
              if (current_frame->periodicity == 2)
                {
                  current_frame->cell_lengths[2] = 1.0;
                  current_frame->cell_angles[0] = 90.0;
                  current_frame->cell_angles[1] = 90.0;
                }
            }
          else if (token_count >= 4)
            {
              current_frame->cell_lengths[0] = g_ascii_strtod(tokens[1], NULL);
              current_frame->cell_lengths[1] = g_ascii_strtod(tokens[2], NULL);
              current_frame->cell_lengths[2] = 1.0;
              current_frame->cell_angles[0] = 90.0;
              current_frame->cell_angles[1] = 90.0;
              current_frame->cell_angles[2] = g_ascii_strtod(tokens[3], NULL);
              current_frame->periodic = TRUE;
              current_frame->periodicity = 2u;
            }
          g_strfreev(tokens);
          continue;
        }

      tokens = gdis_split_simple(line, &token_count);
      if (token_count >= 8)
        {
          gchar *element;
          gchar *ff_type;
          gdouble x;
          gdouble y;
          gdouble z;
          gdouble occupancy;
          gint region;
          gboolean skip_atom;

          if (!gdis_try_parse_double(tokens[1], &x) ||
              !gdis_try_parse_double(tokens[2], &y) ||
              !gdis_try_parse_double(tokens[3], &z))
            {
              g_strfreev(tokens);
              continue;
            }

          skip_atom = FALSE;
          region = -1;
          if (g_ascii_strcasecmp(tokens[4], "SHELL") == 0)
            skip_atom = TRUE;
          else if (gdis_arc_like_type_is_marvin_label(tokens[4]))
            {
              region = gdis_arc_like_region_from_type(tokens[4]);
              if (tokens[4][3] == 'S' || tokens[4][3] == 's')
                skip_atom = TRUE;
            }

          if (skip_atom)
            {
              g_strfreev(tokens);
              continue;
            }

          if (token_count >= 9)
            occupancy = 1.0;
          else
            occupancy = 1.0;

          ff_type = NULL;
          if (token_count >= 7)
            ff_type = gdis_strdup_strip(tokens[6]);

          if (token_count >= 8)
            element = gdis_normalize_element_symbol(tokens[7]);
          else
            element = gdis_normalize_element_symbol(tokens[0]);

          g_ptr_array_add(current_frame->atoms,
                          gdis_atom_new(tokens[0],
                                        element,
                                        (ff_type && ff_type[0] != '\0') ? ff_type : element,
                                        x,
                                        y,
                                        z,
                                        occupancy,
                                        region,
                                        current_frame->atoms->len + 1u));
          g_free(ff_type);
          g_free(element);
        }
      g_strfreev(tokens);
    }

  if (current_frame)
    {
      if (current_frame->atoms->len > 0u)
        g_ptr_array_add(frames, current_frame);
      else
        gdis_model_frame_free(current_frame);
      current_frame = NULL;
    }

  if (frames->len == 0u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No atoms were parsed from any %s frame in '%s'.",
                  gdis_model_format_label(format),
                  model->path);
      goto arc_done;
    }

  gdis_model_frame_apply(model, g_ptr_array_index(frames, 0u));
  model->current_frame_index = 0u;
  if (frames->len > 1u)
    {
      model->frames = frames;
      frames = NULL;
    }
  ok = TRUE;

arc_done:
  g_strfreev(lines);
  g_free(header_title);
  if (frames)
    g_ptr_array_free(frames, TRUE);

  return ok;
}

static gboolean
gdis_model_load_gulp_input_like(GdisModel *model,
                                const gchar *contents,
                                GdisModelFormat format,
                                GError **error)
{
  typedef enum
  {
    GDIS_GULP_COORD_NONE = 0,
    GDIS_GULP_COORD_CARTESIAN,
    GDIS_GULP_COORD_FRACTIONAL,
    GDIS_GULP_COORD_SURFACE_FRACTIONAL
  } GdisGulpCoordMode;

  gchar **lines;
  GdisGulpCoordMode coord_mode;
  gint current_region;
  gboolean seen_name;
  gboolean have_cell_matrix;
  gboolean have_surface_vectors;
  gdouble cell_matrix[9];
  gdouble cell_inverse[9];
  gdouble surface_a[3];
  gdouble surface_b[3];

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(contents != NULL, FALSE);

  lines = g_strsplit(contents, "\n", -1);
  coord_mode = GDIS_GULP_COORD_NONE;
  current_region = -1;
  seen_name = FALSE;
  have_cell_matrix = FALSE;
  have_surface_vectors = FALSE;
  gdis_vec3_set(surface_a, 0.0, 0.0, 0.0);
  gdis_vec3_set(surface_b, 0.0, 0.0, 0.0);

  for (guint i = 0; lines[i] != NULL; i++)
    {
      gchar *line;
      gchar *trimmed;
      g_auto(GStrv) tokens = NULL;
      gint token_count = 0;
      g_autofree gchar *first_lower = NULL;

      line = lines[i];
      g_strchomp(line);
      trimmed = g_strstrip(line);
      if (trimmed[0] == '\0' || trimmed[0] == '#')
        continue;

      tokens = gdis_split_simple(trimmed, &token_count);
      if (token_count <= 0)
        continue;

      first_lower = g_ascii_strdown(tokens[0], -1);

      if (g_str_equal(first_lower, "name"))
        {
          const gchar *name_text;

          if (seen_name && model->atoms->len > 0u)
            break;

          seen_name = TRUE;
          name_text = trimmed + strlen(tokens[0]);
          while (*name_text && g_ascii_isspace(*name_text))
            name_text++;
          if (*name_text)
            {
              g_clear_pointer(&model->title, g_free);
              model->title = g_strdup(name_text);
            }
          continue;
        }

      if (g_str_equal(first_lower, "title"))
        {
          const gchar *title_text;

          title_text = trimmed + strlen(tokens[0]);
          while (*title_text && g_ascii_isspace(*title_text))
            title_text++;

          if (*title_text)
            {
              g_clear_pointer(&model->title, g_free);
              model->title = g_strdup(title_text);
            }
          else
            {
              guint j;

              for (j = i + 1; lines[j] != NULL; j++)
                {
                  gchar *title_line;
                  gchar *title_trimmed;

                  title_line = lines[j];
                  g_strchomp(title_line);
                  title_trimmed = g_strstrip(title_line);
                  if (title_trimmed[0] == '\0')
                    continue;
                  if (g_ascii_strcasecmp(title_trimmed, "end") == 0)
                    break;

                  g_clear_pointer(&model->title, g_free);
                  model->title = g_strdup(title_trimmed);
                  break;
                }
            }
          continue;
        }

      if (g_str_equal(first_lower, "space"))
        {
          if (token_count > 1)
            {
              const gchar *space_text;

              space_text = trimmed + strlen(tokens[0]);
              while (*space_text && g_ascii_isspace(*space_text))
                space_text++;
              if (*space_text)
                {
                  g_clear_pointer(&model->space_group, g_free);
                  model->space_group = g_strdup(space_text);
                }
            }
          else if (lines[i + 1] != NULL)
            {
              gchar *space_line;
              gchar *space_trimmed;

              space_line = lines[i + 1];
              g_strchomp(space_line);
              space_trimmed = g_strstrip(space_line);
              if (space_trimmed[0] != '\0')
                {
                  g_clear_pointer(&model->space_group, g_free);
                  model->space_group = g_strdup(space_trimmed);
                }
            }
          continue;
        }

      if (g_str_equal(first_lower, "cell") || g_str_equal(first_lower, "scell"))
        {
          gdouble values[6] = {0.0, 0.0, 0.0, 90.0, 90.0, 90.0};
          const gboolean is_surface_cell = g_str_equal(first_lower, "scell");
          const guint required = is_surface_cell ? 3u : 6u;
          guint found;
          guint line_index;

          found = gdis_collect_doubles_from_line(trimmed, values, required);
          line_index = i;
          while (found < required && lines[line_index + 1] != NULL)
            {
              guint added;
              gchar *next_line;
              gchar *next_trimmed;

              line_index++;
              next_line = lines[line_index];
              g_strchomp(next_line);
              next_trimmed = g_strstrip(next_line);
              if (next_trimmed[0] == '\0' || next_trimmed[0] == '#')
                continue;

              added = gdis_collect_doubles_from_line(next_trimmed,
                                                     values + found,
                                                     required - found);
              if (added == 0u)
                {
                  line_index--;
                  break;
                }
              found += added;
            }
          i = line_index;

          if (found >= required)
            {
              if (is_surface_cell)
                {
                  model->cell_lengths[0] = values[0];
                  model->cell_lengths[1] = values[1];
                  model->cell_lengths[2] = 1.0;
                  model->cell_angles[0] = 90.0;
                  model->cell_angles[1] = 90.0;
                  model->cell_angles[2] = values[2];
                  model->periodic = TRUE;
                  model->periodicity = 2u;
                }
              else
                {
                  if (values[2] <= 1.0e-8)
                    {
                      model->cell_lengths[0] = values[0];
                      model->cell_lengths[1] = values[1];
                      model->cell_lengths[2] = 1.0;
                      model->cell_angles[0] = 90.0;
                      model->cell_angles[1] = 90.0;
                      model->cell_angles[2] = values[5];
                      model->periodic = TRUE;
                      model->periodicity = 2u;
                    }
                  else
                    {
                      memcpy(model->cell_lengths, values, 3u * sizeof(gdouble));
                      memcpy(model->cell_angles, values + 3, 3u * sizeof(gdouble));
                      model->periodic = TRUE;
                      model->periodicity = 3u;
                    }
                }

              if (g_str_equal(first_lower, "scell"))
                {
                  model->periodic = TRUE;
                  model->periodicity = 2u;
                }
              else if (values[2] > 1.0e-8)
                {
                  model->periodic = TRUE;
                  model->periodicity = 3u;
                }

              have_cell_matrix = gdis_model_build_cell_matrix(model,
                                                              cell_matrix,
                                                              cell_inverse);
            }
          continue;
        }

      if (g_str_equal(first_lower, "vectors"))
        {
          gdouble a_vec[3];
          gdouble b_vec[3];
          gdouble c_vec[3];
          gboolean have_vectors;
          guint parsed = 0u;

          have_vectors = TRUE;
          for (guint j = 0; j < 3u; j++)
            {
              gdouble *target;
              gdouble values[3];

              if (lines[i + 1] == NULL)
                {
                  have_vectors = FALSE;
                  break;
                }

              i++;
              target = (j == 0u) ? a_vec : (j == 1u) ? b_vec : c_vec;
              parsed = gdis_collect_doubles_from_line(lines[i], values, 3u);
              if (parsed < 3u)
                {
                  have_vectors = FALSE;
                  break;
                }
              target[0] = values[0];
              target[1] = values[1];
              target[2] = values[2];
            }

          if (have_vectors &&
              gdis_surface_build_cell_from_vectors(a_vec,
                                                   b_vec,
                                                   c_vec,
                                                   model->cell_lengths,
                                                   model->cell_angles))
            {
              model->periodic = TRUE;
              model->periodicity = 3u;
              have_cell_matrix = gdis_model_build_cell_matrix(model,
                                                              cell_matrix,
                                                              cell_inverse);
            }
          continue;
        }

      if (g_str_equal(first_lower, "svectors"))
        {
          gdouble values[3];
          gdouble c_vec[3];

          if (lines[i + 1] != NULL && lines[i + 2] != NULL)
            {
              guint first_count;
              guint second_count;

              i++;
              first_count = gdis_collect_doubles_from_line(lines[i], values, 3u);
              if (first_count >= 2u)
                {
                  surface_a[0] = values[0];
                  surface_a[1] = values[1];
                  surface_a[2] = (first_count >= 3u) ? values[2] : 0.0;
                }

              i++;
              second_count = gdis_collect_doubles_from_line(lines[i], values, 3u);
              if (second_count >= 2u)
                {
                  surface_b[0] = values[0];
                  surface_b[1] = values[1];
                  surface_b[2] = (second_count >= 3u) ? values[2] : 0.0;
                }

              if (first_count >= 2u && second_count >= 2u)
                {
                  have_surface_vectors = TRUE;
                  c_vec[0] = 0.0;
                  c_vec[1] = 0.0;
                  c_vec[2] = 1.0;
                  if (gdis_surface_build_cell_from_vectors(surface_a,
                                                           surface_b,
                                                           c_vec,
                                                           model->cell_lengths,
                                                           model->cell_angles))
                    {
                      model->periodic = TRUE;
                      model->periodicity = 2u;
                      have_cell_matrix = gdis_model_build_cell_matrix(model,
                                                                      cell_matrix,
                                                                      cell_inverse);
                    }
                }
            }
          continue;
        }

      if (g_str_equal(first_lower, "region"))
        {
          guint parsed_region;

          if (token_count > 1 && gdis_try_parse_uint(tokens[1], &parsed_region))
            current_region = (gint) parsed_region;
          continue;
        }

      if (g_str_equal(first_lower, "cartesian") || g_str_equal(first_lower, "cart"))
        {
          coord_mode = GDIS_GULP_COORD_CARTESIAN;
          continue;
        }

      if (g_str_equal(first_lower, "fractional") || g_str_equal(first_lower, "frac"))
        {
          coord_mode = GDIS_GULP_COORD_FRACTIONAL;
          continue;
        }

      if (g_str_equal(first_lower, "sfractional") || g_str_equal(first_lower, "sfrac"))
        {
          coord_mode = GDIS_GULP_COORD_SURFACE_FRACTIONAL;
          for (gint j = 1; j + 1 < token_count; j++)
            {
              g_autofree gchar *keyword = g_ascii_strdown(tokens[j], -1);
              guint parsed_region;

              if (!g_str_equal(keyword, "region"))
                continue;
              if (gdis_try_parse_uint(tokens[j + 1], &parsed_region))
                current_region = (gint) parsed_region;
            }
          continue;
        }

      if (coord_mode != GDIS_GULP_COORD_NONE && token_count >= 4)
        {
          g_autofree gchar *type = NULL;
          gdouble values[3];
          gdouble cart[3];
          gboolean is_core = FALSE;
          gboolean is_shell = FALSE;
          gboolean has_explicit_type = FALSE;
          gboolean exact_coords = FALSE;
          gint coord_start = 1;
          g_autofree gchar *element = NULL;
          gint region;

          if (g_str_equal(first_lower, "library") ||
              g_str_equal(first_lower, "species") ||
              g_str_equal(first_lower, "print") ||
              g_str_equal(first_lower, "output") ||
              g_str_equal(first_lower, "switch") ||
              g_str_equal(first_lower, "maxcyc") ||
              g_str_equal(first_lower, "title") ||
              g_str_equal(first_lower, "name") ||
              g_str_equal(first_lower, "harm") ||
              g_str_equal(first_lower, "harmonic") ||
              g_str_equal(first_lower, "three") ||
              g_str_equal(first_lower, "torsion") ||
              g_str_equal(first_lower, "morse") ||
              g_str_equal(first_lower, "buck") ||
              g_str_equal(first_lower, "spring"))
            {
              coord_mode = GDIS_GULP_COORD_NONE;
              continue;
            }

          if (token_count >= 5)
            {
              type = g_ascii_strdown(tokens[1], -1);
              is_core = g_str_equal(type, "core") ||
                        g_str_equal(tokens[1], "c") ||
                        g_str_has_prefix(type, "bcor");
              is_shell = g_str_has_prefix(type, "shel") ||
                         g_str_equal(tokens[1], "s") ||
                         g_str_has_prefix(type, "bshe");
              has_explicit_type = is_core || is_shell;
              if (has_explicit_type)
                coord_start = 2;
            }

          if (is_shell)
            continue;

          if (coord_start + 2 < token_count &&
              gdis_try_parse_double_relaxed(tokens[coord_start], &values[0]) &&
              gdis_try_parse_double_relaxed(tokens[coord_start + 1], &values[1]) &&
              gdis_try_parse_double_relaxed(tokens[coord_start + 2], &values[2]))
            {
              exact_coords = TRUE;
            }
          else if (!has_explicit_type &&
                   token_count >= 4 &&
                   gdis_try_parse_double_relaxed(tokens[1], &values[0]) &&
                   gdis_try_parse_double_relaxed(tokens[2], &values[1]) &&
                   gdis_try_parse_double_relaxed(tokens[3], &values[2]))
            {
              coord_start = 1;
              exact_coords = TRUE;
            }

          if (!exact_coords)
            continue;

          if (coord_mode == GDIS_GULP_COORD_CARTESIAN)
            {
              cart[0] = values[0];
              cart[1] = values[1];
              cart[2] = values[2];
            }
          else if (coord_mode == GDIS_GULP_COORD_FRACTIONAL)
            {
              if (!have_cell_matrix)
                continue;
              gdis_frac_to_cart(cell_matrix, values, cart);
            }
          else
            {
              if (have_surface_vectors)
                {
                  cart[0] = surface_a[0] * values[0] + surface_b[0] * values[1];
                  cart[1] = surface_a[1] * values[0] + surface_b[1] * values[1];
                  cart[2] = surface_a[2] * values[0] + surface_b[2] * values[1] + values[2];
                }
              else if (have_cell_matrix)
                {
                  cart[0] = cell_matrix[0] * values[0] + cell_matrix[1] * values[1];
                  cart[1] = cell_matrix[3] * values[0] + cell_matrix[4] * values[1];
                  cart[2] = values[2];
                }
              else
                {
                  continue;
                }
            }

          element = gdis_normalize_element_symbol(tokens[0]);
          region = (coord_mode == GDIS_GULP_COORD_SURFACE_FRACTIONAL) ? current_region : -1;
          g_ptr_array_add(model->atoms,
                          gdis_atom_new(tokens[0],
                                        element,
                                        element,
                                        cart[0],
                                        cart[1],
                                        cart[2],
                                        1.0,
                                        region,
                                        model->atoms->len + 1u));
        }
    }

  g_strfreev(lines);

  if (model->atoms->len == 0u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No atoms were parsed from %s '%s'.",
                  gdis_model_format_label(format),
                  model->path);
      return FALSE;
    }

  return TRUE;
}

static gboolean
gdis_model_load_gulp_output(GdisModel *model, const gchar *contents, GError **error)
{
  typedef enum
  {
    GDIS_GULP_BLOCK_NONE = 0,
    GDIS_GULP_BLOCK_CART,
    GDIS_GULP_BLOCK_FRACTIONAL,
    GDIS_GULP_BLOCK_SURFACE_MIXED
  } GdisGulpBlockMode;

  gchar **lines;
  GPtrArray *best_atoms;
  gboolean best_periodic;
  guint best_periodicity;
  gdouble best_lengths[3];
  gdouble best_angles[3];

  gboolean context_periodic;
  guint context_periodicity;
  gdouble context_lengths[3];
  gdouble context_angles[3];
  gboolean context_have_surface_vectors;
  gdouble context_surface_a[3];
  gdouble context_surface_b[3];

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(contents != NULL, FALSE);

  lines = g_strsplit(contents, "\n", -1);
  best_atoms = NULL;
  best_periodic = FALSE;
  best_periodicity = 0u;
  best_lengths[0] = best_lengths[1] = best_lengths[2] = 0.0;
  best_angles[0] = best_angles[1] = best_angles[2] = 90.0;

  context_periodic = FALSE;
  context_periodicity = 0u;
  context_lengths[0] = context_lengths[1] = context_lengths[2] = 0.0;
  context_angles[0] = context_angles[1] = context_angles[2] = 90.0;
  context_have_surface_vectors = FALSE;
  gdis_vec3_set(context_surface_a, 0.0, 0.0, 0.0);
  gdis_vec3_set(context_surface_b, 0.0, 0.0, 0.0);

  for (guint i = 0; lines[i] != NULL; i++)
    {
      gchar *line;
      gchar *trimmed;
      g_autofree gchar *lower = NULL;
      GdisGulpBlockMode mode;

      line = lines[i];
      g_strchomp(line);
      trimmed = g_strstrip(line);
      if (trimmed[0] == '\0')
        continue;

      lower = g_ascii_strdown(trimmed, -1);

      if (g_strstr_len(lower, -1, "cartesian lattice vectors (angstroms)") != NULL)
        {
          gdouble a_vec[3];
          gdouble b_vec[3];
          gdouble c_vec[3];
          gdouble values[3];
          guint scan_index;
          guint parsed_vectors;
          guint count_a;
          guint count_b;
          guint count_c;

          count_a = 0u;
          count_b = 0u;
          count_c = 0u;
          scan_index = i + 1u;
          parsed_vectors = 0u;
          while (lines[scan_index] != NULL && parsed_vectors < 3u)
            {
              gchar *vector_line;
              guint count;

              vector_line = g_strstrip(lines[scan_index]);
              if (vector_line[0] == '\0')
                {
                  scan_index++;
                  continue;
                }

              count = gdis_collect_doubles_from_line(vector_line, values, 3u);
              if (count < 3u)
                break;

              if (parsed_vectors == 0u)
                {
                  a_vec[0] = values[0];
                  a_vec[1] = values[1];
                  a_vec[2] = values[2];
                  count_a = count;
                }
              else if (parsed_vectors == 1u)
                {
                  b_vec[0] = values[0];
                  b_vec[1] = values[1];
                  b_vec[2] = values[2];
                  count_b = count;
                }
              else
                {
                  c_vec[0] = values[0];
                  c_vec[1] = values[1];
                  c_vec[2] = values[2];
                  count_c = count;
                }

              parsed_vectors++;
              scan_index++;
            }

          if (count_a >= 3u &&
              count_b >= 3u &&
              count_c >= 3u &&
              gdis_surface_build_cell_from_vectors(a_vec,
                                                   b_vec,
                                                   c_vec,
                                                   context_lengths,
                                                   context_angles))
            {
              context_periodic = TRUE;
              context_periodicity = 3u;
            }

          if (scan_index > i + 1u)
            i = scan_index - 1u;
          continue;
        }

      if (g_strstr_len(lower, -1, "surface cartesian vectors (angstroms)") != NULL)
        {
          gdouble values[3];
          guint scan_index;
          guint parsed_vectors;
          guint count_a;
          guint count_b;
          gdouble c_vec[3];

          count_a = 0u;
          count_b = 0u;
          scan_index = i + 1u;
          parsed_vectors = 0u;
          while (lines[scan_index] != NULL && parsed_vectors < 2u)
            {
              gchar *vector_line;
              guint count;

              vector_line = g_strstrip(lines[scan_index]);
              if (vector_line[0] == '\0')
                {
                  scan_index++;
                  continue;
                }

              count = gdis_collect_doubles_from_line(vector_line, values, 3u);
              if (count < 2u)
                break;

              if (parsed_vectors == 0u)
                {
                  context_surface_a[0] = values[0];
                  context_surface_a[1] = values[1];
                  context_surface_a[2] = (count >= 3u) ? values[2] : 0.0;
                  count_a = count;
                }
              else
                {
                  context_surface_b[0] = values[0];
                  context_surface_b[1] = values[1];
                  context_surface_b[2] = (count >= 3u) ? values[2] : 0.0;
                  count_b = count;
                }

              parsed_vectors++;
              scan_index++;
            }

          if (count_a >= 2u && count_b >= 2u)
            {
              context_have_surface_vectors = TRUE;
              c_vec[0] = 0.0;
              c_vec[1] = 0.0;
              c_vec[2] = 1.0;
              if (gdis_surface_build_cell_from_vectors(context_surface_a,
                                                       context_surface_b,
                                                       c_vec,
                                                       context_lengths,
                                                       context_angles))
                {
                  context_periodic = TRUE;
                  context_periodicity = 2u;
                }
            }

          if (scan_index > i + 1u)
            i = scan_index - 1u;
          continue;
        }

      mode = GDIS_GULP_BLOCK_NONE;
      if (g_strstr_len(lower, -1, "fractional coordinates of asymmetric unit") != NULL ||
          g_strstr_len(lower, -1, "final asymmetric unit coordinates") != NULL)
        mode = GDIS_GULP_BLOCK_FRACTIONAL;
      else if (g_strstr_len(lower, -1, "mixed fractional/cartesian coordinates of surface") != NULL)
        mode = GDIS_GULP_BLOCK_SURFACE_MIXED;
      else if (g_strstr_len(lower, -1, "cartesian coordinates of cluster") != NULL ||
               g_strstr_len(lower, -1, "final cartesian coordinates of atoms") != NULL)
        mode = GDIS_GULP_BLOCK_CART;

      if (mode != GDIS_GULP_BLOCK_NONE)
        {
          GPtrArray *block_atoms;
          gboolean block_periodic;
          guint block_periodicity;
          gdouble block_lengths[3];
          gdouble block_angles[3];
          gboolean block_have_surface_vectors;
          gdouble block_surface_a[3];
          gdouble block_surface_b[3];
          gboolean block_have_cell_matrix;
          gdouble block_matrix[9];
          gdouble block_inverse[9];
          gint block_region;
          gboolean block_had_rows;
          GdisModel block_reference = {0};

          block_atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
          block_periodic = context_periodic;
          block_periodicity = context_periodicity;
          memcpy(block_lengths, context_lengths, sizeof(block_lengths));
          memcpy(block_angles, context_angles, sizeof(block_angles));
          block_have_surface_vectors = context_have_surface_vectors;
          gdis_vec3_copy(block_surface_a, context_surface_a);
          gdis_vec3_copy(block_surface_b, context_surface_b);
          if (mode == GDIS_GULP_BLOCK_SURFACE_MIXED && block_periodicity == 0u)
            {
              block_periodic = TRUE;
              block_periodicity = 2u;
            }

          block_reference.periodic = block_periodic;
          block_reference.periodicity = block_periodicity;
          memcpy(block_reference.cell_lengths, block_lengths, sizeof(block_lengths));
          memcpy(block_reference.cell_angles, block_angles, sizeof(block_angles));
          block_have_cell_matrix = gdis_model_build_cell_matrix(&block_reference,
                                                                block_matrix,
                                                                block_inverse);
          block_region = -1;
          block_had_rows = FALSE;

          for (i = i + 1; lines[i] != NULL; i++)
            {
              gchar *block_line;
              gchar *block_trimmed;
              g_autofree gchar *block_lower = NULL;
              g_auto(GStrv) tokens = NULL;
              gint token_count = 0;
              guint serial = 0u;
              gdouble values[3];
              gdouble cart[3];
              g_autofree gchar *element = NULL;
              g_autofree gchar *type = NULL;

              block_line = lines[i];
              g_strchomp(block_line);
              block_trimmed = g_strstrip(block_line);
              if (block_trimmed[0] == '\0')
                {
                  if (block_had_rows)
                    break;
                  continue;
                }

              block_lower = g_ascii_strdown(block_trimmed, -1);
              if (g_strstr_len(block_lower, -1, "coordinates of") != NULL && block_had_rows)
                {
                  i--;
                  break;
                }

              if (g_str_has_prefix(block_lower, "region"))
                {
                  gdouble region_value[1];
                  if (gdis_collect_doubles_from_line(block_trimmed, region_value, 1u) >= 1u)
                    block_region = (gint) llround(region_value[0]);
                  continue;
                }

              if (g_str_has_prefix(block_trimmed, "-") || g_str_has_prefix(block_trimmed, "*"))
                continue;

              tokens = gdis_split_simple(block_trimmed, &token_count);
              if (token_count < 4)
                {
                  if (block_had_rows)
                    {
                      i--;
                      break;
                    }
                  continue;
                }

              if (!gdis_try_parse_uint(tokens[0], &serial))
                {
                  if (block_had_rows)
                    {
                      i--;
                      break;
                    }
                  continue;
                }

              if (!gdis_collect_three_doubles_from_tokens(tokens, token_count, 3, values) &&
                  !gdis_collect_three_doubles_from_tokens(tokens, token_count, 2, values))
                continue;

              type = (token_count > 2) ? g_ascii_strdown(tokens[2], -1) : g_strdup("c");
              if (g_str_has_prefix(type, "shel") || g_str_equal(type, "s"))
                {
                  block_had_rows = TRUE;
                  continue;
                }

              if (mode == GDIS_GULP_BLOCK_CART)
                {
                  cart[0] = values[0];
                  cart[1] = values[1];
                  cart[2] = values[2];
                }
              else if (mode == GDIS_GULP_BLOCK_FRACTIONAL)
                {
                  if (!block_have_cell_matrix)
                    continue;
                  gdis_frac_to_cart(block_matrix, values, cart);
                }
              else
                {
                  if (block_have_surface_vectors)
                    {
                      cart[0] = block_surface_a[0] * values[0] + block_surface_b[0] * values[1];
                      cart[1] = block_surface_a[1] * values[0] + block_surface_b[1] * values[1];
                      cart[2] = block_surface_a[2] * values[0] +
                                block_surface_b[2] * values[1] +
                                values[2];
                    }
                  else if (block_have_cell_matrix)
                    {
                      cart[0] = block_matrix[0] * values[0] + block_matrix[1] * values[1];
                      cart[1] = block_matrix[3] * values[0] + block_matrix[4] * values[1];
                      cart[2] = values[2];
                    }
                  else
                    {
                      cart[0] = values[0];
                      cart[1] = values[1];
                      cart[2] = values[2];
                    }
                }

              element = gdis_normalize_element_symbol(tokens[1]);
              g_ptr_array_add(block_atoms,
                              gdis_atom_new(tokens[1],
                                            element,
                                            element,
                                            cart[0],
                                            cart[1],
                                            cart[2],
                                            1.0,
                                            (mode == GDIS_GULP_BLOCK_SURFACE_MIXED)
                                              ? block_region
                                              : -1,
                                            serial));
              block_had_rows = TRUE;
            }

          if (block_atoms->len > 0u)
            {
              if (best_atoms)
                g_ptr_array_free(best_atoms, TRUE);
              best_atoms = block_atoms;
              best_periodic = block_periodic;
              best_periodicity = block_periodicity;
              memcpy(best_lengths, block_lengths, sizeof(best_lengths));
              memcpy(best_angles, block_angles, sizeof(best_angles));
            }
          else
            {
              g_ptr_array_free(block_atoms, TRUE);
            }

          continue;
        }
    }

  g_strfreev(lines);

  if (!best_atoms)
    {
      GPtrArray *fallback_atoms;
      GdisModel fallback_reference = {0};
      gboolean fallback_have_matrix;
      gdouble fallback_matrix[9];
      gdouble fallback_inverse[9];
      gchar **fallback_lines;

      fallback_atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
      fallback_reference.periodic = context_periodic;
      fallback_reference.periodicity = context_periodicity;
      memcpy(fallback_reference.cell_lengths, context_lengths, sizeof(context_lengths));
      memcpy(fallback_reference.cell_angles, context_angles, sizeof(context_angles));
      fallback_have_matrix = gdis_model_build_cell_matrix(&fallback_reference,
                                                          fallback_matrix,
                                                          fallback_inverse);

      fallback_lines = g_strsplit(contents, "\n", -1);
      for (guint i = 0; fallback_lines[i] != NULL; i++)
        {
          g_auto(GStrv) tokens = NULL;
          gint token_count = 0;
          guint serial;
          gdouble values[3];
          gdouble cart[3];
          g_autofree gchar *type = NULL;
          g_autofree gchar *element = NULL;
          gboolean looks_fractional;

          tokens = gdis_split_simple(fallback_lines[i], &token_count);
          if (token_count < 5)
            continue;
          if (!gdis_try_parse_uint(tokens[0], &serial))
            continue;

          type = (token_count > 2) ? g_ascii_strdown(tokens[2], -1) : g_strdup("c");
          if (g_str_has_prefix(type, "shel") || g_str_equal(type, "s"))
            continue;

          if (!gdis_collect_three_doubles_from_tokens(tokens, token_count, 3, values) &&
              !gdis_collect_three_doubles_from_tokens(tokens, token_count, 2, values))
            continue;

          looks_fractional = fabs(values[0]) <= 1.5 && fabs(values[1]) <= 1.5 && fabs(values[2]) <= 1.5;
          if (fallback_have_matrix && looks_fractional)
            gdis_frac_to_cart(fallback_matrix, values, cart);
          else
            {
              cart[0] = values[0];
              cart[1] = values[1];
              cart[2] = values[2];
            }

          element = gdis_normalize_element_symbol(tokens[1]);
          g_ptr_array_add(fallback_atoms,
                          gdis_atom_new(tokens[1],
                                        element,
                                        element,
                                        cart[0],
                                        cart[1],
                                        cart[2],
                                        1.0,
                                        -1,
                                        serial));
        }
      g_strfreev(fallback_lines);

      if (fallback_atoms->len > 0u)
        {
          best_atoms = fallback_atoms;
          best_periodic = context_periodic;
          best_periodicity = context_periodicity;
          memcpy(best_lengths, context_lengths, sizeof(best_lengths));
          memcpy(best_angles, context_angles, sizeof(best_angles));
        }
      else
        {
          g_ptr_array_free(fallback_atoms, TRUE);
        }
    }

  if (!best_atoms || best_atoms->len == 0u)
    {
      if (best_atoms)
        g_ptr_array_free(best_atoms, TRUE);
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No atom coordinate table could be parsed from GULP output '%s'.",
                  model->path);
      return FALSE;
    }

  if (model->atoms)
    g_ptr_array_free(model->atoms, TRUE);
  model->atoms = best_atoms;
  model->periodic = best_periodic;
  model->periodicity = best_periodicity;
  memcpy(model->cell_lengths, best_lengths, sizeof(best_lengths));
  memcpy(model->cell_angles, best_angles, sizeof(best_angles));

  return TRUE;
}

static gboolean
gdis_model_load_xtl(GdisModel *model, const gchar *contents, GError **error)
{
  gchar **lines;
  gboolean in_atoms;
  gboolean have_cell_matrix;
  gdouble cell_matrix[9];
  gdouble cell_inverse[9];

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(contents != NULL, FALSE);

  lines = g_strsplit(contents, "\n", -1);
  in_atoms = FALSE;
  have_cell_matrix = FALSE;

  for (guint i = 0; lines[i] != NULL; i++)
    {
      gchar *line;
      gchar *trimmed;
      g_auto(GStrv) tokens = NULL;
      gint token_count = 0;
      g_autofree gchar *first_lower = NULL;

      line = lines[i];
      g_strchomp(line);
      trimmed = g_strstrip(line);
      if (trimmed[0] == '\0')
        continue;

      tokens = gdis_split_simple(trimmed, &token_count);
      if (token_count == 0)
        continue;

      first_lower = g_ascii_strdown(tokens[0], -1);

      if (g_str_equal(first_lower, "title"))
        {
          const gchar *title_text;

          title_text = trimmed + strlen(tokens[0]);
          while (*title_text && g_ascii_isspace(*title_text))
            title_text++;
          if (*title_text)
            {
              g_clear_pointer(&model->title, g_free);
              model->title = g_strdup(title_text);
            }
          continue;
        }

      if (g_str_equal(first_lower, "cell"))
        {
          gdouble values[6] = {0.0, 0.0, 0.0, 90.0, 90.0, 90.0};
          guint found;
          guint line_index;

          found = gdis_collect_doubles_from_line(trimmed, values, 6u);
          line_index = i;
          while (found < 6u && lines[line_index + 1] != NULL)
            {
              guint added;

              line_index++;
              added = gdis_collect_doubles_from_line(lines[line_index],
                                                     values + found,
                                                     6u - found);
              if (added == 0u)
                {
                  line_index--;
                  break;
                }
              found += added;
            }
          i = line_index;

          if (found >= 6u)
            {
              memcpy(model->cell_lengths, values, 3u * sizeof(gdouble));
              memcpy(model->cell_angles, values + 3, 3u * sizeof(gdouble));
              model->periodic = TRUE;
              model->periodicity = 3u;
              have_cell_matrix = gdis_model_build_cell_matrix(model,
                                                              cell_matrix,
                                                              cell_inverse);
            }
          continue;
        }

      if (g_str_equal(first_lower, "symmetry"))
        {
          for (gint j = 1; j + 1 < token_count; j++)
            {
              g_autofree gchar *keyword = g_ascii_strdown(tokens[j], -1);
              if (g_str_equal(keyword, "label"))
                {
                  GString *space_group;

                  space_group = g_string_new(tokens[j + 1]);
                  for (gint k = j + 2; k < token_count; k++)
                    {
                      g_string_append_c(space_group, ' ');
                      g_string_append(space_group, tokens[k]);
                    }
                  g_clear_pointer(&model->space_group, g_free);
                  model->space_group = g_string_free(space_group, FALSE);
                  break;
                }
            }
          continue;
        }

      if (g_str_equal(first_lower, "atoms"))
        {
          in_atoms = TRUE;
          continue;
        }

      if (in_atoms)
        {
          gdouble frac_or_cart[3];
          gdouble cart[3];
          gdouble occupancy;
          g_autofree gchar *element = NULL;

          if (g_str_equal(first_lower, "eof"))
            break;

          if (token_count < 4)
            continue;

          if (!gdis_collect_three_doubles_from_tokens(tokens, token_count, 1, frac_or_cart))
            continue;

          if (have_cell_matrix)
            gdis_frac_to_cart(cell_matrix, frac_or_cart, cart);
          else
            {
              cart[0] = frac_or_cart[0];
              cart[1] = frac_or_cart[1];
              cart[2] = frac_or_cart[2];
            }

          occupancy = 1.0;
          if (token_count > 6)
            {
              gdouble parsed_occ;
              if (gdis_try_parse_double_relaxed(tokens[6], &parsed_occ))
                occupancy = parsed_occ;
            }

          element = gdis_normalize_element_symbol(tokens[0]);
          g_ptr_array_add(model->atoms,
                          gdis_atom_new(tokens[0],
                                        element,
                                        element,
                                        cart[0],
                                        cart[1],
                                        cart[2],
                                        occupancy,
                                        -1,
                                        model->atoms->len + 1u));
        }
    }

  g_strfreev(lines);

  if (model->atoms->len == 0u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No atoms were parsed from XTL '%s'.",
                  model->path);
      return FALSE;
    }

  return TRUE;
}

static gboolean
gdis_qe_unit_scale_from_header(const gchar *header,
                               gdouble alat_bohr,
                               gdouble *scale_out,
                               gboolean *fractional_out)
{
  g_autofree gchar *lower = NULL;
  gdouble scale = 1.0;
  gboolean fractional = FALSE;
  gdouble parsed[1];

  g_return_val_if_fail(scale_out != NULL, FALSE);
  g_return_val_if_fail(fractional_out != NULL, FALSE);

  if (!header)
    {
      *scale_out = 1.0;
      *fractional_out = FALSE;
      return TRUE;
    }

  lower = g_ascii_strdown(header, -1);
  fractional = (g_strstr_len(lower, -1, "crystal") != NULL);

  if (fractional)
    scale = 1.0;
  else if (g_strstr_len(lower, -1, "angstrom") != NULL)
    scale = 1.0;
  else if (g_strstr_len(lower, -1, "bohr") != NULL)
    scale = 0.529177210903;
  else if (g_strstr_len(lower, -1, "alat") != NULL)
    {
      guint count;

      count = gdis_collect_doubles_from_line(header, parsed, 1u);
      if (count == 0u)
        parsed[0] = alat_bohr;
      if (parsed[0] <= 0.0)
        return FALSE;
      scale = parsed[0] * 0.529177210903;
    }

  *scale_out = scale;
  *fractional_out = fractional;
  return TRUE;
}

static gboolean
gdis_model_load_qe_like(GdisModel *model, const gchar *contents, GError **error)
{
  gchar **lines;
  gdouble alat_bohr;
  gint qe_ibrav;
  gdouble qe_celldm1;
  gdouble qe_celldm3;
  gboolean qe_have_celldm1;
  gboolean qe_have_celldm3;
  gboolean context_periodic;
  guint context_periodicity;
  gdouble context_lengths[3];
  gdouble context_angles[3];
  GPtrArray *best_atoms;
  gboolean best_periodic;
  guint best_periodicity;
  gdouble best_lengths[3];
  gdouble best_angles[3];

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(contents != NULL, FALSE);

  lines = g_strsplit(contents, "\n", -1);
  alat_bohr = 0.0;
  qe_ibrav = 0;
  qe_celldm1 = 0.0;
  qe_celldm3 = 0.0;
  qe_have_celldm1 = FALSE;
  qe_have_celldm3 = FALSE;
  context_periodic = FALSE;
  context_periodicity = 0u;
  context_lengths[0] = context_lengths[1] = context_lengths[2] = 0.0;
  context_angles[0] = context_angles[1] = context_angles[2] = 90.0;
  best_atoms = NULL;
  best_periodic = FALSE;
  best_periodicity = 0u;
  best_lengths[0] = best_lengths[1] = best_lengths[2] = 0.0;
  best_angles[0] = best_angles[1] = best_angles[2] = 90.0;

  for (guint i = 0; lines[i] != NULL; i++)
    {
      gchar *line;
      gchar *trimmed;
      g_autofree gchar *lower = NULL;

      line = lines[i];
      g_strchomp(line);
      trimmed = g_strstrip(line);
      if (trimmed[0] == '\0')
        continue;

      lower = g_ascii_strdown(trimmed, -1);

      if (g_strstr_len(lower, -1, "lattice parameter (alat)") != NULL)
        {
          gdouble values[1];
          if (gdis_collect_doubles_from_line(trimmed, values, 1u) >= 1u)
            alat_bohr = values[0];
          continue;
        }

      if (g_strstr_len(lower, -1, "ibrav") != NULL && strchr(trimmed, '=') != NULL)
        {
          const gchar *eq;
          gdouble parsed;

          eq = strchr(trimmed, '=');
          if (eq && gdis_try_parse_double_relaxed(eq + 1, &parsed))
            qe_ibrav = (gint) llround(parsed);
        }
      if (g_strstr_len(lower, -1, "celldm(1)") != NULL && strchr(trimmed, '=') != NULL)
        {
          const gchar *eq;
          gdouble parsed;

          eq = strchr(trimmed, '=');
          if (eq && gdis_try_parse_double_relaxed(eq + 1, &parsed))
            {
              qe_celldm1 = parsed;
              qe_have_celldm1 = TRUE;
            }
        }
      if (g_strstr_len(lower, -1, "celldm(3)") != NULL && strchr(trimmed, '=') != NULL)
        {
          const gchar *eq;
          gdouble parsed;

          eq = strchr(trimmed, '=');
          if (eq && gdis_try_parse_double_relaxed(eq + 1, &parsed))
            {
              qe_celldm3 = parsed;
              qe_have_celldm3 = TRUE;
            }
        }

      if (!context_periodic && qe_ibrav == 4 && qe_have_celldm1 && qe_have_celldm3)
        {
          const gdouble a_ang = qe_celldm1 * 0.529177210903;
          const gdouble c_ang = qe_celldm3 * a_ang;
          if (a_ang > 0.0 && c_ang > 0.0)
            {
              context_periodic = TRUE;
              context_periodicity = 3u;
              context_lengths[0] = a_ang;
              context_lengths[1] = a_ang;
              context_lengths[2] = c_ang;
              context_angles[0] = 90.0;
              context_angles[1] = 90.0;
              context_angles[2] = 120.0;
            }
        }

      if (!model->title || model->title[0] == '\0')
        {
          if (g_str_has_prefix(lower, "prefix"))
            {
              const gchar *eq;

              eq = strchr(trimmed, '=');
              if (eq)
                {
                  g_autofree gchar *prefix = g_strdup(eq + 1);
                  g_strstrip(prefix);
                  if (prefix[0] == '\'' || prefix[0] == '"')
                    memmove(prefix, prefix + 1, strlen(prefix));
                  g_strchomp(prefix);
                  g_strstrip(prefix);
                  while (prefix[0] != '\0')
                    {
                      gsize len;

                      len = strlen(prefix);
                      if (prefix[len - 1] == ',' ||
                          prefix[len - 1] == '\'' ||
                          prefix[len - 1] == '"')
                        {
                          prefix[len - 1] = '\0';
                          g_strstrip(prefix);
                          continue;
                        }
                      break;
                    }
                  if (prefix[0] != '\0')
                    {
                      g_clear_pointer(&model->title, g_free);
                      model->title = g_strdup(prefix);
                    }
                }
            }
        }

      if (g_str_has_prefix(lower, "cell_parameters"))
        {
          gdouble unit_scale;
          gboolean fractional_units;
          gdouble a_vec[3];
          gdouble b_vec[3];
          gdouble c_vec[3];
          gdouble values[3];
          guint count_a;
          guint count_b;
          guint count_c;

          if (!gdis_qe_unit_scale_from_header(trimmed,
                                              alat_bohr,
                                              &unit_scale,
                                              &fractional_units))
            continue;
          (void) fractional_units;

          if (lines[i + 1] == NULL || lines[i + 2] == NULL || lines[i + 3] == NULL)
            continue;

          count_a = gdis_collect_doubles_from_line(lines[i + 1], values, 3u);
          if (count_a >= 3u)
            {
              a_vec[0] = values[0] * unit_scale;
              a_vec[1] = values[1] * unit_scale;
              a_vec[2] = values[2] * unit_scale;
            }

          count_b = gdis_collect_doubles_from_line(lines[i + 2], values, 3u);
          if (count_b >= 3u)
            {
              b_vec[0] = values[0] * unit_scale;
              b_vec[1] = values[1] * unit_scale;
              b_vec[2] = values[2] * unit_scale;
            }

          count_c = gdis_collect_doubles_from_line(lines[i + 3], values, 3u);
          if (count_c >= 3u)
            {
              c_vec[0] = values[0] * unit_scale;
              c_vec[1] = values[1] * unit_scale;
              c_vec[2] = values[2] * unit_scale;
            }

          if (count_a >= 3u &&
              count_b >= 3u &&
              count_c >= 3u &&
              gdis_surface_build_cell_from_vectors(a_vec,
                                                   b_vec,
                                                   c_vec,
                                                   context_lengths,
                                                   context_angles))
            {
              context_periodic = TRUE;
              context_periodicity = 3u;
            }

          i += 3u;
          continue;
        }

      if (g_str_has_prefix(lower, "atomic_positions"))
        {
          GPtrArray *block_atoms;
          gboolean block_periodic;
          guint block_periodicity;
          gdouble block_lengths[3];
          gdouble block_angles[3];
          gdouble unit_scale;
          gboolean fractional_units;
          GdisModel block_reference = {0};
          gboolean block_have_matrix;
          gdouble block_matrix[9];
          gdouble block_inverse[9];

          if (!gdis_qe_unit_scale_from_header(trimmed,
                                              alat_bohr,
                                              &unit_scale,
                                              &fractional_units))
            continue;

          block_atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
          block_periodic = context_periodic;
          block_periodicity = context_periodicity;
          memcpy(block_lengths, context_lengths, sizeof(block_lengths));
          memcpy(block_angles, context_angles, sizeof(block_angles));

          block_reference.periodic = block_periodic;
          block_reference.periodicity = block_periodicity;
          memcpy(block_reference.cell_lengths, block_lengths, sizeof(block_lengths));
          memcpy(block_reference.cell_angles, block_angles, sizeof(block_angles));
          block_have_matrix = gdis_model_build_cell_matrix(&block_reference,
                                                           block_matrix,
                                                           block_inverse);

          for (i = i + 1; lines[i] != NULL; i++)
            {
              gchar *atom_line;
              gchar *atom_trimmed;
              g_autofree gchar *atom_lower = NULL;
              g_auto(GStrv) tokens = NULL;
              gint token_count = 0;
              gdouble values[3];
              gdouble cart[3];
              gint species_index;
              gint coord_start;
              guint serial_probe;
              const gchar *species;
              g_autofree gchar *element = NULL;

              atom_line = lines[i];
              g_strchomp(atom_line);
              atom_trimmed = g_strstrip(atom_line);
              if (atom_trimmed[0] == '\0')
                {
                  if (block_atoms->len > 0u)
                    break;
                  continue;
                }

              atom_lower = g_ascii_strdown(atom_trimmed, -1);
              if (g_str_has_prefix(atom_lower, "k_points") ||
                  g_str_has_prefix(atom_lower, "cell_parameters") ||
                  g_str_has_prefix(atom_lower, "end") ||
                  g_str_has_prefix(atom_lower, "begin") ||
                  g_str_has_prefix(atom_lower, "&"))
                {
                  if (block_atoms->len > 0u)
                    {
                      i--;
                      break;
                    }
                  continue;
                }

              tokens = gdis_split_simple(atom_trimmed, &token_count);
              if (token_count < 4)
                {
                  if (block_atoms->len > 0u)
                    {
                      i--;
                      break;
                    }
                  continue;
                }

              species_index = 0;
              coord_start = 1;
              if (gdis_try_parse_uint(tokens[0], &serial_probe) && token_count >= 5)
                {
                  species_index = 1;
                  coord_start = 2;
                }

              if (!gdis_collect_three_doubles_from_tokens(tokens,
                                                          token_count,
                                                          coord_start,
                                                          values))
                {
                  if (block_atoms->len > 0u)
                    {
                      i--;
                      break;
                    }
                  continue;
                }

              if (fractional_units)
                {
                  if (!block_have_matrix)
                    continue;
                  gdis_frac_to_cart(block_matrix, values, cart);
                }
              else
                {
                  cart[0] = values[0] * unit_scale;
                  cart[1] = values[1] * unit_scale;
                  cart[2] = values[2] * unit_scale;
                }

              species = (species_index < token_count) ? tokens[species_index] : "X";
              element = gdis_normalize_element_symbol(species);
              g_ptr_array_add(block_atoms,
                              gdis_atom_new(species,
                                            element,
                                            element,
                                            cart[0],
                                            cart[1],
                                            cart[2],
                                            1.0,
                                            -1,
                                            block_atoms->len + 1u));
            }

          if (block_atoms->len > 0u)
            {
              if (best_atoms)
                g_ptr_array_free(best_atoms, TRUE);
              best_atoms = block_atoms;
              best_periodic = block_periodic;
              best_periodicity = block_periodicity;
              memcpy(best_lengths, block_lengths, sizeof(best_lengths));
              memcpy(best_angles, block_angles, sizeof(best_angles));
            }
          else
            {
              g_ptr_array_free(block_atoms, TRUE);
            }

          continue;
        }
    }

  g_strfreev(lines);

  if (!best_atoms || best_atoms->len == 0u)
    {
      if (best_atoms)
        g_ptr_array_free(best_atoms, TRUE);
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No QE atomic positions were parsed from '%s'.",
                  model->path);
      return FALSE;
    }

  if (model->atoms)
    g_ptr_array_free(model->atoms, TRUE);
  model->atoms = best_atoms;
  model->periodic = best_periodic;
  model->periodicity = best_periodicity;
  memcpy(model->cell_lengths, best_lengths, sizeof(best_lengths));
  memcpy(model->cell_angles, best_angles, sizeof(best_angles));

  return TRUE;
}

static gboolean
gdis_model_load_qe_input(GdisModel *model, const gchar *contents, GError **error)
{
  return gdis_model_load_qe_like(model, contents, error);
}

static gboolean
gdis_model_load_qe_output(GdisModel *model, const gchar *contents, GError **error)
{
  return gdis_model_load_qe_like(model, contents, error);
}

static gboolean
gdis_model_load_gamess_input(GdisModel *model, const gchar *contents, GError **error)
{
  gchar **lines;
  gboolean in_data_block;
  guint data_header_lines;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(contents != NULL, FALSE);

  lines = g_strsplit(contents, "\n", -1);
  in_data_block = FALSE;
  data_header_lines = 0u;

  for (guint i = 0; lines[i] != NULL; i++)
    {
      gchar *line;
      gchar *trimmed;
      g_autofree gchar *lower = NULL;
      g_auto(GStrv) tokens = NULL;
      gint token_count = 0;
      gdouble coords[3];
      g_autofree gchar *element = NULL;

      line = lines[i];
      g_strchomp(line);
      trimmed = g_strstrip(line);
      if (trimmed[0] == '\0')
        continue;

      lower = g_ascii_strdown(trimmed, -1);
      if (!in_data_block)
        {
          if (g_str_has_prefix(lower, "$data"))
            {
              in_data_block = TRUE;
              data_header_lines = 0u;
            }
          continue;
        }

      if (g_str_has_prefix(lower, "$end"))
        break;

      if (data_header_lines == 0u)
        {
          g_clear_pointer(&model->title, g_free);
          model->title = g_strdup(trimmed);
          data_header_lines++;
          continue;
        }
      if (data_header_lines == 1u)
        {
          data_header_lines++;
          continue;
        }

      tokens = gdis_split_simple(trimmed, &token_count);
      if (token_count < 5)
        continue;
      if (!gdis_collect_three_doubles_from_tokens(tokens, token_count, 2, coords))
        continue;

      element = gdis_normalize_element_symbol(tokens[0]);
      g_ptr_array_add(model->atoms,
                      gdis_atom_new(tokens[0],
                                    element,
                                    element,
                                    coords[0],
                                    coords[1],
                                    coords[2],
                                    1.0,
                                    -1,
                                    model->atoms->len + 1u));
    }

  g_strfreev(lines);

  if (model->atoms->len == 0u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No atoms were parsed from GAMESS input '%s'.",
                  model->path);
      return FALSE;
    }

  return TRUE;
}

static gboolean
gdis_model_load_cif(GdisModel *model, const gchar *contents, GError **error)
{
  gchar **lines;
  GPtrArray *records;
  gboolean in_multiline;
  guint i;

  lines = g_strsplit(contents, "\n", -1);
  records = g_ptr_array_new_with_free_func(gdis_cif_atom_record_free);
  in_multiline = FALSE;

  for (i = 0; lines[i] != NULL; )
    {
      gchar *line;
      gchar *trimmed;

      line = lines[i];
      g_strchomp(line);
      trimmed = g_strstrip(line);

      if (trimmed[0] == ';')
        {
          in_multiline = !in_multiline;
          i++;
          continue;
        }

      if (in_multiline || trimmed[0] == '\0' || trimmed[0] == '#')
        {
          i++;
          continue;
        }

      if (g_ascii_strncasecmp(trimmed, "data_", 5) == 0)
        {
          if (!model->title || model->title[0] == '\0')
            {
              g_clear_pointer(&model->title, g_free);
              model->title = gdis_strdup_strip(trimmed + 5);
            }
          i++;
          continue;
        }

      if (g_ascii_strcasecmp(trimmed, "loop_") == 0)
        {
          GPtrArray *headers;
          gboolean atom_loop;

          headers = g_ptr_array_new_with_free_func(g_free);
          atom_loop = FALSE;
          i++;

          while (lines[i] != NULL)
            {
              gchar *header_line;
              gchar *header_trimmed;

              header_line = lines[i];
              g_strchomp(header_line);
              header_trimmed = g_strstrip(header_line);

              if (header_trimmed[0] != '_')
                break;

              g_ptr_array_add(headers, g_strdup(header_trimmed));
              if (g_str_has_prefix(header_trimmed, "_atom_site_"))
                atom_loop = TRUE;
              i++;
            }

          if (atom_loop && headers->len > 0)
            {
              while (lines[i] != NULL)
                {
                  GPtrArray *tokens;
                  gchar *row_line;
                  gchar *row_trimmed;
                  gint label_pos;
                  gint symbol_pos;
                  gint occ_pos;
                  gint frac_x_pos;
                  gint frac_y_pos;
                  gint frac_z_pos;
                  gint cart_x_pos;
                  gint cart_y_pos;
                  gint cart_z_pos;
                  guint h;

                  row_line = lines[i];
                  g_strchomp(row_line);
                  row_trimmed = g_strstrip(row_line);

                  if (row_trimmed[0] == ';')
                    {
                      in_multiline = !in_multiline;
                      i++;
                      continue;
                    }

                  if (in_multiline)
                    {
                      i++;
                      continue;
                    }

                  if (gdis_cif_is_control_line(row_trimmed))
                    break;

                  tokens = gdis_tokenize_cif(row_trimmed);
                  if (tokens->len >= headers->len)
                    {
                      GdisCifAtomRecord *record;
                      gchar *label;
                      gchar *symbol;
                      gdouble x;
                      gdouble y;
                      gdouble z;
                      gboolean fractional;

                      label_pos = -1;
                      symbol_pos = -1;
                      occ_pos = -1;
                      frac_x_pos = -1;
                      frac_y_pos = -1;
                      frac_z_pos = -1;
                      cart_x_pos = -1;
                      cart_y_pos = -1;
                      cart_z_pos = -1;

                      for (h = 0; h < headers->len; h++)
                        {
                          const gchar *header;
                          gchar *lower_header;

                          header = g_ptr_array_index(headers, h);
                          lower_header = g_ascii_strdown(header, -1);

                          if (g_str_equal(lower_header, "_atom_site_label"))
                            label_pos = (gint) h;
                          else if (g_str_equal(lower_header, "_atom_site_type_symbol"))
                            symbol_pos = (gint) h;
                          else if (g_str_equal(lower_header, "_atom_site_occupancy"))
                            occ_pos = (gint) h;
                          else if (g_str_equal(lower_header, "_atom_site_fract_x"))
                            frac_x_pos = (gint) h;
                          else if (g_str_equal(lower_header, "_atom_site_fract_y"))
                            frac_y_pos = (gint) h;
                          else if (g_str_equal(lower_header, "_atom_site_fract_z"))
                            frac_z_pos = (gint) h;
                          else if (g_str_equal(lower_header, "_atom_site_cartn_x") ||
                                   g_str_equal(lower_header, "_atom_site_cart_x"))
                            cart_x_pos = (gint) h;
                          else if (g_str_equal(lower_header, "_atom_site_cartn_y") ||
                                   g_str_equal(lower_header, "_atom_site_cart_y"))
                            cart_y_pos = (gint) h;
                          else if (g_str_equal(lower_header, "_atom_site_cartn_z") ||
                                   g_str_equal(lower_header, "_atom_site_cart_z"))
                            cart_z_pos = (gint) h;

                          g_free(lower_header);
                        }

                      fractional = (frac_x_pos >= 0 && frac_y_pos >= 0 && frac_z_pos >= 0);
                      if (fractional || (cart_x_pos >= 0 && cart_y_pos >= 0 && cart_z_pos >= 0))
                        {
                          record = g_new0(GdisCifAtomRecord, 1);

                          label = label_pos >= 0 ? g_ptr_array_index(tokens, label_pos) : NULL;
                          symbol = symbol_pos >= 0 ? g_ptr_array_index(tokens, symbol_pos) : NULL;

                          record->label = g_strdup(label ? label : (symbol ? symbol : "X"));
                          record->element = symbol ?
                                            gdis_normalize_element_symbol(symbol) :
                                            gdis_normalize_element_symbol(record->label);
                          record->occupancy = occ_pos >= 0 ?
                                              gdis_parse_cif_number(g_ptr_array_index(tokens, occ_pos)) :
                                              1.0;
                          record->fractional = fractional;

                          if (fractional)
                            {
                              x = gdis_parse_cif_number(g_ptr_array_index(tokens, frac_x_pos));
                              y = gdis_parse_cif_number(g_ptr_array_index(tokens, frac_y_pos));
                              z = gdis_parse_cif_number(g_ptr_array_index(tokens, frac_z_pos));
                            }
                          else
                            {
                              x = gdis_parse_cif_number(g_ptr_array_index(tokens, cart_x_pos));
                              y = gdis_parse_cif_number(g_ptr_array_index(tokens, cart_y_pos));
                              z = gdis_parse_cif_number(g_ptr_array_index(tokens, cart_z_pos));
                            }

                          record->coords[0] = x;
                          record->coords[1] = y;
                          record->coords[2] = z;
                          g_ptr_array_add(records, record);
                        }
                    }

                  g_ptr_array_free(tokens, TRUE);
                  i++;
                }
            }
          else
            {
              while (lines[i] != NULL)
                {
                  gchar *row_line;
                  gchar *row_trimmed;

                  row_line = lines[i];
                  g_strchomp(row_line);
                  row_trimmed = g_strstrip(row_line);

                  if (gdis_cif_is_control_line(row_trimmed))
                    break;
                  i++;
                }
            }

          g_ptr_array_free(headers, TRUE);
          continue;
        }

      if (trimmed[0] == '_')
        {
          gchar *tag;
          gchar *value;
          gchar *lower_tag;

          tag = NULL;
          value = NULL;
          if (gdis_cif_parse_value_line(trimmed, &tag, &value))
            {
              lower_tag = g_ascii_strdown(tag, -1);

              if (g_str_equal(lower_tag, "_cell_length_a"))
                model->cell_lengths[0] = gdis_parse_cif_number(value);
              else if (g_str_equal(lower_tag, "_cell_length_b"))
                model->cell_lengths[1] = gdis_parse_cif_number(value);
              else if (g_str_equal(lower_tag, "_cell_length_c"))
                model->cell_lengths[2] = gdis_parse_cif_number(value);
              else if (g_str_equal(lower_tag, "_cell_angle_alpha"))
                model->cell_angles[0] = gdis_parse_cif_number(value);
              else if (g_str_equal(lower_tag, "_cell_angle_beta"))
                model->cell_angles[1] = gdis_parse_cif_number(value);
              else if (g_str_equal(lower_tag, "_cell_angle_gamma"))
                model->cell_angles[2] = gdis_parse_cif_number(value);
              else if (g_str_equal(lower_tag, "_symmetry_space_group_name_h-m") ||
                       g_str_equal(lower_tag, "_space_group_name_h-m_alt") ||
                       g_str_equal(lower_tag, "_space_group_name_h-m"))
                {
                  g_clear_pointer(&model->space_group, g_free);
                  model->space_group = gdis_strdup_strip(value);
                }
              else if (g_str_equal(lower_tag, "_symmetry_int_tables_number") ||
                       g_str_equal(lower_tag, "_space_group_it_number"))
                {
                  if (!model->space_group || model->space_group[0] == '\0')
                    {
                      g_clear_pointer(&model->space_group, g_free);
                      model->space_group = g_strdup(value);
                    }
                }

              g_free(lower_tag);
            }

          g_free(tag);
          g_free(value);
        }

      i++;
    }

  if (records->len > 0)
    {
      guint idx;
      gdouble matrix[9];
      gdouble inverse[9];
      gboolean has_cell_matrix;

      if (model->cell_lengths[0] > 0.0 &&
          model->cell_lengths[1] > 0.0 &&
          model->cell_lengths[2] > 0.0)
        {
          model->periodic = TRUE;
          model->periodicity = 3;
        }

      has_cell_matrix = gdis_model_build_cell_matrix(model, matrix, inverse);

      for (idx = 0; idx < records->len; idx++)
        {
          GdisCifAtomRecord *record;
          gdouble cart[3];

          record = g_ptr_array_index(records, idx);
          if (record->fractional)
            {
              if (!has_cell_matrix)
                {
                  g_set_error(error,
                              GDIS_MODEL_ERROR,
                              GDIS_MODEL_ERROR_PARSE,
                              "CIF atom loop in '%s' uses fractional coordinates without a valid cell.",
                              model->path);
                  g_ptr_array_free(records, TRUE);
                  g_strfreev(lines);
                  return FALSE;
                }

              gdis_frac_to_cart(matrix, record->coords, cart);
            }
          else
            {
              cart[0] = record->coords[0];
              cart[1] = record->coords[1];
              cart[2] = record->coords[2];
            }

          g_ptr_array_add(model->atoms,
                          gdis_atom_new(record->label,
                                        record->element,
                                        record->element,
                                        cart[0],
                                        cart[1],
                                        cart[2],
                                        record->occupancy,
                                        -1,
                                        idx + 1));
        }
    }

  g_ptr_array_free(records, TRUE);
  g_strfreev(lines);

  if (model->atoms->len == 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No supported atom loop was parsed from CIF '%s'.",
                  model->path);
      return FALSE;
    }

  return TRUE;
}

static gboolean
gdis_model_load_qbox_xml(GdisModel *model, const gchar *contents, GError **error)
{
  GdisQboxXmlParseState state = {0};
  GMarkupParser parser = {
    gdis_qbox_xml_start_element,
    gdis_qbox_xml_end_element,
    gdis_qbox_xml_text,
    NULL,
    NULL
  };
  GMarkupParseContext *context;
  gboolean ok;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(contents != NULL, FALSE);

  state.model = model;
  state.species_symbols = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
  state.text = g_string_new("");

  context = g_markup_parse_context_new(&parser, G_MARKUP_TREAT_CDATA_AS_TEXT, &state, NULL);
  ok = g_markup_parse_context_parse(context, contents, strlen(contents), error);
  if (ok)
    ok = g_markup_parse_context_end_parse(context, error);
  g_markup_parse_context_free(context);

  g_clear_pointer(&state.current_species_name, g_free);
  g_clear_pointer(&state.current_atom_name, g_free);
  g_clear_pointer(&state.current_atom_species, g_free);
  if (state.text)
    g_string_free(state.text, TRUE);
  g_clear_pointer(&state.species_symbols, g_hash_table_unref);

  if (!ok)
    return FALSE;

  if (!state.looks_like_qbox)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_UNSUPPORTED_FORMAT,
                  "'%s' is not a supported Qbox XML result file.",
                  model->path);
      return FALSE;
    }

  if (state.have_cell)
    {
      gdouble lengths[3];
      gdouble angles[3];

      if (!gdis_surface_build_cell_from_vectors(state.a_vec, state.b_vec, state.c_vec, lengths, angles))
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_PARSE,
                      "Qbox XML cell vectors in '%s' are invalid.",
                      model->path);
          return FALSE;
        }

      model->periodic = TRUE;
      model->periodicity = 3;
      memcpy(model->cell_lengths, lengths, sizeof(lengths));
      memcpy(model->cell_angles, angles, sizeof(angles));
    }

  if (model->atoms->len == 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No atoms were parsed from Qbox XML '%s'.",
                  model->path);
      return FALSE;
    }

  return TRUE;
}

static void
gdis_qbox_xml_start_element(GMarkupParseContext *context,
                            const gchar         *element_name,
                            const gchar        **attribute_names,
                            const gchar        **attribute_values,
                            gpointer             user_data,
                            GError             **parse_error)
{
  GdisQboxXmlParseState *state;
  const gchar *local_name;

  (void) context;

  state = user_data;
  local_name = gdis_xml_local_name(element_name);
  if (!local_name)
    return;

  if (g_str_equal(local_name, "sample") || g_str_equal(local_name, "simulation"))
    state->looks_like_qbox = TRUE;

  if (g_str_equal(local_name, "unit_cell"))
    {
      const gchar *attr_a = NULL;
      const gchar *attr_b = NULL;
      const gchar *attr_c = NULL;

      for (guint i = 0; attribute_names && attribute_names[i] != NULL; i++)
        {
          const gchar *attr_local_name;

          attr_local_name = gdis_xml_local_name(attribute_names[i]);
          if (g_str_equal(attr_local_name, "a"))
            attr_a = attribute_values[i];
          else if (g_str_equal(attr_local_name, "b"))
            attr_b = attribute_values[i];
          else if (g_str_equal(attr_local_name, "c"))
            attr_c = attribute_values[i];
        }

      if (!attr_a || !attr_b || !attr_c ||
          !gdis_parse_vector3_text(attr_a, state->a_vec) ||
          !gdis_parse_vector3_text(attr_b, state->b_vec) ||
          !gdis_parse_vector3_text(attr_c, state->c_vec))
        {
          g_set_error(parse_error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_PARSE,
                      "Qbox XML unit_cell in '%s' is missing or invalid.",
                      state->model->path);
          return;
        }

      for (guint axis = 0; axis < 3; axis++)
        {
          state->a_vec[axis] *= 0.529177210903;
          state->b_vec[axis] *= 0.529177210903;
          state->c_vec[axis] *= 0.529177210903;
        }
      state->have_cell = TRUE;
      return;
    }

  if (g_str_equal(local_name, "species"))
    {
      g_clear_pointer(&state->current_species_name, g_free);
      for (guint i = 0; attribute_names && attribute_names[i] != NULL; i++)
        {
          if (g_str_equal(gdis_xml_local_name(attribute_names[i]), "name"))
            {
              state->current_species_name = g_strdup(attribute_values[i]);
              break;
            }
        }
      return;
    }

  if (g_str_equal(local_name, "symbol"))
    {
      state->inside_species_symbol = TRUE;
      g_string_truncate(state->text, 0);
      return;
    }

  if (g_str_equal(local_name, "atom"))
    {
      g_clear_pointer(&state->current_atom_name, g_free);
      g_clear_pointer(&state->current_atom_species, g_free);
      state->have_current_atom_position = FALSE;

      for (guint i = 0; attribute_names && attribute_names[i] != NULL; i++)
        {
          const gchar *attr_local_name;

          attr_local_name = gdis_xml_local_name(attribute_names[i]);
          if (g_str_equal(attr_local_name, "name"))
            state->current_atom_name = g_strdup(attribute_values[i]);
          else if (g_str_equal(attr_local_name, "species"))
            state->current_atom_species = g_strdup(attribute_values[i]);
        }
      return;
    }

  if (g_str_equal(local_name, "position"))
    {
      state->inside_atom_position = TRUE;
      g_string_truncate(state->text, 0);
    }
}

static void
gdis_qbox_xml_end_element(GMarkupParseContext *context,
                          const gchar         *element_name,
                          gpointer             user_data,
                          GError             **parse_error)
{
  GdisQboxXmlParseState *state;
  const gchar *local_name;

  (void) context;

  state = user_data;
  local_name = gdis_xml_local_name(element_name);
  if (!local_name)
    return;

  if (g_str_equal(local_name, "symbol"))
    {
      g_autofree gchar *symbol = NULL;

      state->inside_species_symbol = FALSE;
      symbol = gdis_normalize_element_symbol(state->text->str);
      if (state->current_species_name && symbol && symbol[0])
        g_hash_table_replace(state->species_symbols,
                             g_strdup(state->current_species_name),
                             g_steal_pointer(&symbol));
      return;
    }

  if (g_str_equal(local_name, "position"))
    {
      state->inside_atom_position = FALSE;
      if (!gdis_parse_vector3_text(state->text->str, state->current_atom_position))
        {
          g_set_error(parse_error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_PARSE,
                      "Qbox XML atom position in '%s' is invalid.",
                      state->model->path);
          return;
        }
      for (guint axis = 0; axis < 3; axis++)
        state->current_atom_position[axis] *= 0.529177210903;
      state->have_current_atom_position = TRUE;
      return;
    }

  if (g_str_equal(local_name, "atom"))
    {
      g_autofree gchar *species_symbol = NULL;
      const gchar *mapped_symbol;
      const gchar *element_text;
      const gchar *label_text;

      if (!state->have_current_atom_position)
        {
          g_set_error(parse_error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_PARSE,
                      "Qbox XML atom in '%s' is missing a position.",
                      state->model->path);
          return;
        }

      mapped_symbol = state->current_atom_species ?
        g_hash_table_lookup(state->species_symbols, state->current_atom_species) :
        NULL;
      species_symbol = gdis_normalize_element_symbol(mapped_symbol ? mapped_symbol : state->current_atom_species);
      element_text = (species_symbol && species_symbol[0]) ? species_symbol : "X";
      label_text = (state->current_atom_name && state->current_atom_name[0]) ?
        state->current_atom_name : element_text;

      g_ptr_array_add(state->model->atoms,
                      gdis_atom_new(label_text,
                                    element_text,
                                    element_text,
                                    state->current_atom_position[0],
                                    state->current_atom_position[1],
                                    state->current_atom_position[2],
                                    1.0,
                                    -1,
                                    state->model->atoms->len + 1));

      g_clear_pointer(&state->current_atom_name, g_free);
      g_clear_pointer(&state->current_atom_species, g_free);
      state->have_current_atom_position = FALSE;
      return;
    }

  if (g_str_equal(local_name, "species"))
    g_clear_pointer(&state->current_species_name, g_free);
}

static void
gdis_qbox_xml_text(GMarkupParseContext *context,
                   const gchar         *text,
                   gsize                text_len,
                   gpointer             user_data,
                   GError             **parse_error)
{
  GdisQboxXmlParseState *state;

  (void) context;
  (void) parse_error;

  state = user_data;
  if (state->inside_species_symbol || state->inside_atom_position)
    g_string_append_len(state->text, text, text_len);
}

static const gchar *
gdis_xml_local_name(const gchar *name)
{
  const gchar *colon;

  if (!name)
    return NULL;

  colon = strrchr(name, ':');
  if (colon && colon[1] != '\0')
    return colon + 1;

  return name;
}

static gboolean
gdis_parse_vector3_text(const gchar *text, gdouble vector[3])
{
  g_auto(GStrv) tokens = NULL;
  guint token_count = 0u;

  g_return_val_if_fail(vector != NULL, FALSE);

  if (!text)
    return FALSE;

  tokens = g_strsplit_set(text, " \t\r\n,", -1);
  for (guint i = 0; tokens && tokens[i] != NULL; i++)
    {
      if (!tokens[i][0])
        continue;
      if (token_count >= 3u || !gdis_try_parse_double(tokens[i], &vector[token_count]))
        return FALSE;
      token_count++;
    }

  return token_count == 3u;
}

static gboolean
gdis_model_add_bond(GdisModel *model,
                    guint atom_index_a,
                    guint atom_index_b,
                    guint8 order,
                    gboolean inferred)
{
  GdisBond bond;

  if (!model || !model->bonds || atom_index_a == atom_index_b)
    return FALSE;

  bond.atom_index_a = MIN(atom_index_a, atom_index_b);
  bond.atom_index_b = MAX(atom_index_a, atom_index_b);
  bond.order = order > 0 ? order : 1;
  bond.inferred = inferred;
  g_array_append_val(model->bonds, bond);
  gdis_model_refresh_counts(model);
  return TRUE;
}

static gboolean
gdis_model_infer_bonds(GdisModel *model)
{
  GHashTable *seen;
  gdouble matrix[9];
  gdouble inverse[9];
  gboolean use_periodic_image;
  guint i;
  guint j;

  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->atoms || model->atoms->len < 2)
    {
      gdis_model_refresh_counts(model);
      return TRUE;
    }

  if (model->bonds->len > model->explicit_bond_count)
    g_array_set_size(model->bonds, model->explicit_bond_count);

  seen = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
  use_periodic_image = gdis_model_build_cell_matrix(model, matrix, inverse) && model->periodic;

  for (i = 0; i + 1 < model->atoms->len; i++)
    {
      GdisAtom *atom_i;
      gdouble radius_i;

      atom_i = g_ptr_array_index(model->atoms, i);
      radius_i = gdis_lookup_covalent_radius(atom_i->element);

      for (j = i + 1; j < model->atoms->len; j++)
        {
          GdisAtom *atom_j;
          gdouble radius_j;
          gdouble max_distance;
          gdouble distance2;
          gchar *key;

          atom_j = g_ptr_array_index(model->atoms, j);
          radius_j = gdis_lookup_covalent_radius(atom_j->element);
          max_distance = MIN(radius_i + radius_j + 0.45, 3.2);

          distance2 = gdis_model_distance2(model,
                                           atom_i->position,
                                           atom_j->position,
                                           matrix,
                                           inverse,
                                           use_periodic_image);
          if (distance2 < 0.16 || distance2 > max_distance * max_distance)
            continue;

          key = g_strdup_printf("%u:%u", i, j);
          if (g_hash_table_add(seen, key))
            gdis_model_add_bond(model, i, j, 1, TRUE);
        }
    }

  g_hash_table_destroy(seen);
  gdis_model_refresh_counts(model);
  return TRUE;
}

static void
gdis_model_refresh_counts(GdisModel *model)
{
  if (!model)
    return;

  model->atom_count = model->atoms ? model->atoms->len : 0;
  model->bond_count = model->bonds ? model->bonds->len : 0;
}

static void
gdis_model_update_identity(GdisModel *model,
                           const char *path,
                           GdisModelFormat format)
{
  g_return_if_fail(model != NULL);
  g_return_if_fail(path != NULL);

  g_free(model->path);
  g_free(model->basename);
  g_free(model->format_label);

  model->path = g_strdup(path);
  model->basename = g_path_get_basename(path);
  model->format = format;
  model->format_label = g_strdup(gdis_model_format_label(format));
}

static void
gdis_model_make_all_bonds_explicit(GdisModel *model)
{
  guint i;

  g_return_if_fail(model != NULL);
  g_return_if_fail(model->bonds != NULL);

  for (i = 0; i < model->bonds->len; i++)
    {
      GdisBond *bond;

      bond = &g_array_index(model->bonds, GdisBond, i);
      bond->inferred = FALSE;
    }
  model->explicit_bond_count = model->bonds->len;
}

static gint
gdis_model_find_bond_index(const GdisModel *model,
                           guint atom_index_a,
                           guint atom_index_b)
{
  guint i;
  guint min_index;
  guint max_index;

  g_return_val_if_fail(model != NULL, -1);
  g_return_val_if_fail(model->bonds != NULL, -1);

  min_index = MIN(atom_index_a, atom_index_b);
  max_index = MAX(atom_index_a, atom_index_b);

  for (i = 0; i < model->bonds->len; i++)
    {
      const GdisBond *bond;
      guint bond_min;
      guint bond_max;

      bond = &g_array_index(model->bonds, GdisBond, i);
      bond_min = MIN(bond->atom_index_a, bond->atom_index_b);
      bond_max = MAX(bond->atom_index_a, bond->atom_index_b);
      if (bond_min == min_index && bond_max == max_index)
        return (gint) i;
    }

  return -1;
}

static void
gdis_model_finalize_metadata(GdisModel *model)
{
  if (!model)
    return;

  if (!model->title || model->title[0] == '\0')
    {
      g_clear_pointer(&model->title, g_free);
      model->title = g_strdup(model->basename ? model->basename : "");
    }

  if ((!model->space_group || model->space_group[0] == '\0') && model->periodic)
    {
      g_clear_pointer(&model->space_group, g_free);
      model->space_group = g_strdup("P 1");
    }

  gdis_model_refresh_counts(model);
}

static GdisAtom *
gdis_atom_new(const gchar *label,
              const gchar *element,
              const gchar *ff_type,
              gdouble x,
              gdouble y,
              gdouble z,
              gdouble occupancy,
              gint region,
              guint serial)
{
  GdisAtom *atom;

  atom = g_new0(GdisAtom, 1);
  atom->serial = serial;
  atom->label = g_strdup(label ? label : "");
  atom->element = g_strdup(element ? element : "X");
  atom->ff_type = g_strdup((ff_type && ff_type[0] != '\0') ? ff_type :
                           (element && element[0] != '\0') ? element : "X");
  atom->position[0] = x;
  atom->position[1] = y;
  atom->position[2] = z;
  atom->occupancy = occupancy;
  atom->region = region;
  return atom;
}

static void
gdis_atom_free(gpointer data)
{
  GdisAtom *atom;

  atom = data;
  if (!atom)
    return;

  g_free(atom->label);
  g_free(atom->element);
  g_free(atom->ff_type);
  g_free(atom);
}

static void
gdis_cif_atom_record_free(gpointer data)
{
  GdisCifAtomRecord *record;

  record = data;
  if (!record)
    return;

  g_free(record->label);
  g_free(record->element);
  g_free(record);
}

static gchar *
gdis_model_write_xyz(const GdisModel *model)
{
  GString *out;
  guint i;

  out = g_string_new("");
  g_string_append_printf(out, "%u\n", model->atom_count);
  g_string_append_printf(out, "%s\n", model->title ? model->title : model->basename);

  for (i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;

      atom = g_ptr_array_index(model->atoms, i);
      g_string_append_printf(out,
                             "%-2s % .8f % .8f % .8f\n",
                             atom->element ? atom->element : "X",
                             atom->position[0],
                             atom->position[1],
                             atom->position[2]);
    }

  return g_string_free(out, FALSE);
}

static gchar *
gdis_model_write_pdb(const GdisModel *model)
{
  GString *out;
  guint i;

  out = g_string_new("");
  if (model->title && model->title[0] != '\0')
    g_string_append_printf(out, "TITLE     %s\n", model->title);

  if (model->periodic &&
      model->cell_lengths[0] > 0.0 &&
      model->cell_lengths[1] > 0.0 &&
      model->cell_lengths[2] > 0.0)
    {
      g_string_append_printf(out,
                             "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s\n",
                             model->cell_lengths[0],
                             model->cell_lengths[1],
                             model->cell_lengths[2],
                             model->cell_angles[0],
                             model->cell_angles[1],
                             model->cell_angles[2],
                             model->space_group ? model->space_group : "P 1");
    }

  for (i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      gchar atom_name[5];
      gchar residue_name[4];

      atom = g_ptr_array_index(model->atoms, i);
      g_snprintf(atom_name, sizeof(atom_name), "%-4.4s",
                 (atom->label && atom->label[0] != '\0') ? atom->label : atom->element);
      g_snprintf(residue_name, sizeof(residue_name), "%-3.3s", "MOL");

      g_string_append_printf(out,
                             "ATOM  %5u %-4s %-3s %1s%4u    %8.3f%8.3f%8.3f%6.2f%6.2f          %-2s\n",
                             atom->serial,
                             atom_name,
                             residue_name,
                             "A",
                             1u,
                             atom->position[0],
                             atom->position[1],
                             atom->position[2],
                             atom->occupancy,
                             0.0,
                             atom->element ? atom->element : "X");
    }

  for (i = 0; i < model->bonds->len; i++)
    {
      const GdisBond *bond;

      bond = &g_array_index(model->bonds, GdisBond, i);
      g_string_append_printf(out,
                             "CONECT%5u%5u\n",
                             bond->atom_index_a + 1,
                             bond->atom_index_b + 1);
    }

  g_string_append(out, "END\n");
  return g_string_free(out, FALSE);
}

static gchar *
gdis_model_write_arc_like(const GdisModel *model, GdisModelFormat format)
{
  GString *out;
  guint i;
  gint archive_version;

  out = g_string_new("");
  archive_version = (format == GDIS_MODEL_FORMAT_CAR) ? 3 : 2;

  g_string_append_printf(out, "!BIOSYM archive %d\n", archive_version);
  if (model->periodic && model->periodicity == 2)
    g_string_append(out, "PBC=2D\n");
  else
    g_string_append_printf(out, "PBC=%s\n", model->periodic ? "ON" : "OFF");
  g_string_append_printf(out, "%s\n", model->title ? model->title : model->basename);
  g_string_append(out, "!DATE\n");

  if (model->periodic)
    {
      if (model->periodicity == 2)
        {
          g_string_append_printf(out,
                                 "PBC %10.4f %10.4f %9.4f\n",
                                 model->cell_lengths[0],
                                 model->cell_lengths[1],
                                 model->cell_angles[2]);
        }
      else
        {
          g_string_append_printf(out,
                                 "PBC %10.4f %10.4f %10.4f %9.4f %9.4f %9.4f\n",
                                 model->cell_lengths[0],
                                 model->cell_lengths[1],
                                 model->cell_lengths[2],
                                 model->cell_angles[0],
                                 model->cell_angles[1],
                                 model->cell_angles[2]);
        }
    }

  for (i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      const char *type_token;
      gchar region_token[8];

      atom = g_ptr_array_index(model->atoms, i);
      type_token = "CORE";
      if (atom->region >= 0 && atom->region < 4)
        {
          g_snprintf(region_token, sizeof(region_token), "R%dAC", atom->region + 1);
          type_token = region_token;
        }
      g_string_append_printf(out,
                             "%-8s %14.9f %14.9f %14.9f %-4s %5u %-2s %-2s %8.4f %5u\n",
                             atom->label && atom->label[0] != '\0' ? atom->label : atom->element,
                             atom->position[0],
                             atom->position[1],
                             atom->position[2],
                             type_token,
                             atom->serial,
                             atom->ff_type ? atom->ff_type : (atom->element ? atom->element : "X"),
                             atom->element ? atom->element : "X",
                             0.0,
                             atom->serial);
    }

  g_string_append(out, "end\nend\n");
  return g_string_free(out, FALSE);
}

static gchar *
gdis_model_write_cif(const GdisModel *model, GError **error)
{
  GString *out;
  gdouble matrix[9];
  gdouble inverse[9];
  guint i;
  gchar *name;

  if (!model->periodic ||
      model->cell_lengths[0] <= 0.0 ||
      model->cell_lengths[1] <= 0.0 ||
      model->cell_lengths[2] <= 0.0 ||
      !gdis_model_build_cell_matrix(model, matrix, inverse))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "CIF export currently requires a periodic model with a valid unit cell.");
      return NULL;
    }

  out = g_string_new("");
  name = g_strdup(model->basename ? model->basename : "gdis_model");
  if (strrchr(name, '.'))
    *strrchr(name, '.') = '\0';

  g_string_append_printf(out, "data_%s\n", name);
  g_free(name);
  g_string_append_printf(out, "_symmetry_space_group_name_H-M '%s'\n",
                         model->space_group ? model->space_group : "P 1");
  g_string_append_printf(out, "_cell_length_a %.6f\n", model->cell_lengths[0]);
  g_string_append_printf(out, "_cell_length_b %.6f\n", model->cell_lengths[1]);
  g_string_append_printf(out, "_cell_length_c %.6f\n", model->cell_lengths[2]);
  g_string_append_printf(out, "_cell_angle_alpha %.6f\n", model->cell_angles[0]);
  g_string_append_printf(out, "_cell_angle_beta %.6f\n", model->cell_angles[1]);
  g_string_append_printf(out, "_cell_angle_gamma %.6f\n", model->cell_angles[2]);
  g_string_append(out,
                  "\nloop_\n"
                  "_atom_site_label\n"
                  "_atom_site_type_symbol\n"
                  "_atom_site_fract_x\n"
                  "_atom_site_fract_y\n"
                  "_atom_site_fract_z\n");

  for (i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      gdouble frac[3];

      atom = g_ptr_array_index(model->atoms, i);
      gdis_cart_to_frac(inverse, atom->position, frac);
      g_string_append_printf(out,
                             "%s %s %.8f %.8f %.8f\n",
                             atom->label && atom->label[0] != '\0' ? atom->label : atom->element,
                             atom->element ? atom->element : "X",
                             frac[0],
                             frac[1],
                             frac[2]);
    }

  return g_string_free(out, FALSE);
}

static gboolean
gdis_model_has_valid_cell(const GdisModel *model)
{
  return model &&
         model->cell_lengths[0] > 0.0 &&
         model->cell_lengths[1] > 0.0 &&
         model->cell_lengths[2] > 0.0;
}

static gboolean
gdis_model_invert_matrix3(const gdouble matrix[9], gdouble inverse[9])
{
  gdouble det;

  g_return_val_if_fail(matrix != NULL, FALSE);
  g_return_val_if_fail(inverse != NULL, FALSE);

  det = matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7]) -
        matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6]) +
        matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
  if (fabs(det) < 1.0e-12)
    return FALSE;

  inverse[0] =  (matrix[4] * matrix[8] - matrix[5] * matrix[7]) / det;
  inverse[1] = -(matrix[1] * matrix[8] - matrix[2] * matrix[7]) / det;
  inverse[2] =  (matrix[1] * matrix[5] - matrix[2] * matrix[4]) / det;
  inverse[3] = -(matrix[3] * matrix[8] - matrix[5] * matrix[6]) / det;
  inverse[4] =  (matrix[0] * matrix[8] - matrix[2] * matrix[6]) / det;
  inverse[5] = -(matrix[0] * matrix[5] - matrix[2] * matrix[3]) / det;
  inverse[6] =  (matrix[3] * matrix[7] - matrix[4] * matrix[6]) / det;
  inverse[7] = -(matrix[0] * matrix[7] - matrix[1] * matrix[6]) / det;
  inverse[8] =  (matrix[0] * matrix[4] - matrix[1] * matrix[3]) / det;
  return TRUE;
}

static void
gdis_vec3_set(gdouble vector[3], gdouble x, gdouble y, gdouble z)
{
  vector[0] = x;
  vector[1] = y;
  vector[2] = z;
}

static void
gdis_vec3_copy(gdouble dest[3], const gdouble src[3])
{
  dest[0] = src[0];
  dest[1] = src[1];
  dest[2] = src[2];
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
  if (length <= 1.0e-12)
    return FALSE;

  vector[0] /= length;
  vector[1] /= length;
  vector[2] /= length;
  return TRUE;
}

static gint
gdis_int_gcd(gint left, gint right)
{
  left = ABS(left);
  right = ABS(right);

  while (right != 0)
    {
      gint remainder;

      remainder = left % right;
      left = right;
      right = remainder;
    }

  return MAX(left, 1);
}

static void
gdis_int_vector_reduce(gint vector[3])
{
  gint divisor;

  divisor = gdis_int_gcd(vector[0], vector[1]);
  divisor = gdis_int_gcd(divisor, vector[2]);
  if (divisor <= 1)
    return;

  vector[0] /= divisor;
  vector[1] /= divisor;
  vector[2] /= divisor;
}

static gboolean
gdis_surface_build_inplane_indices(gint h,
                                   gint k,
                                   gint l,
                                   gint first[3],
                                   gint second[3])
{
  gint hkl[3];

  g_return_val_if_fail(first != NULL, FALSE);
  g_return_val_if_fail(second != NULL, FALSE);

  hkl[0] = h;
  hkl[1] = k;
  hkl[2] = l;

  if (h != 0 || k != 0)
    {
      first[0] = k;
      first[1] = -h;
      first[2] = 0;
    }
  else if (l != 0)
    {
      first[0] = 1;
      first[1] = 0;
      first[2] = 0;
    }
  else
    {
      return FALSE;
    }

  second[0] = hkl[1] * first[2] - hkl[2] * first[1];
  second[1] = hkl[2] * first[0] - hkl[0] * first[2];
  second[2] = hkl[0] * first[1] - hkl[1] * first[0];

  gdis_int_vector_reduce(first);
  gdis_int_vector_reduce(second);

  return !(first[0] == 0 && first[1] == 0 && first[2] == 0) &&
         !(second[0] == 0 && second[1] == 0 && second[2] == 0);
}

static gboolean
gdis_surface_build_cell_from_vectors(const gdouble a_vec[3],
                                     const gdouble b_vec[3],
                                     const gdouble c_vec[3],
                                     gdouble lengths[3],
                                     gdouble angles[3])
{
  gdouble a_length;
  gdouble b_length;
  gdouble c_length;

  g_return_val_if_fail(a_vec != NULL, FALSE);
  g_return_val_if_fail(b_vec != NULL, FALSE);
  g_return_val_if_fail(c_vec != NULL, FALSE);
  g_return_val_if_fail(lengths != NULL, FALSE);
  g_return_val_if_fail(angles != NULL, FALSE);

  a_length = gdis_vec3_length(a_vec);
  b_length = gdis_vec3_length(b_vec);
  c_length = gdis_vec3_length(c_vec);
  if (a_length <= 1.0e-8 || b_length <= 1.0e-8 || c_length <= 1.0e-8)
    return FALSE;

  lengths[0] = a_length;
  lengths[1] = b_length;
  lengths[2] = c_length;
  angles[0] = acos(CLAMP(gdis_vec3_dot(b_vec, c_vec) / (b_length * c_length), -1.0, 1.0)) * (180.0 / G_PI);
  angles[1] = acos(CLAMP(gdis_vec3_dot(a_vec, c_vec) / (a_length * c_length), -1.0, 1.0)) * (180.0 / G_PI);
  angles[2] = acos(CLAMP(gdis_vec3_dot(a_vec, b_vec) / (a_length * b_length), -1.0, 1.0)) * (180.0 / G_PI);
  return TRUE;
}

static gchar *
gdis_model_surface_default_path(const GdisModel *source,
                                gint h,
                                gint k,
                                gint l)
{
  g_autofree gchar *base = NULL;
  gchar *dot;

  g_return_val_if_fail(source != NULL, NULL);

  base = g_strdup(source->basename ? source->basename : "surface");
  dot = strrchr(base, '.');
  if (dot)
    *dot = '\0';

  return g_strdup_printf("%s-surface-%d_%d_%d.cif", base, h, k, l);
}

static gboolean
gdis_build_cell_matrix_from_parameters(const gdouble lengths[3],
                                       const gdouble angles[3],
                                       gdouble matrix[9],
                                       gdouble inverse[9])
{
  gdouble alpha;
  gdouble beta;
  gdouble gamma;
  gdouble ax;
  gdouble bx;
  gdouble by;
  gdouble cx;
  gdouble cy;
  gdouble cz_sq;

  if (!lengths || !angles || lengths[0] <= 0.0 || lengths[1] <= 0.0 || lengths[2] <= 0.0)
    return FALSE;

  alpha = angles[0] * (G_PI / 180.0);
  beta = angles[1] * (G_PI / 180.0);
  gamma = angles[2] * (G_PI / 180.0);

  ax = lengths[0];
  bx = lengths[1] * cos(gamma);
  by = lengths[1] * sin(gamma);
  cx = lengths[2] * cos(beta);

  if (fabs(by) < 1.0e-10)
    return FALSE;

  cy = lengths[2] * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma);
  cz_sq = lengths[2] * lengths[2] - cx * cx - cy * cy;
  if (cz_sq < 0.0 && fabs(cz_sq) < 1.0e-10)
    cz_sq = 0.0;
  if (cz_sq < 0.0)
    return FALSE;

  matrix[0] = ax;
  matrix[1] = bx;
  matrix[2] = cx;
  matrix[3] = 0.0;
  matrix[4] = by;
  matrix[5] = cy;
  matrix[6] = 0.0;
  matrix[7] = 0.0;
  matrix[8] = sqrt(cz_sq);

  if (!inverse)
    return TRUE;

  return gdis_model_invert_matrix3(matrix, inverse);
}

static gboolean
gdis_model_build_cell_matrix(const GdisModel *model, gdouble matrix[9], gdouble inverse[9])
{
  if (!gdis_model_has_valid_cell(model))
    return FALSE;

  return gdis_build_cell_matrix_from_parameters(model->cell_lengths,
                                                model->cell_angles,
                                                matrix,
                                                inverse);
}

static void
gdis_frac_to_cart(const gdouble matrix[9], const gdouble frac[3], gdouble cart[3])
{
  cart[0] = matrix[0] * frac[0] + matrix[1] * frac[1] + matrix[2] * frac[2];
  cart[1] = matrix[3] * frac[0] + matrix[4] * frac[1] + matrix[5] * frac[2];
  cart[2] = matrix[6] * frac[0] + matrix[7] * frac[1] + matrix[8] * frac[2];
}

static void
gdis_cart_to_frac(const gdouble inverse[9], const gdouble cart[3], gdouble frac[3])
{
  frac[0] = inverse[0] * cart[0] + inverse[1] * cart[1] + inverse[2] * cart[2];
  frac[1] = inverse[3] * cart[0] + inverse[4] * cart[1] + inverse[5] * cart[2];
  frac[2] = inverse[6] * cart[0] + inverse[7] * cart[1] + inverse[8] * cart[2];
}

static gdouble
gdis_model_distance2(const GdisModel *model,
                     const gdouble a[3],
                     const gdouble b[3],
                     const gdouble matrix[9],
                     const gdouble inverse[9],
                     gboolean use_periodic_image)
{
  gdouble diff[3];

  if (use_periodic_image)
    {
      gdouble frac_a[3];
      gdouble frac_b[3];
      gdouble frac_diff[3];
      guint dims;

      gdis_cart_to_frac(inverse, a, frac_a);
      gdis_cart_to_frac(inverse, b, frac_b);

      frac_diff[0] = frac_b[0] - frac_a[0];
      frac_diff[1] = frac_b[1] - frac_a[1];
      frac_diff[2] = frac_b[2] - frac_a[2];

      dims = MAX(0u, MIN(model->periodicity, 3u));
      if (dims > 0)
        frac_diff[0] -= floor(frac_diff[0] + 0.5);
      if (dims > 1)
        frac_diff[1] -= floor(frac_diff[1] + 0.5);
      if (dims > 2)
        frac_diff[2] -= floor(frac_diff[2] + 0.5);

      gdis_frac_to_cart(matrix, frac_diff, diff);
    }
  else
    {
      diff[0] = b[0] - a[0];
      diff[1] = b[1] - a[1];
      diff[2] = b[2] - a[2];
    }

  return diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
}

static gchar *
gdis_strdup_strip(const gchar *text)
{
  gchar *copy;

  if (!text)
    return g_strdup("");

  copy = g_strdup(text);
  g_strstrip(copy);
  return copy;
}

static gchar *
gdis_copy_field(const gchar *line, gsize start, gsize width)
{
  gsize len;
  gsize available;

  if (!line)
    return g_strdup("");

  len = strlen(line);
  if (start >= len)
    return g_strdup("");

  available = MIN(width, len - start);
  return g_strndup(line + start, available);
}

static gchar *
gdis_normalize_element_symbol(const gchar *text)
{
  gchar symbol[3];
  guint out;
  const gchar *p;

  symbol[0] = 'X';
  symbol[1] = '\0';
  symbol[2] = '\0';

  if (!text)
    return g_strdup(symbol);

  p = text;
  while (*p && !g_ascii_isalpha(*p))
    p++;

  if (!*p)
    return g_strdup(symbol);

  out = 0;
  symbol[out++] = g_ascii_toupper(*p++);

  while (*p && out < 2)
    {
      if (g_ascii_isalpha(*p))
        {
          symbol[out++] = g_ascii_tolower(*p);
          break;
        }
      if (g_ascii_isdigit(*p))
        break;
      p++;
    }

  symbol[out] = '\0';
  return g_strdup(symbol);
}

static gchar *
gdis_pdb_guess_element(const gchar *atom_name)
{
  return gdis_normalize_element_symbol(atom_name);
}

static gboolean
gdis_arc_like_type_is_marvin_label(const gchar *text)
{
  g_return_val_if_fail(text != NULL, FALSE);

  if (strlen(text) < 4)
    return FALSE;

  return (text[0] == 'R' || text[0] == 'r') &&
         g_ascii_isdigit(text[1]) &&
         (text[2] == 'A' || text[2] == 'a' ||
          text[2] == 'B' || text[2] == 'b') &&
         (text[3] == 'C' || text[3] == 'c' ||
          text[3] == 'S' || text[3] == 's');
}

static gint
gdis_arc_like_region_from_type(const gchar *text)
{
  g_return_val_if_fail(text != NULL, -1);

  if (!gdis_arc_like_type_is_marvin_label(text))
    return -1;

  return (gint) (text[1] - '1');
}

static gboolean
gdis_try_parse_double(const gchar *text, gdouble *value)
{
  gchar *endptr;
  gdouble parsed;

  g_return_val_if_fail(value != NULL, FALSE);

  if (!text)
    return FALSE;

  parsed = g_ascii_strtod(text, &endptr);
  if (endptr == text)
    return FALSE;

  while (*endptr != '\0')
    {
      if (!g_ascii_isspace(*endptr))
        return FALSE;
      endptr++;
    }

  *value = parsed;
  return TRUE;
}

static gboolean
gdis_try_parse_double_relaxed(const gchar *text, gdouble *value)
{
  g_autofree gchar *clean = NULL;
  gchar *endptr;
  gdouble parsed;
  gsize len;

  g_return_val_if_fail(value != NULL, FALSE);

  if (gdis_try_parse_double(text, value))
    return TRUE;

  if (!text)
    return FALSE;

  clean = g_strdup(text);
  g_strstrip(clean);
  if (clean[0] == '\0' || g_str_equal(clean, "*"))
    return FALSE;

  while (clean[0] == '(' || clean[0] == '[' || clean[0] == '{' || clean[0] == '<')
    memmove(clean, clean + 1, strlen(clean));

  for (gchar *p = clean; *p != '\0'; p++)
    {
      if (*p == 'd' || *p == 'D')
        *p = 'E';
    }

  len = strlen(clean);
  while (len > 0)
    {
      gchar tail;

      tail = clean[len - 1];
      if (g_ascii_isspace(tail) ||
          tail == '*' ||
          tail == ',' ||
          tail == ';' ||
          tail == ':' ||
          tail == ')' ||
          tail == ']' ||
          tail == '}' ||
          tail == '>')
        {
          clean[len - 1] = '\0';
          len--;
          continue;
        }
      break;
    }

  if (clean[0] == '\0')
    return FALSE;

  parsed = g_ascii_strtod(clean, &endptr);
  if (endptr == clean)
    return FALSE;

  while (*endptr != '\0')
    {
      if (!g_ascii_isspace(*endptr))
        return FALSE;
      endptr++;
    }

  *value = parsed;
  return TRUE;
}

static gboolean
gdis_try_parse_uint(const gchar *text, guint *value)
{
  gchar *endptr;
  guint64 parsed;

  g_return_val_if_fail(value != NULL, FALSE);

  if (!text)
    return FALSE;

  parsed = g_ascii_strtoull(text, &endptr, 10);
  if (endptr == text)
    return FALSE;

  while (*endptr != '\0')
    {
      if (!g_ascii_isspace(*endptr))
        return FALSE;
      endptr++;
    }

  if (parsed > G_MAXUINT)
    return FALSE;

  *value = (guint) parsed;
  return TRUE;
}

static guint
gdis_collect_doubles_from_line(const gchar *line, gdouble *values, guint max_values)
{
  g_auto(GStrv) tokens = NULL;
  guint found = 0u;
  gint token_count = 0;

  if (!line || !values || max_values == 0u)
    return 0u;

  tokens = gdis_split_simple(line, &token_count);
  for (gint i = 0; i < token_count && found < max_values; i++)
    {
      gdouble parsed;

      if (gdis_try_parse_double_relaxed(tokens[i], &parsed))
        values[found++] = parsed;
    }

  return found;
}

static gboolean
gdis_collect_three_doubles_from_tokens(gchar **tokens,
                                       gint token_count,
                                       gint start_index,
                                       gdouble out[3])
{
  gint i;
  guint found = 0u;

  g_return_val_if_fail(out != NULL, FALSE);

  if (!tokens || token_count <= 0 || start_index >= token_count)
    return FALSE;

  for (i = MAX(start_index, 0); i < token_count && found < 3u; i++)
    {
      gdouble parsed;

      if (gdis_try_parse_double_relaxed(tokens[i], &parsed))
        out[found++] = parsed;
    }

  return found == 3u;
}

static gchar **
gdis_split_simple(const gchar *line, gint *count)
{
  gchar **raw_tokens;
  GPtrArray *tokens;
  gint i;
  gchar **result;

  raw_tokens = g_strsplit_set(line ? line : "", " \t", -1);
  tokens = g_ptr_array_new();

  for (i = 0; raw_tokens[i] != NULL; i++)
    {
      if (raw_tokens[i][0] != '\0')
        g_ptr_array_add(tokens, raw_tokens[i]);
      else
        g_free(raw_tokens[i]);
    }

  g_free(raw_tokens);

  g_ptr_array_add(tokens, NULL);
  result = (gchar **) g_ptr_array_free(tokens, FALSE);

  if (count)
    *count = g_strv_length(result);

  return result;
}

static GPtrArray *
gdis_tokenize_cif(const gchar *line)
{
  GPtrArray *tokens;
  const gchar *p;

  tokens = g_ptr_array_new_with_free_func(g_free);
  p = line;

  while (p && *p)
    {
      GString *token;
      gchar quote;

      while (*p && g_ascii_isspace(*p))
        p++;

      if (!*p || *p == '#')
        break;

      token = g_string_new(NULL);
      quote = '\0';

      if (*p == '\'' || *p == '"')
        {
          quote = *p;
          p++;
        }

      while (*p)
        {
          if (quote != '\0')
            {
              if (*p == quote)
                {
                  p++;
                  break;
                }
            }
          else
            {
              if (g_ascii_isspace(*p) || *p == '#')
                break;
            }

          g_string_append_c(token, *p);
          p++;
        }

      g_ptr_array_add(tokens, g_string_free(token, FALSE));

      if (*p == '#')
        break;
    }

  return tokens;
}

static gboolean
gdis_cif_parse_value_line(const gchar *line, gchar **tag_out, gchar **value_out)
{
  GPtrArray *tokens;
  gboolean ok;

  g_return_val_if_fail(tag_out != NULL, FALSE);
  g_return_val_if_fail(value_out != NULL, FALSE);

  *tag_out = NULL;
  *value_out = NULL;

  tokens = gdis_tokenize_cif(line);
  ok = tokens->len >= 2;
  if (ok)
    {
      *tag_out = g_strdup(g_ptr_array_index(tokens, 0));
      *value_out = g_strdup(g_ptr_array_index(tokens, 1));
    }
  g_ptr_array_free(tokens, TRUE);

  return ok;
}

static gdouble
gdis_parse_cif_number(const gchar *text)
{
  gchar *copy;
  gchar *paren;
  gdouble value;

  if (!text)
    return 0.0;

  if (g_str_equal(text, ".") || g_str_equal(text, "?"))
    return 0.0;

  copy = g_strdup(text);
  g_strstrip(copy);

  paren = strchr(copy, '(');
  if (paren)
    *paren = '\0';

  value = g_ascii_strtod(copy, NULL);
  g_free(copy);
  return value;
}

static gboolean
gdis_cif_is_control_line(const gchar *line)
{
  if (!line || line[0] == '\0')
    return TRUE;

  if (line[0] == '#')
    return TRUE;

  return g_str_has_prefix(line, "loop_") ||
         g_str_has_prefix(line, "data_") ||
         line[0] == '_';
}

static gdouble
gdis_lookup_covalent_radius(const gchar *element)
{
  static const GdisRadiusEntry radii[] = {
    { "H", 0.31 },  { "B", 0.85 },  { "C", 0.76 },  { "N", 0.71 },
    { "O", 0.66 },  { "F", 0.57 },  { "Na", 1.66 }, { "Mg", 1.41 },
    { "Al", 1.21 }, { "Si", 1.11 }, { "P", 1.07 },  { "S", 1.05 },
    { "Cl", 1.02 }, { "K", 2.03 },  { "Ca", 1.76 }, { "Ti", 1.60 },
    { "Cr", 1.39 }, { "Mn", 1.39 }, { "Fe", 1.32 }, { "Co", 1.26 },
    { "Ni", 1.24 }, { "Cu", 1.32 }, { "Zn", 1.22 }, { "Br", 1.20 },
    { "Ag", 1.45 }, { "I", 1.39 },  { "Ba", 1.98 }, { "Pt", 1.36 },
    { "Au", 1.36 }, { "Pb", 1.46 }, { "Se", 1.20 }, { "As", 1.19 },
    { "Li", 1.28 }, { "Sn", 1.39 }
  };
  guint i;

  if (!element || element[0] == '\0')
    return 0.77;

  for (i = 0; i < G_N_ELEMENTS(radii); i++)
    {
      if (g_ascii_strcasecmp(radii[i].symbol, element) == 0)
        return radii[i].radius;
    }

  return 0.77;
}
