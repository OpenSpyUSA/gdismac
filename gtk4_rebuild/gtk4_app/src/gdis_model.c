#include "gdis_model.h"
#include "gdis_elements.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../../legacy_snapshot/src/sginfo.h"

typedef struct
{
  gchar *label;
  gchar *element;
  gdouble coords[3];
  gdouble occupancy;
  gboolean fractional;
  guint multiplicity;
  gboolean have_multiplicity;
} GdisCifAtomRecord;

typedef struct
{
  gchar *title;
  gchar *space_group;
  gboolean periodic;
  guint periodicity;
  gdouble cell_lengths[3];
  gdouble cell_angles[3];
  gboolean has_energy;
  gdouble energy_ev;
  gboolean has_force_rms;
  gdouble force_rms_ev_ang;
  GPtrArray *atoms;
} GdisModelFrame;

typedef struct
{
  GdisModel *model;
  GHashTable *species_symbols;
  GPtrArray *frames;
  GdisModelFrame *current_frame;
  gchar *current_species_name;
  gchar *current_atom_name;
  gchar *current_atom_species;
  GString *text;
  gdouble current_atom_position[3];
  gdouble default_a_vec[3];
  gdouble default_b_vec[3];
  gdouble default_c_vec[3];
  guint atomset_count;
  guint current_iteration_count;
  gboolean inside_iteration;
  gboolean inside_species_symbol;
  gboolean inside_atom_position;
  gboolean inside_etotal;
  gboolean inside_atom_force;
  gboolean have_current_atom_position;
  gboolean have_current_atom_force;
  gboolean current_frame_has_cell;
  gboolean have_default_cell;
  gboolean have_current_iteration_energy;
  gdouble current_iteration_energy_ev;
  gdouble current_atom_force[3];
  gdouble current_frame_force_sum_sq;
  guint current_frame_force_count;
  gboolean looks_like_qbox;
} GdisQboxXmlParseState;

typedef struct
{
  const char *symbol;
  gdouble radius;
} GdisRadiusEntry;

typedef struct
{
  gdouble matrix[9];
  gdouble offset[3];
} GdisSymmetryOp;

typedef struct
{
  gdouble values[3];
} GdisFractionalPosition;

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
static gboolean gdis_model_load_vasp_xml(GdisModel *model, const gchar *contents, GError **error);
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
static guint gdis_cif_sum_reported_multiplicities(const GPtrArray *records);
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
static gboolean gdis_model_atoms_inside_primary_cell(GPtrArray *atoms,
                                                     guint periodicity,
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
static gboolean gdis_gulp_type_is_core(const gchar *text);
static gboolean gdis_gulp_type_is_shell(const gchar *text);
static gchar *gdis_gulp_make_species_key(const gchar *name, const gchar *type_text);
static const gdouble *gdis_gulp_lookup_species_charge(GHashTable *species_charges,
                                                      const gchar *label,
                                                      const gchar *element,
                                                      const gchar *type_text);
static guint gdis_collect_doubles_from_line(const gchar *line,
                                            gdouble *values,
                                            guint max_values);
static gchar *gdis_xtl_compose_space_group(const gchar *label,
                                           guint number,
                                           const gchar *qualifier);
static gboolean gdis_collect_three_doubles_from_tokens(gchar **tokens,
                                                       gint token_count,
                                                       gint start_index,
                                                       gdouble out[3]);
static gchar **gdis_split_simple(const gchar *line, gint *count);
static GPtrArray *gdis_tokenize_cif(const gchar *line);
static gboolean gdis_cif_parse_value_line(const gchar *line, gchar **tag_out, gchar **value_out);
static gdouble gdis_parse_cif_number(const gchar *text);
static gboolean gdis_cif_is_control_line(const gchar *line);
static gboolean gdis_parse_cif_symmetry_operation(const gchar *text, GdisSymmetryOp *op_out);
gboolean gdis_model_set_element_covalent_override(GdisModel *model,
                                                  const gchar *symbol,
                                                  gdouble radius);
void gdis_model_copy_element_overrides(GdisModel *dest, const GdisModel *src);
gdouble gdis_model_lookup_covalent_radius(const GdisModel *model, const GdisAtom *atom);
static gdouble gdis_lookup_covalent_radius(const gchar *element);
gboolean gdis_model_should_infer_bond(const GdisModel *model,
                                      const GdisAtom *atom_i,
                                      const GdisAtom *atom_j,
                                      gdouble radius_i,
                                      gdouble radius_j,
                                      gdouble *max_distance_out);
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
static GArray *gdis_double_array_new(void);
static GArray *gdis_double_array_clone(const GArray *source);
static GPtrArray *gdis_double_array_ptr_array_clone(const GPtrArray *source);
static void gdis_model_clear_plot_arrays(GdisModel *model);
static void gdis_model_build_stick_spectrum_arrays(const GArray *positions,
                                                   const GArray *intensities,
                                                   GArray **x_values_out,
                                                   GArray **y_values_out);
static guint gdis_collect_doubles_from_text_scan(const gchar *text,
                                                 gdouble *values,
                                                 guint max_values);
static guint gdis_model_parse_vasp_atom_symbols(const gchar *contents,
                                                GPtrArray *symbols_out);
static gboolean gdis_model_parse_vasp_structure_block(GdisModel *model,
                                                      gchar **lines,
                                                      guint start_index,
                                                      const GPtrArray *atom_symbols,
                                                      GPtrArray **atoms_out,
                                                      gdouble lengths_out[3],
                                                      gdouble angles_out[3],
                                                      guint *end_index_out);
static gboolean gdis_model_build_band_arrays_from_triplets(const GArray *kpoint_triplets,
                                                           GPtrArray *kpoint_bands,
                                                           gdouble fermi_energy_ev,
                                                           GArray **x_values_out,
                                                           GArray **y_values_out,
                                                           guint *path_count_out,
                                                           guint *series_count_out);
static gboolean gdis_model_parse_vasp_dos_block(gchar **lines,
                                                guint start_index,
                                                GArray **dos_x_values_out,
                                                GArray **dos_y_values_out,
                                                gboolean *have_fermi_energy_out,
                                                gdouble *fermi_energy_ev_out,
                                                guint *end_index_out);
static gboolean gdis_model_parse_vasp_eigenvalues_block(gchar **lines,
                                                        guint start_index,
                                                        const GArray *kpoint_triplets,
                                                        gboolean have_fermi_energy,
                                                        gdouble fermi_energy_ev,
                                                        GArray **band_x_values_out,
                                                        GArray **band_y_values_out,
                                                        guint *path_count_out,
                                                        guint *series_count_out,
                                                        guint *end_index_out);
static gint gdis_str_ibegin(const char *s1, const char *s2);
static gint gdis_build_sginfo_compat(T_SgInfo *sg_info, const gchar *space_group_text);
static void gdis_free_sginfo_compat(T_SgInfo *sg_info);
static gchar *gdis_format_space_group_label(const T_SgInfo *sg_info);
static gboolean gdis_model_uses_rhombohedral_setting(const GdisModel *model,
                                                     const gchar *space_group_text);
static gboolean gdis_build_space_group_operations(const GdisModel *model,
                                                  const gchar *space_group_text,
                                                  gint origin_choice,
                                                  GArray **ops_out,
                                                  gchar **canonical_label_out);
static gdouble gdis_wrap_fractional_coordinate(gdouble value);
static void gdis_wrap_fractional_vector(gdouble frac[3], guint periodicity);
static gboolean gdis_fractional_positions_equal(const gdouble left[3],
                                                const gdouble right[3],
                                                guint periodicity,
                                                gdouble tolerance);
static GdisAtom *gdis_atom_clone_with_position(const GdisAtom *source,
                                               const gdouble position[3],
                                               guint serial);
static gboolean gdis_model_expand_atoms_with_symmetry_ops(GdisModel *model,
                                                          const GArray *ops);
gboolean gdis_model_expand_space_group_atoms(GdisModel *model,
                                             gint origin_choice);
static gboolean gdis_model_expand_space_group_atoms_best_match(GdisModel *model,
                                                               guint expected_atom_count);
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
gboolean gdis_qbox_apply_cell_vectors_to_frame(GdisQboxXmlParseState *state,
                                               GdisModelFrame *frame,
                                               const gdouble a_vec[3],
                                               const gdouble b_vec[3],
                                               const gdouble c_vec[3],
                                               GError **error);
gboolean gdis_qbox_begin_frame(GdisQboxXmlParseState *state,
                               GError **error);
gboolean gdis_qbox_finish_current_frame(GdisQboxXmlParseState *state,
                                        GError **error);
static gchar *gdis_model_surface_default_path(const GdisModel *source,
                                              gint h,
                                              gint k,
                                              gint l);
static gboolean gdis_model_sum_charge_atoms(const GPtrArray *atoms,
                                            gboolean periodic,
                                            guint periodicity,
                                            const gdouble lengths[3],
                                            const gdouble angles[3],
                                            const gchar *space_group,
                                            gint origin_choice,
                                            gboolean expand_space_group,
                                            gdouble *total_charge_out);
static gboolean gdis_model_sum_known_charge_atoms(const GPtrArray *atoms,
                                                  gboolean periodic,
                                                  guint periodicity,
                                                  const gdouble lengths[3],
                                                  const gdouble angles[3],
                                                  const gchar *space_group,
                                                  gint origin_choice,
                                                  gboolean expand_space_group,
                                                  gdouble *total_charge_out);
static gboolean gdis_model_sum_legacy_visible_charge(const GPtrArray *atoms,
                                                     gdouble *total_charge_out);

static const gdouble GDIS_QBOX_BOHR_TO_ANGSTROM = 0.529177210903;
static const gdouble GDIS_QBOX_HARTREE_TO_EV = 27.211386245988;
static const gdouble GDIS_QBOX_HARTREE_PER_BOHR_TO_EV_PER_ANG =
  GDIS_QBOX_HARTREE_TO_EV / GDIS_QBOX_BOHR_TO_ANGSTROM;

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
    case GDIS_MODEL_FORMAT_VASP_XML:
      ok = gdis_model_load_vasp_xml(model, contents, error);
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
  copy->has_total_charge = model->has_total_charge;
  copy->total_charge_e = model->total_charge_e;
  copy->has_density = model->has_density;
  copy->density_g_cm3 = model->density_g_cm3;
  copy->has_energy = model->has_energy;
  copy->energy_ev = model->energy_ev;
  copy->has_force_rms = model->has_force_rms;
  copy->force_rms_ev_ang = model->force_rms_ev_ang;
  copy->has_pressure = model->has_pressure;
  copy->pressure_gpa = model->pressure_gpa;
  copy->has_fermi_energy = model->has_fermi_energy;
  copy->fermi_energy_ev = model->fermi_energy_ev;
  copy->pressure_x_values = gdis_double_array_clone(model->pressure_x_values);
  copy->pressure_y_values_gpa = gdis_double_array_clone(model->pressure_y_values_gpa);
  copy->dos_x_values_ev = gdis_double_array_clone(model->dos_x_values_ev);
  copy->dos_y_values = gdis_double_array_clone(model->dos_y_values);
  copy->band_x_values = gdis_double_array_clone(model->band_x_values);
  copy->band_y_values_ev = gdis_double_array_clone(model->band_y_values_ev);
  copy->band_path_count = model->band_path_count;
  copy->band_series_count = model->band_series_count;
  copy->frequency_x_values_cm1 = gdis_double_array_clone(model->frequency_x_values_cm1);
  copy->frequency_y_values = gdis_double_array_clone(model->frequency_y_values);
  copy->raman_x_values_cm1 = gdis_double_array_clone(model->raman_x_values_cm1);
  copy->raman_y_values = gdis_double_array_clone(model->raman_y_values);
  copy->atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  copy->bonds = g_array_sized_new(FALSE,
                                  FALSE,
                                  sizeof(GdisBond),
                                  model->bonds ? model->bonds->len : 0u);
  copy->element_covalent_overrides = NULL;

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
          ((GdisAtom *) g_ptr_array_index(copy->atoms, copy->atoms->len - 1u))->charge = atom->charge;
          ((GdisAtom *) g_ptr_array_index(copy->atoms, copy->atoms->len - 1u))->has_charge = atom->has_charge;
        }
    }

  if (model->bonds && model->bonds->len > 0)
    g_array_append_vals(copy->bonds, model->bonds->data, model->bonds->len);

  gdis_model_copy_element_overrides(copy, model);

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
  gchar *preserved_path;
  gchar *preserved_basename;
  gchar *preserved_format_label;
  GdisModelFormat preserved_format;

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

  preserved_path = g_strdup(dest->path);
  preserved_basename = g_strdup(dest->basename);
  preserved_format_label = g_strdup(dest->format_label);
  preserved_format = dest->format;

  gdis_model_clear_contents(dest);
  *dest = *copy;
  g_free(copy);
  g_free(dest->path);
  g_free(dest->basename);
  g_free(dest->format_label);
  dest->path = preserved_path;
  dest->basename = preserved_basename;
  dest->format_label = preserved_format_label;
  dest->format = preserved_format;
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
  g_clear_pointer(&model->element_covalent_overrides, g_hash_table_unref);
  gdis_model_clear_plot_arrays(model);

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

static GArray *
gdis_double_array_new(void)
{
  return g_array_new(FALSE, FALSE, sizeof(gdouble));
}

static GArray *
gdis_double_array_clone(const GArray *source)
{
  GArray *copy;

  if (!source)
    return NULL;

  copy = gdis_double_array_new();
  if (source->len > 0u)
    g_array_append_vals(copy, source->data, source->len);
  return copy;
}

static GPtrArray *
gdis_double_array_ptr_array_clone(const GPtrArray *source)
{
  GPtrArray *copy;

  if (!source)
    return NULL;

  copy = g_ptr_array_new_with_free_func((GDestroyNotify) g_array_unref);
  for (guint i = 0; i < source->len; i++)
    {
      const GArray *row;

      row = g_ptr_array_index((GPtrArray *) source, i);
      if (!row)
        continue;
      g_ptr_array_add(copy, gdis_double_array_clone(row));
    }
  return copy;
}

static void
gdis_model_clear_plot_arrays(GdisModel *model)
{
  g_return_if_fail(model != NULL);

  g_clear_pointer(&model->pressure_x_values, g_array_unref);
  g_clear_pointer(&model->pressure_y_values_gpa, g_array_unref);
  g_clear_pointer(&model->dos_x_values_ev, g_array_unref);
  g_clear_pointer(&model->dos_y_values, g_array_unref);
  g_clear_pointer(&model->band_x_values, g_array_unref);
  g_clear_pointer(&model->band_y_values_ev, g_array_unref);
  g_clear_pointer(&model->frequency_x_values_cm1, g_array_unref);
  g_clear_pointer(&model->frequency_y_values, g_array_unref);
  g_clear_pointer(&model->raman_x_values_cm1, g_array_unref);
  g_clear_pointer(&model->raman_y_values, g_array_unref);
  model->has_pressure = FALSE;
  model->pressure_gpa = 0.0;
  model->has_fermi_energy = FALSE;
  model->fermi_energy_ev = 0.0;
  model->band_path_count = 0u;
  model->band_series_count = 0u;
}

static void
gdis_model_build_stick_spectrum_arrays(const GArray *positions,
                                       const GArray *intensities,
                                       GArray **x_values_out,
                                       GArray **y_values_out)
{
  GArray *x_values;
  GArray *y_values;
  guint count;

  g_return_if_fail(x_values_out != NULL);
  g_return_if_fail(y_values_out != NULL);

  *x_values_out = NULL;
  *y_values_out = NULL;

  if (!positions || positions->len == 0u)
    return;

  count = positions->len;
  if (intensities && intensities->len > 0u)
    count = MIN(count, intensities->len);

  if (count == 0u)
    return;

  x_values = gdis_double_array_new();
  y_values = gdis_double_array_new();

  for (guint i = 0; i < count; i++)
    {
      const gdouble x_value = g_array_index(positions, gdouble, i);
      const gdouble y_value = intensities
                                ? fabs(g_array_index(intensities, gdouble, i))
                                : 1.0;
      const gdouble baseline = 0.0;

      g_array_append_val(x_values, x_value);
      g_array_append_val(y_values, baseline);
      g_array_append_val(x_values, x_value);
      g_array_append_val(y_values, y_value);
      g_array_append_val(x_values, x_value);
      g_array_append_val(y_values, baseline);
    }

  *x_values_out = x_values;
  *y_values_out = y_values;
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
    case GDIS_MODEL_FORMAT_VASP_XML:
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_UNSUPPORTED_FORMAT,
                  "Saving VASP XML is not available from this app.");
      break;
    case GDIS_MODEL_FORMAT_QBOX_XML:
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_UNSUPPORTED_FORMAT,
                  "Saving Qbox XML is not available from this app.");
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
                  "Saving %s is not available from this app.",
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
                  "Surface construction requires a fully 3D periodic crystal model.");
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
  gdis_model_copy_element_overrides(surface, source);

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
        format = GDIS_MODEL_FORMAT_VASP_XML;
      else
        format = GDIS_MODEL_FORMAT_QBOX_XML;
    }
  else if (g_str_equal(lower, "r"))
    format = GDIS_MODEL_FORMAT_QBOX_XML;
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
    case GDIS_MODEL_FORMAT_VASP_XML:
      return "VASP XML";
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

guint
gdis_model_get_component_count(const GdisModel *model)
{
  guint *component_ids;
  guint component_count;

  g_return_val_if_fail(model != NULL, 0u);

  component_ids = gdis_model_build_component_ids(model, &component_count);
  g_free(component_ids);
  return component_count;
}

gboolean
gdis_model_get_cell_volume(const GdisModel *model, gdouble *volume_out)
{
  gdouble matrix[9];
  gdouble inverse[9];
  gdouble a_vec[3];
  gdouble b_vec[3];
  gdouble c_vec[3];
  gdouble cross_ab[3];

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(volume_out != NULL, FALSE);

  *volume_out = 0.0;
  if (!gdis_model_build_cell_matrix(model, matrix, inverse))
    return FALSE;

  gdis_vec3_set(a_vec, matrix[0], matrix[3], matrix[6]);
  gdis_vec3_set(b_vec, matrix[1], matrix[4], matrix[7]);
  gdis_vec3_set(c_vec, matrix[2], matrix[5], matrix[8]);
  gdis_vec3_cross(a_vec, b_vec, cross_ab);
  *volume_out = fabs(gdis_vec3_dot(cross_ab, c_vec));
  return TRUE;
}

gboolean
gdis_model_get_cell_matrix(const GdisModel *model,
                           gdouble matrix[9],
                           gdouble inverse[9])
{
  return gdis_model_build_cell_matrix(model, matrix, inverse);
}

void
gdis_model_compute_minimum_image_delta(const GdisModel *model,
                                       const gdouble matrix[9],
                                       const gdouble inverse[9],
                                       const gdouble from[3],
                                       const gdouble to[3],
                                       gdouble delta[3])
{
  g_return_if_fail(model != NULL);
  g_return_if_fail(matrix != NULL);
  g_return_if_fail(inverse != NULL);
  g_return_if_fail(from != NULL);
  g_return_if_fail(to != NULL);
  g_return_if_fail(delta != NULL);

  if (!model->periodic)
    {
      delta[0] = to[0] - from[0];
      delta[1] = to[1] - from[1];
      delta[2] = to[2] - from[2];
      return;
    }

  gdis_model_minimum_image_delta(model, matrix, inverse, from, to, delta);
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
  model->element_covalent_overrides = NULL;
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
  copy->has_energy = frame->has_energy;
  copy->energy_ev = frame->energy_ev;
  copy->has_force_rms = frame->has_force_rms;
  copy->force_rms_ev_ang = frame->force_rms_ev_ang;

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
          ((GdisAtom *) g_ptr_array_index(copy->atoms, copy->atoms->len - 1u))->charge = atom->charge;
          ((GdisAtom *) g_ptr_array_index(copy->atoms, copy->atoms->len - 1u))->has_charge = atom->has_charge;
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
          ((GdisAtom *) g_ptr_array_index(atoms, atoms->len - 1u))->charge = atom->charge;
          ((GdisAtom *) g_ptr_array_index(atoms, atoms->len - 1u))->has_charge = atom->has_charge;
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
  model->has_energy = frame->has_energy;
  model->energy_ev = frame->energy_ev;
  model->has_force_rms = frame->has_force_rms;
  model->force_rms_ev_ang = frame->force_rms_ev_ang;

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
gdis_model_atoms_inside_primary_cell(GPtrArray *atoms,
                                     guint periodicity,
                                     const gdouble inverse[9])
{
  static const gdouble tolerance = 1.0e-5;
  guint dims;

  g_return_val_if_fail(atoms != NULL, FALSE);
  g_return_val_if_fail(inverse != NULL, FALSE);

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
        {
          if (frac[axis] < -tolerance || frac[axis] > 1.0 + tolerance)
            return FALSE;
        }
    }

  return TRUE;
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
  gboolean prefer_atom_wrap;

  g_return_if_fail(model != NULL);

  if (!model->periodic ||
      !model->atoms ||
      model->atoms->len == 0u ||
      !gdis_model_build_cell_matrix(model, matrix, inverse))
    return;

  /* Legacy BIOSYM ARC/CAR import wrapped each atom into the primary cell
   * while reading the Cartesian coordinates. Keep that default presentation
   * in the GTK4 rebuild so crystal networks do not open with part of the
   * bonded framework sitting outside the drawn unit cell. */
  prefer_atom_wrap = (model->format == GDIS_MODEL_FORMAT_ARC ||
                      model->format == GDIS_MODEL_FORMAT_CAR ||
                      model->format == GDIS_MODEL_FORMAT_GULP_INPUT);

  if (prefer_atom_wrap)
    gdis_model_wrap_atoms_to_cell(model->atoms,
                                  model->periodicity,
                                  matrix,
                                  inverse);
  else
    {
      if (!gdis_model_repack_periodic_components(model, matrix, inverse) ||
          !gdis_model_atoms_inside_primary_cell(model->atoms,
                                                model->periodicity,
                                                inverse))
        {
          gdis_model_wrap_atoms_to_cell(model->atoms,
                                        model->periodicity,
                                        matrix,
                                        inverse);
        }
    }

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

          if (prefer_atom_wrap)
            gdis_model_wrap_atoms_to_cell(frame->atoms,
                                          frame->periodicity,
                                          matrix,
                                          inverse);
          else
            {
              GdisModel snapshot;

              memset(&snapshot, 0, sizeof(snapshot));
              snapshot.atoms = frame->atoms;
              snapshot.bonds = model->bonds;
              snapshot.periodic = frame->periodic;
              snapshot.periodicity = frame->periodicity;
              memcpy(snapshot.cell_lengths, frame->cell_lengths, sizeof(snapshot.cell_lengths));
              memcpy(snapshot.cell_angles, frame->cell_angles, sizeof(snapshot.cell_angles));
              if (!gdis_model_repack_periodic_components(&snapshot, matrix, inverse) ||
                  !gdis_model_atoms_inside_primary_cell(frame->atoms,
                                                        frame->periodicity,
                                                        inverse))
                gdis_model_wrap_atoms_to_cell(frame->atoms,
                                              frame->periodicity,
                                              matrix,
                                              inverse);
            }
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
          GdisAtom *atom;
          gchar *element;
          gchar *ff_type;
          gdouble x;
          gdouble y;
          gdouble z;
          gdouble occupancy;
          gdouble charge;
          gboolean have_charge;
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

          charge = 0.0;
          have_charge = token_count >= 9 &&
                        gdis_try_parse_double_relaxed(tokens[8], &charge);

          ff_type = NULL;
          if (token_count >= 7)
            ff_type = gdis_strdup_strip(tokens[6]);

          if (token_count >= 8)
            element = gdis_normalize_element_symbol(tokens[7]);
          else
            element = gdis_normalize_element_symbol(tokens[0]);

          atom = gdis_atom_new(tokens[0],
                               element,
                               (ff_type && ff_type[0] != '\0') ? ff_type : element,
                               x,
                               y,
                               z,
                               occupancy,
                               region,
                               current_frame->atoms->len + 1u);
          if (have_charge)
            {
              atom->charge = charge;
              atom->has_charge = TRUE;
            }
          g_ptr_array_add(current_frame->atoms, atom);
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
  if (gdis_model_sum_known_charge_atoms(model->atoms,
                                        model->periodic,
                                        model->periodicity,
                                        model->cell_lengths,
                                        model->cell_angles,
                                        model->space_group,
                                        0,
                                        FALSE,
                                        &model->total_charge_e))
    model->has_total_charge = TRUE;
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
  gboolean saw_fractional_atoms;
  gboolean saw_non_fractional_atoms;
  gboolean ignore_additional_geometry;
  gboolean inside_element_block;
  GHashTable *species_charges;
  GPtrArray *charge_atoms;
  gint origin_choice;
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
  saw_fractional_atoms = FALSE;
  saw_non_fractional_atoms = FALSE;
  ignore_additional_geometry = FALSE;
  inside_element_block = FALSE;
  species_charges = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
  charge_atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  origin_choice = 0;
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

      if (inside_element_block ||
          g_str_equal(first_lower, "element") ||
          g_str_equal(first_lower, "elements"))
        {
          const gint property_index = inside_element_block ? 0 : 1;

          if (!inside_element_block)
            {
              inside_element_block = TRUE;
              coord_mode = GDIS_GULP_COORD_NONE;
            }

          if (token_count > property_index &&
              g_ascii_strcasecmp(tokens[property_index], "end") == 0)
            {
              inside_element_block = FALSE;
              continue;
            }

          if (token_count >= property_index + 3)
            {
              g_autofree gchar *property_lower = NULL;
              const gchar *symbol_text;
              const GdisElementInfo *element_info = NULL;
              gdouble radius;
              gchar *symbol = NULL;
              guint atomic_number;

              property_lower = g_ascii_strdown(tokens[property_index], -1);
              if (!g_str_equal(property_lower, "cova") &&
                  !g_str_equal(property_lower, "covalent"))
                continue;

              symbol_text = tokens[property_index + 1];
              if (gdis_try_parse_uint(symbol_text, &atomic_number))
                element_info = gdis_element_lookup_atomic_number(atomic_number);

              if (element_info && element_info->symbol && element_info->symbol[0] != '\0')
                symbol = g_strdup(element_info->symbol);
              else
                symbol = gdis_normalize_element_symbol(symbol_text);

              if (symbol &&
                  symbol[0] != '\0' &&
                  gdis_try_parse_double_relaxed(tokens[property_index + 2], &radius))
                gdis_model_set_element_covalent_override(model, symbol, radius);

              g_free(symbol);
            }
          continue;
        }

      if (g_str_equal(first_lower, "species"))
        {
          gboolean species_count_limited = FALSE;
          guint species_remaining = 0u;

          if (token_count > 1 && gdis_try_parse_uint(tokens[1], &species_remaining))
            species_count_limited = TRUE;

          for (;;)
            {
              guint next_index;
              gchar *species_line;
              gchar *species_trimmed;
              g_auto(GStrv) species_tokens = NULL;
              gint species_token_count = 0;

              if (species_count_limited && species_remaining == 0u)
                break;

              next_index = i + 1u;
              if (lines[next_index] == NULL)
                break;

              species_line = lines[next_index];
              g_strchomp(species_line);
              species_trimmed = g_strstrip(species_line);
              if (species_trimmed[0] == '\0' || species_trimmed[0] == '#')
                {
                  i = next_index;
                  continue;
                }

              species_tokens = gdis_split_simple(species_trimmed, &species_token_count);
              if (species_token_count <= 0)
                {
                  i = next_index;
                  continue;
                }

              if (g_ascii_strcasecmp(species_tokens[0], "end") == 0)
                {
                  i = next_index;
                  break;
                }

              {
                const gchar *type_name;
                gint charge_index;
                gdouble parsed_charge;

                if (species_token_count >= 3 &&
                    (gdis_gulp_type_is_core(species_tokens[1]) ||
                     gdis_gulp_type_is_shell(species_tokens[1])))
                  {
                    type_name = species_tokens[1];
                    charge_index = 2;
                  }
                else if (species_token_count >= 2)
                  {
                    type_name = "core";
                    charge_index = 1;
                  }
                else
                  {
                    break;
                  }

                if (gdis_try_parse_double_relaxed(species_tokens[charge_index], &parsed_charge))
                  {
                    gdouble *stored_charge;
                    g_autofree gchar *exact_key = NULL;
                    g_autofree gchar *element_symbol = NULL;
                    g_autofree gchar *element_key = NULL;

                    stored_charge = g_new(gdouble, 1);
                    *stored_charge = parsed_charge;
                    exact_key = gdis_gulp_make_species_key(species_tokens[0], type_name);
                    g_hash_table_replace(species_charges,
                                         g_strdup(exact_key),
                                         stored_charge);

                    element_symbol = gdis_normalize_element_symbol(species_tokens[0]);
                    if (element_symbol &&
                        element_symbol[0] != '\0' &&
                        g_ascii_strcasecmp(element_symbol, species_tokens[0]) != 0)
                      {
                        gdouble *element_charge;

                        element_charge = g_new(gdouble, 1);
                        *element_charge = parsed_charge;
                        element_key = gdis_gulp_make_species_key(element_symbol, type_name);
                        g_hash_table_replace(species_charges,
                                             g_strdup(element_key),
                                             element_charge);
                      }
                  }
                else
                  {
                    break;
                  }
              }

              i = next_index;
              if (species_count_limited && species_remaining > 0u)
                species_remaining--;
            }
          continue;
        }

      if (g_str_equal(first_lower, "name"))
        {
          const gchar *name_text;

          if (seen_name && model->atoms->len > 0u)
            {
              ignore_additional_geometry = TRUE;
              coord_mode = GDIS_GULP_COORD_NONE;
              current_region = -1;
              continue;
            }

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
          if (ignore_additional_geometry)
            continue;

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
                  i++;
                }
            }
          continue;
        }

      if (g_str_equal(first_lower, "cell") || g_str_equal(first_lower, "scell"))
        {
          if (ignore_additional_geometry)
            continue;

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
          if (ignore_additional_geometry)
            continue;

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
          if (ignore_additional_geometry)
            continue;

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
          if (ignore_additional_geometry)
            continue;

          guint parsed_region;

          if (token_count > 1 && gdis_try_parse_uint(tokens[1], &parsed_region))
            current_region = (gint) parsed_region;
          continue;
        }

      if (g_str_equal(first_lower, "origin"))
        {
          if (ignore_additional_geometry)
            continue;

          guint parsed_origin;

          if (token_count > 1 && gdis_try_parse_uint(tokens[1], &parsed_origin))
            origin_choice = (gint) parsed_origin;
          continue;
        }

      if (g_str_equal(first_lower, "cartesian") || g_str_equal(first_lower, "cart"))
        {
          if (ignore_additional_geometry)
            continue;

          coord_mode = GDIS_GULP_COORD_CARTESIAN;
          continue;
        }

      if (g_str_equal(first_lower, "fractional") || g_str_equal(first_lower, "frac"))
        {
          if (ignore_additional_geometry)
            continue;

          coord_mode = GDIS_GULP_COORD_FRACTIONAL;
          continue;
        }

      if (g_str_equal(first_lower, "sfractional") || g_str_equal(first_lower, "sfrac"))
        {
          if (ignore_additional_geometry)
            continue;

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
          if (ignore_additional_geometry)
            continue;

          GdisAtom *atom;
          g_autofree gchar *type = NULL;
          gdouble values[3];
          gdouble cart[3];
          gdouble charge = 0.0;
          const gdouble *mapped_charge = NULL;
          gdouble occupancy;
          gboolean have_charge = FALSE;
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

          if (!g_ascii_isalpha(tokens[0][0]))
            continue;

          if (token_count >= 5)
            {
              type = g_ascii_strdown(tokens[1], -1);
              is_core = gdis_gulp_type_is_core(type);
              is_shell = gdis_gulp_type_is_shell(type);
              has_explicit_type = is_core || is_shell;
              if (has_explicit_type)
                coord_start = 2;
            }

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

          if (coord_mode == GDIS_GULP_COORD_FRACTIONAL)
            saw_fractional_atoms = TRUE;
          else
            saw_non_fractional_atoms = TRUE;

          if (coord_start + 3 < token_count &&
              token_count == coord_start + 4 &&
              gdis_try_parse_double_relaxed(tokens[coord_start + 3], &charge))
            {
              have_charge = TRUE;
            }
          else if (coord_start + 4 < token_count &&
              gdis_try_parse_double_relaxed(tokens[coord_start + 3], &charge))
            {
              have_charge = TRUE;
            }
          else if (coord_start + 3 < token_count &&
                   gdis_try_parse_double_relaxed(tokens[coord_start + 3], &charge) &&
                   fabs(charge) > 1.05)
            {
              have_charge = TRUE;
            }

          occupancy = 1.0;
          if (coord_start + 4 < token_count)
            {
              gdouble parsed_occupancy;

              if (gdis_try_parse_double_relaxed(tokens[coord_start + 4], &parsed_occupancy) &&
                  parsed_occupancy >= 0.0 &&
                  parsed_occupancy <= 1.05)
                occupancy = parsed_occupancy;
            }
          else if (!have_charge &&
                   coord_start + 3 < token_count)
            {
              gdouble parsed_occupancy;

              if (gdis_try_parse_double_relaxed(tokens[coord_start + 3], &parsed_occupancy) &&
                  parsed_occupancy >= 0.0 &&
                  parsed_occupancy <= 1.05)
                occupancy = parsed_occupancy;
            }

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
          if (!have_charge)
            {
              mapped_charge = gdis_gulp_lookup_species_charge(species_charges,
                                                              tokens[0],
                                                              element,
                                                              is_shell ? "shell" : "core");
              if (mapped_charge)
                {
                  charge = *mapped_charge;
                  have_charge = TRUE;
                }
            }

            {
              GdisAtom *charge_atom;

              charge_atom = gdis_atom_new(tokens[0],
                                          element,
                                          has_explicit_type ? tokens[1] : "core",
                                          cart[0],
                                          cart[1],
                                          cart[2],
                                          occupancy,
                                          (coord_mode == GDIS_GULP_COORD_SURFACE_FRACTIONAL)
                                            ? current_region
                                            : -1,
                                          charge_atoms->len + 1u);
              if (have_charge)
                {
                  charge_atom->charge = charge;
                  charge_atom->has_charge = TRUE;
                }
              g_ptr_array_add(charge_atoms, charge_atom);
            }

          if (is_shell)
            continue;

          region = (coord_mode == GDIS_GULP_COORD_SURFACE_FRACTIONAL) ? current_region : -1;
          atom = gdis_atom_new(tokens[0],
                               element,
                               tokens[0],
                               cart[0],
                               cart[1],
                               cart[2],
                               occupancy,
                               region,
                               model->atoms->len + 1u);
          if (have_charge)
            {
              atom->charge = charge;
              atom->has_charge = TRUE;
            }
          g_ptr_array_add(model->atoms, atom);
        }
    }

  g_strfreev(lines);

  for (guint i = 0; i < charge_atoms->len; i++)
    {
      GdisAtom *charge_atom;
      const gdouble *resolved_charge;

      charge_atom = g_ptr_array_index(charge_atoms, i);
      if (!charge_atom || charge_atom->has_charge)
        continue;

      resolved_charge = gdis_gulp_lookup_species_charge(species_charges,
                                                        charge_atom->label,
                                                        charge_atom->element,
                                                        charge_atom->ff_type);
      if (!resolved_charge)
        continue;

      charge_atom->charge = *resolved_charge;
      charge_atom->has_charge = TRUE;
    }

  for (guint i = 0; i < model->atoms->len; i++)
    {
      GdisAtom *atom_iter;
      const gdouble *resolved_charge;

      atom_iter = g_ptr_array_index(model->atoms, i);
      if (!atom_iter || atom_iter->has_charge)
        continue;

      resolved_charge = gdis_gulp_lookup_species_charge(species_charges,
                                                        atom_iter->label,
                                                        atom_iter->element,
                                                        "core");
      if (!resolved_charge)
        continue;

      atom_iter->charge = *resolved_charge;
      atom_iter->has_charge = TRUE;
    }

  if (charge_atoms->len > 0u &&
      gdis_model_sum_known_charge_atoms(charge_atoms,
                                        model->periodic,
                                        model->periodicity,
                                        model->cell_lengths,
                                        model->cell_angles,
                                        model->space_group,
                                        origin_choice,
                                        saw_fractional_atoms &&
                                          !saw_non_fractional_atoms &&
                                          model->periodic &&
                                          model->periodicity == 3u &&
                                          model->space_group &&
                                          model->space_group[0] != '\0',
                                        &model->total_charge_e))
    model->has_total_charge = TRUE;

  g_ptr_array_free(charge_atoms, TRUE);
  g_hash_table_unref(species_charges);

  if (saw_fractional_atoms &&
      !saw_non_fractional_atoms &&
      model->periodic &&
      model->periodicity == 3u &&
      model->space_group &&
      model->space_group[0] != '\0')
    gdis_model_expand_space_group_atoms(model, origin_choice);

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
  g_autofree gchar *context_space_group = NULL;
  gdouble context_surface_a[3];
  gdouble context_surface_b[3];
  g_autofree gchar *best_space_group = NULL;
  gboolean best_from_fractional_block;
  guint best_priority;
  guint best_config_index;
  guint selected_config_index;
  guint current_config_index;
  gboolean inside_output_config;
  gboolean selected_context_valid;
  gboolean selected_context_periodic;
  guint selected_context_periodicity;
  gdouble selected_context_lengths[3];
  gdouble selected_context_angles[3];
  gboolean selected_context_have_surface_vectors;
  guint best_charge_priority;
  guint best_charge_config_index;
  gboolean best_has_total_charge;
  gboolean have_pressure;
  gdouble last_pressure_gpa;
  GArray *pressure_x_values;
  GArray *pressure_y_values;
  GArray *frequency_values;
  GArray *frequency_intensities;
  GArray *raman_values;
  gboolean selected_config_has_shell_species;
  g_autofree gchar *selected_context_space_group = NULL;
  gdouble selected_context_surface_a[3];
  gdouble selected_context_surface_b[3];

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
  best_from_fractional_block = FALSE;
  best_priority = 0u;
  best_config_index = 0u;
  selected_config_index = 0u;
  current_config_index = 0u;
  inside_output_config = FALSE;
  selected_context_valid = FALSE;
  selected_context_periodic = FALSE;
  selected_context_periodicity = 0u;
  selected_context_lengths[0] = selected_context_lengths[1] = selected_context_lengths[2] = 0.0;
  selected_context_angles[0] = selected_context_angles[1] = selected_context_angles[2] = 90.0;
  selected_context_have_surface_vectors = FALSE;
  best_charge_priority = 0u;
  best_charge_config_index = 0u;
  best_has_total_charge = FALSE;
  have_pressure = FALSE;
  last_pressure_gpa = 0.0;
  pressure_x_values = gdis_double_array_new();
  pressure_y_values = gdis_double_array_new();
  frequency_values = NULL;
  frequency_intensities = NULL;
  raman_values = NULL;
  selected_config_has_shell_species = FALSE;
  gdis_vec3_set(context_surface_a, 0.0, 0.0, 0.0);
  gdis_vec3_set(context_surface_b, 0.0, 0.0, 0.0);
  gdis_vec3_set(selected_context_surface_a, 0.0, 0.0, 0.0);
  gdis_vec3_set(selected_context_surface_b, 0.0, 0.0, 0.0);

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

      if (current_config_index == selected_config_index &&
          g_strstr_len(lower, -1, "final energy") != NULL)
        {
          gdouble values[1];

          if (gdis_collect_doubles_from_line(trimmed, values, 1u) >= 1u)
            {
              model->energy_ev = values[0];
              model->has_energy = TRUE;
            }
          continue;
        }

      if (current_config_index == selected_config_index &&
          !model->has_energy &&
          g_strstr_len(lower, -1, "total lattice energy") != NULL &&
          g_strstr_len(lower, -1, "ev") != NULL)
        {
          gdouble values[1];

          if (gdis_collect_doubles_from_line(trimmed, values, 1u) >= 1u)
            {
              model->energy_ev = values[0];
              model->has_energy = TRUE;
            }
          continue;
        }

      if (current_config_index == selected_config_index &&
          g_strstr_len(lower, -1, "density of cell") != NULL)
        {
          gdouble values[1];

          if (gdis_collect_doubles_from_line(trimmed, values, 1u) >= 1u)
            {
              model->density_g_cm3 = values[0];
              model->has_density = TRUE;
            }
          continue;
        }

      if (current_config_index == selected_config_index &&
          g_strstr_len(lower, -1, "pressure of configuration") != NULL &&
          g_strstr_len(lower, -1, "gpa") != NULL)
        {
          gdouble values[1];

          if (gdis_collect_doubles_from_text_scan(trimmed, values, 1u) >= 1u)
            {
              const gdouble x_value = (gdouble) pressure_y_values->len + 1.0;

              g_array_append_val(pressure_x_values, x_value);
              g_array_append_val(pressure_y_values, values[0]);
              last_pressure_gpa = values[0];
              have_pressure = TRUE;
            }
          continue;
        }

      if (current_config_index == selected_config_index &&
          g_strstr_len(lower, -1, "frequencies (cm-1)") != NULL)
        {
          GArray *parsed_frequencies;
          guint scan_index;

          parsed_frequencies = gdis_double_array_new();
          scan_index = i + 1u;
          for (; lines[scan_index] != NULL; scan_index++)
            {
              gchar *freq_line;
              g_autofree gchar *freq_lower = NULL;
              gdouble values[64];
              guint count;

              freq_line = g_strstrip(lines[scan_index]);
              if (freq_line[0] == '\0')
                {
                  if (parsed_frequencies->len > 0u)
                    break;
                  continue;
                }

              freq_lower = g_ascii_strdown(freq_line, -1);
              if (freq_line[0] == '-')
                {
                  if (parsed_frequencies->len > 0u)
                    break;
                  continue;
                }
              if (g_strstr_len(freq_lower, -1, "phonon properties") != NULL ||
                  g_strstr_len(freq_lower, -1, "phonon density of states") != NULL)
                {
                  if (parsed_frequencies->len > 0u)
                    {
                      scan_index--;
                      break;
                    }
                  continue;
                }

              count = gdis_collect_doubles_from_text_scan(freq_line,
                                                          values,
                                                          G_N_ELEMENTS(values));
              if (count == 0u)
                {
                  if (parsed_frequencies->len > 0u)
                    {
                      scan_index--;
                      break;
                    }
                  continue;
                }
              g_array_append_vals(parsed_frequencies, values, count);
            }

          if (parsed_frequencies->len > 0u)
            {
              g_clear_pointer(&frequency_values, g_array_unref);
              frequency_values = parsed_frequencies;
            }
          else
            {
              g_array_unref(parsed_frequencies);
            }

          if (lines[scan_index] != NULL)
            i = scan_index;
          else
            i = scan_index - 1u;
          continue;
        }

      if (current_config_index == selected_config_index &&
          g_str_has_prefix(lower, "raman"))
        {
          gdouble values[64];
          guint count;

          count = gdis_collect_doubles_from_text_scan(trimmed,
                                                      values,
                                                      G_N_ELEMENTS(values));
          if (count > 0u)
            {
              if (!raman_values)
                raman_values = gdis_double_array_new();
              g_array_append_vals(raman_values, values, count);
            }
          continue;
        }

      if (g_strstr_len(lower, -1, "input for configuration") != NULL)
        {
          gdouble config_value[1];
          const gchar *title_text;

          if (gdis_collect_doubles_from_line(trimmed, config_value, 1u) >= 1u)
            current_config_index = (guint) MAX(0.0, llround(config_value[0]));
          else
            current_config_index = 0u;

          if (selected_config_index == 0u)
            selected_config_index = current_config_index;

          inside_output_config = FALSE;
          context_periodic = FALSE;
          context_periodicity = 0u;
          context_lengths[0] = context_lengths[1] = context_lengths[2] = 0.0;
          context_angles[0] = context_angles[1] = context_angles[2] = 90.0;
          context_have_surface_vectors = FALSE;
          g_clear_pointer(&context_space_group, g_free);
          gdis_vec3_set(context_surface_a, 0.0, 0.0, 0.0);
          gdis_vec3_set(context_surface_b, 0.0, 0.0, 0.0);

          if (current_config_index == selected_config_index)
            {
              title_text = strchr(trimmed, ':');
              if (title_text)
                {
                  g_autofree gchar *parsed_title = NULL;

                  parsed_title = gdis_strdup_strip(title_text + 1);
                  while (parsed_title && parsed_title[0] != '\0')
                    {
                      gsize len;

                      len = strlen(parsed_title);
                      if (len == 0u || parsed_title[len - 1] != '*')
                        break;
                      parsed_title[len - 1] = '\0';
                      g_strstrip(parsed_title);
                    }

                  if (parsed_title && parsed_title[0] != '\0')
                    {
                      g_clear_pointer(&model->title, g_free);
                      model->title = g_strdup(parsed_title);
                    }
                }
            }
          continue;
        }

      if (g_strstr_len(lower, -1, "output for configuration") != NULL)
        {
          gdouble config_value[1];
          const gchar *title_text;

          if (gdis_collect_doubles_from_line(trimmed, config_value, 1u) >= 1u)
            current_config_index = (guint) MAX(0.0, llround(config_value[0]));
          else
            current_config_index = 0u;

          if (selected_config_index == 0u)
            selected_config_index = current_config_index;

          inside_output_config = TRUE;

          if (current_config_index == selected_config_index && selected_context_valid)
            {
              context_periodic = selected_context_periodic;
              context_periodicity = selected_context_periodicity;
              memcpy(context_lengths, selected_context_lengths, sizeof(context_lengths));
              memcpy(context_angles, selected_context_angles, sizeof(context_angles));
              context_have_surface_vectors = selected_context_have_surface_vectors;
              g_clear_pointer(&context_space_group, g_free);
              context_space_group = selected_context_space_group
                                      ? g_strdup(selected_context_space_group)
                                      : NULL;
              gdis_vec3_copy(context_surface_a, selected_context_surface_a);
              gdis_vec3_copy(context_surface_b, selected_context_surface_b);
            }

          if (current_config_index == selected_config_index)
            {
              title_text = strchr(trimmed, ':');
              if (title_text)
                {
                  g_autofree gchar *parsed_title = NULL;

                  parsed_title = gdis_strdup_strip(title_text + 1);
                  while (parsed_title && parsed_title[0] != '\0')
                    {
                      gsize len;

                      len = strlen(parsed_title);
                      if (len == 0u || parsed_title[len - 1] != '*')
                        break;
                      parsed_title[len - 1] = '\0';
                      g_strstrip(parsed_title);
                    }

                  if (parsed_title && parsed_title[0] != '\0')
                    {
                      g_clear_pointer(&model->title, g_free);
                      model->title = g_strdup(parsed_title);
                    }
                }
            }
          continue;
        }

      if (g_str_has_prefix(lower, "space group"))
        {
          const gchar *value_text;

          value_text = strchr(trimmed, ':');
          if (value_text)
            {
              g_autofree gchar *parsed_group = NULL;

              value_text++;
              while (*value_text && g_ascii_isspace(*value_text))
                value_text++;

              parsed_group = gdis_strdup_strip(value_text);
              if (parsed_group && parsed_group[0] != '\0')
                {
                  g_clear_pointer(&context_space_group, g_free);
                  context_space_group = g_steal_pointer(&parsed_group);
                }
            }

          if (!inside_output_config &&
              selected_config_index != 0u &&
              current_config_index == selected_config_index)
            {
              selected_context_valid = TRUE;
              selected_context_periodic = context_periodic;
              selected_context_periodicity = context_periodicity;
              memcpy(selected_context_lengths, context_lengths, sizeof(context_lengths));
              memcpy(selected_context_angles, context_angles, sizeof(context_angles));
              selected_context_have_surface_vectors = context_have_surface_vectors;
              g_clear_pointer(&selected_context_space_group, g_free);
              selected_context_space_group = context_space_group
                                               ? g_strdup(context_space_group)
                                               : NULL;
              gdis_vec3_copy(selected_context_surface_a, context_surface_a);
              gdis_vec3_copy(selected_context_surface_b, context_surface_b);
            }
          continue;
        }

      if (g_str_has_prefix(lower, "non-standard setting of group"))
        {
          const gchar *value_text;

          value_text = strchr(trimmed, ':');
          if (value_text)
            {
              g_autofree gchar *parsed_group = NULL;

              value_text++;
              while (*value_text && g_ascii_isspace(*value_text))
                value_text++;

              parsed_group = gdis_strdup_strip(value_text);
              if (parsed_group && parsed_group[0] != '\0')
                {
                  g_clear_pointer(&context_space_group, g_free);
                  context_space_group = g_steal_pointer(&parsed_group);
                }
            }

          if (!inside_output_config &&
              selected_config_index != 0u &&
              current_config_index == selected_config_index)
            {
              selected_context_valid = TRUE;
              selected_context_periodic = context_periodic;
              selected_context_periodicity = context_periodicity;
              memcpy(selected_context_lengths, context_lengths, sizeof(context_lengths));
              memcpy(selected_context_angles, context_angles, sizeof(context_angles));
              selected_context_have_surface_vectors = context_have_surface_vectors;
              g_clear_pointer(&selected_context_space_group, g_free);
              selected_context_space_group = context_space_group
                                               ? g_strdup(context_space_group)
                                               : NULL;
              gdis_vec3_copy(selected_context_surface_a, context_surface_a);
              gdis_vec3_copy(selected_context_surface_b, context_surface_b);
            }
          continue;
        }

      if (g_strstr_len(lower, -1, "species output for") != NULL)
        {
          gboolean parsed_any;
          guint scan_index;

          parsed_any = FALSE;
          scan_index = i + 1u;
          while (lines[scan_index] != NULL)
            {
              gchar *species_line;
              gchar *species_trimmed;
              g_auto(GStrv) species_tokens = NULL;
              gint species_token_count = 0;
              gdouble radius;
              g_autofree gchar *symbol = NULL;

              species_line = lines[scan_index];
              g_strchomp(species_line);
              species_trimmed = g_strstrip(species_line);
              if (species_trimmed[0] == '\0')
                {
                  scan_index++;
                  continue;
                }
              if (species_trimmed[0] == '-' || species_trimmed[0] == '*')
                {
                  scan_index++;
                  continue;
                }

              species_tokens = gdis_split_simple(species_trimmed, &species_token_count);
              if (species_token_count >= 6 &&
                  g_ascii_isalpha(species_tokens[0][0]) &&
                  gdis_try_parse_double_relaxed(species_tokens[5], &radius))
                {
                  if (current_config_index == selected_config_index &&
                      species_token_count > 1 &&
                      g_ascii_strncasecmp(species_tokens[1], "shell", 5) == 0)
                    selected_config_has_shell_species = TRUE;
                  symbol = gdis_normalize_element_symbol(species_tokens[0]);
                  if (symbol && symbol[0] != '\0')
                    gdis_model_set_element_covalent_override(model, symbol, radius);
                  parsed_any = TRUE;
                  scan_index++;
                  continue;
                }

              if (parsed_any)
                break;
              scan_index++;
            }

          if (scan_index > i + 1u)
            i = scan_index - 1u;
          continue;
        }

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

          if (!inside_output_config &&
              selected_config_index != 0u &&
              current_config_index == selected_config_index)
            {
              selected_context_valid = TRUE;
              selected_context_periodic = context_periodic;
              selected_context_periodicity = context_periodicity;
              memcpy(selected_context_lengths, context_lengths, sizeof(context_lengths));
              memcpy(selected_context_angles, context_angles, sizeof(context_angles));
              selected_context_have_surface_vectors = context_have_surface_vectors;
              g_clear_pointer(&selected_context_space_group, g_free);
              selected_context_space_group = context_space_group
                                               ? g_strdup(context_space_group)
                                               : NULL;
              gdis_vec3_copy(selected_context_surface_a, context_surface_a);
              gdis_vec3_copy(selected_context_surface_b, context_surface_b);
            }
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

          if (!inside_output_config &&
              selected_config_index != 0u &&
              current_config_index == selected_config_index)
            {
              selected_context_valid = TRUE;
              selected_context_periodic = context_periodic;
              selected_context_periodicity = context_periodicity;
              memcpy(selected_context_lengths, context_lengths, sizeof(context_lengths));
              memcpy(selected_context_angles, context_angles, sizeof(context_angles));
              selected_context_have_surface_vectors = context_have_surface_vectors;
              g_clear_pointer(&selected_context_space_group, g_free);
              selected_context_space_group = context_space_group
                                               ? g_strdup(context_space_group)
                                               : NULL;
              gdis_vec3_copy(selected_context_surface_a, context_surface_a);
              gdis_vec3_copy(selected_context_surface_b, context_surface_b);
            }
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
          GPtrArray *block_charge_atoms;
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
          block_charge_atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
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
              gdouble charge = 0.0;
              gboolean have_charge = FALSE;
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

              if (mode == GDIS_GULP_BLOCK_FRACTIONAL)
                {
                  gdouble numeric_values[8];
                  guint numeric_count;

                  numeric_count = 0u;
                  for (gint token_index = 2; token_index < token_count && numeric_count < G_N_ELEMENTS(numeric_values); token_index++)
                    {
                      gdouble parsed_value;

                      if (!gdis_try_parse_double_relaxed(tokens[token_index], &parsed_value))
                        continue;
                      numeric_values[numeric_count++] = parsed_value;
                    }
                  if (numeric_count >= 5u)
                    {
                      charge = numeric_values[3];
                      have_charge = TRUE;
                    }
                }

              element = gdis_normalize_element_symbol(tokens[1]);
              if (have_charge)
                {
                  GdisAtom *charge_atom;

                  charge_atom = gdis_atom_new(tokens[1],
                                              element,
                                              type ? type : element,
                                              cart[0],
                                              cart[1],
                                              cart[2],
                                              1.0,
                                              (mode == GDIS_GULP_BLOCK_SURFACE_MIXED)
                                                ? block_region
                                                : -1,
                                              serial);
                  charge_atom->charge = charge;
                  charge_atom->has_charge = TRUE;
                  g_ptr_array_add(block_charge_atoms, charge_atom);
                }

              if (g_str_has_prefix(type, "shel") || g_str_equal(type, "s"))
                {
                  block_had_rows = TRUE;
                  continue;
                }

              {
                GdisAtom *atom;

                atom = gdis_atom_new(tokens[1],
                                     element,
                                     element,
                                     cart[0],
                                     cart[1],
                                     cart[2],
                                     1.0,
                                     (mode == GDIS_GULP_BLOCK_SURFACE_MIXED)
                                       ? block_region
                                       : -1,
                                     serial);
                if (have_charge)
                  {
                    atom->charge = charge;
                    atom->has_charge = TRUE;
                  }
                g_ptr_array_add(block_atoms, atom);
              }
              block_had_rows = TRUE;
            }

          if (block_atoms->len > 0u)
            {
              const gboolean allow_block = (selected_config_index == 0u ||
                                            current_config_index == 0u ||
                                            current_config_index == selected_config_index);
              const guint block_priority = inside_output_config ? 2u : 1u;
              gdouble block_total_charge;

              if (allow_block &&
                  (!best_atoms ||
                   block_priority > best_priority ||
                   (block_priority == best_priority &&
                    current_config_index == best_config_index)))
                {
                  if (best_atoms)
                    g_ptr_array_free(best_atoms, TRUE);
                  best_atoms = block_atoms;
                  g_clear_pointer(&best_space_group, g_free);
                  best_space_group = context_space_group ? g_strdup(context_space_group) : NULL;
                  best_from_fractional_block = (mode == GDIS_GULP_BLOCK_FRACTIONAL);
                  best_periodic = block_periodic;
                  best_periodicity = block_periodicity;
                  best_priority = block_priority;
                  best_config_index = current_config_index;
                  memcpy(best_lengths, block_lengths, sizeof(best_lengths));
                  memcpy(best_angles, block_angles, sizeof(best_angles));
                }
              else
                {
                  g_ptr_array_free(block_atoms, TRUE);
                }

              if (allow_block &&
                  block_charge_atoms->len > 0u &&
                  gdis_model_sum_charge_atoms(block_charge_atoms,
                                              block_periodic,
                                              block_periodicity,
                                              block_lengths,
                                              block_angles,
                                              context_space_group,
                                              0,
                                              mode == GDIS_GULP_BLOCK_FRACTIONAL,
                                              &block_total_charge) &&
                  (!selected_config_has_shell_species ||
                   fabs(block_total_charge) <= 1.0e-3) &&
                  (!best_has_total_charge ||
                   block_priority > best_charge_priority ||
                   (block_priority == best_charge_priority &&
                    current_config_index == best_charge_config_index)))
                {
                  model->has_total_charge = TRUE;
                  model->total_charge_e = block_total_charge;
                  best_has_total_charge = TRUE;
                  best_charge_priority = block_priority;
                  best_charge_config_index = current_config_index;
                }
            }
          else
            {
              g_ptr_array_free(block_atoms, TRUE);
            }
          g_ptr_array_free(block_charge_atoms, TRUE);

          continue;
        }
    }

  g_strfreev(lines);
  gdis_model_clear_plot_arrays(model);

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
      g_clear_pointer(&pressure_x_values, g_array_unref);
      g_clear_pointer(&pressure_y_values, g_array_unref);
      g_clear_pointer(&frequency_values, g_array_unref);
      g_clear_pointer(&frequency_intensities, g_array_unref);
      g_clear_pointer(&raman_values, g_array_unref);
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
  model->has_pressure = have_pressure;
  model->pressure_gpa = last_pressure_gpa;
  if (best_space_group && best_space_group[0] != '\0')
    {
      g_clear_pointer(&model->space_group, g_free);
      model->space_group = g_strdup(best_space_group);
    }

  if (have_pressure && pressure_y_values->len > 0u)
    {
      model->pressure_x_values = pressure_x_values;
      model->pressure_y_values_gpa = pressure_y_values;
      pressure_x_values = NULL;
      pressure_y_values = NULL;
    }

  if (frequency_values && frequency_values->len > 0u)
    {
      frequency_intensities = gdis_double_array_new();
      for (guint i = 0; i < frequency_values->len; i++)
        {
          const gdouble intensity = 1.0;

          g_array_append_val(frequency_intensities, intensity);
        }
      gdis_model_build_stick_spectrum_arrays(frequency_values,
                                             frequency_intensities,
                                             &model->frequency_x_values_cm1,
                                             &model->frequency_y_values);
    }

  if (raman_values && raman_values->len > 0u)
    {
      GArray *raman_positions;

      if (frequency_values && frequency_values->len > 0u)
        raman_positions = gdis_double_array_clone(frequency_values);
      else
        {
          raman_positions = gdis_double_array_new();
          for (guint i = 0; i < raman_values->len; i++)
            {
              const gdouble mode_index = (gdouble) i + 1.0;

              g_array_append_val(raman_positions, mode_index);
            }
        }

      gdis_model_build_stick_spectrum_arrays(raman_positions,
                                             raman_values,
                                             &model->raman_x_values_cm1,
                                             &model->raman_y_values);
      g_array_unref(raman_positions);
    }

  if (best_from_fractional_block &&
      model->periodic &&
      model->periodicity == 3u &&
      model->space_group &&
      model->space_group[0] != '\0')
    gdis_model_expand_space_group_atoms(model, 0);

  g_clear_pointer(&pressure_x_values, g_array_unref);
  g_clear_pointer(&pressure_y_values, g_array_unref);
  g_clear_pointer(&frequency_values, g_array_unref);
  g_clear_pointer(&frequency_intensities, g_array_unref);
  g_clear_pointer(&raman_values, g_array_unref);

  return TRUE;
}

static gchar *
gdis_xtl_compose_space_group(const gchar *label,
                             guint number,
                             const gchar *qualifier)
{
  g_autofree gchar *base = NULL;
  g_autofree gchar *qualifier_lower = NULL;

  if (label && label[0] != '\0')
    base = gdis_strdup_strip(label);
  else if (number > 0u)
    base = g_strdup_printf("%u", number);

  if (!base || base[0] == '\0')
    return NULL;

  if (!qualifier || qualifier[0] == '\0' || strchr(base, ':') != NULL)
    return g_steal_pointer(&base);

  qualifier_lower = g_ascii_strdown(qualifier, -1);
  g_strstrip(qualifier_lower);

  if (g_str_has_prefix(qualifier_lower, "origin_"))
    {
      const gchar *origin_text;

      origin_text = qualifier_lower + strlen("origin_");
      if (origin_text[0] != '\0')
        return g_strdup_printf("%s:%s", base, origin_text);
    }

  if (g_str_equal(qualifier_lower, "hexagonal"))
    return g_strdup_printf("%s:H", base);
  if (g_str_equal(qualifier_lower, "rhombohedral"))
    return g_strdup_printf("%s:R", base);
  if (g_str_equal(qualifier_lower, "a_unique"))
    return g_strdup_printf("%s:a", base);
  if (g_str_equal(qualifier_lower, "b_unique"))
    return g_strdup_printf("%s:b", base);
  if (g_str_equal(qualifier_lower, "c_unique"))
    return g_strdup_printf("%s:c", base);

  return g_steal_pointer(&base);
}

static gboolean
gdis_model_load_xtl(GdisModel *model, const gchar *contents, GError **error)
{
  gchar **lines;
  gboolean in_atoms;
  gboolean have_cell_matrix;
  gboolean atoms_are_fractional;
  gdouble cell_matrix[9];
  gdouble cell_inverse[9];
  guint symmetry_number;
  g_autofree gchar *symmetry_label = NULL;
  g_autofree gchar *symmetry_qualifier = NULL;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(contents != NULL, FALSE);

  lines = g_strsplit(contents, "\n", -1);
  in_atoms = FALSE;
  have_cell_matrix = FALSE;
  atoms_are_fractional = TRUE;
  symmetry_number = 0u;
  model->periodic = TRUE;
  model->periodicity = 3u;

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

      if (g_str_equal(first_lower, "dimension"))
        {
          guint parsed_dimension;

          if (token_count > 1 && gdis_try_parse_uint(tokens[1], &parsed_dimension))
            {
              if (parsed_dimension >= 1u && parsed_dimension <= 3u)
                {
                  model->periodic = TRUE;
                  model->periodicity = parsed_dimension;
                }
            }
          continue;
        }

      if (g_str_equal(first_lower, "cell"))
        {
          gdouble values[6] = {0.0, 0.0, 1.0, 90.0, 90.0, 90.0};
          guint required;
          guint found;
          guint line_index;

          required = (model->periodicity == 2u) ? 3u : 6u;
          found = gdis_collect_doubles_from_line(trimmed, values, required);
          line_index = i;
          while (found < required && lines[line_index + 1] != NULL)
            {
              guint added;

              line_index++;
              added = gdis_collect_doubles_from_line(lines[line_index],
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
              if (model->periodicity == 2u)
                {
                  model->cell_lengths[0] = values[0];
                  model->cell_lengths[1] = values[1];
                  model->cell_lengths[2] = 1.0;
                  model->cell_angles[0] = 90.0;
                  model->cell_angles[1] = 90.0;
                  model->cell_angles[2] = values[2];
                }
              else
                {
                  memcpy(model->cell_lengths, values, 3u * sizeof(gdouble));
                  memcpy(model->cell_angles, values + 3, 3u * sizeof(gdouble));
                  model->periodic = TRUE;
                  model->periodicity = 3u;
                }
              have_cell_matrix = gdis_model_build_cell_matrix(model,
                                                              cell_matrix,
                                                              cell_inverse);
            }
          continue;
        }

      if (g_str_equal(first_lower, "symmetry"))
        {
          for (gint j = 1; j < token_count; )
            {
              g_autofree gchar *keyword = g_ascii_strdown(tokens[j], -1);

              if (g_str_equal(keyword, "number") && j + 1 < token_count)
                {
                  guint parsed_number;

                  if (gdis_try_parse_uint(tokens[j + 1], &parsed_number))
                    symmetry_number = parsed_number;
                  j += 2;
                  continue;
                }

              if ((g_str_equal(keyword, "label") ||
                   g_str_equal(keyword, "qualifier")) &&
                  j + 1 < token_count)
                {
                  GString *value;

                  value = g_string_new("");
                  for (j = j + 1; j < token_count; j++)
                    {
                      g_autofree gchar *next_keyword = NULL;

                      next_keyword = g_ascii_strdown(tokens[j], -1);
                      if (g_str_equal(next_keyword, "number") ||
                          g_str_equal(next_keyword, "label") ||
                          g_str_equal(next_keyword, "qualifier"))
                        break;

                      if (value->len > 0)
                        g_string_append_c(value, ' ');
                      g_string_append(value, tokens[j]);
                    }

                  if (g_str_equal(keyword, "label"))
                    {
                      g_clear_pointer(&symmetry_label, g_free);
                      symmetry_label = g_string_free(value, FALSE);
                    }
                  else
                    {
                      g_clear_pointer(&symmetry_qualifier, g_free);
                      symmetry_qualifier = g_string_free(value, FALSE);
                    }
                  continue;
                }

              j++;
            }

          g_clear_pointer(&model->space_group, g_free);
          model->space_group = gdis_xtl_compose_space_group(symmetry_label,
                                                            symmetry_number,
                                                            symmetry_qualifier);
          continue;
        }

      if (g_str_equal(first_lower, "atoms"))
        {
          in_atoms = TRUE;
          atoms_are_fractional = TRUE;
          if (lines[i + 1] != NULL)
            {
              g_auto(GStrv) header_tokens = NULL;
              gint header_count = 0;

              header_tokens = gdis_split_simple(lines[i + 1], &header_count);
              for (gint h = 0; h < header_count; h++)
                {
                  g_autofree gchar *header_lower = g_ascii_strdown(header_tokens[h], -1);

                  if (g_str_equal(header_lower, "carx") ||
                      g_str_equal(header_lower, "cary") ||
                      g_str_equal(header_lower, "carz"))
                    {
                      atoms_are_fractional = FALSE;
                      break;
                    }
                }
            }
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

          if (atoms_are_fractional && have_cell_matrix)
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

  if (have_cell_matrix &&
      model->periodic &&
      model->periodicity == 3u &&
      model->space_group &&
      model->space_group[0] != '\0')
    gdis_model_expand_space_group_atoms(model, 0);

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
                               gboolean default_bohr_when_unspecified,
                               gdouble *scale_out,
                               gboolean *fractional_out)
{
  g_autofree gchar *lower = NULL;
  gdouble scale = default_bohr_when_unspecified ? 0.529177210903 : 1.0;
  gboolean fractional = FALSE;
  gdouble parsed[1];

  g_return_val_if_fail(scale_out != NULL, FALSE);
  g_return_val_if_fail(fractional_out != NULL, FALSE);

  if (!header)
    {
      *scale_out = scale;
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
gdis_model_load_qe_like(GdisModel *model,
                        const gchar *contents,
                        gboolean default_bohr_when_unspecified,
                        GError **error)
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
  GArray *pressure_x_values;
  GArray *pressure_y_values;
  GPtrArray *current_band_rows;
  GPtrArray *best_band_rows;
  GArray *current_band_kpoint_triplets;
  GArray *best_band_kpoint_triplets;
  GArray *band_x_values;
  GArray *band_y_values;
  guint expected_kpoint_count;
  guint pressure_step;
  guint band_path_count;
  guint band_series_count;
  gboolean have_pressure;
  gdouble last_pressure_gpa;
  gboolean have_fermi_energy;
  gdouble fermi_energy_ev;
  gboolean have_energy;
  gdouble last_energy_ev;
  const gdouble ry_to_ev = 13.605693122994;

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
  pressure_x_values = gdis_double_array_new();
  pressure_y_values = gdis_double_array_new();
  current_band_rows = g_ptr_array_new_with_free_func((GDestroyNotify) g_array_unref);
  best_band_rows = NULL;
  current_band_kpoint_triplets = gdis_double_array_new();
  best_band_kpoint_triplets = NULL;
  band_x_values = NULL;
  band_y_values = NULL;
  expected_kpoint_count = 0u;
  pressure_step = 0u;
  band_path_count = 0u;
  band_series_count = 0u;
  have_pressure = FALSE;
  last_pressure_gpa = 0.0;
  have_fermi_energy = FALSE;
  fermi_energy_ev = 0.0;
  have_energy = FALSE;
  last_energy_ev = 0.0;

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

      if (g_strstr_len(lower, -1, "number of k points") != NULL && strchr(trimmed, '=') != NULL)
        {
          gdouble values[1];

          if (gdis_collect_doubles_from_line(trimmed, values, 1u) >= 1u &&
              values[0] > 0.0)
            expected_kpoint_count = (guint) llround(values[0]);
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

      if (g_strstr_len(lower, -1, "total energy") != NULL &&
          g_strstr_len(lower, -1, "ry") != NULL)
        {
          gdouble values[1];

          if (gdis_collect_doubles_from_text_scan(trimmed, values, 1u) >= 1u)
            {
              last_energy_ev = values[0] * ry_to_ev;
              have_energy = TRUE;
            }
          continue;
        }

      if (g_strstr_len(lower, -1, "fermi energy") != NULL ||
          g_strstr_len(lower, -1, "highest occupied level") != NULL)
        {
          gdouble values[2];

          if (gdis_collect_doubles_from_text_scan(trimmed, values, 2u) >= 1u)
            {
              fermi_energy_ev = values[0];
              have_fermi_energy = TRUE;
            }
          continue;
        }

      if (g_strstr_len(lower, -1, "total   stress") != NULL)
        {
          const gchar *pressure_text;
          gdouble values[1];

          pressure_text = strstr(trimmed, "P=");
          if (!pressure_text)
            pressure_text = strstr(trimmed, "p=");
          if (pressure_text &&
              gdis_collect_doubles_from_text_scan(pressure_text, values, 1u) >= 1u)
            {
              const gdouble pressure_gpa = values[0] * 0.1;
              const gdouble x_value = (gdouble) (++pressure_step);

              g_array_append_val(pressure_x_values, x_value);
              g_array_append_val(pressure_y_values, pressure_gpa);
              last_pressure_gpa = pressure_gpa;
              have_pressure = TRUE;
            }
          continue;
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
                                              default_bohr_when_unspecified,
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
                                              default_bohr_when_unspecified,
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

      if (g_strstr_len(lower, -1, "bands (ev):") != NULL &&
          strchr(trimmed, '=') != NULL)
        {
          gdouble header_values[16];
          guint header_count;
          GArray *band_values;
          guint scan_index;

          header_count = gdis_collect_doubles_from_text_scan(trimmed,
                                                             header_values,
                                                             G_N_ELEMENTS(header_values));
          if (header_count < 3u)
            continue;

          if (expected_kpoint_count > 0u &&
              current_band_rows->len >= expected_kpoint_count)
            {
              if (best_band_rows)
                g_ptr_array_free(best_band_rows, TRUE);
              best_band_rows = gdis_double_array_ptr_array_clone(current_band_rows);
              g_clear_pointer(&best_band_kpoint_triplets, g_array_unref);
              best_band_kpoint_triplets = gdis_double_array_clone(current_band_kpoint_triplets);
              g_ptr_array_set_size(current_band_rows, 0u);
              g_array_set_size(current_band_kpoint_triplets, 0u);
            }

          band_values = gdis_double_array_new();
          scan_index = i + 1u;
          for (; lines[scan_index] != NULL; scan_index++)
            {
              gchar *row_line;
              g_autofree gchar *row_lower = NULL;
              gdouble row_values[64];
              guint row_count;

              row_line = g_strstrip(lines[scan_index]);
              if (row_line[0] == '\0')
                {
                  if (band_values->len > 0u)
                    break;
                  continue;
                }

              row_lower = g_ascii_strdown(row_line, -1);
              if (g_strstr_len(row_lower, -1, "bands (ev):") != NULL ||
                  g_strstr_len(row_lower, -1, "highest occupied level") != NULL ||
                  g_strstr_len(row_lower, -1, "fermi energy") != NULL ||
                  g_strstr_len(row_lower, -1, "total   stress") != NULL ||
                  g_str_has_prefix(row_lower, "atomic_positions") ||
                  g_str_has_prefix(row_lower, "cell_parameters") ||
                  g_str_has_prefix(row_lower, "end") ||
                  g_str_has_prefix(row_lower, "begin") ||
                  g_str_has_prefix(row_lower, "&"))
                {
                  if (band_values->len > 0u)
                    {
                      scan_index--;
                      break;
                    }
                  continue;
                }

              row_count = gdis_collect_doubles_from_text_scan(row_line,
                                                              row_values,
                                                              G_N_ELEMENTS(row_values));
              if (row_count == 0u)
                {
                  if (band_values->len > 0u)
                    {
                      scan_index--;
                      break;
                    }
                  continue;
                }

              g_array_append_vals(band_values, row_values, row_count);
            }

          if (band_values->len > 0u)
            {
              const gdouble coords[3] = {
                header_values[0],
                header_values[1],
                header_values[2]
              };

              g_array_append_vals(current_band_kpoint_triplets, coords, 3u);
              g_ptr_array_add(current_band_rows, band_values);
            }
          else
            {
              g_array_unref(band_values);
            }

          if (expected_kpoint_count > 0u &&
              current_band_rows->len == expected_kpoint_count)
            {
              if (best_band_rows)
                g_ptr_array_free(best_band_rows, TRUE);
              best_band_rows = gdis_double_array_ptr_array_clone(current_band_rows);
              g_clear_pointer(&best_band_kpoint_triplets, g_array_unref);
              best_band_kpoint_triplets = gdis_double_array_clone(current_band_kpoint_triplets);
              g_ptr_array_set_size(current_band_rows, 0u);
              g_array_set_size(current_band_kpoint_triplets, 0u);
            }

          if (lines[scan_index] != NULL)
            i = scan_index;
          else
            i = scan_index - 1u;
          continue;
        }
    }

  g_strfreev(lines);
  if (current_band_rows->len > 0u &&
      (!best_band_rows || current_band_rows->len >= best_band_rows->len))
    {
      if (best_band_rows)
        g_ptr_array_free(best_band_rows, TRUE);
      best_band_rows = gdis_double_array_ptr_array_clone(current_band_rows);
      g_clear_pointer(&best_band_kpoint_triplets, g_array_unref);
      best_band_kpoint_triplets = gdis_double_array_clone(current_band_kpoint_triplets);
    }
  g_ptr_array_free(current_band_rows, TRUE);
  g_array_unref(current_band_kpoint_triplets);
  gdis_model_clear_plot_arrays(model);

  if (!best_atoms || best_atoms->len == 0u)
    {
      if (best_atoms)
        g_ptr_array_free(best_atoms, TRUE);
      g_clear_pointer(&pressure_x_values, g_array_unref);
      g_clear_pointer(&pressure_y_values, g_array_unref);
      g_clear_pointer(&best_band_kpoint_triplets, g_array_unref);
      g_clear_pointer(&band_x_values, g_array_unref);
      g_clear_pointer(&band_y_values, g_array_unref);
      if (best_band_rows)
        g_ptr_array_free(best_band_rows, TRUE);
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
  model->has_energy = have_energy;
  model->energy_ev = last_energy_ev;
  model->has_pressure = have_pressure;
  model->pressure_gpa = last_pressure_gpa;
  model->has_fermi_energy = have_fermi_energy;
  model->fermi_energy_ev = fermi_energy_ev;

  if (pressure_y_values->len > 0u)
    {
      model->pressure_x_values = pressure_x_values;
      model->pressure_y_values_gpa = pressure_y_values;
    }
  else
    {
      g_array_unref(pressure_x_values);
      g_array_unref(pressure_y_values);
    }

  if (best_band_rows && best_band_rows->len > 0u)
    {
      gdis_model_build_band_arrays_from_triplets(best_band_kpoint_triplets,
                                                 best_band_rows,
                                                 have_fermi_energy ? fermi_energy_ev : 0.0,
                                                 &band_x_values,
                                                 &band_y_values,
                                                 &band_path_count,
                                                 &band_series_count);
    }
  model->band_x_values = band_x_values;
  model->band_y_values_ev = band_y_values;
  model->band_path_count = band_path_count;
  model->band_series_count = band_series_count;

  g_clear_pointer(&best_band_kpoint_triplets, g_array_unref);
  if (best_band_rows)
    g_ptr_array_free(best_band_rows, TRUE);

  return TRUE;
}

static gboolean
gdis_model_load_vasp_xml(GdisModel *model, const gchar *contents, GError **error)
{
  gchar **lines;
  GPtrArray *atom_symbols;
  GPtrArray *best_atoms;
  gdouble best_lengths[3];
  gdouble best_angles[3];
  GArray *kpoint_triplets;
  GArray *pressure_x_values;
  GArray *pressure_y_values;
  GArray *dos_x_values;
  GArray *dos_y_values;
  GArray *band_x_values;
  GArray *band_y_values;
  guint band_path_count;
  guint band_series_count;
  guint pressure_step;
  gdouble last_pressure_gpa;
  gboolean have_pressure;
  gboolean have_energy;
  gdouble last_energy_ev;
  gboolean have_fermi_energy;
  gdouble fermi_energy_ev;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(contents != NULL, FALSE);

  lines = g_strsplit(contents, "\n", -1);
  atom_symbols = g_ptr_array_new_with_free_func(g_free);
  best_atoms = NULL;
  best_lengths[0] = best_lengths[1] = best_lengths[2] = 0.0;
  best_angles[0] = best_angles[1] = best_angles[2] = 90.0;
  kpoint_triplets = gdis_double_array_new();
  pressure_x_values = gdis_double_array_new();
  pressure_y_values = gdis_double_array_new();
  dos_x_values = NULL;
  dos_y_values = NULL;
  band_x_values = NULL;
  band_y_values = NULL;
  band_path_count = 0u;
  band_series_count = 0u;
  pressure_step = 0u;
  last_pressure_gpa = 0.0;
  have_pressure = FALSE;
  have_energy = FALSE;
  last_energy_ev = 0.0;
  have_fermi_energy = FALSE;
  fermi_energy_ev = 0.0;

  gdis_model_parse_vasp_atom_symbols(contents, atom_symbols);

  for (guint i = 0; lines[i] != NULL; i++)
    {
      gchar *trimmed;

      trimmed = g_strstrip(lines[i]);
      if (trimmed[0] == '\0')
        continue;

      if ((!model->title || model->title[0] == '\0') &&
          g_strstr_len(trimmed, -1, "name=\"SYSTEM\"") != NULL)
        {
          const gchar *start;
          const gchar *end;

          start = strchr(trimmed, '>');
          end = start ? strstr(start + 1, "</i>") : NULL;
          if (start && end && end > start + 1)
            {
              g_clear_pointer(&model->title, g_free);
              model->title = g_strndup(start + 1, (gsize) (end - (start + 1)));
              g_strstrip(model->title);
            }
          continue;
        }

      if (g_str_has_prefix(trimmed, "<varray name=\"kpointlist\""))
        {
          for (i = i + 1u; lines[i] != NULL; i++)
            {
              gchar *k_line;
              gdouble values[3];

              k_line = g_strstrip(lines[i]);
              if (g_str_has_prefix(k_line, "</varray>"))
                break;
              if (gdis_collect_doubles_from_text_scan(k_line, values, 3u) < 3u)
                continue;
              g_array_append_vals(kpoint_triplets, values, 3u);
            }
          continue;
        }

      if (g_str_has_prefix(trimmed, "<calculation>"))
        {
          gboolean calc_have_stress;
          gboolean calc_have_energy;
          gdouble calc_pressure_gpa;
          gdouble calc_energy_ev;

          calc_have_stress = FALSE;
          calc_have_energy = FALSE;
          calc_pressure_gpa = 0.0;
          calc_energy_ev = 0.0;

          for (i = i + 1u; lines[i] != NULL; i++)
            {
              gchar *calc_line;

              calc_line = g_strstrip(lines[i]);
              if (calc_line[0] == '\0')
                continue;
              if (g_str_has_prefix(calc_line, "</calculation>"))
                break;

              if (g_strstr_len(calc_line, -1, "name=\"e_fr_energy\"") != NULL)
                {
                  gdouble values[1];

                  if (gdis_collect_doubles_from_text_scan(calc_line, values, 1u) >= 1u)
                    {
                      calc_energy_ev = values[0];
                      calc_have_energy = TRUE;
                    }
                  continue;
                }

              if (g_str_has_prefix(calc_line, "<varray name=\"stress\""))
                {
                  gdouble diag_sum;
                  guint row_count;

                  diag_sum = 0.0;
                  row_count = 0u;
                  for (i = i + 1u; lines[i] != NULL && row_count < 3u; i++)
                    {
                      gchar *stress_line;
                      gdouble values[3];

                      stress_line = g_strstrip(lines[i]);
                      if (g_str_has_prefix(stress_line, "</varray>"))
                        break;
                      if (gdis_collect_doubles_from_text_scan(stress_line, values, 3u) < 3u)
                        continue;
                      diag_sum += values[row_count];
                      row_count++;
                    }
                  if (row_count == 3u)
                    {
                      calc_pressure_gpa = -(diag_sum / 3.0) * 0.1;
                      calc_have_stress = TRUE;
                    }
                  continue;
                }

              if (g_str_has_prefix(calc_line, "<structure>"))
                {
                  GPtrArray *structure_atoms;
                  gdouble structure_lengths[3];
                  gdouble structure_angles[3];
                  guint structure_end = i;

                  structure_atoms = NULL;
                  if (gdis_model_parse_vasp_structure_block(model,
                                                            lines,
                                                            i,
                                                            atom_symbols,
                                                            &structure_atoms,
                                                            structure_lengths,
                                                            structure_angles,
                                                            &structure_end))
                    {
                      if (best_atoms)
                        g_ptr_array_free(best_atoms, TRUE);
                      best_atoms = structure_atoms;
                      memcpy(best_lengths, structure_lengths, sizeof(best_lengths));
                      memcpy(best_angles, structure_angles, sizeof(best_angles));
                    }
                  i = structure_end;
                  continue;
                }
            }

          if (calc_have_stress)
            {
              gdouble x_value;

              pressure_step++;
              x_value = (gdouble) pressure_step;
              g_array_append_val(pressure_x_values, x_value);
              g_array_append_val(pressure_y_values, calc_pressure_gpa);
              last_pressure_gpa = calc_pressure_gpa;
              have_pressure = TRUE;
            }
          if (calc_have_energy)
            {
              last_energy_ev = calc_energy_ev;
              have_energy = TRUE;
            }
          continue;
        }

      if (g_str_equal(trimmed, "<dos>"))
        {
          GArray *energies;
          GArray *totals;
          guint spin_index;

          energies = gdis_double_array_new();
          totals = gdis_double_array_new();
          spin_index = 0u;

          for (i = i + 1u; lines[i] != NULL; i++)
            {
              gchar *dos_line;

              dos_line = g_strstrip(lines[i]);
              if (g_str_has_prefix(dos_line, "</dos>"))
                break;
              if (g_strstr_len(dos_line, -1, "name=\"efermi\"") != NULL)
                {
                  gdouble values[1];

                  if (gdis_collect_doubles_from_text_scan(dos_line, values, 1u) >= 1u)
                    {
                      fermi_energy_ev = values[0];
                      have_fermi_energy = TRUE;
                    }
                  continue;
                }
              if (g_str_has_prefix(dos_line, "<set comment=\"spin "))
                {
                  spin_index++;
                  guint dos_index = 0u;

                  for (i = i + 1u; lines[i] != NULL; i++)
                    {
                      gchar *row_line;
                      gdouble values[3];

                      row_line = g_strstrip(lines[i]);
                      if (g_str_has_prefix(row_line, "</set>"))
                        break;
                      if (gdis_collect_doubles_from_text_scan(row_line, values, 3u) < 2u)
                        continue;

                      if (spin_index == 1u)
                        {
                          g_array_append_val(energies, values[0]);
                          g_array_append_val(totals, values[1]);
                        }
                      else if (dos_index < totals->len)
                        {
                          gdouble sum_value;

                          sum_value = g_array_index(totals, gdouble, dos_index) + values[1];
                          g_array_index(totals, gdouble, dos_index) = sum_value;
                        }
                      dos_index++;
                    }
                }
            }

          if (energies->len > 0u && totals->len == energies->len)
            {
              for (guint idx = 0; idx < energies->len; idx++)
                {
                  gdouble shifted;

                  shifted = g_array_index(energies, gdouble, idx) -
                            (have_fermi_energy ? fermi_energy_ev : 0.0);
                  g_array_index(energies, gdouble, idx) = shifted;
                }
              g_clear_pointer(&dos_x_values, g_array_unref);
              g_clear_pointer(&dos_y_values, g_array_unref);
              dos_x_values = energies;
              dos_y_values = totals;
              energies = NULL;
              totals = NULL;
            }

          g_clear_pointer(&energies, g_array_unref);
          g_clear_pointer(&totals, g_array_unref);
          continue;
        }

      if (g_str_equal(trimmed, "<eigenvalues>"))
        {
          GPtrArray *kpoint_bands;
          gboolean in_spin_one;

          kpoint_bands = g_ptr_array_new_with_free_func((GDestroyNotify) g_array_unref);
          in_spin_one = FALSE;

          for (i = i + 1u; lines[i] != NULL; i++)
            {
              gchar *band_line;

              band_line = g_strstrip(lines[i]);
              if (g_str_has_prefix(band_line, "<set comment=\"spin 1\""))
                {
                  in_spin_one = TRUE;
                  continue;
                }
              if (g_str_has_prefix(band_line, "<set comment=\"spin 2\"") ||
                  g_str_has_prefix(band_line, "</eigenvalues>"))
                break;
              if (!in_spin_one)
                continue;

              if (g_str_has_prefix(band_line, "<set comment=\"kpoint "))
                {
                  GArray *band_values;

                  band_values = gdis_double_array_new();
                  for (i = i + 1u; lines[i] != NULL; i++)
                    {
                      gchar *row_line;
                      gdouble values[2];

                      row_line = g_strstrip(lines[i]);
                      if (g_str_has_prefix(row_line, "</set>"))
                        break;
                      if (gdis_collect_doubles_from_text_scan(row_line, values, 2u) < 1u)
                        continue;
                      g_array_append_val(band_values, values[0]);
                    }

                  if (band_values->len > 0u)
                    g_ptr_array_add(kpoint_bands, band_values);
                  else
                    g_array_unref(band_values);
                }
            }

          if (kpoint_bands->len > 0u)
            {
              g_clear_pointer(&band_x_values, g_array_unref);
              g_clear_pointer(&band_y_values, g_array_unref);
              gdis_model_build_band_arrays_from_triplets(kpoint_triplets,
                                                         kpoint_bands,
                                                         have_fermi_energy ? fermi_energy_ev : 0.0,
                                                         &band_x_values,
                                                         &band_y_values,
                                                         &band_path_count,
                                                         &band_series_count);
            }

          g_ptr_array_free(kpoint_bands, TRUE);
          continue;
        }
    }

  for (guint i = 0; lines[i] != NULL; i++)
    {
      gchar *trimmed;
      GArray *parsed_dos_x;
      GArray *parsed_dos_y;
      gboolean parsed_have_fermi;
      gdouble parsed_fermi_energy_ev;
      guint end_index;

      trimmed = g_strstrip(lines[i]);
      if (!g_str_equal(trimmed, "<dos>"))
        continue;

      parsed_dos_x = NULL;
      parsed_dos_y = NULL;
      parsed_have_fermi = have_fermi_energy;
      parsed_fermi_energy_ev = fermi_energy_ev;
      end_index = i;
      if (gdis_model_parse_vasp_dos_block(lines,
                                          i,
                                          &parsed_dos_x,
                                          &parsed_dos_y,
                                          &parsed_have_fermi,
                                          &parsed_fermi_energy_ev,
                                          &end_index))
        {
          g_clear_pointer(&dos_x_values, g_array_unref);
          g_clear_pointer(&dos_y_values, g_array_unref);
          dos_x_values = parsed_dos_x;
          dos_y_values = parsed_dos_y;
          have_fermi_energy = parsed_have_fermi;
          fermi_energy_ev = parsed_fermi_energy_ev;
          i = end_index;
        }
    }

  for (guint i = 0; lines[i] != NULL; i++)
    {
      gchar *trimmed;
      GArray *parsed_band_x;
      GArray *parsed_band_y;
      guint parsed_path_count;
      guint parsed_series_count;
      guint end_index;

      trimmed = g_strstrip(lines[i]);
      if (!g_str_equal(trimmed, "<eigenvalues>"))
        continue;

      parsed_band_x = NULL;
      parsed_band_y = NULL;
      parsed_path_count = 0u;
      parsed_series_count = 0u;
      end_index = i;
      if (gdis_model_parse_vasp_eigenvalues_block(lines,
                                                  i,
                                                  kpoint_triplets,
                                                  have_fermi_energy,
                                                  fermi_energy_ev,
                                                  &parsed_band_x,
                                                  &parsed_band_y,
                                                  &parsed_path_count,
                                                  &parsed_series_count,
                                                  &end_index))
        {
          g_clear_pointer(&band_x_values, g_array_unref);
          g_clear_pointer(&band_y_values, g_array_unref);
          band_x_values = parsed_band_x;
          band_y_values = parsed_band_y;
          band_path_count = parsed_path_count;
          band_series_count = parsed_series_count;
          i = end_index;
        }
    }

  g_strfreev(lines);
  g_ptr_array_free(atom_symbols, TRUE);
  g_array_unref(kpoint_triplets);
  gdis_model_clear_plot_arrays(model);

  if (!best_atoms || best_atoms->len == 0u)
    {
      if (best_atoms)
        g_ptr_array_free(best_atoms, TRUE);
      g_clear_pointer(&pressure_x_values, g_array_unref);
      g_clear_pointer(&pressure_y_values, g_array_unref);
      g_clear_pointer(&dos_x_values, g_array_unref);
      g_clear_pointer(&dos_y_values, g_array_unref);
      g_clear_pointer(&band_x_values, g_array_unref);
      g_clear_pointer(&band_y_values, g_array_unref);
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No final structure was parsed from VASP XML '%s'.",
                  model->path);
      return FALSE;
    }

  if (model->atoms)
    g_ptr_array_free(model->atoms, TRUE);
  model->atoms = best_atoms;
  model->periodic = TRUE;
  model->periodicity = 3u;
  memcpy(model->cell_lengths, best_lengths, sizeof(best_lengths));
  memcpy(model->cell_angles, best_angles, sizeof(best_angles));
  model->has_energy = have_energy;
  model->energy_ev = last_energy_ev;
  model->has_pressure = have_pressure;
  model->pressure_gpa = last_pressure_gpa;
  model->has_fermi_energy = have_fermi_energy;
  model->fermi_energy_ev = fermi_energy_ev;
  if (pressure_y_values->len > 0u)
    {
      model->pressure_x_values = pressure_x_values;
      model->pressure_y_values_gpa = pressure_y_values;
    }
  else
    {
      g_array_unref(pressure_x_values);
      g_array_unref(pressure_y_values);
    }
  model->dos_x_values_ev = dos_x_values;
  model->dos_y_values = dos_y_values;
  model->band_x_values = band_x_values;
  model->band_y_values_ev = band_y_values;
  model->band_path_count = band_path_count;
  model->band_series_count = band_series_count;

  return TRUE;
}

static gboolean
gdis_model_load_qe_input(GdisModel *model, const gchar *contents, GError **error)
{
  return gdis_model_load_qe_like(model, contents, FALSE, error);
}

static gboolean
gdis_model_load_qe_output(GdisModel *model, const gchar *contents, GError **error)
{
  return gdis_model_load_qe_like(model, contents, TRUE, error);
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
  GArray *symmetry_ops;
  gboolean in_multiline;
  gboolean can_expand_symmetry;
  guint expected_symmetry_atom_count;
  guint i;

  lines = g_strsplit(contents, "\n", -1);
  records = g_ptr_array_new_with_free_func(gdis_cif_atom_record_free);
  symmetry_ops = g_array_new(FALSE, FALSE, sizeof(GdisSymmetryOp));
  in_multiline = FALSE;
  can_expand_symmetry = FALSE;
  expected_symmetry_atom_count = 0u;

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
          if (records->len > 0u)
            break;

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
          gint symop_pos;

          headers = g_ptr_array_new_with_free_func(g_free);
          atom_loop = FALSE;
          symop_pos = -1;
          i++;

          while (lines[i] != NULL)
            {
              gchar *header_line;
              gchar *header_trimmed;
              g_autofree gchar *lower_header = NULL;

              header_line = lines[i];
              g_strchomp(header_line);
              header_trimmed = g_strstrip(header_line);

              if (header_trimmed[0] != '_')
                break;

              lower_header = g_ascii_strdown(header_trimmed, -1);
              g_ptr_array_add(headers, g_strdup(header_trimmed));
              if (g_str_has_prefix(header_trimmed, "_atom_site_"))
                atom_loop = TRUE;
              else if (g_str_equal(lower_header, "_symmetry_equiv_pos_as_xyz") ||
                       g_str_equal(lower_header, "_equiv_pos_as_xyz") ||
                       g_str_equal(lower_header, "_space_group_symop_operation_xyz") ||
                       g_str_equal(lower_header, "_space_group_symop.operation_xyz"))
                symop_pos = (gint) headers->len - 1;
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
                  gint mult_pos;
                  gint calc_flag_pos;
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
                      mult_pos = -1;
                      calc_flag_pos = -1;

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
                          else if (g_str_equal(lower_header, "_atom_site_symmetry_multiplicity"))
                            mult_pos = (gint) h;
                          else if (g_str_equal(lower_header, "_atom_site_calc_flag"))
                            calc_flag_pos = (gint) h;

                          g_free(lower_header);
                        }

                      fractional = (frac_x_pos >= 0 && frac_y_pos >= 0 && frac_z_pos >= 0);
                      if (fractional || (cart_x_pos >= 0 && cart_y_pos >= 0 && cart_z_pos >= 0))
                        {
                          if (calc_flag_pos >= 0 && tokens->len > (guint) calc_flag_pos)
                            {
                              const gchar *calc_flag;

                              calc_flag = g_ptr_array_index(tokens, calc_flag_pos);
                              if (g_ascii_strcasecmp(calc_flag, "dum") == 0 ||
                                  g_ascii_strcasecmp(calc_flag, "dummy") == 0)
                                {
                                  g_ptr_array_free(tokens, TRUE);
                                  i++;
                                  continue;
                                }
                            }

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
                          record->multiplicity = 0u;
                          record->have_multiplicity = FALSE;
                          if (mult_pos >= 0 && tokens->len > (guint) mult_pos)
                            {
                              const gchar *multiplicity_text;
                              guint multiplicity;

                              multiplicity_text = g_ptr_array_index(tokens, mult_pos);
                              multiplicity = 0u;
                              if (gdis_try_parse_uint(multiplicity_text, &multiplicity) &&
                                  multiplicity > 0u)
                                {
                                  record->multiplicity = multiplicity;
                                  record->have_multiplicity = TRUE;
                                }
                              else
                                {
                                  gdouble multiplicity_value;
                                  glong rounded;

                                  multiplicity_value = gdis_parse_cif_number(multiplicity_text);
                                  rounded = lround(multiplicity_value);
                                  if (multiplicity_value > 0.0 &&
                                      fabs(multiplicity_value - (gdouble) rounded) < 1.0e-4)
                                    {
                                      record->multiplicity = (guint) rounded;
                                      record->have_multiplicity = TRUE;
                                    }
                                }
                            }
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
          else if (symop_pos >= 0)
            {
              while (lines[i] != NULL)
                {
                  GPtrArray *tokens;
                  gchar *row_line;
                  gchar *row_trimmed;

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
                  if (symop_pos >= 0 && tokens->len > (guint) symop_pos)
                    {
                      GdisSymmetryOp op;

                      if (gdis_parse_cif_symmetry_operation(g_ptr_array_index(tokens, symop_pos), &op))
                        g_array_append_val(symmetry_ops, op);
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
      can_expand_symmetry = (has_cell_matrix &&
                             model->periodic &&
                             model->periodicity == 3u);

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
                  g_array_unref(symmetry_ops);
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

  expected_symmetry_atom_count = gdis_cif_sum_reported_multiplicities(records);

  if (can_expand_symmetry)
    {
      if (symmetry_ops->len > 0u)
        gdis_model_expand_atoms_with_symmetry_ops(model, symmetry_ops);
      else if (model->space_group && model->space_group[0] != '\0')
        {
          if (!gdis_model_expand_space_group_atoms_best_match(model,
                                                              expected_symmetry_atom_count))
            gdis_model_expand_space_group_atoms(model, 0);
        }
    }

  g_array_unref(symmetry_ops);
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
  state.frames = g_ptr_array_new_with_free_func(gdis_model_frame_free);
  state.text = g_string_new("");

  context = g_markup_parse_context_new(&parser, G_MARKUP_TREAT_CDATA_AS_TEXT, &state, NULL);
  ok = g_markup_parse_context_parse(context, contents, strlen(contents), error);
  if (ok)
    ok = g_markup_parse_context_end_parse(context, error);
  if (ok)
    ok = gdis_qbox_finish_current_frame(&state, error);
  g_markup_parse_context_free(context);

  if (state.current_frame)
    gdis_model_frame_free(state.current_frame);
  g_clear_pointer(&state.current_species_name, g_free);
  g_clear_pointer(&state.current_atom_name, g_free);
  g_clear_pointer(&state.current_atom_species, g_free);
  if (state.text)
    g_string_free(state.text, TRUE);
  g_clear_pointer(&state.species_symbols, g_hash_table_unref);

  if (!ok)
    {
      g_clear_pointer(&state.frames, g_ptr_array_unref);
      return FALSE;
    }

  if (!state.looks_like_qbox)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_UNSUPPORTED_FORMAT,
                  "'%s' is not a supported Qbox XML result file.",
                  model->path);
      return FALSE;
    }

  if (!state.frames || state.frames->len == 0u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "No atoms were parsed from Qbox result '%s'.",
                  model->path);
      return FALSE;
    }

  gdis_model_frame_apply(model, g_ptr_array_index(state.frames, 0u));
  model->current_frame_index = 0u;
  if (state.frames->len > 1u)
    model->frames = g_steal_pointer(&state.frames);
  else
    g_clear_pointer(&state.frames, g_ptr_array_unref);

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

  if (g_str_equal(local_name, "iteration"))
    {
      state->inside_iteration = TRUE;
      state->current_iteration_count = 0u;
      state->have_current_iteration_energy = FALSE;

      for (guint i = 0; attribute_names && attribute_names[i] != NULL; i++)
        {
          guint parsed_count;

          if (!g_str_equal(gdis_xml_local_name(attribute_names[i]), "count"))
            continue;
          if (gdis_try_parse_uint(attribute_values[i], &parsed_count))
            state->current_iteration_count = parsed_count;
          break;
        }
      return;
    }

  if (g_str_equal(local_name, "atomset"))
    {
      if (state->current_frame && !gdis_qbox_finish_current_frame(state, parse_error))
        return;
      if (!gdis_qbox_begin_frame(state, parse_error))
        return;
      return;
    }

  if (g_str_equal(local_name, "etotal"))
    {
      state->inside_etotal = TRUE;
      g_string_truncate(state->text, 0);
      return;
    }

  if (g_str_equal(local_name, "unit_cell"))
    {
      const gchar *attr_a = NULL;
      const gchar *attr_b = NULL;
      const gchar *attr_c = NULL;
      gdouble a_vec[3];
      gdouble b_vec[3];
      gdouble c_vec[3];

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
          !gdis_parse_vector3_text(attr_a, a_vec) ||
          !gdis_parse_vector3_text(attr_b, b_vec) ||
          !gdis_parse_vector3_text(attr_c, c_vec))
        {
          g_set_error(parse_error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_PARSE,
                      "Qbox unit_cell in '%s' is missing or invalid.",
                      state->model->path);
          return;
        }

      for (guint axis = 0; axis < 3; axis++)
        {
          a_vec[axis] *= 0.529177210903;
          b_vec[axis] *= 0.529177210903;
          c_vec[axis] *= 0.529177210903;
        }

      if (state->current_frame)
        {
          if (!gdis_qbox_apply_cell_vectors_to_frame(state,
                                                     state->current_frame,
                                                     a_vec,
                                                     b_vec,
                                                     c_vec,
                                                     parse_error))
            return;
          state->current_frame_has_cell = TRUE;
        }
      else
        {
          gdis_vec3_copy(state->default_a_vec, a_vec);
          gdis_vec3_copy(state->default_b_vec, b_vec);
          gdis_vec3_copy(state->default_c_vec, c_vec);
          state->have_default_cell = TRUE;
        }
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
      state->have_current_atom_force = FALSE;

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
      return;
    }

  if (g_str_equal(local_name, "force"))
    {
      state->inside_atom_force = TRUE;
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

  if (g_str_equal(local_name, "force"))
    {
      state->inside_atom_force = FALSE;
      state->have_current_atom_force = gdis_parse_vector3_text(state->text->str,
                                                               state->current_atom_force);
      return;
    }

  if (g_str_equal(local_name, "etotal"))
    {
      gdouble parsed_energy = 0.0;

      state->inside_etotal = FALSE;
      if (!gdis_try_parse_double_relaxed(state->text->str, &parsed_energy))
        return;

      state->have_current_iteration_energy = TRUE;
      state->current_iteration_energy_ev = parsed_energy * GDIS_QBOX_HARTREE_TO_EV;
      if (state->current_frame)
        {
          state->current_frame->has_energy = TRUE;
          state->current_frame->energy_ev = state->current_iteration_energy_ev;
        }
      return;
    }

  if (g_str_equal(local_name, "atom"))
    {
      g_autofree gchar *species_symbol = NULL;
      const gchar *mapped_symbol;
      const gchar *element_text;
      const gchar *label_text;
      GPtrArray *target_atoms;

      if (!state->current_frame && !gdis_qbox_begin_frame(state, parse_error))
        return;

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

      target_atoms = state->current_frame ? state->current_frame->atoms : state->model->atoms;
      g_ptr_array_add(target_atoms,
                      gdis_atom_new(label_text,
                                    element_text,
                                    element_text,
                                    state->current_atom_position[0],
                                    state->current_atom_position[1],
                                    state->current_atom_position[2],
                                    1.0,
                                    -1,
                                    target_atoms->len + 1u));

      if (state->have_current_atom_force && state->current_frame)
        {
          const gdouble fx = state->current_atom_force[0] * GDIS_QBOX_HARTREE_PER_BOHR_TO_EV_PER_ANG;
          const gdouble fy = state->current_atom_force[1] * GDIS_QBOX_HARTREE_PER_BOHR_TO_EV_PER_ANG;
          const gdouble fz = state->current_atom_force[2] * GDIS_QBOX_HARTREE_PER_BOHR_TO_EV_PER_ANG;

          state->current_frame_force_sum_sq += fx * fx + fy * fy + fz * fz;
          state->current_frame_force_count++;
        }

      g_clear_pointer(&state->current_atom_name, g_free);
      g_clear_pointer(&state->current_atom_species, g_free);
      state->have_current_atom_position = FALSE;
      state->have_current_atom_force = FALSE;
      return;
    }

  if (g_str_equal(local_name, "atomset"))
    {
      if (!gdis_qbox_finish_current_frame(state, parse_error))
        return;
      return;
    }

  if (g_str_equal(local_name, "species"))
    g_clear_pointer(&state->current_species_name, g_free);

  if (g_str_equal(local_name, "iteration"))
    {
      state->inside_iteration = FALSE;
      state->current_iteration_count = 0u;
      state->have_current_iteration_energy = FALSE;
    }
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
  if (state->inside_species_symbol ||
      state->inside_atom_position ||
      state->inside_etotal ||
      state->inside_atom_force)
    g_string_append_len(state->text, text, text_len);
}

gboolean
gdis_qbox_apply_cell_vectors_to_frame(GdisQboxXmlParseState *state,
                                      GdisModelFrame *frame,
                                      const gdouble a_vec[3],
                                      const gdouble b_vec[3],
                                      const gdouble c_vec[3],
                                      GError **error)
{
  gdouble lengths[3];
  gdouble angles[3];

  g_return_val_if_fail(state != NULL, FALSE);
  g_return_val_if_fail(frame != NULL, FALSE);

  if (!gdis_surface_build_cell_from_vectors(a_vec, b_vec, c_vec, lengths, angles))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_PARSE,
                  "Qbox cell vectors in '%s' are invalid.",
                  state->model->path);
      return FALSE;
    }

  frame->periodic = TRUE;
  frame->periodicity = 3u;
  memcpy(frame->cell_lengths, lengths, sizeof(lengths));
  memcpy(frame->cell_angles, angles, sizeof(angles));
  return TRUE;
}

gboolean
gdis_qbox_begin_frame(GdisQboxXmlParseState *state,
                      GError **error)
{
  const char *base_title;

  g_return_val_if_fail(state != NULL, FALSE);

  if (state->current_frame)
    return TRUE;

  state->atomset_count++;
  state->current_frame = gdis_model_frame_new();
  state->current_frame_has_cell = FALSE;
  state->current_frame_force_sum_sq = 0.0;
  state->current_frame_force_count = 0u;

  base_title = state->model->basename ? state->model->basename : "qbox";
  if (state->inside_iteration)
    {
      guint iteration;

      iteration = state->current_iteration_count > 0u
        ? state->current_iteration_count
        : state->atomset_count;
      state->current_frame->title = g_strdup_printf("%s iteration %u",
                                                    base_title,
                                                    iteration);
    }
  else if (state->atomset_count == 1u)
    {
      state->current_frame->title = g_strdup(base_title);
    }
  else
    {
      state->current_frame->title = g_strdup_printf("%s frame %u",
                                                    base_title,
                                                    state->atomset_count);
    }

  if (state->have_default_cell)
    {
      if (!gdis_qbox_apply_cell_vectors_to_frame(state,
                                                 state->current_frame,
                                                 state->default_a_vec,
                                                 state->default_b_vec,
                                                 state->default_c_vec,
                                                 error))
        return FALSE;
      state->current_frame_has_cell = TRUE;
    }

  if (state->have_current_iteration_energy)
    {
      state->current_frame->has_energy = TRUE;
      state->current_frame->energy_ev = state->current_iteration_energy_ev;
    }

  return TRUE;
}

gboolean
gdis_qbox_finish_current_frame(GdisQboxXmlParseState *state,
                               GError **error)
{
  g_return_val_if_fail(state != NULL, FALSE);

  if (!state->current_frame)
    return TRUE;

  if (state->current_frame->atoms->len == 0u)
    {
      gdis_model_frame_free(state->current_frame);
      state->current_frame = NULL;
      state->current_frame_has_cell = FALSE;
      return TRUE;
    }

  if (!state->current_frame_has_cell && state->have_default_cell)
    {
      if (!gdis_qbox_apply_cell_vectors_to_frame(state,
                                                 state->current_frame,
                                                 state->default_a_vec,
                                                 state->default_b_vec,
                                                 state->default_c_vec,
                                                 error))
        return FALSE;
      state->current_frame_has_cell = TRUE;
    }

  if (state->current_frame_force_count > 0u)
    {
      state->current_frame->has_force_rms = TRUE;
      state->current_frame->force_rms_ev_ang =
        sqrt(state->current_frame_force_sum_sq / (gdouble) state->current_frame_force_count);
    }

  g_ptr_array_add(state->frames, state->current_frame);
  state->current_frame = NULL;
  state->current_frame_has_cell = FALSE;
  state->current_frame_force_sum_sq = 0.0;
  state->current_frame_force_count = 0u;
  return TRUE;
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

gboolean
gdis_model_set_element_covalent_override(GdisModel *model,
                                         const gchar *symbol,
                                         gdouble radius)
{
  g_autofree gchar *normalized_symbol = NULL;
  gdouble *radius_copy;

  g_return_val_if_fail(model != NULL, FALSE);

  normalized_symbol = gdis_normalize_element_symbol(symbol);
  if (!normalized_symbol || normalized_symbol[0] == '\0')
    return FALSE;

  if (!model->element_covalent_overrides)
    {
      model->element_covalent_overrides = g_hash_table_new_full(g_str_hash,
                                                                g_str_equal,
                                                                g_free,
                                                                g_free);
    }

  radius_copy = g_new(gdouble, 1);
  *radius_copy = MAX(radius, 0.0);
  g_hash_table_replace(model->element_covalent_overrides,
                       g_strdup(normalized_symbol),
                       radius_copy);
  return TRUE;
}

void
gdis_model_copy_element_overrides(GdisModel *dest, const GdisModel *src)
{
  GHashTableIter iter;
  gpointer key;
  gpointer value;

  g_return_if_fail(dest != NULL);

  g_clear_pointer(&dest->element_covalent_overrides, g_hash_table_unref);
  if (!src || !src->element_covalent_overrides)
    return;

  g_hash_table_iter_init(&iter, src->element_covalent_overrides);
  while (g_hash_table_iter_next(&iter, &key, &value))
    gdis_model_set_element_covalent_override(dest, key, *((const gdouble *) value));
}

gdouble
gdis_model_lookup_covalent_radius(const GdisModel *model, const GdisAtom *atom)
{
  const gdouble *override_radius;

  if (!atom)
    return gdis_lookup_covalent_radius(NULL);

  if (model &&
      model->element_covalent_overrides &&
      atom->element &&
      atom->element[0] != '\0')
    {
      override_radius = g_hash_table_lookup(model->element_covalent_overrides, atom->element);
      if (override_radius)
        return MAX(*override_radius, 0.0);
    }

  return gdis_lookup_covalent_radius(atom->element);
}

gboolean
gdis_model_should_infer_bond(const GdisModel *model,
                             const GdisAtom *atom_i,
                             const GdisAtom *atom_j,
                             gdouble radius_i,
                             gdouble radius_j,
                             gdouble *max_distance_out)
{
  gdouble max_distance;

  g_return_val_if_fail(atom_i != NULL, FALSE);
  g_return_val_if_fail(atom_j != NULL, FALSE);

  max_distance = (MAX(radius_i, 0.0) + MAX(radius_j, 0.0)) * 1.10;
  if (model && model->periodic)
    max_distance = MIN(max_distance, 3.2);

  if (max_distance_out)
    *max_distance_out = max_distance;

  return max_distance > 0.0;
}

static gboolean
gdis_model_infer_bonds(GdisModel *model)
{
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

  use_periodic_image = gdis_model_build_cell_matrix(model, matrix, inverse) && model->periodic;

  for (i = 0; i + 1 < model->atoms->len; i++)
    {
      GdisAtom *atom_i;
      gdouble radius_i;

      atom_i = g_ptr_array_index(model->atoms, i);
      radius_i = gdis_model_lookup_covalent_radius(model, atom_i);

      for (j = i + 1; j < model->atoms->len; j++)
        {
          GdisAtom *atom_j;
          gdouble radius_j;
          gdouble max_distance;
          gdouble distance2;

          atom_j = g_ptr_array_index(model->atoms, j);
          radius_j = gdis_model_lookup_covalent_radius(model, atom_j);
          if (!gdis_model_should_infer_bond(model,
                                            atom_i,
                                            atom_j,
                                            radius_i,
                                            radius_j,
                                            &max_distance))
            continue;

          distance2 = gdis_model_distance2(model,
                                           atom_i->position,
                                           atom_j->position,
                                           matrix,
                                           inverse,
                                           use_periodic_image);
          if (distance2 < 0.16 || distance2 > max_distance * max_distance)
            continue;

          gdis_model_add_bond(model, i, j, 1, TRUE);
        }
    }

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
  gdouble visible_charge_total;

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

  if (!model->has_total_charge)
    {
      if (gdis_model_sum_legacy_visible_charge(model->atoms, &visible_charge_total))
        model->total_charge_e = visible_charge_total;
      else
        model->total_charge_e = 0.0;
      model->has_total_charge = TRUE;
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

static guint
gdis_cif_sum_reported_multiplicities(const GPtrArray *records)
{
  guint64 total;

  if (!records || records->len == 0u)
    return 0u;

  total = 0u;
  for (guint i = 0; i < records->len; i++)
    {
      const GdisCifAtomRecord *record;

      record = g_ptr_array_index(records, i);
      if (!record || !record->have_multiplicity || record->multiplicity == 0u)
        return 0u;

      total += record->multiplicity;
      if (total > G_MAXUINT)
        return 0u;
    }

  return (guint) total;
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

static gint
gdis_str_ibegin(const char *s1, const char *s2)
{
  char u1;
  char u2;

  while (*s1 && *s2)
    {
      u1 = (char) g_ascii_toupper((guchar) *s1++);
      u2 = (char) g_ascii_toupper((guchar) *s2++);
      if (u1 < u2)
        return -1;
      if (u1 > u2)
        return 1;
    }

  if (*s2)
    return -1;

  return 0;
}

static gint
gdis_build_sginfo_compat(T_SgInfo *sg_info, const gchar *space_group_text)
{
  gint volume_letter;
  const T_TabSgName *table_entry;

  g_return_val_if_fail(sg_info != NULL, -1);
  g_return_val_if_fail(space_group_text != NULL, -1);

  while (*space_group_text && g_ascii_isspace(*space_group_text))
    space_group_text++;

  volume_letter = -1;
  if (g_ascii_isdigit(*space_group_text))
    volume_letter = 'A';
  else if (gdis_str_ibegin(space_group_text, "VolA") == 0)
    {
      volume_letter = 'A';
      space_group_text += 4;
    }
  else if (gdis_str_ibegin(space_group_text, "VolI") == 0 ||
           gdis_str_ibegin(space_group_text, "Vol1") == 0)
    {
      volume_letter = 'I';
      space_group_text += 4;
    }
  else if (gdis_str_ibegin(space_group_text, "Hall") == 0)
    {
      volume_letter = 0;
      space_group_text += 4;
    }

  while (*space_group_text && g_ascii_isspace(*space_group_text))
    space_group_text++;

  if (volume_letter == -1)
    volume_letter = 'A';

  table_entry = NULL;
  if (volume_letter != 0)
    {
      table_entry = FindTabSgNameEntry(space_group_text, volume_letter);
      if (!table_entry)
        return -1;
      space_group_text = table_entry->HallSymbol;
    }

  SgError = NULL;
  sg_info->MaxList = 192;
  sg_info->ListSeitzMx = g_malloc(sg_info->MaxList * sizeof(*sg_info->ListSeitzMx));
  sg_info->ListRotMxInfo = g_malloc(sg_info->MaxList * sizeof(*sg_info->ListRotMxInfo));

  InitSgInfo(sg_info);
  sg_info->TabSgName = table_entry;

  ParseHallSymbol(space_group_text, sg_info);
  if (SgError != NULL)
    return -1;

  return CompleteSgInfo(sg_info);
}

static void
gdis_free_sginfo_compat(T_SgInfo *sg_info)
{
  if (!sg_info)
    return;

  g_free(sg_info->ListSeitzMx);
  g_free(sg_info->ListRotMxInfo);
  sg_info->ListSeitzMx = NULL;
  sg_info->ListRotMxInfo = NULL;
}

static gchar *
gdis_format_space_group_label(const T_SgInfo *sg_info)
{
  gchar *label;
  gchar *equals;
  gchar *stripped;

  if (!sg_info || !sg_info->TabSgName || !sg_info->TabSgName->SgLabels)
    return NULL;

  label = g_strdup(sg_info->TabSgName->SgLabels);
  equals = strstr(label, "=");
  if (equals)
    *equals = '\0';

  for (gchar *iter = label; *iter != '\0'; iter++)
    {
      if (*iter == '_')
        *iter = ' ';
    }

  stripped = gdis_strdup_strip(label);
  g_free(label);
  return stripped;
}

static gboolean
gdis_model_uses_rhombohedral_setting(const GdisModel *model,
                                     const gchar *space_group_text)
{
  const gchar *iter;
  gchar *endptr;
  long sg_number;
  gdouble dx;
  gdouble x2;

  g_return_val_if_fail(space_group_text != NULL, FALSE);

  iter = space_group_text;
  while (*iter && g_ascii_isspace(*iter))
    iter++;

  if (strchr(iter, ':') != NULL)
    return g_strrstr(iter, ":R") != NULL;

  if (g_ascii_toupper(*iter) == 'R')
    goto cell_check;

  sg_number = strtol(iter, &endptr, 10);
  if (endptr == iter ||
      (sg_number != 146 && sg_number != 148 && sg_number != 155 &&
       sg_number != 160 && sg_number != 161 && sg_number != 166 &&
       sg_number != 167))
    return FALSE;

cell_check:
  if (!model)
    return FALSE;

  dx = model->cell_lengths[0] - model->cell_lengths[1];
  x2 = dx * dx;
  dx = model->cell_lengths[0] - model->cell_lengths[2];
  x2 += dx * dx;
  dx = model->cell_lengths[1] - model->cell_lengths[2];
  x2 += dx * dx;

  return x2 < 0.001;
}

static gboolean
gdis_space_group_is_r_lattice(const gchar *space_group_text)
{
  const gchar *iter;
  gchar *endptr;
  long sg_number;

  g_return_val_if_fail(space_group_text != NULL, FALSE);

  iter = space_group_text;
  while (*iter && g_ascii_isspace(*iter))
    iter++;

  if (g_ascii_toupper(*iter) == 'R')
    return TRUE;

  sg_number = strtol(iter, &endptr, 10);
  if (endptr == iter)
    return FALSE;

  return sg_number == 146 || sg_number == 148 || sg_number == 155 ||
         sg_number == 160 || sg_number == 161 || sg_number == 166 ||
         sg_number == 167;
}

static gboolean
gdis_build_space_group_operations(const GdisModel *model,
                                  const gchar *space_group_text,
                                  gint origin_choice,
                                  GArray **ops_out,
                                  gchar **canonical_label_out)
{
  T_SgInfo sg_info = {0};
  g_autofree gchar *lookup_name = NULL;
  GArray *ops;
  const gint *tr_vector;
  gint translation_count;
  gint inversion_loops;

  g_return_val_if_fail(space_group_text != NULL, FALSE);
  g_return_val_if_fail(ops_out != NULL, FALSE);

  if (!strchr(space_group_text, ':') &&
      gdis_model_uses_rhombohedral_setting(model, space_group_text))
    lookup_name = g_strdup_printf("%s:R", space_group_text);
  else if (!strchr(space_group_text, ':') &&
           gdis_space_group_is_r_lattice(space_group_text))
    lookup_name = g_strdup_printf("%s:H", space_group_text);
  else
    lookup_name = g_strdup(space_group_text);

  if (origin_choice > 0 && !strchr(lookup_name, ':'))
    {
      g_autofree gchar *with_origin = g_strdup_printf("%s:%d", lookup_name, origin_choice);
      g_free(lookup_name);
      lookup_name = g_steal_pointer(&with_origin);
    }

  if (gdis_build_sginfo_compat(&sg_info, lookup_name) != 0)
    return FALSE;

  ops = g_array_new(FALSE, FALSE, sizeof(GdisSymmetryOp));
  translation_count = sg_info.LatticeInfo ? sg_info.LatticeInfo->nTrVector : 0;
  tr_vector = sg_info.LatticeInfo ? sg_info.LatticeInfo->TrVector : NULL;
  inversion_loops = Sg_nLoopInv(&sg_info);

  if (translation_count < 1 || tr_vector == NULL || sg_info.nList < 1)
    {
      g_array_unref(ops);
      gdis_free_sginfo_compat(&sg_info);
      return FALSE;
    }

  for (gint i_translation = 0; i_translation < translation_count; i_translation++, tr_vector += 3)
    {
      for (gint i_inversion = 0; i_inversion < inversion_loops; i_inversion++)
        {
          gint sign;

          sign = (i_inversion == 0) ? 1 : -1;
          for (gint i_list = 0; i_list < sg_info.nList; i_list++)
            {
              const T_RTMx *seitz;
              GdisSymmetryOp op;

              seitz = &sg_info.ListSeitzMx[i_list];
              memset(&op, 0, sizeof(op));

              for (gint row = 0; row < 3; row++)
                {
                  gint translated;

                  for (gint col = 0; col < 3; col++)
                    op.matrix[row * 3 + col] = (gdouble) (sign * seitz->s.R[row * 3 + col]);

                  translated = iModPositive(sign * seitz->s.T[row] + tr_vector[row], STBF);
                  if (translated > STBF / 2)
                    translated -= STBF;
                  op.offset[row] = (gdouble) translated / (gdouble) STBF;
                }

              g_array_append_val(ops, op);
            }
        }
    }

  if (canonical_label_out)
    *canonical_label_out = gdis_format_space_group_label(&sg_info);

  gdis_free_sginfo_compat(&sg_info);
  *ops_out = ops;
  return TRUE;
}

static gdouble
gdis_wrap_fractional_coordinate(gdouble value)
{
  value -= floor(value);
  if (value < 0.0)
    value += 1.0;
  if (value >= 1.0)
    value -= floor(value);
  if (fabs(value - 1.0) < 1.0e-10)
    value = 0.0;
  return value;
}

static void
gdis_wrap_fractional_vector(gdouble frac[3], guint periodicity)
{
  guint dims;

  dims = MIN(periodicity, 3u);
  for (guint i = 0; i < dims; i++)
    frac[i] = gdis_wrap_fractional_coordinate(frac[i]);
}

static gboolean
gdis_fractional_positions_equal(const gdouble left[3],
                                const gdouble right[3],
                                guint periodicity,
                                gdouble tolerance)
{
  guint dims;

  dims = MIN(periodicity, 3u);
  for (guint i = 0; i < 3u; i++)
    {
      gdouble delta;

      delta = fabs(left[i] - right[i]);
      if (i < dims)
        delta = MIN(delta, fabs(1.0 - delta));
      if (delta > tolerance)
        return FALSE;
    }

  return TRUE;
}

static GdisAtom *
gdis_atom_clone_with_position(const GdisAtom *source,
                              const gdouble position[3],
                              guint serial)
{
  GdisAtom *copy;

  g_return_val_if_fail(source != NULL, NULL);
  g_return_val_if_fail(position != NULL, NULL);

  copy = g_new0(GdisAtom, 1);
  copy->serial = serial;
  copy->label = g_strdup(source->label);
  copy->element = g_strdup(source->element);
  copy->ff_type = g_strdup(source->ff_type);
  copy->position[0] = position[0];
  copy->position[1] = position[1];
  copy->position[2] = position[2];
  copy->occupancy = source->occupancy;
  copy->charge = source->charge;
  copy->has_charge = source->has_charge;
  copy->region = source->region;
  return copy;
}

static gboolean
gdis_model_sum_charge_atoms_impl(const GPtrArray *atoms,
                                 gboolean periodic,
                                 guint periodicity,
                                 const gdouble lengths[3],
                                 const gdouble angles[3],
                                 const gchar *space_group,
                                 gint origin_choice,
                                 gboolean expand_space_group,
                                 gboolean require_mixed_signs,
                                 gdouble *total_charge_out)
{
  GdisModel charge_model = {0};
  gboolean saw_positive;
  gboolean saw_negative;
  gdouble total_charge;

  g_return_val_if_fail(total_charge_out != NULL, FALSE);

  *total_charge_out = 0.0;
  if (!atoms || atoms->len == 0u)
    return FALSE;

  saw_positive = FALSE;
  saw_negative = FALSE;
  total_charge = 0.0;
  charge_model.atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  for (guint i = 0; i < atoms->len; i++)
    {
      const GdisAtom *source_atom;

      source_atom = g_ptr_array_index(atoms, i);
      if (!source_atom || !source_atom->has_charge)
        continue;
      if (source_atom->charge > 1.0e-8)
        saw_positive = TRUE;
      if (source_atom->charge < -1.0e-8)
        saw_negative = TRUE;

      g_ptr_array_add(charge_model.atoms,
                      gdis_atom_clone_with_position(source_atom,
                                                    source_atom->position,
                                                    source_atom->serial));
    }

  if (charge_model.atoms->len == 0u)
    {
      g_ptr_array_free(charge_model.atoms, TRUE);
      return FALSE;
    }

  if (require_mixed_signs && (!saw_positive || !saw_negative))
    {
      g_ptr_array_free(charge_model.atoms, TRUE);
      return FALSE;
    }

  charge_model.periodic = periodic;
  charge_model.periodicity = periodicity;
  memcpy(charge_model.cell_lengths, lengths, sizeof(charge_model.cell_lengths));
  memcpy(charge_model.cell_angles, angles, sizeof(charge_model.cell_angles));
  charge_model.space_group = g_strdup(space_group);

  if (expand_space_group &&
      charge_model.periodic &&
      charge_model.periodicity == 3u &&
      charge_model.space_group &&
      charge_model.space_group[0] != '\0')
    gdis_model_expand_space_group_atoms(&charge_model, origin_choice);

  for (guint i = 0; i < charge_model.atoms->len; i++)
    {
      const GdisAtom *atom;

      atom = g_ptr_array_index(charge_model.atoms, i);
      if (!atom || !atom->has_charge)
        continue;
      total_charge += atom->charge;
    }

  g_ptr_array_free(charge_model.atoms, TRUE);
  g_free(charge_model.space_group);
  *total_charge_out = total_charge;
  return TRUE;
}

static gboolean
gdis_model_sum_charge_atoms(const GPtrArray *atoms,
                            gboolean periodic,
                            guint periodicity,
                            const gdouble lengths[3],
                            const gdouble angles[3],
                            const gchar *space_group,
                            gint origin_choice,
                            gboolean expand_space_group,
                            gdouble *total_charge_out)
{
  return gdis_model_sum_charge_atoms_impl(atoms,
                                          periodic,
                                          periodicity,
                                          lengths,
                                          angles,
                                          space_group,
                                          origin_choice,
                                          expand_space_group,
                                          TRUE,
                                          total_charge_out);
}

static gboolean
gdis_model_sum_known_charge_atoms(const GPtrArray *atoms,
                                  gboolean periodic,
                                  guint periodicity,
                                  const gdouble lengths[3],
                                  const gdouble angles[3],
                                  const gchar *space_group,
                                  gint origin_choice,
                                  gboolean expand_space_group,
                                  gdouble *total_charge_out)
{
  return gdis_model_sum_charge_atoms_impl(atoms,
                                          periodic,
                                          periodicity,
                                          lengths,
                                          angles,
                                          space_group,
                                          origin_choice,
                                          expand_space_group,
                                          FALSE,
                                          total_charge_out);
}

static gboolean
gdis_model_sum_legacy_visible_charge(const GPtrArray *atoms,
                                     gdouble *total_charge_out)
{
  gboolean saw_atom;
  gdouble total_charge;

  g_return_val_if_fail(total_charge_out != NULL, FALSE);

  *total_charge_out = 0.0;
  if (!atoms || atoms->len == 0u)
    return FALSE;

  saw_atom = FALSE;
  total_charge = 0.0;
  for (guint i = 0; i < atoms->len; i++)
    {
      const GdisAtom *atom;

      atom = g_ptr_array_index(atoms, i);
      if (!atom)
        continue;

      saw_atom = TRUE;
      if (atom->has_charge)
        {
          total_charge += atom->charge;
        }
      else
        {
          gdouble legacy_charge;

          if (gdis_element_lookup_legacy_charge(atom->element, &legacy_charge))
            total_charge += legacy_charge;
        }
    }

  if (!saw_atom)
    return FALSE;

  if (fabs(total_charge) < 1.0e-8)
    total_charge = 0.0;

  *total_charge_out = total_charge;
  return TRUE;
}

static gboolean
gdis_model_expand_atoms_with_symmetry_ops(GdisModel *model,
                                          const GArray *ops)
{
  GArray *fractional_positions;
  GPtrArray *expanded_atoms;
  gdouble matrix[9];
  gdouble inverse[9];
  gboolean expanded;
  guint original_count;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(ops != NULL, FALSE);

  if (!model->periodic ||
      model->periodicity < 3u ||
      !model->atoms ||
      model->atoms->len == 0u ||
      ops->len == 0u ||
      !gdis_model_build_cell_matrix(model, matrix, inverse))
    return FALSE;

  expanded_atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  fractional_positions = g_array_new(FALSE, FALSE, sizeof(GdisFractionalPosition));
  expanded = FALSE;
  original_count = model->atoms->len;

  for (guint atom_index = 0; atom_index < original_count; atom_index++)
    {
      const GdisAtom *source_atom;
      gdouble source_frac[3];

      source_atom = g_ptr_array_index(model->atoms, atom_index);
      gdis_cart_to_frac(inverse, source_atom->position, source_frac);

      for (guint op_index = 0; op_index < ops->len; op_index++)
        {
          const GdisSymmetryOp *op;
          gdouble frac[3];
          gdouble cart[3];
          gboolean duplicate;

          op = &g_array_index(ops, GdisSymmetryOp, op_index);
          frac[0] = op->matrix[0] * source_frac[0] +
                    op->matrix[1] * source_frac[1] +
                    op->matrix[2] * source_frac[2] +
                    op->offset[0];
          frac[1] = op->matrix[3] * source_frac[0] +
                    op->matrix[4] * source_frac[1] +
                    op->matrix[5] * source_frac[2] +
                    op->offset[1];
          frac[2] = op->matrix[6] * source_frac[0] +
                    op->matrix[7] * source_frac[1] +
                    op->matrix[8] * source_frac[2] +
                    op->offset[2];
          gdis_wrap_fractional_vector(frac, model->periodicity);

          duplicate = FALSE;
          for (guint existing_index = 0; existing_index < expanded_atoms->len; existing_index++)
            {
              const GdisAtom *existing_atom;
              const GdisFractionalPosition *existing_frac;

              existing_atom = g_ptr_array_index(expanded_atoms, existing_index);
              existing_frac = &g_array_index(fractional_positions, GdisFractionalPosition, existing_index);

              if (existing_atom->region != source_atom->region)
                continue;
              if (g_strcmp0(existing_atom->label, source_atom->label) != 0)
                continue;
              if (g_strcmp0(existing_atom->ff_type, source_atom->ff_type) != 0)
                continue;
              if (!gdis_fractional_positions_equal(existing_frac->values,
                                                   frac,
                                                   model->periodicity,
                                                   5.0e-5))
                continue;

              duplicate = TRUE;
              break;
            }

          if (duplicate)
            continue;

          gdis_frac_to_cart(matrix, frac, cart);
          g_ptr_array_add(expanded_atoms,
                          gdis_atom_clone_with_position(source_atom,
                                                        cart,
                                                        expanded_atoms->len + 1u));

          {
            GdisFractionalPosition stored_frac;

            stored_frac.values[0] = frac[0];
            stored_frac.values[1] = frac[1];
            stored_frac.values[2] = frac[2];
            g_array_append_val(fractional_positions, stored_frac);
          }
        }
    }

  if (expanded_atoms->len >= original_count && expanded_atoms->len > 0u)
    {
      g_ptr_array_unref(model->atoms);
      model->atoms = expanded_atoms;
      expanded = TRUE;
      expanded_atoms = NULL;
    }

  if (expanded_atoms)
    g_ptr_array_unref(expanded_atoms);
  if (fractional_positions)
    g_array_unref(fractional_positions);

  return expanded;
}

gboolean
gdis_model_expand_space_group_atoms(GdisModel *model,
                                    gint origin_choice)
{
  GArray *ops;
  gchar *canonical_label;
  gboolean expanded;

  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->periodic ||
      model->periodicity < 3u ||
      !model->space_group ||
      model->space_group[0] == '\0')
    return FALSE;

  ops = NULL;
  canonical_label = NULL;
  if (!gdis_build_space_group_operations(model,
                                         model->space_group,
                                         origin_choice,
                                         &ops,
                                         &canonical_label))
    return FALSE;

  expanded = gdis_model_expand_atoms_with_symmetry_ops(model, ops);

  if (canonical_label && canonical_label[0] != '\0')
    {
      g_clear_pointer(&model->space_group, g_free);
      model->space_group = canonical_label;
      canonical_label = NULL;
    }

  if (ops)
    g_array_unref(ops);
  g_free(canonical_label);

  return expanded;
}

static gboolean
gdis_model_expand_space_group_atoms_best_match(GdisModel *model,
                                               guint expected_atom_count)
{
  static const gint choices[] = { 0, 1, 2 };
  GdisModel *best_model;
  guint best_delta;

  g_return_val_if_fail(model != NULL, FALSE);

  if (expected_atom_count == 0u)
    return FALSE;

  best_model = NULL;
  best_delta = G_MAXUINT;

  for (guint i = 0; i < G_N_ELEMENTS(choices); i++)
    {
      GdisModel *trial;
      guint trial_atom_count;
      guint delta;

      trial = gdis_model_clone(model);
      if (!trial)
        continue;

      if (!gdis_model_expand_space_group_atoms(trial, choices[i]))
        {
          gdis_model_free(trial);
          continue;
        }

      trial_atom_count = trial->atoms ? trial->atoms->len : 0u;
      if (trial_atom_count > expected_atom_count)
        delta = trial_atom_count - expected_atom_count;
      else
        delta = expected_atom_count - trial_atom_count;

      if (!best_model ||
          delta < best_delta ||
          (delta == best_delta && choices[i] == 0))
        {
          if (best_model)
            gdis_model_free(best_model);
          best_model = trial;
          best_delta = delta;
        }
      else
        {
          gdis_model_free(trial);
        }
    }

  if (!best_model)
    return FALSE;

  {
    GError *copy_error;
    gboolean copied;

    copy_error = NULL;
    copied = gdis_model_copy_from(model, best_model, &copy_error);
    g_clear_error(&copy_error);
    gdis_model_free(best_model);
    return copied;
  }
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

  {
    gchar *slash;

    slash = strchr(clean, '/');
    if (slash &&
        slash != clean &&
        slash[1] != '\0' &&
        strchr(slash + 1, '/') == NULL)
      {
        gdouble numerator;
        gdouble denominator;

        *slash = '\0';
        if (gdis_try_parse_double(clean, &numerator) &&
            gdis_try_parse_double(slash + 1, &denominator) &&
            fabs(denominator) > 1.0e-12)
          {
            *value = numerator / denominator;
            return TRUE;
          }
        *slash = '/';
      }
  }

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

static gboolean
gdis_gulp_type_is_core(const gchar *text)
{
  g_autofree gchar *lower = NULL;

  if (!text || text[0] == '\0')
    return FALSE;

  lower = g_ascii_strdown(text, -1);
  return g_str_equal(lower, "core") ||
         g_str_equal(lower, "c") ||
         g_str_has_prefix(lower, "bcor");
}

static gboolean
gdis_gulp_type_is_shell(const gchar *text)
{
  g_autofree gchar *lower = NULL;

  if (!text || text[0] == '\0')
    return FALSE;

  lower = g_ascii_strdown(text, -1);
  return g_str_has_prefix(lower, "shel") ||
         g_str_equal(lower, "s") ||
         g_str_has_prefix(lower, "bshe");
}

static gchar *
gdis_gulp_make_species_key(const gchar *name, const gchar *type_text)
{
  g_autofree gchar *name_lower = NULL;

  name_lower = g_ascii_strdown(name ? name : "", -1);
  return g_strdup_printf("%s|%s",
                         name_lower,
                         gdis_gulp_type_is_shell(type_text) ? "shell" : "core");
}

static const gdouble *
gdis_gulp_lookup_species_charge(GHashTable *species_charges,
                                const gchar *label,
                                const gchar *element,
                                const gchar *type_text)
{
  const gdouble *mapped_charge = NULL;
  g_autofree gchar *exact_key = NULL;
  g_autofree gchar *element_key = NULL;

  if (!species_charges)
    return NULL;

  if (label && label[0] != '\0')
    {
      exact_key = gdis_gulp_make_species_key(label, type_text);
      mapped_charge = g_hash_table_lookup(species_charges, exact_key);
    }

  if (!mapped_charge && element && element[0] != '\0')
    {
      element_key = gdis_gulp_make_species_key(element, type_text);
      mapped_charge = g_hash_table_lookup(species_charges, element_key);
    }

  return mapped_charge;
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

static guint
gdis_collect_doubles_from_text_scan(const gchar *text, gdouble *values, guint max_values)
{
  const gchar *cursor;
  guint found;

  if (!text || !values || max_values == 0u)
    return 0u;

  cursor = text;
  found = 0u;
  while (*cursor != '\0' && found < max_values)
    {
      gchar *endptr = NULL;
      gdouble parsed;

      if (!(g_ascii_isdigit(*cursor) ||
            *cursor == '+' ||
            *cursor == '-' ||
            *cursor == '.'))
        {
          cursor++;
          continue;
        }

      parsed = g_ascii_strtod(cursor, &endptr);
      if (endptr != cursor)
        {
          values[found++] = parsed;
          cursor = endptr;
          continue;
        }

      cursor++;
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

static guint
gdis_model_parse_vasp_atom_symbols(const gchar *contents, GPtrArray *symbols_out)
{
  gchar **lines;
  guint parsed_count;

  g_return_val_if_fail(symbols_out != NULL, 0u);

  parsed_count = 0u;
  lines = g_strsplit(contents ? contents : "", "\n", -1);
  for (guint i = 0; lines[i] != NULL; i++)
    {
      gchar *trimmed;

      trimmed = g_strstrip(lines[i]);
      if (!g_str_has_prefix(trimmed, "<array name=\"atoms\""))
        continue;

      for (i = i + 1u; lines[i] != NULL; i++)
        {
          gchar *row;
          const gchar *cell_start;
          const gchar *cell_end;
          g_autofree gchar *symbol = NULL;

          row = g_strstrip(lines[i]);
          if (g_str_has_prefix(row, "</array>"))
            break;
          if (!g_strstr_len(row, -1, "<rc><c>"))
            continue;

          cell_start = strstr(row, "<c>");
          cell_end = cell_start ? strstr(cell_start + 3, "</c>") : NULL;
          if (!cell_start || !cell_end || cell_end <= cell_start + 3)
            continue;

          symbol = g_strndup(cell_start + 3, (gsize) (cell_end - (cell_start + 3)));
          g_strstrip(symbol);
          if (!symbol[0])
            continue;

          symbol = gdis_normalize_element_symbol(symbol);
          if (symbol && symbol[0])
            {
              g_ptr_array_add(symbols_out, g_steal_pointer(&symbol));
              parsed_count++;
            }
        }
      break;
    }

  g_strfreev(lines);
  return parsed_count;
}

static gboolean
gdis_model_parse_vasp_structure_block(GdisModel *model,
                                      gchar **lines,
                                      guint start_index,
                                      const GPtrArray *atom_symbols,
                                      GPtrArray **atoms_out,
                                      gdouble lengths_out[3],
                                      gdouble angles_out[3],
                                      guint *end_index_out)
{
  GPtrArray *atoms;
  gdouble basis[9];
  gboolean have_basis;
  gboolean in_positions;
  guint atom_index;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(lines != NULL, FALSE);
  g_return_val_if_fail(atoms_out != NULL, FALSE);
  g_return_val_if_fail(lengths_out != NULL, FALSE);
  g_return_val_if_fail(angles_out != NULL, FALSE);

  atoms = g_ptr_array_new_with_free_func(gdis_atom_free);
  memset(basis, 0, sizeof(basis));
  lengths_out[0] = lengths_out[1] = lengths_out[2] = 0.0;
  angles_out[0] = angles_out[1] = angles_out[2] = 90.0;
  if (end_index_out)
    *end_index_out = start_index;
  have_basis = FALSE;
  in_positions = FALSE;
  atom_index = 0u;

  for (guint i = start_index + 1u; lines[i] != NULL; i++)
    {
      gchar *trimmed;

      trimmed = g_strstrip(lines[i]);
      if (trimmed[0] == '\0')
        continue;

      if (g_str_has_prefix(trimmed, "</structure>"))
        {
          if (end_index_out)
            *end_index_out = i;
          break;
        }

      if (g_str_has_prefix(trimmed, "<varray name=\"basis\""))
        {
          gdouble values[3];
          guint row;

          row = 0u;
          for (i = i + 1u; lines[i] != NULL && row < 3u; i++)
            {
              gchar *vec_line;

              vec_line = g_strstrip(lines[i]);
              if (g_str_has_prefix(vec_line, "</varray>"))
                break;
              if (gdis_collect_doubles_from_text_scan(vec_line, values, 3u) < 3u)
                continue;

              basis[row] = values[0];
              basis[row + 3u] = values[1];
              basis[row + 6u] = values[2];
              row++;
            }

          if (row == 3u &&
              gdis_surface_build_cell_from_vectors((gdouble[3]) {basis[0], basis[3], basis[6]},
                                                   (gdouble[3]) {basis[1], basis[4], basis[7]},
                                                   (gdouble[3]) {basis[2], basis[5], basis[8]},
                                                   lengths_out,
                                                   angles_out))
            have_basis = TRUE;
          if (lines[i] == NULL)
            break;
          continue;
        }

      if (g_str_has_prefix(trimmed, "<varray name=\"positions\""))
        {
          in_positions = TRUE;
          continue;
        }

      if (in_positions)
        {
          gdouble frac[3];
          gdouble cart[3];
          g_autofree gchar *element = NULL;
          const gchar *label_text;

          if (g_str_has_prefix(trimmed, "</varray>"))
            {
              in_positions = FALSE;
              continue;
            }

          if (!have_basis || gdis_collect_doubles_from_text_scan(trimmed, frac, 3u) < 3u)
            continue;

          gdis_frac_to_cart(basis, frac, cart);
          if (atom_symbols && atom_index < atom_symbols->len)
            label_text = g_ptr_array_index((GPtrArray *) atom_symbols, atom_index);
          else
            label_text = "X";
          element = gdis_normalize_element_symbol(label_text);
          if (!element || !element[0])
            {
              g_clear_pointer(&element, g_free);
              element = g_strdup("X");
            }

          g_ptr_array_add(atoms,
                          gdis_atom_new(label_text,
                                        element,
                                        element,
                                        cart[0],
                                        cart[1],
                                        cart[2],
                                        1.0,
                                        -1,
                                        atom_index + 1u));
          atom_index++;
        }
    }

  if (!have_basis || atoms->len == 0u)
    {
      g_ptr_array_free(atoms, TRUE);
      return FALSE;
    }

  *atoms_out = atoms;
  return TRUE;
}

static gboolean
gdis_model_build_band_arrays_from_triplets(const GArray *kpoint_triplets,
                                           GPtrArray *kpoint_bands,
                                           gdouble fermi_energy_ev,
                                           GArray **x_values_out,
                                           GArray **y_values_out,
                                           guint *path_count_out,
                                           guint *series_count_out)
{
  GArray *kpath;
  GArray *x_values;
  GArray *y_values;
  guint kpoint_count;
  guint band_count;
  gdouble cumulative;

  g_return_val_if_fail(x_values_out != NULL, FALSE);
  g_return_val_if_fail(y_values_out != NULL, FALSE);

  *x_values_out = NULL;
  *y_values_out = NULL;
  if (path_count_out)
    *path_count_out = 0u;
  if (series_count_out)
    *series_count_out = 0u;

  if (!kpoint_triplets || !kpoint_bands || kpoint_bands->len == 0u)
    return FALSE;

  kpoint_count = kpoint_bands->len;
  band_count = G_MAXUINT;
  for (guint i = 0; i < kpoint_count; i++)
    {
      GArray *band_row;

      band_row = g_ptr_array_index(kpoint_bands, i);
      if (!band_row || band_row->len == 0u)
        return FALSE;
      band_count = MIN(band_count, band_row->len);
    }
  if (band_count == 0u || band_count == G_MAXUINT)
    return FALSE;

  kpath = gdis_double_array_new();
  cumulative = 0.0;
  for (guint i = 0; i < kpoint_count; i++)
    {
      gdouble path_value;

      if (i == 0u || !kpoint_triplets || kpoint_triplets->len < (i + 1u) * 3u)
        {
          path_value = (gdouble) i;
        }
      else
        {
          const gdouble *coords = &g_array_index(kpoint_triplets, gdouble, i * 3u);
          const gdouble *prev = &g_array_index(kpoint_triplets, gdouble, (i - 1u) * 3u);
          const gdouble dx = coords[0] - prev[0];
          const gdouble dy = coords[1] - prev[1];
          const gdouble dz = coords[2] - prev[2];

          cumulative += sqrt(dx * dx + dy * dy + dz * dz);
          path_value = cumulative;
        }
      g_array_append_val(kpath, path_value);
    }

  x_values = gdis_double_array_new();
  y_values = gdis_double_array_new();
  for (guint band_index = 0; band_index < band_count; band_index++)
    {
      for (guint kpoint_index = 0; kpoint_index < kpoint_count; kpoint_index++)
        {
          GArray *band_row;
          gdouble x_value;
          gdouble y_value;

          band_row = g_ptr_array_index(kpoint_bands, kpoint_index);
          x_value = g_array_index(kpath, gdouble, kpoint_index);
          y_value = g_array_index(band_row, gdouble, band_index) - fermi_energy_ev;
          g_array_append_val(x_values, x_value);
          g_array_append_val(y_values, y_value);
        }
      if (band_index + 1u < band_count)
        {
          gdouble gap_x;
          gdouble gap_y;

          gap_x = g_array_index(kpath, gdouble, kpoint_count - 1u);
          gap_y = NAN;
          g_array_append_val(x_values, gap_x);
          g_array_append_val(y_values, gap_y);
        }
    }

  *x_values_out = x_values;
  *y_values_out = y_values;
  if (path_count_out)
    *path_count_out = kpoint_count;
  if (series_count_out)
    *series_count_out = band_count;
  g_array_unref(kpath);
  return TRUE;
}

static gboolean
gdis_model_parse_vasp_dos_block(gchar **lines,
                                guint start_index,
                                GArray **dos_x_values_out,
                                GArray **dos_y_values_out,
                                gboolean *have_fermi_energy_out,
                                gdouble *fermi_energy_ev_out,
                                guint *end_index_out)
{
  GArray *energies;
  GArray *totals;
  gboolean have_fermi_energy;
  gdouble fermi_energy_ev;
  guint spin_index;

  g_return_val_if_fail(lines != NULL, FALSE);
  g_return_val_if_fail(dos_x_values_out != NULL, FALSE);
  g_return_val_if_fail(dos_y_values_out != NULL, FALSE);

  *dos_x_values_out = NULL;
  *dos_y_values_out = NULL;
  if (end_index_out)
    *end_index_out = start_index;

  energies = gdis_double_array_new();
  totals = gdis_double_array_new();
  have_fermi_energy = have_fermi_energy_out ? *have_fermi_energy_out : FALSE;
  fermi_energy_ev = fermi_energy_ev_out ? *fermi_energy_ev_out : 0.0;
  spin_index = 0u;

  for (guint i = start_index + 1u; lines[i] != NULL; i++)
    {
      gchar *dos_line;

      dos_line = g_strstrip(lines[i]);
      if (g_str_has_prefix(dos_line, "</dos>"))
        {
          if (end_index_out)
            *end_index_out = i;
          break;
        }
      if (g_strstr_len(dos_line, -1, "name=\"efermi\"") != NULL)
        {
          gdouble values[1];

          if (gdis_collect_doubles_from_text_scan(dos_line, values, 1u) >= 1u)
            {
              fermi_energy_ev = values[0];
              have_fermi_energy = TRUE;
            }
          continue;
        }
      if (g_str_has_prefix(dos_line, "<set comment=\"spin "))
        {
          guint dos_index = 0u;

          spin_index++;
          for (i = i + 1u; lines[i] != NULL; i++)
            {
              gchar *row_line;
              gdouble values[3];

              row_line = g_strstrip(lines[i]);
              if (g_str_has_prefix(row_line, "</set>"))
                break;
              if (gdis_collect_doubles_from_text_scan(row_line, values, 3u) < 2u)
                continue;

              if (spin_index == 1u)
                {
                  g_array_append_val(energies, values[0]);
                  g_array_append_val(totals, values[1]);
                }
              else if (dos_index < totals->len)
                {
                  gdouble sum_value;

                  sum_value = g_array_index(totals, gdouble, dos_index) + values[1];
                  g_array_index(totals, gdouble, dos_index) = sum_value;
                }
              dos_index++;
            }
        }
    }

  if (energies->len == 0u || totals->len != energies->len)
    {
      g_array_unref(energies);
      g_array_unref(totals);
      return FALSE;
    }

  for (guint idx = 0; idx < energies->len; idx++)
    {
      gdouble shifted;

      shifted = g_array_index(energies, gdouble, idx) -
                (have_fermi_energy ? fermi_energy_ev : 0.0);
      g_array_index(energies, gdouble, idx) = shifted;
    }

  if (have_fermi_energy_out)
    *have_fermi_energy_out = have_fermi_energy;
  if (fermi_energy_ev_out)
    *fermi_energy_ev_out = fermi_energy_ev;
  *dos_x_values_out = energies;
  *dos_y_values_out = totals;
  return TRUE;
}

static gboolean
gdis_model_parse_vasp_eigenvalues_block(gchar **lines,
                                        guint start_index,
                                        const GArray *kpoint_triplets,
                                        gboolean have_fermi_energy,
                                        gdouble fermi_energy_ev,
                                        GArray **band_x_values_out,
                                        GArray **band_y_values_out,
                                        guint *path_count_out,
                                        guint *series_count_out,
                                        guint *end_index_out)
{
  GPtrArray *kpoint_bands;
  gboolean in_spin_one;
  gboolean ok;

  g_return_val_if_fail(lines != NULL, FALSE);
  g_return_val_if_fail(band_x_values_out != NULL, FALSE);
  g_return_val_if_fail(band_y_values_out != NULL, FALSE);

  *band_x_values_out = NULL;
  *band_y_values_out = NULL;
  if (path_count_out)
    *path_count_out = 0u;
  if (series_count_out)
    *series_count_out = 0u;
  if (end_index_out)
    *end_index_out = start_index;

  kpoint_bands = g_ptr_array_new_with_free_func((GDestroyNotify) g_array_unref);
  in_spin_one = FALSE;
  ok = FALSE;

  for (guint i = start_index + 1u; lines[i] != NULL; i++)
    {
      gchar *band_line;

      band_line = g_strstrip(lines[i]);
      if (g_str_has_prefix(band_line, "<set comment=\"spin 1\""))
        {
          in_spin_one = TRUE;
          continue;
        }
      if (g_str_has_prefix(band_line, "<set comment=\"spin 2\"") ||
          g_str_has_prefix(band_line, "</eigenvalues>"))
        {
          if (end_index_out)
            *end_index_out = i;
          break;
        }
      if (!in_spin_one)
        continue;

      if (g_str_has_prefix(band_line, "<set comment=\"kpoint "))
        {
          GArray *band_values;

          band_values = gdis_double_array_new();
          for (i = i + 1u; lines[i] != NULL; i++)
            {
              gchar *row_line;
              gdouble values[2];

              row_line = g_strstrip(lines[i]);
              if (g_str_has_prefix(row_line, "</set>"))
                break;
              if (gdis_collect_doubles_from_text_scan(row_line, values, 2u) < 1u)
                continue;
              g_array_append_val(band_values, values[0]);
            }

          if (band_values->len > 0u)
            g_ptr_array_add(kpoint_bands, band_values);
          else
            g_array_unref(band_values);
        }
    }

  if (kpoint_bands->len > 0u)
    {
      ok = gdis_model_build_band_arrays_from_triplets(kpoint_triplets,
                                                      kpoint_bands,
                                                      have_fermi_energy ? fermi_energy_ev : 0.0,
                                                      band_x_values_out,
                                                      band_y_values_out,
                                                      path_count_out,
                                                      series_count_out);
    }

  g_ptr_array_free(kpoint_bands, TRUE);
  return ok;
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

static gboolean
gdis_parse_cif_symmetry_operation(const gchar *text, GdisSymmetryOp *op_out)
{
  g_auto(GStrv) components = NULL;
  g_autofree gchar *normalized = NULL;
  GdisSymmetryOp op;
  gsize len;

  g_return_val_if_fail(op_out != NULL, FALSE);

  if (!text)
    return FALSE;

  normalized = gdis_strdup_strip(text);
  if (!normalized || normalized[0] == '\0')
    return FALSE;

  len = strlen(normalized);
  if (len >= 2 &&
      ((normalized[0] == '\'' && normalized[len - 1] == '\'') ||
       (normalized[0] == '"' && normalized[len - 1] == '"')))
    {
      normalized[len - 1] = '\0';
      memmove(normalized, normalized + 1, len - 1);
    }

  components = g_strsplit(normalized, ",", -1);
  if (!components || !components[0] || !components[1] || !components[2] || components[3] != NULL)
    return FALSE;

  memset(&op, 0, sizeof(op));
  for (guint row = 0; row < 3u; row++)
    {
      gchar *component;
      const gchar *cursor;
      gint sign;

      component = g_strstrip(components[row]);
      if (!component || component[0] == '\0')
        return FALSE;

      cursor = component;
      sign = 1;
      while (*cursor)
        {
          if (g_ascii_isspace(*cursor))
            {
              cursor++;
              continue;
            }

          if (*cursor == '+')
            {
              sign = 1;
              cursor++;
              continue;
            }
          if (*cursor == '-')
            {
              sign = -1;
              cursor++;
              continue;
            }

          switch (g_ascii_tolower(*cursor))
            {
            case 'x':
              op.matrix[row * 3u] += (gdouble) sign;
              sign = 1;
              cursor++;
              continue;
            case 'y':
              op.matrix[row * 3u + 1u] += (gdouble) sign;
              sign = 1;
              cursor++;
              continue;
            case 'z':
              op.matrix[row * 3u + 2u] += (gdouble) sign;
              sign = 1;
              cursor++;
              continue;
            default:
              break;
            }

          if (g_ascii_isdigit(*cursor) || *cursor == '.')
            {
              gchar *endptr;
              gdouble numerator;
              gdouble value;

              numerator = g_ascii_strtod(cursor, &endptr);
              if (endptr == cursor)
                return FALSE;

              value = numerator;
              if (*endptr == '/')
                {
                  gchar *den_endptr;
                  gdouble denominator;

                  denominator = g_ascii_strtod(endptr + 1, &den_endptr);
                  if (den_endptr == endptr + 1 || fabs(denominator) <= 1.0e-12)
                    return FALSE;
                  value = numerator / denominator;
                  endptr = den_endptr;
                }

              op.offset[row] += (gdouble) sign * value;
              sign = 1;
              cursor = endptr;
              continue;
            }

          return FALSE;
        }
    }

  *op_out = op;
  return TRUE;
}

static gdouble
gdis_lookup_covalent_radius(const gchar *element)
{
  static const GdisRadiusEntry radii[] = {
    { "X", 0.100000 },
    { "H", 0.330000 },
    { "He", 0.460000 },
    { "Li", 0.680000 },
    { "Be", 0.350000 },
    { "B", 0.820000 },
    { "C", 0.770000 },
    { "N", 0.750000 },
    { "O", 0.730000 },
    { "F", 0.720000 },
    { "Ne", 0.670000 },
    { "Na", 0.970000 },
    { "Mg", 1.100000 },
    { "Al", 1.180000 },
    { "Si", 1.110000 },
    { "P", 1.060000 },
    { "S", 1.020000 },
    { "Cl", 0.990000 },
    { "Ar", 1.060000 },
    { "K", 1.330000 },
    { "Ca", 0.990000 },
    { "Sc", 1.440000 },
    { "Ti", 1.470000 },
    { "V", 1.330000 },
    { "Cr", 1.350000 },
    { "Mn", 1.350000 },
    { "Fe", 1.340000 },
    { "Co", 1.330000 },
    { "Ni", 1.500000 },
    { "Cu", 1.520000 },
    { "Zn", 1.450000 },
    { "Ga", 1.220000 },
    { "Ge", 1.170000 },
    { "As", 1.210000 },
    { "Se", 1.220000 },
    { "Br", 1.210000 },
    { "Kr", 1.890000 },
    { "Rb", 1.470000 },
    { "Sr", 1.120000 },
    { "Y", 1.780000 },
    { "Zr", 1.560000 },
    { "Nb", 1.480000 },
    { "Mo", 1.470000 },
    { "Tc", 1.350000 },
    { "Ru", 1.400000 },
    { "Rh", 1.450000 },
    { "Pd", 1.500000 },
    { "Ag", 1.590000 },
    { "Cd", 1.690000 },
    { "In", 1.630000 },
    { "Sn", 1.460000 },
    { "Sb", 1.460000 },
    { "Te", 1.470000 },
    { "I", 1.330000 },
    { "Xe", 1.310000 },
    { "Cs", 1.670000 },
    { "Ba", 1.340000 },
    { "La", 1.870000 },
    { "Ce", 1.830000 },
    { "Pr", 1.820000 },
    { "Nd", 1.810000 },
    { "Pm", 1.800000 },
    { "Sm", 1.800000 },
    { "Eu", 1.990000 },
    { "Gd", 1.790000 },
    { "Tb", 1.760000 },
    { "Dy", 1.750000 },
    { "Ho", 1.740000 },
    { "Er", 1.730000 },
    { "Tm", 1.720000 },
    { "Yb", 1.940000 },
    { "Lu", 1.720000 },
    { "Hf", 1.570000 },
    { "Ta", 1.430000 },
    { "W", 1.370000 },
    { "Re", 1.350000 },
    { "Os", 1.370000 },
    { "Ir", 1.320000 },
    { "Pt", 1.500000 },
    { "Au", 1.500000 },
    { "Hg", 1.700000 },
    { "Tl", 1.550000 },
    { "Pb", 1.540000 },
    { "Bi", 1.540000 },
    { "Po", 1.680000 },
    { "At", 1.470000 },
    { "Rn", 1.420000 },
    { "Fr", 2.600000 },
    { "Ra", 1.900000 },
    { "Ac", 1.880000 },
    { "Th", 1.790000 },
    { "Pa", 1.610000 },
    { "U", 1.580000 },
    { "Np", 1.550000 },
    { "Pu", 1.530000 },
    { "Am", 1.510000 },
    { "Cm", 1.690000 },
    { "Bk", 2.600000 },
    { "Cf", 2.600000 },
    { "Es", 2.600000 },
    { "Fm", 2.600000 },
    { "Md", 2.600000 },
    { "No", 2.600000 },
    { "Lr", 2.600000 },
    { "Rf", 2.600000 },
    { "Db", 2.600000 },
    { "Sg", 2.600000 },
    { "Bh", 2.600000 },
    { "Hs", 2.600000 },
    { "Mt", 2.600000 }
  };
  guint i;

  if (!element || element[0] == '\0')
    return 0.10;

  for (i = 0; i < G_N_ELEMENTS(radii); i++)
    {
      if (g_ascii_strcasecmp(radii[i].symbol, element) == 0)
        return radii[i].radius;
    }

  return 0.77;
}
