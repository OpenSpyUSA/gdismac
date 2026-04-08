#pragma once

#include <glib.h>

G_BEGIN_DECLS

typedef enum
{
  GDIS_MODEL_FORMAT_UNKNOWN = 0,
  GDIS_MODEL_FORMAT_XYZ,
  GDIS_MODEL_FORMAT_PDB,
  GDIS_MODEL_FORMAT_ARC,
  GDIS_MODEL_FORMAT_CAR,
  GDIS_MODEL_FORMAT_CIF,
  GDIS_MODEL_FORMAT_QBOX_XML,
  GDIS_MODEL_FORMAT_GULP_INPUT,
  GDIS_MODEL_FORMAT_GULP_OUTPUT,
  GDIS_MODEL_FORMAT_XTL,
  GDIS_MODEL_FORMAT_QE_INPUT,
  GDIS_MODEL_FORMAT_QE_OUTPUT,
  GDIS_MODEL_FORMAT_GAMESS_INPUT
} GdisModelFormat;

typedef enum
{
  GDIS_MODEL_ERROR_FAILED,
  GDIS_MODEL_ERROR_IO,
  GDIS_MODEL_ERROR_UNSUPPORTED_FORMAT,
  GDIS_MODEL_ERROR_PARSE
} GdisModelError;

typedef struct
{
  guint serial;
  gchar *label;
  gchar *element;
  gchar *ff_type;
  gdouble position[3];
  gdouble occupancy;
  gint region;
} GdisAtom;

typedef struct
{
  guint atom_index_a;
  guint atom_index_b;
  guint8 order;
  gboolean inferred;
} GdisBond;

typedef struct
{
  gchar *path;
  gchar *basename;
  gchar *title;
  gchar *format_label;
  gchar *space_group;

  GdisModelFormat format;

  gboolean periodic;
  guint periodicity;

  gdouble cell_lengths[3];
  gdouble cell_angles[3];
  gint image_limits[6];

  guint atom_count;
  guint bond_count;
  guint explicit_bond_count;

  GPtrArray *atoms;
  GArray *bonds;
  GPtrArray *frames;
  guint current_frame_index;
} GdisModel;

GQuark gdis_model_error_quark(void);

#define GDIS_MODEL_ERROR (gdis_model_error_quark())

GdisModel *gdis_model_create_empty(const char *path, GdisModelFormat format);
GdisModel *gdis_model_load(const char *path, GError **error);
GdisModel *gdis_model_clone(const GdisModel *model);
gboolean gdis_model_copy_from(GdisModel *dest,
                              const GdisModel *src,
                              GError **error);
void gdis_model_free(GdisModel *model);
gboolean gdis_model_save(GdisModel *model, const char *path, GError **error);
gboolean gdis_model_delete_atoms(GdisModel *model,
                                 const guint *atom_indices,
                                 guint count,
                                 GError **error);
gboolean gdis_model_add_atom(GdisModel *model,
                             const char *label,
                             const char *element,
                             const char *ff_type,
                             gint region,
                             gdouble x,
                             gdouble y,
                             gdouble z,
                             GError **error);
gboolean gdis_model_update_atom(GdisModel *model,
                                guint atom_index,
                                const char *label,
                                const char *element,
                                const char *ff_type,
                                gint region,
                                gdouble x,
                                gdouble y,
                                gdouble z,
                                GError **error);
gboolean gdis_model_add_explicit_bond(GdisModel *model,
                                      guint atom_index_a,
                                      guint atom_index_b,
                                      guint8 order,
                                      GError **error);
gboolean gdis_model_remove_bond(GdisModel *model,
                                guint atom_index_a,
                                guint atom_index_b,
                                GError **error);
void gdis_model_reset_image_limits(GdisModel *model);
void gdis_model_get_image_limits(const GdisModel *model, gint limits_out[6]);
gboolean gdis_model_set_image_limits(GdisModel *model,
                                     const gint limits[6],
                                     GError **error);
gboolean gdis_model_confine_atoms_to_cell(GdisModel *model, GError **error);
gboolean gdis_model_confine_molecules_to_cell(GdisModel *model, GError **error);
gboolean gdis_model_force_p1(GdisModel *model, GError **error);
gboolean gdis_model_make_supercell(GdisModel *model,
                                   guint repeat_a,
                                   guint repeat_b,
                                   guint repeat_c,
                                   GError **error);
gboolean gdis_model_make_supercell_from_image_limits(GdisModel *model,
                                                     GError **error);
gboolean gdis_model_build_surface_slab(const GdisModel *source,
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
                                       GError **error);

GdisModelFormat gdis_model_format_from_path(const char *path);
const char *gdis_model_format_label(GdisModelFormat format);

gboolean gdis_model_reset_inferred_bonds(GdisModel *model);
guint gdis_model_get_frame_count(const GdisModel *model);
guint gdis_model_get_current_frame_index(const GdisModel *model);
const char *gdis_model_get_frame_title(const GdisModel *model,
                                       guint frame_index);
gboolean gdis_model_set_frame_index(GdisModel *model,
                                    guint frame_index,
                                    GError **error);
void gdis_model_discard_frames(GdisModel *model);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GdisModel, gdis_model_free)

G_END_DECLS
