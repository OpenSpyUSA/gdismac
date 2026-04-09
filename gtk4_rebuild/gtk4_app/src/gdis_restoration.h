#pragma once

#include <glib.h>

#include "gdis_model.h"

G_BEGIN_DECLS

typedef enum
{
  GDIS_DISLOCATION_TYPE_SCREW = 0,
  GDIS_DISLOCATION_TYPE_EDGE
} GdisDislocationType;

typedef struct
{
  GdisDislocationType type;
  gdouble line_direction[3];
  gdouble burgers_vector[3];
  gdouble origin[3];
  gdouble radius;
  gdouble poisson_ratio;
  const GArray *selected_atoms;
} GdisDislocationSettings;

typedef struct
{
  const gchar *output_dir;
  const gchar *project_name;
  guint grid_counts[3];
  guint rotation_counts[3];
  gdouble translation_span[3];
  guint preview_limit;
  const GArray *selected_atoms;
} GdisDockingSettings;

typedef struct
{
  guint box_dim;
  guint solvent_index;
  const guint *component_counts;
  guint component_count;
  gboolean random_rotate;
} GdisMdiSettings;

typedef struct
{
  guint atom_index;
  guint distance_ref;
  guint angle_ref;
  guint torsion_ref;
  gdouble distance;
  gdouble angle;
  gdouble torsion;
} GdisZmatrixRow;

void gdis_dislocation_settings_init(GdisDislocationSettings *settings);
const char *gdis_dislocation_type_label(GdisDislocationType type);
gboolean gdis_model_apply_dislocation(GdisModel *model,
                                      const GdisDislocationSettings *settings,
                                      gchar **summary_out,
                                      GError **error);

void gdis_docking_settings_init(GdisDockingSettings *settings);
gboolean gdis_generate_docking_project(const GdisModel *model,
                                       const GdisDockingSettings *settings,
                                       gchar **summary_out,
                                       GError **error);

void gdis_mdi_settings_init(GdisMdiSettings *settings);
gboolean gdis_generate_mdi_model(const GPtrArray *source_models,
                                 const GdisMdiSettings *settings,
                                 GdisModel **model_out,
                                 gchar **summary_out,
                                 GError **error);

gboolean gdis_build_zmatrix_rows(const GdisModel *model,
                                 const GArray *selected_atoms,
                                 gboolean use_selection,
                                 GArray **scope_out,
                                 GArray **rows_out,
                                 GError **error);
char *gdis_format_zmatrix_rows(const GdisModel *model,
                               const GArray *scope,
                               const GArray *rows,
                               gboolean use_selection);
gboolean gdis_apply_zmatrix_rows(GdisModel *model,
                                 const GArray *scope,
                                 const GArray *rows,
                                 gchar **summary_out,
                                 GError **error);
char *gdis_build_zmatrix_report(const GdisModel *model,
                                const GArray *selected_atoms,
                                gboolean use_selection,
                                GError **error);

G_END_DECLS
