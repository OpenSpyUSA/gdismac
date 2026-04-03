#pragma once

#include <glib.h>

#include "gdis_model.h"

G_BEGIN_DECLS

typedef enum
{
  GDIS_MEASURE_MODE_AUTO = 0,
  GDIS_MEASURE_MODE_DISTANCE,
  GDIS_MEASURE_MODE_ANGLE,
  GDIS_MEASURE_MODE_TORSION
} GdisMeasureMode;

typedef enum
{
  GDIS_DIFFRACTION_RADIATION_XRAY = 0,
  GDIS_DIFFRACTION_RADIATION_NEUTRON,
  GDIS_DIFFRACTION_RADIATION_ELECTRON
} GdisDiffractionRadiation;

typedef enum
{
  GDIS_DIFFRACTION_BROADENING_GAUSSIAN = 0,
  GDIS_DIFFRACTION_BROADENING_LORENTZIAN,
  GDIS_DIFFRACTION_BROADENING_PSEUDO_VOIGT
} GdisDiffractionBroadening;

typedef struct
{
  GdisDiffractionRadiation radiation;
  GdisDiffractionBroadening broadening;
  gdouble wavelength;
  gdouble theta_min;
  gdouble theta_max;
  gdouble theta_step;
  gdouble asym;
  gdouble u;
  gdouble v;
  gdouble w;
} GdisDiffractionSettings;

typedef struct
{
  gint h;
  gint k;
  gint l;
  guint reflections_merged;
  gdouble d_spacing;
  gdouble two_theta;
  gdouble intensity;
  gdouble relative_intensity;
} GdisDiffractionPeak;

typedef struct
{
  GdisDiffractionSettings settings;
  GArray *x_values;
  GArray *y_values;
  GArray *peaks;
  gdouble max_intensity;
} GdisDiffractionPattern;

char *gdis_report_measurements(const GdisModel *model,
                               const guint *atom_indices,
                               guint count,
                               GdisMeasureMode mode);
char *gdis_report_surface(const GdisModel *model);
char *gdis_report_diffraction(const GdisModel *model);
char *gdis_report_isosurface(const GdisModel *model);
void gdis_diffraction_settings_init(GdisDiffractionSettings *settings);
const char *gdis_diffraction_radiation_label(GdisDiffractionRadiation radiation);
const char *gdis_diffraction_broadening_label(GdisDiffractionBroadening broadening);
gboolean gdis_diffraction_calculate(const GdisModel *model,
                                    const GdisDiffractionSettings *settings,
                                    GdisDiffractionPattern **pattern_out,
                                    GError **error);
gboolean gdis_diffraction_export(const GdisDiffractionPattern *pattern,
                                 const char *path,
                                 gboolean append,
                                 GError **error);
void gdis_diffraction_pattern_free(GdisDiffractionPattern *pattern);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GdisDiffractionPattern, gdis_diffraction_pattern_free)

G_END_DECLS
