#include "gdis_reports.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct
{
  gint h;
  gint k;
  gint l;
  gdouble d_spacing;
} GdisPlaneEntry;

typedef struct
{
  gdouble values[11];
} GdisScatteringFactor;

typedef struct
{
  const char *symbol;
  guint atomic_number;
} GdisAtomicNumber;

static gboolean gdis_build_cell_vectors(const GdisModel *model,
                                        gdouble a_vec[3],
                                        gdouble b_vec[3],
                                        gdouble c_vec[3]);
static gboolean gdis_build_cell_matrix(const GdisModel *model,
                                       gdouble matrix[9],
                                       gdouble inverse[9]);
static gboolean gdis_invert_matrix3(const gdouble matrix[9], gdouble inverse[9]);
static guint gdis_atomic_number_for_element(const char *element);
static void gdis_normalize_element_symbol(const char *element, gchar normalized[4]);
static gboolean gdis_measure_prepare_periodic(const GdisModel *model,
                                              gdouble matrix[9],
                                              gdouble inverse[9]);
static void gdis_cart_to_frac(const gdouble inverse[9],
                              const gdouble cart[3],
                              gdouble frac[3]);
static void gdis_measure_minimum_image_delta(const gdouble matrix[9],
                                             const gdouble inverse[9],
                                             const gdouble from[3],
                                             const gdouble to[3],
                                             gdouble delta[3]);
static gdouble gdis_d_spacing(const gdouble inverse[9], gint h, gint k, gint l);
static gdouble gdis_diffraction_structure_factor_intensity(const GdisModel *model,
                                                           const gdouble inverse[9],
                                                           GdisDiffractionRadiation radiation,
                                                           gint h,
                                                           gint k,
                                                           gint l,
                                                           gdouble sol);
static gdouble gdis_measure_distance(const GdisModel *model,
                                     const gdouble a[3],
                                     const gdouble b[3]);
static gdouble gdis_measure_angle(const GdisModel *model,
                                  const gdouble a[3],
                                  const gdouble b[3],
                                  const gdouble c[3]);
static gdouble gdis_measure_torsion(const GdisModel *model,
                                    const gdouble a[3],
                                    const gdouble b[3],
                                    const gdouble c[3],
                                    const gdouble d[3]);
static gchar *gdis_find_support_file(const char *relative_path);
static gboolean gdis_diffraction_load_scattering_data(void);
static const GdisScatteringFactor *gdis_diffraction_lookup_scattering_factor(const char *element);
static gdouble gdis_diffraction_atomic_scattering_factor(const char *element,
                                                         GdisDiffractionRadiation radiation,
                                                         gdouble sol);
static gdouble gdis_diffraction_lorentz_polarization_factor(GdisDiffractionRadiation radiation,
                                                            gdouble two_theta);
static gdouble gdis_diffraction_peak_profile(const GdisDiffractionSettings *settings,
                                             gdouble fwhm,
                                             gdouble delta_two_theta);
static gboolean gdis_diffraction_validate_settings(const GdisDiffractionSettings *settings,
                                                   GError **error);
static gchar *gdis_diffraction_pattern_summary(const GdisModel *model,
                                               const GdisDiffractionPattern *pattern);
static gint gdis_plane_compare(gconstpointer left, gconstpointer right);
static gint gdis_peak_compare(gconstpointer left, gconstpointer right);
static void gdis_canonicalize_hkl(gint *h, gint *k, gint *l);

static GHashTable *gdis_scattering_factor_table;
static gboolean gdis_scattering_factor_table_loaded;

static const GdisAtomicNumber gdis_atomic_numbers[] = {
  {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},    {"C", 6},
  {"N", 7},   {"O", 8},   {"F", 9},   {"Ne", 10}, {"Na", 11},  {"Mg", 12},
  {"Al", 13}, {"Si", 14}, {"P", 15},  {"S", 16},  {"Cl", 17},  {"Ar", 18},
  {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23},   {"Cr", 24},
  {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29},  {"Zn", 30},
  {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35},  {"Kr", 36},
  {"Rb", 37}, {"Sr", 38}, {"Y", 39},  {"Zr", 40}, {"Nb", 41},  {"Mo", 42},
  {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50}, {"Sb", 51},  {"Te", 52},
  {"I", 53},  {"Xe", 54}, {"Cs", 55}, {"Ba", 56}, {"La", 57},  {"Ce", 58},
  {"Pr", 59}, {"Nd", 60}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64},  {"Tb", 65},
  {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70},  {"Lu", 71},
  {"Hf", 72}, {"Ta", 73}, {"W", 74},  {"Re", 75}, {"Os", 76},  {"Ir", 77},
  {"Pt", 78}, {"Au", 79}, {"Hg", 80}, {"Pb", 82}, {"Bi", 83}
};

void
gdis_diffraction_settings_init(GdisDiffractionSettings *settings)
{
  g_return_if_fail(settings != NULL);

  settings->radiation = GDIS_DIFFRACTION_RADIATION_XRAY;
  settings->broadening = GDIS_DIFFRACTION_BROADENING_GAUSSIAN;
  settings->wavelength = 1.5418;
  settings->theta_min = 0.0;
  settings->theta_max = 90.0;
  settings->theta_step = 0.1;
  settings->asym = 0.18;
  settings->u = 0.0;
  settings->v = 0.0;
  settings->w = 0.0;
}

const char *
gdis_diffraction_radiation_label(GdisDiffractionRadiation radiation)
{
  switch (radiation)
    {
    case GDIS_DIFFRACTION_RADIATION_NEUTRON:
      return "Neutrons";
    case GDIS_DIFFRACTION_RADIATION_ELECTRON:
      return "Electrons";
    case GDIS_DIFFRACTION_RADIATION_XRAY:
    default:
      return "X-Rays";
    }
}

const char *
gdis_diffraction_broadening_label(GdisDiffractionBroadening broadening)
{
  switch (broadening)
    {
    case GDIS_DIFFRACTION_BROADENING_LORENTZIAN:
      return "Lorentzian";
    case GDIS_DIFFRACTION_BROADENING_PSEUDO_VOIGT:
      return "Pseudo-Voigt";
    case GDIS_DIFFRACTION_BROADENING_GAUSSIAN:
    default:
      return "Gaussian";
    }
}

void
gdis_diffraction_pattern_free(GdisDiffractionPattern *pattern)
{
  if (!pattern)
    return;

  if (pattern->x_values)
    g_array_free(pattern->x_values, TRUE);
  if (pattern->y_values)
    g_array_free(pattern->y_values, TRUE);
  if (pattern->peaks)
    g_array_free(pattern->peaks, TRUE);
  g_free(pattern);
}

gboolean
gdis_diffraction_export(const GdisDiffractionPattern *pattern,
                        const char *path,
                        gboolean append,
                        GError **error)
{
  FILE *fp;
  guint i;

  g_return_val_if_fail(pattern != NULL, FALSE);
  g_return_val_if_fail(path != NULL, FALSE);

  if (!pattern->x_values || !pattern->y_values || pattern->x_values->len != pattern->y_values->len)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "No diffraction spectrum is available to export.");
      return FALSE;
    }

  fp = fopen(path, append ? "at" : "wt");
  if (!fp)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not open %s for diffraction export.",
                  path);
      return FALSE;
    }

  for (i = 0; i < pattern->x_values->len; i++)
    {
      const gdouble x = g_array_index(pattern->x_values, gdouble, i);
      const gdouble y = g_array_index(pattern->y_values, gdouble, i);

      fprintf(fp, "%.6f %.10f\n", x, y);
    }

  fclose(fp);
  return TRUE;
}

static gboolean
gdis_diffraction_validate_settings(const GdisDiffractionSettings *settings,
                                   GError **error)
{
  g_return_val_if_fail(settings != NULL, FALSE);

  if (settings->wavelength <= 0.0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Wavelength must be greater than zero.");
      return FALSE;
    }

  if (settings->theta_step <= 0.0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "2Theta step must be greater than zero.");
      return FALSE;
    }

  if (settings->theta_min < 0.0 || settings->theta_max <= settings->theta_min)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "2Theta minimum must be non-negative and strictly smaller than the maximum.");
      return FALSE;
    }

  if (settings->asym < 0.0 || settings->asym > 1.0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "The mixing parameter must lie between 0.0 and 1.0.");
      return FALSE;
    }

  return TRUE;
}

static gdouble
gdis_diffraction_lorentz_polarization_factor(GdisDiffractionRadiation radiation,
                                             gdouble two_theta)
{
  const gdouble theta = 0.5 * two_theta * (G_PI / 180.0);
  const gdouble sin_theta = sin(theta);
  const gdouble sin_two_theta = sin(2.0 * theta);
  gdouble factor;

  if (fabs(sin_theta) < 1.0e-8 || fabs(sin_two_theta) < 1.0e-8)
    return 0.0;

  if (radiation == GDIS_DIFFRACTION_RADIATION_XRAY)
    {
      const gdouble cos_two_theta = cos(2.0 * theta);

      factor = 0.5 * (1.0 + cos_two_theta * cos_two_theta);
    }
  else
    {
      factor = 1.0;
    }

  return factor / (sin_theta * sin_two_theta);
}

static gdouble
gdis_diffraction_peak_profile(const GdisDiffractionSettings *settings,
                              gdouble fwhm,
                              gdouble delta_two_theta)
{
  const gdouble c1 = 4.0 * log(2.0);
  const gdouble sqrt_c1 = sqrt(c1);
  const gdouble sqrt_pi = sqrt(G_PI);

  g_return_val_if_fail(settings != NULL, 0.0);

  if (fabs(fwhm) <= 1.0e-5)
    {
      return fabs(delta_two_theta) <= (0.5 * settings->theta_step)
               ? (1.0 / settings->theta_step)
               : 0.0;
    }

  switch (settings->broadening)
    {
    case GDIS_DIFFRACTION_BROADENING_LORENTZIAN:
      return fwhm / (2.0 * G_PI * (delta_two_theta * delta_two_theta + 0.25 * fwhm * fwhm));

    case GDIS_DIFFRACTION_BROADENING_PSEUDO_VOIGT:
      {
        const gdouble lorentz =
          2.0 / (fwhm * G_PI * (1.0 + 4.0 * (delta_two_theta * delta_two_theta / (fwhm * fwhm))));
        const gdouble gaussian =
          sqrt_c1 * exp(-c1 * delta_two_theta * delta_two_theta / (fwhm * fwhm)) /
          (fwhm * G_PI);

        return settings->asym * lorentz + (1.0 - settings->asym) * gaussian;
      }

    case GDIS_DIFFRACTION_BROADENING_GAUSSIAN:
    default:
      return sqrt_c1 * exp(-c1 * delta_two_theta * delta_two_theta / (fwhm * fwhm)) /
             (fwhm * sqrt_pi);
    }
}

static gdouble
gdis_diffraction_structure_factor_intensity(const GdisModel *model,
                                            const gdouble inverse[9],
                                            GdisDiffractionRadiation radiation,
                                            gint h,
                                            gint k,
                                            gint l,
                                            gdouble sol)
{
  gdouble real_part;
  gdouble imag_part;
  guint i;

  (void) sol;

  real_part = 0.0;
  imag_part = 0.0;

  for (i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      gdouble frac[3];
      gdouble phase;
      gdouble weight;

      atom = g_ptr_array_index(model->atoms, i);
      gdis_cart_to_frac(inverse, atom->position, frac);
      frac[0] -= floor(frac[0]);
      frac[1] -= floor(frac[1]);
      frac[2] -= floor(frac[2]);

      phase = 2.0 * G_PI * ((gdouble) h * frac[0] +
                            (gdouble) k * frac[1] +
                            (gdouble) l * frac[2]);
      weight = (gdouble) gdis_atomic_number_for_element(atom->element);

      if (radiation == GDIS_DIFFRACTION_RADIATION_NEUTRON)
        weight = sqrt(weight);
      else if (radiation == GDIS_DIFFRACTION_RADIATION_ELECTRON)
        weight *= 1.0 + (0.15 * weight);

      real_part += weight * cos(phase);
      imag_part += weight * sin(phase);
    }

  return real_part * real_part + imag_part * imag_part;
}

static gchar *
gdis_diffraction_pattern_summary(const GdisModel *model,
                                 const GdisDiffractionPattern *pattern)
{
  GString *summary;
  guint i;

  g_return_val_if_fail(pattern != NULL, NULL);

  summary = g_string_new("");
  g_string_append_printf(summary,
                         "Diffraction Preview\nModel: %s\nRadiation: %s\nWavelength: %.6f\n"
                         "2theta range: %.2f to %.2f step %.2f\n"
                         "Broadening: %s\n\n",
                         model ? model->basename : "none",
                         gdis_diffraction_radiation_label(pattern->settings.radiation),
                         pattern->settings.wavelength,
                         pattern->settings.theta_min,
                         pattern->settings.theta_max,
                         pattern->settings.theta_step,
                         gdis_diffraction_broadening_label(pattern->settings.broadening));
  g_string_append(summary, "  2theta    rel I     d(A)     h   k   l\n");

  for (i = 0; i < pattern->peaks->len && i < 16; i++)
    {
      const GdisDiffractionPeak *peak;

      peak = &g_array_index(pattern->peaks, GdisDiffractionPeak, i);
      g_string_append_printf(summary,
                             " %7.3f  %6.2f  %7.4f   %3d %3d %3d\n",
                             peak->two_theta,
                             peak->relative_intensity,
                             peak->d_spacing,
                             peak->h,
                             peak->k,
                             peak->l);
    }

  g_string_append(summary,
                  "\nThis summary is backed by the GTK4 diffraction engine and matches the plot window inputs.");
  return g_string_free(summary, FALSE);
}

char *
gdis_report_measurements(const GdisModel *model,
                         const guint *atom_indices,
                         guint count,
                         GdisMeasureMode mode)
{
  GString *report;
  guint valid_count;
  guint i;
  const char *mode_label;

  report = g_string_new("");
  if (!model)
    {
      g_string_append(report, "No active model loaded.");
      return g_string_free(report, FALSE);
    }

  g_string_append_printf(report,
                         "Model: %s\nFormat: %s\n",
                         model->basename,
                         model->format_label);

  switch (mode)
    {
    case GDIS_MEASURE_MODE_DISTANCE:
      mode_label = "Distance";
      break;
    case GDIS_MEASURE_MODE_ANGLE:
      mode_label = "Angle";
      break;
    case GDIS_MEASURE_MODE_TORSION:
      mode_label = "Torsion";
      break;
    case GDIS_MEASURE_MODE_AUTO:
    default:
      mode_label = "Auto";
      break;
    }

  g_string_append_printf(report, "Mode: %s\n\n", mode_label);

  if (!atom_indices || count == 0)
    {
      g_string_append(report,
                      "No atoms have been picked yet.\n"
                      "Click atoms in the viewer to build a measurement set.\n"
                      "The last 2 picked atoms give a distance,\n"
                      "the last 3 give an angle,\n"
                      "and the last 4 give a torsion.\n");
      return g_string_free(report, FALSE);
    }

  if (mode == GDIS_MEASURE_MODE_DISTANCE && count < 2)
    {
      g_string_append(report, "Pick 2 atoms to calculate a distance.");
      return g_string_free(report, FALSE);
    }
  if (mode == GDIS_MEASURE_MODE_ANGLE && count < 3)
    {
      g_string_append(report, "Pick 3 atoms to calculate an angle.");
      return g_string_free(report, FALSE);
    }
  if (mode == GDIS_MEASURE_MODE_TORSION && count < 4)
    {
      g_string_append(report, "Pick 4 atoms to calculate a torsion.");
      return g_string_free(report, FALSE);
    }

  valid_count = 0;
  for (i = 0; i < count; i++)
    {
      const GdisAtom *atom;

      if (atom_indices[i] >= model->atoms->len)
        continue;

      atom = g_ptr_array_index(model->atoms, atom_indices[i]);
      valid_count++;
      g_string_append_printf(report,
                             "Pick %u: %s [%s] #%u  %.5f  %.5f  %.5f\n",
                             valid_count,
                             atom->label,
                             atom->element,
                             atom->serial,
                             atom->position[0],
                             atom->position[1],
                             atom->position[2]);
    }

  if (valid_count < 2 && mode == GDIS_MEASURE_MODE_AUTO)
    {
      g_string_append(report, "\nPick one more atom to measure a distance.");
      return g_string_free(report, FALSE);
    }

  if (mode == GDIS_MEASURE_MODE_DISTANCE ||
      mode == GDIS_MEASURE_MODE_AUTO)
    {
      if (valid_count >= 2)
        {
          const GdisAtom *a1;
          const GdisAtom *a2;
          gdouble distance;

          a1 = g_ptr_array_index(model->atoms, atom_indices[valid_count - 2]);
          a2 = g_ptr_array_index(model->atoms, atom_indices[valid_count - 1]);
          distance = gdis_measure_distance(model, a1->position, a2->position);
          g_string_append_printf(report,
                                 "\nDistance (%s -> %s): %.5f A\n",
                                 a1->label,
                                 a2->label,
                                 distance);
        }
      if (mode == GDIS_MEASURE_MODE_DISTANCE)
        return g_string_free(report, FALSE);
    }

  if (mode == GDIS_MEASURE_MODE_ANGLE ||
      mode == GDIS_MEASURE_MODE_AUTO)
    {
      if (valid_count >= 3)
        {
          const GdisAtom *a1;
          const GdisAtom *a2;
          const GdisAtom *a3;
          gdouble angle;

          a1 = g_ptr_array_index(model->atoms, atom_indices[valid_count - 3]);
          a2 = g_ptr_array_index(model->atoms, atom_indices[valid_count - 2]);
          a3 = g_ptr_array_index(model->atoms, atom_indices[valid_count - 1]);
          angle = gdis_measure_angle(model, a1->position, a2->position, a3->position);
          g_string_append_printf(report,
                                 "Angle (%s-%s-%s): %.3f deg\n",
                                 a1->label,
                                 a2->label,
                                 a3->label,
                                 angle);
        }
      if (mode == GDIS_MEASURE_MODE_ANGLE)
        return g_string_free(report, FALSE);
    }

  if (mode == GDIS_MEASURE_MODE_TORSION ||
      mode == GDIS_MEASURE_MODE_AUTO)
    {
      if (valid_count >= 4)
        {
          const GdisAtom *a1;
          const GdisAtom *a2;
          const GdisAtom *a3;
          const GdisAtom *a4;
          gdouble torsion;

          a1 = g_ptr_array_index(model->atoms, atom_indices[valid_count - 4]);
          a2 = g_ptr_array_index(model->atoms, atom_indices[valid_count - 3]);
          a3 = g_ptr_array_index(model->atoms, atom_indices[valid_count - 2]);
          a4 = g_ptr_array_index(model->atoms, atom_indices[valid_count - 1]);
          torsion = gdis_measure_torsion(model,
                                         a1->position,
                                         a2->position,
                                         a3->position,
                                         a4->position);
          g_string_append_printf(report,
                                 "Torsion (%s-%s-%s-%s): %.3f deg\n",
                                 a1->label,
                                 a2->label,
                                 a3->label,
                                 a4->label,
                                 torsion);
        }
      if (mode == GDIS_MEASURE_MODE_TORSION)
        return g_string_free(report, FALSE);
    }

  if (valid_count > 4 && mode == GDIS_MEASURE_MODE_AUTO)
    g_string_append(report, "\nOnly the last 4 picked atoms are used for torsion.");

  return g_string_free(report, FALSE);
}

char *
gdis_report_surface(const GdisModel *model)
{
  GString *report;
  GArray *planes;
  gdouble matrix[9];
  gdouble inverse[9];
  gint h;
  gint k;
  gint l;
  guint i;

  report = g_string_new("");
  if (!model)
    {
      g_string_append(report, "No active model loaded.");
      return g_string_free(report, FALSE);
    }

  g_string_append_printf(report,
                         "Surface Explorer\nModel: %s\n\n",
                         model->basename);

  if (!model->periodic || !gdis_build_cell_matrix(model, matrix, inverse))
    {
      g_string_append(report,
                      "Surface analysis needs a periodic crystal model with a valid unit cell.\n"
                      "Open a CIF, ARC/CAR, or another periodic structure first.");
      return g_string_free(report, FALSE);
    }

  planes = g_array_new(FALSE, FALSE, sizeof(GdisPlaneEntry));

  for (h = -4; h <= 4; h++)
    {
      for (k = -4; k <= 4; k++)
        {
          for (l = -4; l <= 4; l++)
            {
              GdisPlaneEntry plane;

              if (h == 0 && k == 0 && l == 0)
                continue;

              plane.h = h;
              plane.k = k;
              plane.l = l;
              gdis_canonicalize_hkl(&plane.h, &plane.k, &plane.l);
              if (plane.h == 0 && plane.k == 0 && plane.l == 0)
                continue;

              plane.d_spacing = gdis_d_spacing(inverse, plane.h, plane.k, plane.l);
              if (plane.d_spacing <= 1.0e-6)
                continue;

              {
                gboolean duplicate;
                guint j;

                duplicate = FALSE;
                for (j = 0; j < planes->len; j++)
                  {
                    GdisPlaneEntry *existing;

                    existing = &g_array_index(planes, GdisPlaneEntry, j);
                    if (existing->h == plane.h &&
                        existing->k == plane.k &&
                        existing->l == plane.l)
                      {
                        duplicate = TRUE;
                        break;
                      }
                  }
                if (!duplicate)
                  g_array_append_val(planes, plane);
              }
            }
        }
    }

  g_array_sort(planes, gdis_plane_compare);
  g_string_append_printf(report,
                         "Cell: %.4f %.4f %.4f / %.2f %.2f %.2f\n"
                         "Space group: %s\n\n"
                         "Low-index planes ranked by d-spacing:\n",
                         model->cell_lengths[0],
                         model->cell_lengths[1],
                         model->cell_lengths[2],
                         model->cell_angles[0],
                         model->cell_angles[1],
                         model->cell_angles[2],
                         model->space_group ? model->space_group : "unknown");
  g_string_append(report, "  h   k   l     d(A)\n");

  for (i = 0; i < planes->len && i < 24; i++)
    {
      const GdisPlaneEntry *plane;

      plane = &g_array_index(planes, GdisPlaneEntry, i);
      g_string_append_printf(report,
                             "%3d %3d %3d   %8.4f\n",
                             plane->h,
                             plane->k,
                             plane->l,
                             plane->d_spacing);
    }

  g_string_append(report,
                  "\nThis GTK4 bridge now provides plane ranking and d-spacing lookup.\n"
                  "The full legacy slab-construction workflow is still a deeper port.");
  g_array_free(planes, TRUE);
  return g_string_free(report, FALSE);
}

char *
gdis_report_diffraction(const GdisModel *model)
{
  GdisDiffractionSettings settings;
  GdisDiffractionPattern *pattern;
  GError *error;
  gchar *summary;

  gdis_diffraction_settings_init(&settings);

  pattern = NULL;
  error = NULL;
  if (!gdis_diffraction_calculate(model, &settings, &pattern, &error))
    {
      summary = g_strdup(error ? error->message : "Diffraction calculation failed.");
      g_clear_error(&error);
      return summary;
    }

  summary = gdis_diffraction_pattern_summary(model, pattern);
  gdis_diffraction_pattern_free(pattern);
  return summary;
}

gboolean
gdis_diffraction_calculate(const GdisModel *model,
                           const GdisDiffractionSettings *settings,
                           GdisDiffractionPattern **pattern_out,
                           GError **error)
{
  GdisDiffractionPattern *pattern;
  gdouble matrix[9];
  gdouble inverse[9];
  gdouble theta_limit;
  gdouble dhkl_min;
  gdouble max_peak_intensity;
  gdouble max_spectrum_intensity;
  gdouble merge_tolerance;
  guint sample_count;
  gint h;
  gint k;
  gint l;
  guint i;

  g_return_val_if_fail(pattern_out != NULL, FALSE);
  *pattern_out = NULL;

  if (!model)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "No active model loaded.");
      return FALSE;
    }

  if (!gdis_diffraction_validate_settings(settings, error))
    return FALSE;

  if (!model->periodic || !gdis_build_cell_matrix(model, matrix, inverse))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Diffraction needs a periodic crystal model with a valid unit cell.");
      return FALSE;
    }

  if (model->periodicity > 0 && model->periodicity < 3)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Diffraction currently targets fully 3D periodic crystals.");
      return FALSE;
    }

  theta_limit = settings->theta_max;
  if (settings->wavelength < 0.5)
    {
      const gdouble limit = asin(CLAMP(2.0 * settings->wavelength, -1.0, 1.0)) * (180.0 / G_PI);

      theta_limit = MIN(theta_limit, limit);
      if (settings->theta_min >= theta_limit)
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_FAILED,
                      "The selected wavelength restricts 2Theta to %.3f degrees or less.",
                      limit);
          return FALSE;
        }
    }

  dhkl_min = 0.5 * settings->wavelength /
             sin(0.5 * theta_limit * (G_PI / 180.0));

  pattern = g_new0(GdisDiffractionPattern, 1);
  pattern->settings = *settings;
  pattern->settings.theta_max = theta_limit;
  pattern->x_values = g_array_new(FALSE, FALSE, sizeof(gdouble));
  pattern->y_values = g_array_new(FALSE, FALSE, sizeof(gdouble));
  pattern->peaks = g_array_new(FALSE, FALSE, sizeof(GdisDiffractionPeak));

  merge_tolerance = MAX(0.02, settings->theta_step * 0.2);

  for (h = -19; h <= 19; h++)
    {
      for (k = -19; k <= 19; k++)
        {
          for (l = -19; l <= 19; l++)
            {
              GdisDiffractionPeak peak;
              gdouble sin_theta;
              gdouble sol;
              gdouble intensity;
              gboolean merged;
              guint j;

              if (h == 0 && k == 0 && l == 0)
                continue;

              peak.h = h;
              peak.k = k;
              peak.l = l;
              peak.reflections_merged = 1;
              peak.d_spacing = gdis_d_spacing(inverse, h, k, l);
              if (peak.d_spacing <= dhkl_min || peak.d_spacing <= 1.0e-6)
                continue;

              sin_theta = settings->wavelength / (2.0 * peak.d_spacing);
              if (sin_theta <= 0.0 || sin_theta >= 1.0)
                continue;

              peak.two_theta = 2.0 * asin(sin_theta) * (180.0 / G_PI);
              if (peak.two_theta < settings->theta_min - merge_tolerance ||
                  peak.two_theta > theta_limit + merge_tolerance)
                continue;

              sol = 0.5 / peak.d_spacing;
              intensity = gdis_diffraction_structure_factor_intensity(model,
                                                                     inverse,
                                                                     settings->radiation,
                                                                     h,
                                                                     k,
                                                                     l,
                                                                     sol);
              if (intensity < 1.0e-10)
                continue;

              peak.intensity = intensity;
              peak.relative_intensity = 0.0;

              merged = FALSE;
              for (j = 0; j < pattern->peaks->len; j++)
                {
                  GdisDiffractionPeak *existing;

                  existing = &g_array_index(pattern->peaks, GdisDiffractionPeak, j);
                  if (fabs(existing->two_theta - peak.two_theta) < merge_tolerance)
                    {
                      existing->intensity += peak.intensity;
                      existing->reflections_merged++;
                      merged = TRUE;
                      break;
                    }
                }

              if (!merged)
                g_array_append_val(pattern->peaks, peak);
            }
        }
    }

  if (pattern->peaks->len == 0)
    {
      gdis_diffraction_pattern_free(pattern);
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "No Bragg peaks were generated for the current unit cell and 2Theta range.");
      return FALSE;
    }

  g_array_sort(pattern->peaks, gdis_peak_compare);

  max_peak_intensity = 0.0;
  for (i = 0; i < pattern->peaks->len; i++)
    {
      GdisDiffractionPeak *peak;

      peak = &g_array_index(pattern->peaks, GdisDiffractionPeak, i);
      peak->intensity *= gdis_diffraction_lorentz_polarization_factor(settings->radiation,
                                                                      peak->two_theta);
      max_peak_intensity = MAX(max_peak_intensity, peak->intensity);
    }

  for (i = 0; i < pattern->peaks->len; i++)
    {
      GdisDiffractionPeak *peak;

      peak = &g_array_index(pattern->peaks, GdisDiffractionPeak, i);
      peak->relative_intensity = (max_peak_intensity > 0.0)
                                   ? (100.0 * peak->intensity / max_peak_intensity)
                                   : 0.0;
    }

  sample_count = 1 + (guint) floor((theta_limit - settings->theta_min) / settings->theta_step);
  for (i = 0; i < sample_count; i++)
    {
      const gdouble x = settings->theta_min + ((gdouble) i * settings->theta_step);
      const gdouble y = 0.0;

      g_array_append_val(pattern->x_values, x);
      g_array_append_val(pattern->y_values, y);
    }

  max_spectrum_intensity = 0.0;
  for (i = 0; i < pattern->peaks->len; i++)
    {
      const GdisDiffractionPeak *peak;
      gdouble theta;
      gdouble tan_theta;
      gdouble fwhm_sq;
      gdouble fwhm;
      guint sample_index;

      peak = &g_array_index(pattern->peaks, GdisDiffractionPeak, i);
      theta = 0.5 * peak->two_theta * (G_PI / 180.0);
      tan_theta = tan(theta);
      fwhm_sq = settings->w + settings->v * tan_theta + settings->u * tan_theta * tan_theta;
      fwhm = (fwhm_sq > 1.0e-10) ? sqrt(fwhm_sq) : 0.0;

      for (sample_index = 0; sample_index < pattern->x_values->len; sample_index++)
        {
          const gdouble x = g_array_index(pattern->x_values, gdouble, sample_index);
          const gdouble contribution = peak->intensity *
                                       gdis_diffraction_peak_profile(&pattern->settings,
                                                                     fwhm,
                                                                     x - peak->two_theta);
          gdouble *sample;

          sample = &g_array_index(pattern->y_values, gdouble, sample_index);
          *sample += contribution;
          max_spectrum_intensity = MAX(max_spectrum_intensity, *sample);
        }
    }

  pattern->max_intensity = MAX(max_peak_intensity, max_spectrum_intensity);
  *pattern_out = pattern;
  return TRUE;
}

char *
gdis_report_isosurface(const GdisModel *model)
{
  GString *report;

  report = g_string_new("");
  g_string_append(report, "Iso-surface Status\n\n");
  if (model)
    g_string_append_printf(report, "Active model: %s\n\n", model->basename);

  g_string_append(report,
                  "Iso-surfaces are not just a missing button in this GTK4 rebuild.\n"
                  "They need volumetric scalar-field data, such as electron density or potential grids.\n\n"
                  "The current bridge loads structure models only:\n"
                  "  XYZ\n"
                  "  PDB\n"
                  "  ARC/CAR\n"
                  "  CIF\n\n"
                  "To make iso-surfaces real here, the next port would need:\n"
                  "  1. grid-data loaders\n"
                  "  2. scalar-field storage in the GTK4 model bridge\n"
                  "  3. marching-cubes style surface extraction\n"
                  "  4. shaded mesh rendering in the GTK4 viewer\n\n"
                  "So this one is possible, but it is a separate engine port, not a quick UI hookup.");
  return g_string_free(report, FALSE);
}

static gboolean
gdis_build_cell_vectors(const GdisModel *model,
                        gdouble a_vec[3],
                        gdouble b_vec[3],
                        gdouble c_vec[3])
{
  gdouble alpha;
  gdouble beta;
  gdouble gamma;
  gdouble sin_gamma;
  gdouble cx;
  gdouble cy;
  gdouble cz2;

  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->periodic)
    return FALSE;

  if (model->cell_lengths[0] <= 1.0e-6 ||
      model->cell_lengths[1] <= 1.0e-6 ||
      model->cell_lengths[2] <= 1.0e-6)
    return FALSE;

  alpha = model->cell_angles[0] * (G_PI / 180.0);
  beta = model->cell_angles[1] * (G_PI / 180.0);
  gamma = model->cell_angles[2] * (G_PI / 180.0);
  sin_gamma = sin(gamma);
  if (fabs(sin_gamma) < 1.0e-6)
    return FALSE;

  a_vec[0] = model->cell_lengths[0];
  a_vec[1] = 0.0;
  a_vec[2] = 0.0;

  b_vec[0] = model->cell_lengths[1] * cos(gamma);
  b_vec[1] = model->cell_lengths[1] * sin_gamma;
  b_vec[2] = 0.0;

  cx = model->cell_lengths[2] * cos(beta);
  cy = model->cell_lengths[2] * (cos(alpha) - cos(beta) * cos(gamma)) / sin_gamma;
  cz2 = model->cell_lengths[2] * model->cell_lengths[2] - cx * cx - cy * cy;
  if (cz2 < 0.0)
    cz2 = 0.0;

  c_vec[0] = cx;
  c_vec[1] = cy;
  c_vec[2] = sqrt(cz2);
  return TRUE;
}

static gboolean
gdis_build_cell_matrix(const GdisModel *model,
                       gdouble matrix[9],
                       gdouble inverse[9])
{
  gdouble a_vec[3];
  gdouble b_vec[3];
  gdouble c_vec[3];

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(matrix != NULL, FALSE);
  g_return_val_if_fail(inverse != NULL, FALSE);

  if (!gdis_build_cell_vectors(model, a_vec, b_vec, c_vec))
    return FALSE;

  matrix[0] = a_vec[0];
  matrix[1] = b_vec[0];
  matrix[2] = c_vec[0];
  matrix[3] = a_vec[1];
  matrix[4] = b_vec[1];
  matrix[5] = c_vec[1];
  matrix[6] = a_vec[2];
  matrix[7] = b_vec[2];
  matrix[8] = c_vec[2];

  return gdis_invert_matrix3(matrix, inverse);
}

static gboolean
gdis_invert_matrix3(const gdouble matrix[9], gdouble inverse[9])
{
  gdouble det;

  det = matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7]) -
        matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6]) +
        matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
  if (fabs(det) < 1.0e-12)
    return FALSE;

  inverse[0] = (matrix[4] * matrix[8] - matrix[5] * matrix[7]) / det;
  inverse[1] = (matrix[2] * matrix[7] - matrix[1] * matrix[8]) / det;
  inverse[2] = (matrix[1] * matrix[5] - matrix[2] * matrix[4]) / det;
  inverse[3] = (matrix[5] * matrix[6] - matrix[3] * matrix[8]) / det;
  inverse[4] = (matrix[0] * matrix[8] - matrix[2] * matrix[6]) / det;
  inverse[5] = (matrix[2] * matrix[3] - matrix[0] * matrix[5]) / det;
  inverse[6] = (matrix[3] * matrix[7] - matrix[4] * matrix[6]) / det;
  inverse[7] = (matrix[1] * matrix[6] - matrix[0] * matrix[7]) / det;
  inverse[8] = (matrix[0] * matrix[4] - matrix[1] * matrix[3]) / det;
  return TRUE;
}

static guint
gdis_atomic_number_for_element(const char *element)
{
  gchar normalized[4];
  guint i;

  if (!element || !*element)
    return 6;

  memset(normalized, 0, sizeof(normalized));
  normalized[0] = g_ascii_toupper(element[0]);
  if (element[1] && g_ascii_isalpha(element[1]))
    normalized[1] = g_ascii_tolower(element[1]);

  for (i = 0; i < G_N_ELEMENTS(gdis_atomic_numbers); i++)
    {
      if (g_strcmp0(normalized, gdis_atomic_numbers[i].symbol) == 0)
        return gdis_atomic_numbers[i].atomic_number;
    }

  return 6;
}

static void
gdis_measure_minimum_image_delta(const gdouble matrix[9],
                                 const gdouble inverse[9],
                                 const gdouble from[3],
                                 const gdouble to[3],
                                 gdouble delta[3])
{
  gdouble from_frac[3];
  gdouble to_frac[3];
  gdouble diff_frac[3];
  guint axis;

  gdis_cart_to_frac(inverse, from, from_frac);
  gdis_cart_to_frac(inverse, to, to_frac);

  for (axis = 0; axis < 3; axis++)
    diff_frac[axis] = to_frac[axis] - from_frac[axis];

  for (axis = 0; axis < 3; axis++)
    diff_frac[axis] -= floor(diff_frac[axis] + 0.5);

  delta[0] = matrix[0] * diff_frac[0] + matrix[1] * diff_frac[1] + matrix[2] * diff_frac[2];
  delta[1] = matrix[3] * diff_frac[0] + matrix[4] * diff_frac[1] + matrix[5] * diff_frac[2];
  delta[2] = matrix[6] * diff_frac[0] + matrix[7] * diff_frac[1] + matrix[8] * diff_frac[2];
}

static gboolean
gdis_measure_prepare_periodic(const GdisModel *model,
                              gdouble matrix[9],
                              gdouble inverse[9])
{
  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->periodic)
    return FALSE;

  return gdis_build_cell_matrix(model, matrix, inverse);
}

static void
gdis_cart_to_frac(const gdouble inverse[9], const gdouble cart[3], gdouble frac[3])
{
  frac[0] = inverse[0] * cart[0] + inverse[1] * cart[1] + inverse[2] * cart[2];
  frac[1] = inverse[3] * cart[0] + inverse[4] * cart[1] + inverse[5] * cart[2];
  frac[2] = inverse[6] * cart[0] + inverse[7] * cart[1] + inverse[8] * cart[2];
}

static gdouble
gdis_d_spacing(const gdouble inverse[9], gint h, gint k, gint l)
{
  gdouble gx;
  gdouble gy;
  gdouble gz;
  gdouble magnitude;

  gx = inverse[0] * h + inverse[3] * k + inverse[6] * l;
  gy = inverse[1] * h + inverse[4] * k + inverse[7] * l;
  gz = inverse[2] * h + inverse[5] * k + inverse[8] * l;
  magnitude = sqrt(gx * gx + gy * gy + gz * gz);
  if (magnitude < 1.0e-12)
    return 0.0;
  return 1.0 / magnitude;
}

static gdouble
gdis_structure_factor_intensity(const GdisModel *model,
                                const gdouble inverse[9],
                                gint h,
                                gint k,
                                gint l)
{
  gdouble real_part;
  gdouble imag_part;
  guint i;

  real_part = 0.0;
  imag_part = 0.0;

  for (i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      gdouble frac[3];
      gdouble phase;
      gdouble weight;

      atom = g_ptr_array_index(model->atoms, i);
      gdis_cart_to_frac(inverse, atom->position, frac);
      frac[0] -= floor(frac[0]);
      frac[1] -= floor(frac[1]);
      frac[2] -= floor(frac[2]);

      phase = 2.0 * G_PI * ((gdouble) h * frac[0] +
                            (gdouble) k * frac[1] +
                            (gdouble) l * frac[2]);
      weight = (gdouble) gdis_atomic_number_for_element(atom->element);
      real_part += weight * cos(phase);
      imag_part += weight * sin(phase);
    }

  return real_part * real_part + imag_part * imag_part;
}

static gdouble
gdis_measure_distance(const GdisModel *model, const gdouble a[3], const gdouble b[3])
{
  gdouble dx;
  gdouble dy;
  gdouble dz;
  gdouble matrix[9];
  gdouble inverse[9];
  gdouble delta[3];

  if (model && gdis_measure_prepare_periodic(model, matrix, inverse))
    {
      gdis_measure_minimum_image_delta(matrix, inverse, a, b, delta);
      dx = delta[0];
      dy = delta[1];
      dz = delta[2];
    }
  else
    {
      dx = b[0] - a[0];
      dy = b[1] - a[1];
      dz = b[2] - a[2];
    }

  return sqrt(dx * dx + dy * dy + dz * dz);
}

static gdouble
gdis_measure_angle(const GdisModel *model,
                   const gdouble a[3],
                   const gdouble b[3],
                   const gdouble c[3])
{
  gdouble ba[3];
  gdouble bc[3];
  gdouble dot;
  gdouble mag_ba;
  gdouble mag_bc;
  gdouble cosine;
  gdouble matrix[9];
  gdouble inverse[9];

  if (model && gdis_measure_prepare_periodic(model, matrix, inverse))
    {
      gdis_measure_minimum_image_delta(matrix, inverse, b, a, ba);
      gdis_measure_minimum_image_delta(matrix, inverse, b, c, bc);
    }
  else
    {
      ba[0] = a[0] - b[0];
      ba[1] = a[1] - b[1];
      ba[2] = a[2] - b[2];
      bc[0] = c[0] - b[0];
      bc[1] = c[1] - b[1];
      bc[2] = c[2] - b[2];
    }

  dot = ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2];
  mag_ba = sqrt(ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]);
  mag_bc = sqrt(bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]);
  if (mag_ba < 1.0e-12 || mag_bc < 1.0e-12)
    return 0.0;

  cosine = dot / (mag_ba * mag_bc);
  cosine = CLAMP(cosine, -1.0, 1.0);
  return acos(cosine) * (180.0 / G_PI);
}

static gdouble
gdis_measure_torsion(const GdisModel *model,
                     const gdouble a[3],
                     const gdouble b[3],
                     const gdouble c[3],
                     const gdouble d[3])
{
  gdouble pos_a[3];
  gdouble pos_b[3];
  gdouble pos_c[3];
  gdouble pos_d[3];
  gdouble b1[3];
  gdouble b2[3];
  gdouble b3[3];
  gdouble n1[3];
  gdouble n2[3];
  gdouble m1[3];
  gdouble n1_mag;
  gdouble n2_mag;
  gdouble b2_mag;
  gdouble x;
  gdouble y;
  gdouble matrix[9];
  gdouble inverse[9];
  gdouble delta[3];

  if (model && gdis_measure_prepare_periodic(model, matrix, inverse))
    {
      memcpy(pos_a, a, sizeof(pos_a));
      gdis_measure_minimum_image_delta(matrix, inverse, a, b, delta);
      pos_b[0] = pos_a[0] + delta[0];
      pos_b[1] = pos_a[1] + delta[1];
      pos_b[2] = pos_a[2] + delta[2];

      gdis_measure_minimum_image_delta(matrix, inverse, b, c, delta);
      pos_c[0] = pos_b[0] + delta[0];
      pos_c[1] = pos_b[1] + delta[1];
      pos_c[2] = pos_b[2] + delta[2];

      gdis_measure_minimum_image_delta(matrix, inverse, c, d, delta);
      pos_d[0] = pos_c[0] + delta[0];
      pos_d[1] = pos_c[1] + delta[1];
      pos_d[2] = pos_c[2] + delta[2];
    }
  else
    {
      memcpy(pos_a, a, sizeof(pos_a));
      memcpy(pos_b, b, sizeof(pos_b));
      memcpy(pos_c, c, sizeof(pos_c));
      memcpy(pos_d, d, sizeof(pos_d));
    }

  b1[0] = pos_b[0] - pos_a[0];
  b1[1] = pos_b[1] - pos_a[1];
  b1[2] = pos_b[2] - pos_a[2];
  b2[0] = pos_c[0] - pos_b[0];
  b2[1] = pos_c[1] - pos_b[1];
  b2[2] = pos_c[2] - pos_b[2];
  b3[0] = pos_d[0] - pos_c[0];
  b3[1] = pos_d[1] - pos_c[1];
  b3[2] = pos_d[2] - pos_c[2];

  n1[0] = b1[1] * b2[2] - b1[2] * b2[1];
  n1[1] = b1[2] * b2[0] - b1[0] * b2[2];
  n1[2] = b1[0] * b2[1] - b1[1] * b2[0];

  n2[0] = b2[1] * b3[2] - b2[2] * b3[1];
  n2[1] = b2[2] * b3[0] - b2[0] * b3[2];
  n2[2] = b2[0] * b3[1] - b2[1] * b3[0];

  n1_mag = sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
  n2_mag = sqrt(n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]);
  b2_mag = sqrt(b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]);
  if (n1_mag < 1.0e-12 || n2_mag < 1.0e-12 || b2_mag < 1.0e-12)
    return 0.0;

  m1[0] = n1[1] * b2[2] - n1[2] * b2[1];
  m1[1] = n1[2] * b2[0] - n1[0] * b2[2];
  m1[2] = n1[0] * b2[1] - n1[1] * b2[0];

  x = (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) / (n1_mag * n2_mag);
  y = (m1[0] * n2[0] + m1[1] * n2[1] + m1[2] * n2[2]) / (n1_mag * n2_mag * b2_mag);

  return atan2(y, x) * (180.0 / G_PI);
}

static gint
gdis_plane_compare(gconstpointer left, gconstpointer right)
{
  const GdisPlaneEntry *a;
  const GdisPlaneEntry *b;

  a = left;
  b = right;

  if (a->d_spacing > b->d_spacing)
    return -1;
  if (a->d_spacing < b->d_spacing)
    return 1;

  if (a->h != b->h)
    return a->h - b->h;
  if (a->k != b->k)
    return a->k - b->k;
  return a->l - b->l;
}

static gint
gdis_peak_compare(gconstpointer left, gconstpointer right)
{
  const GdisDiffractionPeak *a;
  const GdisDiffractionPeak *b;

  a = left;
  b = right;

  if (a->two_theta < b->two_theta)
    return -1;
  if (a->two_theta > b->two_theta)
    return 1;
  return 0;
}

static void
gdis_canonicalize_hkl(gint *h, gint *k, gint *l)
{
  if (*h < 0 ||
      (*h == 0 && *k < 0) ||
      (*h == 0 && *k == 0 && *l < 0))
    {
      *h = -*h;
      *k = -*k;
      *l = -*l;
    }
}
