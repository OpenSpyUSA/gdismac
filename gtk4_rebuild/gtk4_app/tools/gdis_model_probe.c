#include "../src/gdis_model.h"

#include <glib.h>
#include <stdio.h>

static void
print_model_summary(const char *path,
                    const GdisModel *model)
{
  g_autofree gchar *charge_text = NULL;
  g_autofree gchar *density_text = NULL;
  g_autofree gchar *energy_text = NULL;
  g_autofree gchar *force_text = NULL;
  g_autofree gchar *pressure_text = NULL;
  const guint pressure_points = (model->pressure_x_values && model->pressure_y_values_gpa)
                                  ? MIN(model->pressure_x_values->len, model->pressure_y_values_gpa->len)
                                  : 0u;
  const guint dos_points = (model->dos_x_values_ev && model->dos_y_values)
                             ? MIN(model->dos_x_values_ev->len, model->dos_y_values->len)
                             : 0u;
  const guint band_points = (model->band_x_values && model->band_y_values_ev)
                              ? MIN(model->band_x_values->len, model->band_y_values_ev->len)
                              : 0u;
  const guint frequency_points = (model->frequency_x_values_cm1 && model->frequency_y_values)
                                   ? MIN(model->frequency_x_values_cm1->len, model->frequency_y_values->len)
                                   : 0u;
  const guint raman_points = (model->raman_x_values_cm1 && model->raman_y_values)
                               ? MIN(model->raman_x_values_cm1->len, model->raman_y_values->len)
                               : 0u;

  charge_text = model->has_total_charge
    ? g_strdup_printf("%.4f", model->total_charge_e)
    : g_strdup("n/a");
  density_text = model->has_density
    ? g_strdup_printf("%.4f", model->density_g_cm3)
    : g_strdup("n/a");
  energy_text = model->has_energy
    ? g_strdup_printf("%.5f", model->energy_ev)
    : g_strdup("n/a");
  force_text = model->has_force_rms
    ? g_strdup_printf("%.5f", model->force_rms_ev_ang)
    : g_strdup("n/a");
  pressure_text = model->has_pressure
    ? g_strdup_printf("%.5f", model->pressure_gpa)
    : g_strdup("n/a");

  g_print("OK\t%s\tformat=%s\tatoms=%u\tmolecules=%u\tbonds=%u\texplicit=%u\tcharge=%s\tdensity=%s\tenergy=%s\tforce_rms=%s\tpressure=%s\tpressure_points=%u\tdos_points=%u\tband_points=%u\tband_paths=%u\tband_series=%u\tfrequency_points=%u\traman_points=%u\tframes=%u\tperiodic=%s\tdims=%u\tspace_group=%s\ttitle=%s\n",
          path,
          gdis_model_format_label(model->format),
          model->atom_count,
          gdis_model_get_component_count(model),
          model->bond_count,
          model->explicit_bond_count,
          charge_text,
          density_text,
          energy_text,
          force_text,
          pressure_text,
          pressure_points,
          dos_points,
          band_points,
          model->band_path_count,
          model->band_series_count,
          frequency_points,
          raman_points,
          gdis_model_get_frame_count(model),
          model->periodic ? "yes" : "no",
          model->periodicity,
          model->space_group ? model->space_group : "",
          model->title ? model->title : "");
}

int
main(int argc, char **argv)
{
  int exit_code = 0;

  if (argc < 2)
    {
      g_printerr("usage: %s <model-path> [more-paths...]\n", argv[0]);
      return 2;
    }

  for (int i = 1; i < argc; i++)
    {
      const char *path = argv[i];
      g_autoptr(GError) error = NULL;
      g_autoptr(GdisModel) model = gdis_model_load(path, &error);

      if (!model)
        {
          exit_code = 1;
          g_printerr("ERR\t%s\t%s\n",
                     path,
                     error && error->message ? error->message : "load failed");
          continue;
        }

      print_model_summary(path, model);
    }

  return exit_code;
}
