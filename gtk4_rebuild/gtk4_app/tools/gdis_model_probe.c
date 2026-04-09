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

  charge_text = model->has_total_charge
    ? g_strdup_printf("%.4f", model->total_charge_e)
    : g_strdup("n/a");
  density_text = model->has_density
    ? g_strdup_printf("%.4f", model->density_g_cm3)
    : g_strdup("n/a");
  energy_text = model->has_energy
    ? g_strdup_printf("%.5f", model->energy_ev)
    : g_strdup("n/a");

  g_print("OK\t%s\tformat=%s\tatoms=%u\tmolecules=%u\tbonds=%u\texplicit=%u\tcharge=%s\tdensity=%s\tenergy=%s\tframes=%u\tperiodic=%s\tdims=%u\tspace_group=%s\ttitle=%s\n",
          path,
          gdis_model_format_label(model->format),
          model->atom_count,
          gdis_model_get_component_count(model),
          model->bond_count,
          model->explicit_bond_count,
          charge_text,
          density_text,
          energy_text,
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
