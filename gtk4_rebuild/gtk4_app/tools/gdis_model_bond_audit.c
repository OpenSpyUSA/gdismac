#include "../src/gdis_model.h"

#include <glib.h>
#include <stdio.h>

typedef struct
{
  guint total;
  GHashTable *pair_counts;
  GHashTable *coordination_counts;
} BondAudit;

static gchar *
make_pair_key(const char *left,
              const char *right)
{
  const char *a = left && left[0] ? left : "X";
  const char *b = right && right[0] ? right : "X";

  if (g_ascii_strcasecmp(a, b) <= 0)
    return g_strdup_printf("%s-%s", a, b);
  return g_strdup_printf("%s-%s", b, a);
}

static void
increment_counter(GHashTable *table,
                  const char *key)
{
  gpointer value;
  guint count;

  value = g_hash_table_lookup(table, key);
  count = value ? GPOINTER_TO_UINT(value) : 0u;
  g_hash_table_replace(table, g_strdup(key), GUINT_TO_POINTER(count + 1u));
}

static void
print_count_entry(gpointer key,
                  gpointer value,
                  gpointer user_data)
{
  (void) user_data;
  g_print("  %s: %u\n", (const char *) key, GPOINTER_TO_UINT(value));
}

static void
print_sorted_counts(const char *header,
                    GHashTable *table)
{
  g_print("%s\n", header);
  g_hash_table_foreach(table, print_count_entry, NULL);
}

int
main(int argc, char **argv)
{
  if (argc != 2)
    {
      g_printerr("usage: %s <model-path>\n", argv[0]);
      return 2;
    }

  g_autoptr(GError) error = NULL;
  g_autoptr(GdisModel) model = gdis_model_load(argv[1], &error);
  BondAudit audit = {0};
  GArray *coordination = NULL;

  if (!model)
    {
      g_printerr("%s\n", error && error->message ? error->message : "load failed");
      return 1;
    }

  audit.pair_counts = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
  audit.coordination_counts = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
  coordination = g_array_sized_new(FALSE, TRUE, sizeof(guint), model->atom_count);
  g_array_set_size(coordination, model->atom_count);

  for (guint i = 0; i < coordination->len; i++)
    g_array_index(coordination, guint, i) = 0u;

  for (guint i = 0; i < model->bonds->len; i++)
    {
      const GdisBond *bond;
      const GdisAtom *atom_a;
      const GdisAtom *atom_b;
      g_autofree gchar *pair_key = NULL;

      bond = &g_array_index(model->bonds, GdisBond, i);
      if (bond->atom_index_a >= model->atom_count || bond->atom_index_b >= model->atom_count)
        continue;

      atom_a = g_ptr_array_index(model->atoms, bond->atom_index_a);
      atom_b = g_ptr_array_index(model->atoms, bond->atom_index_b);
      pair_key = make_pair_key(atom_a ? atom_a->element : NULL,
                               atom_b ? atom_b->element : NULL);
      increment_counter(audit.pair_counts, pair_key);
      g_array_index(coordination, guint, bond->atom_index_a)++;
      g_array_index(coordination, guint, bond->atom_index_b)++;
      audit.total++;
    }

  for (guint i = 0; i < model->atom_count; i++)
    {
      const GdisAtom *atom = g_ptr_array_index(model->atoms, i);
      guint degree = g_array_index(coordination, guint, i);
      g_autofree gchar *coord_key = NULL;

      coord_key = g_strdup_printf("%s:%u",
                                  atom && atom->element && atom->element[0] ? atom->element : "X",
                                  degree);
      increment_counter(audit.coordination_counts, coord_key);
    }

  g_print("Model: %s\n", argv[1]);
  g_print("Atoms: %u\n", model->atom_count);
  g_print("Bonds: %u\n", audit.total);
  g_print("Space group: %s\n", model->space_group ? model->space_group : "");
  print_sorted_counts("Bond pairs:", audit.pair_counts);
  print_sorted_counts("Coordination:", audit.coordination_counts);

  g_hash_table_unref(audit.pair_counts);
  g_hash_table_unref(audit.coordination_counts);
  g_array_free(coordination, TRUE);
  return 0;
}
