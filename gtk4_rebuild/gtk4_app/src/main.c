#include <gtk/gtk.h>

#include "gdis_gtk4_window.h"

static const char *const WINDOW_DATA_KEY = "gdis-gtk4-window";
static gboolean startup_action_consumed = FALSE;

static void maybe_activate_startup_action(GdisGtk4Window *window)
{
  const gchar *action_name;

  if (startup_action_consumed || !window)
    return;

  action_name = g_getenv("GDIS_GTK4_STARTUP_ACTION");
  if (!action_name || !action_name[0])
    return;

  startup_action_consumed = TRUE;
  gdis_gtk4_window_activate_action(window, action_name);
}

static GdisGtk4Window *lookup_active_window_wrapper(GtkApplication *app)
{
  GtkWindow *window;

  window = gtk_application_get_active_window(app);
  if (!window)
    return NULL;

  return g_object_get_data(G_OBJECT(window), WINDOW_DATA_KEY);
}

static void on_activate(GtkApplication *app, gpointer user_data)
{
  GdisGtk4Window *window;

  (void) user_data;

  window = lookup_active_window_wrapper(app);
  if (!window)
    window = gdis_gtk4_window_new(app);
  gdis_gtk4_window_present(window);
  maybe_activate_startup_action(window);
}

static void on_open(GApplication *app,
                    GFile **files,
                    gint n_files,
                    const gchar *hint,
                    gpointer user_data)
{
  GdisGtk4Window *window;
  GPtrArray *paths;
  gint i;

  (void) hint;
  (void) user_data;

  window = lookup_active_window_wrapper(GTK_APPLICATION(app));
  if (!window)
    window = gdis_gtk4_window_new(GTK_APPLICATION(app));

  paths = g_ptr_array_new_with_free_func(g_free);
  for (i = 0; i < n_files; i++)
    {
      gchar *path;

      path = g_file_get_path(files[i]);
      if (path)
        g_ptr_array_add(paths, path);
    }
  g_ptr_array_add(paths, NULL);

  gdis_gtk4_window_open_startup_paths(window, (const char *const *) paths->pdata);
  g_ptr_array_free(paths, TRUE);
  gdis_gtk4_window_present(window);
  maybe_activate_startup_action(window);
}

int main(int argc, char **argv)
{
  GtkApplication *app;
  int status;

  app = gtk_application_new("org.gdis.Gtk4Rebuild",
                            G_APPLICATION_DEFAULT_FLAGS | G_APPLICATION_HANDLES_OPEN);
  g_signal_connect(app, "activate", G_CALLBACK(on_activate), NULL);
  g_signal_connect(app, "open", G_CALLBACK(on_open), NULL);

  status = g_application_run(G_APPLICATION(app), argc, argv);
  g_object_unref(app);

  return status;
}
