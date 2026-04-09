#include <gtk/gtk.h>

#include "gdis_gtk4_window.h"

static const char *const WINDOW_DATA_KEY = "gdis-gtk4-window";
static gboolean startup_action_consumed = FALSE;
static gboolean startup_export_consumed = FALSE;
static gboolean forced_startup_paths_consumed = FALSE;
static gchar **forced_startup_paths = NULL;

typedef struct
{
  GtkWindow *window;
} StartupExportRequest;

static GdisGtk4Window *lookup_window_wrapper(GtkApplication *app);

static gboolean env_flag_enabled(const gchar *value)
{
  if (!value || !value[0])
    return FALSE;

  return g_ascii_strcasecmp(value, "1") == 0 ||
         g_ascii_strcasecmp(value, "true") == 0 ||
         g_ascii_strcasecmp(value, "yes") == 0 ||
         g_ascii_strcasecmp(value, "on") == 0;
}

static gchar **build_forced_startup_paths(int argc, char **argv)
{
  GPtrArray *paths;

  paths = g_ptr_array_new_with_free_func(g_free);
  for (gint i = 1; i < argc; i++)
    {
      if (!argv[i] || !argv[i][0])
        continue;
      if (argv[i][0] == '-' && argv[i][1] != '\0')
        continue;
      g_ptr_array_add(paths, g_strdup(argv[i]));
    }

  if (paths->len == 0u)
    {
      g_ptr_array_free(paths, TRUE);
      return NULL;
    }

  g_ptr_array_add(paths, NULL);
  return (gchar **) g_ptr_array_free(paths, FALSE);
}

static gboolean maybe_quit_if_no_primary_windows_idle(gpointer user_data)
{
  GtkApplication *app;

  app = GTK_APPLICATION(user_data);
  /*
   * Tool dialogs are associated with the GtkApplication on macOS so they
   * integrate with the app menu and focus handling. They must not control
   * whether the whole app exits when closed.
   */
  if (GTK_IS_APPLICATION(app) && lookup_window_wrapper(app) == NULL)
    g_application_quit(G_APPLICATION(app));

  return G_SOURCE_REMOVE;
}

static GdisGtk4Window *lookup_window_wrapper(GtkApplication *app)
{
  GtkWindow *window;
  GList *windows;

  window = gtk_application_get_active_window(app);
  if (window)
    {
      GdisGtk4Window *wrapper;

      wrapper = g_object_get_data(G_OBJECT(window), WINDOW_DATA_KEY);
      if (wrapper)
        return wrapper;
    }

  windows = gtk_application_get_windows(app);
  for (GList *iter = windows; iter != NULL; iter = iter->next)
    {
      GdisGtk4Window *wrapper;

      wrapper = g_object_get_data(G_OBJECT(iter->data), WINDOW_DATA_KEY);
      if (wrapper)
        return wrapper;
    }

  return NULL;
}

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

static gboolean run_startup_export_idle(gpointer user_data)
{
  StartupExportRequest *request;
  GdisGtk4Window *window;
  GError *error;

  request = user_data;
  window = NULL;
  error = NULL;
  if (request && GTK_IS_WINDOW(request->window))
    window = g_object_get_data(G_OBJECT(request->window), WINDOW_DATA_KEY);

  if (window && !gdis_gtk4_window_run_startup_export(window, &error))
    {
      g_printerr("Startup export failed: %s\n",
                 error ? error->message : "unknown error");
      g_clear_error(&error);
    }

  if (request)
    {
      g_clear_object(&request->window);
      g_free(request);
    }

  return G_SOURCE_REMOVE;
}

static void maybe_queue_startup_export(GdisGtk4Window *window)
{
  const gchar *export_path;
  StartupExportRequest *request;
  GtkWindow *gtk_window;

  if (startup_export_consumed || !window)
    return;

  export_path = g_getenv("GDIS_GTK4_STARTUP_EXPORT_PNG");
  if (!export_path || !export_path[0])
    return;

  gtk_window = gdis_gtk4_window_get_window(window);
  if (!GTK_IS_WINDOW(gtk_window))
    return;

  startup_export_consumed = TRUE;
  request = g_new0(StartupExportRequest, 1);
  request->window = g_object_ref(gtk_window);
  g_idle_add_full(G_PRIORITY_DEFAULT_IDLE,
                  run_startup_export_idle,
                  request,
                  NULL);
}

static void maybe_open_forced_startup_paths(GdisGtk4Window *window)
{
  if (forced_startup_paths_consumed || !window || !forced_startup_paths)
    return;

  forced_startup_paths_consumed = TRUE;
  gdis_gtk4_window_open_startup_paths(window,
                                      (const char *const *) forced_startup_paths);
}

static gboolean ensure_initial_window_idle(gpointer user_data)
{
  GtkApplication *app;
  GdisGtk4Window *window;

  app = GTK_APPLICATION(user_data);
  if (!GTK_IS_APPLICATION(app))
    return G_SOURCE_REMOVE;

  window = lookup_window_wrapper(app);
  if (!window)
    window = gdis_gtk4_window_new(app);
  if (window)
    {
      maybe_open_forced_startup_paths(window);
      gdis_gtk4_window_present(window);
      maybe_activate_startup_action(window);
      maybe_queue_startup_export(window);
    }

  return G_SOURCE_REMOVE;
}

static void on_startup(GApplication *app, gpointer user_data)
{
  (void) user_data;

  g_idle_add_full(G_PRIORITY_DEFAULT_IDLE,
                  ensure_initial_window_idle,
                  g_object_ref(app),
                  (GDestroyNotify) g_object_unref);
}

static void on_window_removed(GtkApplication *app,
                              GtkWindow *window,
                              gpointer user_data)
{
  (void) window;
  (void) user_data;

  g_idle_add_full(G_PRIORITY_DEFAULT_IDLE,
                  maybe_quit_if_no_primary_windows_idle,
                  g_object_ref(app),
                  (GDestroyNotify) g_object_unref);
}

static void on_activate(GtkApplication *app, gpointer user_data)
{
  GdisGtk4Window *window;

  (void) user_data;

  window = lookup_window_wrapper(app);
  if (!window)
    window = gdis_gtk4_window_new(app);
  maybe_open_forced_startup_paths(window);
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

  window = lookup_window_wrapper(GTK_APPLICATION(app));
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
  maybe_queue_startup_export(window);
}

int main(int argc, char **argv)
{
  GtkApplication *app;
  GApplicationFlags flags;
  int status;

  if (g_getenv("GDIS_GTK4_STARTUP_EXPORT_PNG") ||
      env_flag_enabled(g_getenv("GDIS_GTK4_FORCE_STARTUP_PATHS")))
    forced_startup_paths = build_forced_startup_paths(argc, argv);

  flags = G_APPLICATION_DEFAULT_FLAGS | G_APPLICATION_HANDLES_OPEN;
  if (env_flag_enabled(g_getenv("GDIS_GTK4_NON_UNIQUE")))
    flags |= G_APPLICATION_NON_UNIQUE;

  app = gtk_application_new("org.gdis.Gtk4Rebuild", flags);
  g_signal_connect(app, "startup", G_CALLBACK(on_startup), NULL);
  g_signal_connect(app, "window-removed", G_CALLBACK(on_window_removed), NULL);
  g_signal_connect(app, "activate", G_CALLBACK(on_activate), NULL);
  g_signal_connect(app, "open", G_CALLBACK(on_open), NULL);

  status = g_application_run(G_APPLICATION(app), argc, argv);
  g_object_unref(app);
  g_strfreev(forced_startup_paths);

  return status;
}
