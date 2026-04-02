#include <gtk/gtk.h>

#include "gdis_gtk4_window.h"

static void on_activate(GtkApplication *app, gpointer user_data)
{
  GdisGtk4Window *window;

  (void) user_data;

  window = gdis_gtk4_window_new(app);
  gdis_gtk4_window_present(window);
}

int main(int argc, char *argv[])
{
  GtkApplication *app;
  int status;

  app = gtk_application_new("org.gdis.Gtk4", G_APPLICATION_DEFAULT_FLAGS);
  g_signal_connect(app, "activate", G_CALLBACK(on_activate), NULL);

  status = g_application_run(G_APPLICATION(app), argc, argv);
  g_object_unref(app);

  return status;
}
