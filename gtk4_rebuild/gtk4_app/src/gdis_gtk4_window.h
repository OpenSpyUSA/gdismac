#pragma once

#include <gtk/gtk.h>

typedef struct _GdisGtk4Window GdisGtk4Window;

GdisGtk4Window *gdis_gtk4_window_new(GtkApplication *app);
void gdis_gtk4_window_present(GdisGtk4Window *self);
GtkWindow *gdis_gtk4_window_get_window(GdisGtk4Window *self);
void gdis_gtk4_window_activate_action(GdisGtk4Window *self, const char *action_name);
void gdis_gtk4_window_open_startup_paths(GdisGtk4Window *self,
                                         const char *const *paths);
gboolean gdis_gtk4_window_run_startup_export(GdisGtk4Window *self,
                                             GError **error);
