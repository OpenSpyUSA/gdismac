#pragma once

#include <gtk/gtk.h>

G_BEGIN_DECLS

void gdis_macos_menu_install(GtkApplication *app);
gchar *gdis_macos_choose_directory(const char *title, const char *initial_dir);

G_END_DECLS
