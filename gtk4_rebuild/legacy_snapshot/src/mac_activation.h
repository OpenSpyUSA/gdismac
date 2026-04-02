#ifndef GDIS_MAC_ACTIVATION_H
#define GDIS_MAC_ACTIVATION_H

#include <glib.h>

#ifdef __APPLE__
void mac_activate_foreground_app_now(void);
gboolean mac_activate_foreground_app(gpointer data);
#endif

#endif
