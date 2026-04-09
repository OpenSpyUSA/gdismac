#ifdef __APPLE__

#import <AppKit/AppKit.h>

#include "gdis_macos_menu.h"

enum
{
  GDIS_QBOX_MENU_TAG = 46051
};

static const char *const GDIS_MACOS_MENU_INSTALL_KEY = "gdis-macos-menu-install-key";

@interface GdisMacosMenuTarget : NSObject
{
  GtkApplication *_app;
}

- (instancetype)initWithApp:(GtkApplication *)app;
- (void)activateActionNamed:(const char *)action_name;
- (void)openQbox:(id)sender;
- (void)useLastXml:(id)sender;
- (void)continueLast:(id)sender;
- (void)importResult:(id)sender;
- (void)showResults:(id)sender;

@end

@implementation GdisMacosMenuTarget

- (instancetype)initWithApp:(GtkApplication *)app
{
  self = [super init];
  if (self)
    _app = app;
  return self;
}

- (void)activateActionNamed:(const char *)action_name
{
  if (!_app || !action_name)
    return;

  g_action_group_activate_action(G_ACTION_GROUP(_app), action_name, NULL);
}

- (void)openQbox:(id)sender
{
  (void) sender;
  [self activateActionNamed:"qbox"];
}

- (void)useLastXml:(id)sender
{
  (void) sender;
  [self activateActionNamed:"qbox-use-last-xml"];
}

- (void)continueLast:(id)sender
{
  (void) sender;
  [self activateActionNamed:"qbox-continue-last"];
}

- (void)importResult:(id)sender
{
  (void) sender;
  [self activateActionNamed:"qbox-import-result"];
}

- (void)showResults:(id)sender
{
  (void) sender;
  [self activateActionNamed:"qbox-results"];
}

@end

static GdisMacosMenuTarget *gdis_qbox_menu_target = nil;

typedef struct
{
  GtkApplication *app;
  guint attempts_remaining;
} GdisMacosMenuInstallState;

static void
gdis_macos_menu_add_item(NSMenu *menu,
                         NSString *title,
                         SEL action,
                         id target)
{
  NSMenuItem *item;

  item = [[NSMenuItem alloc] initWithTitle:title action:action keyEquivalent:@""];
  [item setTarget:target];
  [menu addItem:item];
}

static NSInteger
gdis_macos_menu_find_insert_index(NSMenu *main_menu)
{
  NSInteger index;

  index = [main_menu indexOfItemWithTitle:@"Window"];
  if (index >= 0)
    return index;

  index = [main_menu indexOfItemWithTitle:@"Help"];
  if (index >= 0)
    return index;

  return [main_menu numberOfItems];
}

static void
gdis_macos_menu_remove_qbox_item_at_index(NSMenu *main_menu,
                                          NSInteger index)
{
  if (!main_menu)
    return;

  if (index < 0 || index >= [main_menu numberOfItems])
    return;

  [main_menu removeItemAtIndex:index];
}

static gboolean
gdis_macos_menu_item_is_qbox(NSMenuItem *item)
{
  NSString *title;

  if (!item)
    return FALSE;

  title = [item title];
  return [item tag] == GDIS_QBOX_MENU_TAG || (title && [title isEqualToString:@"Qbox"]);
}

static gboolean
gdis_macos_menu_item_is_foreign_qbox(NSMenuItem *item)
{
  NSString *title;

  if (!item)
    return FALSE;

  title = [item title];
  return [item tag] != GDIS_QBOX_MENU_TAG && title && [title isEqualToString:@"Qbox"];
}

static void
gdis_macos_menu_remove_existing_qbox(NSMenu *main_menu)
{
  NSInteger index;

  for (index = [main_menu numberOfItems] - 1; index >= 0; index--)
    {
      NSMenuItem *item;

      item = [main_menu itemAtIndex:index];
      if (gdis_macos_menu_item_is_qbox(item))
        gdis_macos_menu_remove_qbox_item_at_index(main_menu, index);
    }
}

static gboolean
gdis_macos_menu_preserve_exported_qbox(NSMenu *main_menu)
{
  NSInteger index;
  NSInteger first_foreign_index;
  gboolean foreign_found;

  foreign_found = FALSE;
  first_foreign_index = -1;

  for (index = [main_menu numberOfItems] - 1; index >= 0; index--)
    {
      NSMenuItem *item;

      item = [main_menu itemAtIndex:index];
      if (!gdis_macos_menu_item_is_qbox(item))
        continue;

      if (gdis_macos_menu_item_is_foreign_qbox(item))
        {
          if (!foreign_found)
            {
              foreign_found = TRUE;
              first_foreign_index = index;
            }
          else
            {
              gdis_macos_menu_remove_qbox_item_at_index(main_menu, index);
            }
          continue;
        }

      if (foreign_found)
        gdis_macos_menu_remove_qbox_item_at_index(main_menu, index);
    }

  if (!foreign_found)
    return FALSE;

  for (index = [main_menu numberOfItems] - 1; index > first_foreign_index; index--)
    {
      NSMenuItem *item;

      item = [main_menu itemAtIndex:index];
      if (gdis_macos_menu_item_is_foreign_qbox(item))
        gdis_macos_menu_remove_qbox_item_at_index(main_menu, index);
    }

  return TRUE;
}

static void
gdis_macos_menu_install_now(GtkApplication *app)
{
  NSMenu *main_menu;
  NSMenu *qbox_menu;
  NSMenuItem *root_item;
  NSInteger insert_index;

  if (!GTK_IS_APPLICATION(app))
    return;

  main_menu = [NSApp mainMenu];
  if (!main_menu)
    return;

  if (!gdis_qbox_menu_target)
    gdis_qbox_menu_target = [[GdisMacosMenuTarget alloc] initWithApp:app];

  if (gdis_macos_menu_preserve_exported_qbox(main_menu))
    return;

  gdis_macos_menu_remove_existing_qbox(main_menu);

  qbox_menu = [[NSMenu alloc] initWithTitle:@"Qbox"];
  gdis_macos_menu_add_item(qbox_menu, @"Open Qbox", @selector(openQbox:), gdis_qbox_menu_target);
  [qbox_menu addItem:[NSMenuItem separatorItem]];
  gdis_macos_menu_add_item(qbox_menu, @"Use Last XML", @selector(useLastXml:), gdis_qbox_menu_target);
  gdis_macos_menu_add_item(qbox_menu, @"Continue Last", @selector(continueLast:), gdis_qbox_menu_target);
  gdis_macos_menu_add_item(qbox_menu, @"Import Result", @selector(importResult:), gdis_qbox_menu_target);
  gdis_macos_menu_add_item(qbox_menu, @"Results", @selector(showResults:), gdis_qbox_menu_target);

  root_item = [[NSMenuItem alloc] initWithTitle:@"Qbox" action:nil keyEquivalent:@""];
  [root_item setTag:GDIS_QBOX_MENU_TAG];
  [root_item setSubmenu:qbox_menu];

  insert_index = gdis_macos_menu_find_insert_index(main_menu);
  [main_menu insertItem:root_item atIndex:insert_index];
}

static void
gdis_macos_menu_install_state_free(gpointer data)
{
  GdisMacosMenuInstallState *state;

  state = (GdisMacosMenuInstallState *) data;
  if (!state)
    return;

  g_clear_object(&state->app);
  g_free(state);
}

static gboolean
gdis_macos_menu_install_tick(gpointer data)
{
  GdisMacosMenuInstallState *state;

  state = (GdisMacosMenuInstallState *) data;
  if (!state || !GTK_IS_APPLICATION(state->app))
    return G_SOURCE_REMOVE;

  gdis_macos_menu_install_now(state->app);

  if (state->attempts_remaining > 0)
    state->attempts_remaining--;

  return state->attempts_remaining > 0 ? G_SOURCE_CONTINUE : G_SOURCE_REMOVE;
}

void
gdis_macos_menu_install(GtkApplication *app)
{
  GdisMacosMenuInstallState *state;

  g_return_if_fail(GTK_IS_APPLICATION(app));

  if (g_object_get_data(G_OBJECT(app), GDIS_MACOS_MENU_INSTALL_KEY))
    return;

  g_object_set_data(G_OBJECT(app),
                    GDIS_MACOS_MENU_INSTALL_KEY,
                    GINT_TO_POINTER(1));

  gdis_macos_menu_install_now(app);

  state = g_new0(GdisMacosMenuInstallState, 1);
  state->app = g_object_ref(app);
  state->attempts_remaining = 20;
  g_timeout_add_full(G_PRIORITY_DEFAULT,
                     200,
                     gdis_macos_menu_install_tick,
                     state,
                     gdis_macos_menu_install_state_free);
}

gchar *
gdis_macos_choose_directory(const char *title, const char *initial_dir)
{
  @autoreleasepool
    {
      NSOpenPanel *panel;
      NSURL *directory_url;
      NSString *panel_title;
      NSURL *selected_url;
      NSString *selected_path;

      panel = [NSOpenPanel openPanel];
      panel_title = (title && title[0]) ? [NSString stringWithUTF8String:title]
                                        : @"Open all models in folder";
      [panel setTitle:panel_title];
      [panel setPrompt:@"Select"];
      [panel setCanChooseDirectories:YES];
      [panel setCanChooseFiles:NO];
      [panel setAllowsMultipleSelection:NO];
      [panel setResolvesAliases:YES];

      if (initial_dir && initial_dir[0])
        {
          NSString *path;

          path = [NSString stringWithUTF8String:initial_dir];
          directory_url = [NSURL fileURLWithPath:path isDirectory:YES];
          if (directory_url)
            [panel setDirectoryURL:directory_url];
        }

      if ([panel runModal] != NSModalResponseOK)
        return NULL;

      selected_url = [[panel URLs] firstObject];
      if (!selected_url)
        return NULL;

      selected_path = [selected_url path];
      if (!selected_path)
        return NULL;

      return g_strdup([selected_path UTF8String]);
    }
}

#endif
