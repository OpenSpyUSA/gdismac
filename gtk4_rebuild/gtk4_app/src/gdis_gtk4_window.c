#include "gdis_gtk4_window.h"

#include <stdarg.h>

#include "gdis_legacy_map.h"

struct _GdisGtk4Window
{
  GtkApplication *app;
  GtkWidget *window;
  GtkTextBuffer *status_buffer;
  GtkWidget *model_list;
  GtkWidget *selection_drop_down;
};

static const char *const WINDOW_DATA_KEY = "gdis-gtk4-window";

static GdisGtk4Window *gdis_gtk4_window_lookup(GtkApplication *app)
{
  GtkWindow *window;

  window = gtk_application_get_active_window(app);
  if (!window)
    return NULL;

  return g_object_get_data(G_OBJECT(window), WINDOW_DATA_KEY);
}

static void gdis_gtk4_window_log(GdisGtk4Window *self, const char *format, ...)
{
  GtkTextIter end;
  va_list args;
  char *message;

  g_return_if_fail(self != NULL);
  g_return_if_fail(self->status_buffer != NULL);

  va_start(args, format);
  message = g_strdup_vprintf(format, args);
  va_end(args);

  gtk_text_buffer_get_end_iter(self->status_buffer, &end);
  gtk_text_buffer_insert(self->status_buffer, &end, message, -1);

  g_free(message);
}

static void viewer_draw(GtkDrawingArea *area,
                        cairo_t *cr,
                        int width,
                        int height,
                        gpointer data)
{
  (void) area;
  (void) data;

  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_paint(cr);

  cairo_set_source_rgb(cr, 0.20, 0.20, 0.20);
  cairo_set_line_width(cr, 1.0);
  cairo_rectangle(cr, 1.0, 1.0, width - 2.0, height - 2.0);
  cairo_stroke(cr);

  cairo_set_source_rgb(cr, 0.10, 0.35, 0.45);
  cairo_set_line_width(cr, 1.0);
  cairo_move_to(cr, width * 0.5, 20.0);
  cairo_line_to(cr, width * 0.5, height - 20.0);
  cairo_move_to(cr, 20.0, height * 0.5);
  cairo_line_to(cr, width - 20.0, height * 0.5);
  cairo_stroke(cr);

  cairo_set_source_rgb(cr, 0.86, 0.88, 0.90);
  cairo_select_font_face(cr, "SF Mono", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, 16.0);
  cairo_move_to(cr, width * 0.5 - 110.0, height * 0.5 - 12.0);
  cairo_show_text(cr, "GTK4 viewer scaffold");

  cairo_set_font_size(cr, 13.0);
  cairo_move_to(cr, width * 0.5 - 170.0, height * 0.5 + 18.0);
  cairo_show_text(cr, "Renderer port target: legacy display_box -> GTK4");
}

static GtkWidget *new_section_label(const char *text)
{
  GtkWidget *label;

  label = gtk_label_new(text);
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_widget_add_css_class(label, "heading");

  return label;
}

static GtkWidget *new_entry_row_grid(const char *const *fields)
{
  GtkWidget *grid;
  int row;

  grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(grid), 12);
  gtk_grid_set_row_spacing(GTK_GRID(grid), 6);

  for (row = 0; fields[row] != NULL; row++)
    {
      GtkWidget *name;
      GtkWidget *entry;

      name = gtk_label_new(fields[row]);
      gtk_widget_set_halign(name, GTK_ALIGN_START);
      gtk_grid_attach(GTK_GRID(grid), name, 0, row, 1, 1);

      entry = gtk_entry_new();
      gtk_editable_set_editable(GTK_EDITABLE(entry), FALSE);
      gtk_entry_set_placeholder_text(GTK_ENTRY(entry), "Port legacy value binding later");
      gtk_grid_attach(GTK_GRID(grid), entry, 1, row, 1, 1);
    }

  return grid;
}

static GtkWidget *build_content_page(void)
{
  GtkWidget *box;
  GtkWidget *label;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);

  label = gtk_label_new("Legacy content pane placeholder.\nThis will host model tree details, metadata, and content inspection.");
  gtk_label_set_wrap(GTK_LABEL(label), TRUE);
  gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
  gtk_box_append(GTK_BOX(box), label);

  return box;
}

static GtkWidget *build_editing_page(void)
{
  static const char *const fields[] = {
    "Element",
    "Label",
    "x",
    "y",
    "z",
    "Charge",
    "Mass",
    "Occupancy",
    NULL
  };

  GtkWidget *box;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_box_append(GTK_BOX(box), new_section_label("Model : Editing"));
  gtk_box_append(GTK_BOX(box), new_entry_row_grid(fields));

  return box;
}

static GtkWidget *build_display_page(void)
{
  GtkWidget *box;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_box_append(GTK_BOX(box), new_section_label("Model : Display"));
  gtk_box_append(GTK_BOX(box), gtk_check_button_new_with_label("Show atoms"));
  gtk_box_append(GTK_BOX(box), gtk_check_button_new_with_label("Show bonds"));
  gtk_box_append(GTK_BOX(box), gtk_check_button_new_with_label("Show cell"));
  gtk_box_append(GTK_BOX(box), gtk_check_button_new_with_label("Show labels"));

  return box;
}

static GtkWidget *build_images_page(void)
{
  GtkWidget *box;
  GtkWidget *label;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_box_append(GTK_BOX(box), new_section_label("Model : Images"));

  label = gtk_label_new("Supercell, image limits, and periodic image controls will be ported here.");
  gtk_label_set_wrap(GTK_LABEL(label), TRUE);
  gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
  gtk_box_append(GTK_BOX(box), label);

  return box;
}

static GtkWidget *build_symmetry_page(void)
{
  GtkWidget *box;
  GtkWidget *label;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_box_append(GTK_BOX(box), new_section_label("Model : Symmetry"));

  label = gtk_label_new("Space group, symmetry operators, and equivalent positions will be ported here.");
  gtk_label_set_wrap(GTK_LABEL(label), TRUE);
  gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
  gtk_box_append(GTK_BOX(box), label);

  return box;
}

static GtkWidget *build_viewing_page(void)
{
  GtkWidget *box;
  GtkWidget *row;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_box_append(GTK_BOX(box), new_section_label("Model : Viewing"));

  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_box_append(GTK_BOX(row), gtk_button_new_with_label("Reset View"));
  gtk_box_append(GTK_BOX(row), gtk_button_new_with_label("View X"));
  gtk_box_append(GTK_BOX(row), gtk_button_new_with_label("View Y"));
  gtk_box_append(GTK_BOX(row), gtk_button_new_with_label("View Z"));
  gtk_box_append(GTK_BOX(box), row);

  return box;
}

static GtkWidget *build_sidebar_stack(void)
{
  GtkWidget *stack;
  GtkWidget *page;

  stack = gtk_stack_new();
  gtk_stack_set_transition_type(GTK_STACK(stack), GTK_STACK_TRANSITION_TYPE_SLIDE_LEFT_RIGHT);
  gtk_widget_set_vexpand(stack, TRUE);

  page = build_content_page();
  gtk_stack_add_titled(GTK_STACK(stack), page, "content", "Content");

  page = build_editing_page();
  gtk_stack_add_titled(GTK_STACK(stack), page, "editing", "Editing");

  page = build_display_page();
  gtk_stack_add_titled(GTK_STACK(stack), page, "display", "Display");

  page = build_images_page();
  gtk_stack_add_titled(GTK_STACK(stack), page, "images", "Images");

  page = build_symmetry_page();
  gtk_stack_add_titled(GTK_STACK(stack), page, "symmetry", "Symmetry");

  page = build_viewing_page();
  gtk_stack_add_titled(GTK_STACK(stack), page, "viewing", "Viewing");

  return stack;
}

static void on_model_selected(GtkListBox *box, GtkListBoxRow *row, gpointer user_data)
{
  GdisGtk4Window *self;
  const char *name;

  (void) box;

  self = user_data;
  if (!row)
    return;

  name = g_object_get_data(G_OBJECT(row), "model-path");
  if (!name)
    name = "unknown";

  gdis_gtk4_window_log(self, "Selected model placeholder: %s\n", name);
}

static GtkWidget *build_model_list(GdisGtk4Window *self)
{
  GtkWidget *list;
  int i;

  list = gtk_list_box_new();
  gtk_list_box_set_selection_mode(GTK_LIST_BOX(list), GTK_SELECTION_SINGLE);

  for (i = 0; gdis_legacy_model_samples[i] != NULL; i++)
    {
      GtkWidget *row;
      GtkWidget *label;

      row = gtk_list_box_row_new();
      label = gtk_label_new(gdis_legacy_model_samples[i]);
      gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
      gtk_widget_set_margin_start(label, 8);
      gtk_widget_set_margin_end(label, 8);
      gtk_widget_set_margin_top(label, 8);
      gtk_widget_set_margin_bottom(label, 8);
      gtk_list_box_row_set_child(GTK_LIST_BOX_ROW(row), label);
      g_object_set_data_full(G_OBJECT(row), "model-path",
                             g_strdup(gdis_legacy_model_samples[i]), g_free);
      gtk_list_box_append(GTK_LIST_BOX(list), row);
    }

  g_signal_connect(list, "row-selected", G_CALLBACK(on_model_selected), self);

  return list;
}

static GtkWidget *build_sidebar(GdisGtk4Window *self)
{
  GtkWidget *box;
  GtkWidget *scroller;
  GtkWidget *switcher;
  GtkWidget *stack;
  GtkStringList *drop_down_model;
  GtkWidget *footer;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12);
  gtk_widget_set_size_request(box, 360, -1);
  gtk_widget_set_margin_start(box, 12);
  gtk_widget_set_margin_end(box, 12);
  gtk_widget_set_margin_top(box, 12);
  gtk_widget_set_margin_bottom(box, 12);

  gtk_box_append(GTK_BOX(box), new_section_label("Models"));

  self->model_list = build_model_list(self);
  scroller = gtk_scrolled_window_new();
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                 GTK_POLICY_AUTOMATIC,
                                 GTK_POLICY_AUTOMATIC);
  gtk_widget_set_vexpand(scroller, TRUE);
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), self->model_list);
  gtk_box_append(GTK_BOX(box), scroller);

  stack = build_sidebar_stack();
  switcher = gtk_stack_switcher_new();
  gtk_stack_switcher_set_stack(GTK_STACK_SWITCHER(switcher), GTK_STACK(stack));
  gtk_box_append(GTK_BOX(box), switcher);
  gtk_box_append(GTK_BOX(box), stack);

  gtk_box_append(GTK_BOX(box), new_section_label("Selection"));
  drop_down_model = gtk_string_list_new(gdis_legacy_selection_modes);
  self->selection_drop_down = gtk_drop_down_new(G_LIST_MODEL(drop_down_model), NULL);
  gtk_box_append(GTK_BOX(box), self->selection_drop_down);
  g_object_unref(drop_down_model);

  footer = gtk_label_new("Legacy reference lives in ../legacy_snapshot/\nTarget: Linux-style single-window GTK4 shell");
  gtk_label_set_wrap(GTK_LABEL(footer), TRUE);
  gtk_label_set_xalign(GTK_LABEL(footer), 0.0f);
  gtk_box_append(GTK_BOX(box), footer);

  return box;
}

static GtkWidget *build_viewer_area(void)
{
  GtkWidget *frame;
  GtkWidget *drawing_area;

  frame = gtk_frame_new(NULL);
  gtk_widget_set_hexpand(frame, TRUE);
  gtk_widget_set_vexpand(frame, TRUE);

  drawing_area = gtk_drawing_area_new();
  gtk_widget_set_hexpand(drawing_area, TRUE);
  gtk_widget_set_vexpand(drawing_area, TRUE);
  gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawing_area), viewer_draw, NULL, NULL);
  gtk_frame_set_child(GTK_FRAME(frame), drawing_area);

  return frame;
}

static GtkWidget *build_status_view(GdisGtk4Window *self)
{
  GtkWidget *scroller;
  GtkWidget *text_view;

  scroller = gtk_scrolled_window_new();
  gtk_widget_set_vexpand(scroller, TRUE);

  text_view = gtk_text_view_new();
  gtk_text_view_set_editable(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_WORD_CHAR);
  gtk_text_view_set_monospace(GTK_TEXT_VIEW(text_view), TRUE);
  self->status_buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));

  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  return scroller;
}

static GtkWidget *build_main_content(GdisGtk4Window *self)
{
  GtkWidget *main_paned;
  GtkWidget *right_paned;
  GtkWidget *viewer;
  GtkWidget *status;

  main_paned = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
  gtk_widget_set_hexpand(main_paned, TRUE);
  gtk_widget_set_vexpand(main_paned, TRUE);
  gtk_paned_set_wide_handle(GTK_PANED(main_paned), TRUE);

  right_paned = gtk_paned_new(GTK_ORIENTATION_VERTICAL);
  gtk_paned_set_wide_handle(GTK_PANED(right_paned), TRUE);

  viewer = build_viewer_area();
  status = build_status_view(self);

  gtk_paned_set_start_child(GTK_PANED(main_paned), build_sidebar(self));
  gtk_paned_set_end_child(GTK_PANED(main_paned), right_paned);
  gtk_paned_set_position(GTK_PANED(main_paned), 360);

  gtk_paned_set_start_child(GTK_PANED(right_paned), viewer);
  gtk_paned_set_end_child(GTK_PANED(right_paned), status);
  gtk_paned_set_position(GTK_PANED(right_paned), 650);

  return main_paned;
}

static void action_activate(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  const char *action_name;

  self = user_data;
  action_name = g_object_get_data(G_OBJECT(button), "action-name");
  if (!action_name)
    return;

  g_action_group_activate_action(G_ACTION_GROUP(self->app), action_name, NULL);
}

static GtkWidget *build_toolbar(GdisGtk4Window *self)
{
  GtkWidget *toolbar;
  int i;

  toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_widget_set_margin_start(toolbar, 8);
  gtk_widget_set_margin_end(toolbar, 8);
  gtk_widget_set_margin_top(toolbar, 8);
  gtk_widget_set_margin_bottom(toolbar, 8);

  for (i = 0; gdis_legacy_toolbar_actions[i][0] != NULL; i++)
    {
      GtkWidget *button;

      button = gtk_button_new_with_label(gdis_legacy_toolbar_actions[i][0]);
      g_object_set_data_full(G_OBJECT(button), "action-name",
                             g_strdup(gdis_legacy_toolbar_actions[i][1]), g_free);
      g_signal_connect(button, "clicked", G_CALLBACK(action_activate), self);
      gtk_box_append(GTK_BOX(toolbar), button);
    }

  return toolbar;
}

static GMenuModel *build_menu_bar_model(void)
{
  GMenu *root;
  GMenu *file_menu;
  GMenu *edit_menu;
  GMenu *tools_menu;
  GMenu *view_menu;
  GMenu *help_menu;

  root = g_menu_new();

  file_menu = g_menu_new();
  g_menu_append(file_menu, "Open", "app.open");
  g_menu_append(file_menu, "Save", "app.save");
  g_menu_append(file_menu, "Quit", "app.quit");
  g_menu_append_submenu(root, "File", G_MENU_MODEL(file_menu));

  edit_menu = g_menu_new();
  g_menu_append(edit_menu, "New Model", "app.new-model");
  g_menu_append(edit_menu, "Edit Structure", "app.edit");
  g_menu_append_submenu(root, "Edit", G_MENU_MODEL(edit_menu));

  tools_menu = g_menu_new();
  g_menu_append(tools_menu, "Render", "app.render");
  g_menu_append(tools_menu, "Measure", "app.measure");
  g_menu_append(tools_menu, "Iso-surfaces", "app.isosurface");
  g_menu_append(tools_menu, "Surface", "app.surface");
  g_menu_append(tools_menu, "Diffraction", "app.diffraction");
  g_menu_append_submenu(root, "Tools", G_MENU_MODEL(tools_menu));

  view_menu = g_menu_new();
  g_menu_append(view_menu, "Reset View", "app.reset-view");
  g_menu_append_submenu(root, "View", G_MENU_MODEL(view_menu));

  help_menu = g_menu_new();
  g_menu_append(help_menu, "About", "app.about");
  g_menu_append_submenu(root, "Help", G_MENU_MODEL(help_menu));

  g_object_unref(file_menu);
  g_object_unref(edit_menu);
  g_object_unref(tools_menu);
  g_object_unref(view_menu);
  g_object_unref(help_menu);

  return G_MENU_MODEL(root);
}

static void action_dispatch(GSimpleAction *action, GVariant *parameter, gpointer user_data)
{
  GtkApplication *app;
  GdisGtk4Window *self;
  const char *name;

  (void) parameter;

  app = GTK_APPLICATION(user_data);
  self = gdis_gtk4_window_lookup(app);
  name = g_action_get_name(G_ACTION(action));

  if (!self || !name)
    return;

  if (g_strcmp0(name, "quit") == 0)
    {
      gdis_gtk4_window_log(self, "Quit requested from GTK4 scaffold.\n");
      g_application_quit(G_APPLICATION(app));
      return;
    }

  if (g_strcmp0(name, "open") == 0)
    {
      gdis_gtk4_window_log(self, "Open requested. File dialog port is next.\n");
      return;
    }

  if (g_strcmp0(name, "save") == 0)
    {
      gdis_gtk4_window_log(self, "Save requested. File export port is next.\n");
      return;
    }

  if (g_strcmp0(name, "about") == 0)
    {
      gdis_gtk4_window_log(self, "GTK4 scaffold for GDIS. Legacy UI reference: gui_main.c.\n");
      return;
    }

  gdis_gtk4_window_log(self, "Activated action: %s\n", name);
}

static void install_actions(GtkApplication *app)
{
  const GActionEntry entries[] = {
    {"open", action_dispatch, NULL, NULL, NULL},
    {"save", action_dispatch, NULL, NULL, NULL},
    {"quit", action_dispatch, NULL, NULL, NULL},
    {"new-model", action_dispatch, NULL, NULL, NULL},
    {"edit", action_dispatch, NULL, NULL, NULL},
    {"render", action_dispatch, NULL, NULL, NULL},
    {"measure", action_dispatch, NULL, NULL, NULL},
    {"isosurface", action_dispatch, NULL, NULL, NULL},
    {"surface", action_dispatch, NULL, NULL, NULL},
    {"diffraction", action_dispatch, NULL, NULL, NULL},
    {"reset-view", action_dispatch, NULL, NULL, NULL},
    {"about", action_dispatch, NULL, NULL, NULL}
  };

  g_action_map_add_action_entries(G_ACTION_MAP(app),
                                  entries,
                                  G_N_ELEMENTS(entries),
                                  app);
}

GdisGtk4Window *gdis_gtk4_window_new(GtkApplication *app)
{
  GdisGtk4Window *self;
  GtkWidget *root;
  GtkWidget *menu_bar;

  g_return_val_if_fail(GTK_IS_APPLICATION(app), NULL);

  self = g_new0(GdisGtk4Window, 1);
  self->app = app;

  install_actions(app);

  self->window = gtk_application_window_new(app);
  gtk_window_set_title(GTK_WINDOW(self->window), "GDIS GTK4");
  gtk_window_set_default_size(GTK_WINDOW(self->window), 1440, 900);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_window_set_child(GTK_WINDOW(self->window), root);

  menu_bar = gtk_popover_menu_bar_new_from_model(build_menu_bar_model());
  gtk_box_append(GTK_BOX(root), menu_bar);
  gtk_box_append(GTK_BOX(root), build_toolbar(self));
  gtk_box_append(GTK_BOX(root), gtk_separator_new(GTK_ORIENTATION_HORIZONTAL));
  gtk_box_append(GTK_BOX(root), build_main_content(self));

  g_object_set_data_full(G_OBJECT(self->window), WINDOW_DATA_KEY, self, g_free);

  gdis_gtk4_window_log(self, "GTK4 scaffold initialized.\n");
  gdis_gtk4_window_log(self, "This shell mirrors the legacy Linux single-window layout.\n");
  gdis_gtk4_window_log(self, "Next step: connect legacy model and renderer code.\n");

  return self;
}

void gdis_gtk4_window_present(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);
  g_return_if_fail(GTK_IS_WINDOW(self->window));

  gtk_window_present(GTK_WINDOW(self->window));
}

GtkWindow *gdis_gtk4_window_get_window(GdisGtk4Window *self)
{
  g_return_val_if_fail(self != NULL, NULL);
  g_return_val_if_fail(GTK_IS_WINDOW(self->window), NULL);

  return GTK_WINDOW(self->window);
}
