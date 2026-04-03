#include "gdis_gtk4_window.h"

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "gdis_legacy_map.h"
#include "gdis_model.h"
#include "gdis_reports.h"

typedef struct
{
  guint atom_index;
  guint cell_index;
  gint image_offset[3];
  gdouble screen_x;
  gdouble screen_y;
  gdouble depth;
  gdouble radius;
} GdisProjectedAtom;

typedef struct
{
  gdouble x;
  gdouble y;
  gdouble z;
} GdisVec3;

typedef struct _GdisMeasureTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkTextBuffer *buffer;
  GdisMeasureMode mode;
} GdisMeasureTool;

typedef struct _GdisEditTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkTextBuffer *buffer;
  GtkWidget *label_entry;
  GtkWidget *element_entry;
  GtkWidget *x_entry;
  GtkWidget *y_entry;
  GtkWidget *z_entry;
  GtkWidget *bond_order_entry;
} GdisEditTool;

typedef struct _GdisDiffractionTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkWidget *plot_window;
  GtkWidget *plot_area;
  GtkWidget *execute_button;
  GtkWidget *radiation_dropdown;
  GtkWidget *broadening_dropdown;
  GtkWidget *wavelength_spin;
  GtkWidget *asym_spin;
  GtkWidget *theta_min_spin;
  GtkWidget *theta_max_spin;
  GtkWidget *theta_step_spin;
  GtkWidget *u_spin;
  GtkWidget *v_spin;
  GtkWidget *w_spin;
  GtkWidget *output_entry;
  GtkTextBuffer *report_buffer;
  GdisDiffractionPattern *pattern;
  gchar *plot_title;
  gchar *last_export_path;
} GdisDiffractionTool;

typedef struct
{
  struct _GdisGtk4Window *owner;
  GdisModel *model;
} GdisSaveDialogContext;

typedef enum
{
  GDIS_SELECTION_MODE_ATOMS = 0,
  GDIS_SELECTION_MODE_LABEL,
  GDIS_SELECTION_MODE_FF_TYPE,
  GDIS_SELECTION_MODE_ELEMENTS,
  GDIS_SELECTION_MODE_ELEMENTS_IN_MOLECULE,
  GDIS_SELECTION_MODE_MOLECULES,
  GDIS_SELECTION_MODE_MOLECULE_FRAGMENTS,
  GDIS_SELECTION_MODE_REGIONS
} GdisSelectionMode;

typedef enum
{
  GDIS_CLICK_MODE_SELECT = 0,
  GDIS_CLICK_MODE_ADD_BOND,
  GDIS_CLICK_MODE_REMOVE_BOND
} GdisClickMode;

struct _GdisGtk4Window
{
  GtkApplication *app;
  GtkWidget *window;
  GtkWidget *main_paned;
  GtkWidget *right_paned;
  GtkTextBuffer *status_buffer;
  GtkWidget *model_list;
  GtkTextBuffer *active_summary_buffer;
  GtkTextBuffer *content_buffer;
  GtkTextBuffer *editing_buffer;
  GtkTextBuffer *images_buffer;
  GtkTextBuffer *symmetry_buffer;
  GtkWidget *image_limit_spin[6];
  GtkWidget *supercell_repeat_spin[3];
  GtkWidget *viewer_area;
  GtkWidget *show_atoms_toggle;
  GtkWidget *show_bonds_toggle;
  GtkWidget *show_cell_toggle;
  GtkWidget *show_labels_toggle;
  GPtrArray *models;
  GdisModel *active_model;
  guint selected_atom_index;
  GArray *selected_atoms;
  GArray *picked_atoms;
  GdisMeasureTool *measure_tool;
  GdisEditTool *edit_tool;
  GdisDiffractionTool *diffraction_tool;
  GdisSelectionMode selection_mode;
  GdisClickMode click_mode;
  guint fragment_anchor_index;
  guint untitled_counter;
  guint layout_restore_tick_id;
  gdouble rotation_x;
  gdouble rotation_y;
  gdouble drag_origin_x;
  gdouble drag_origin_y;
  gdouble press_x;
  gdouble press_y;
  gdouble zoom;
};

static const char *const WINDOW_DATA_KEY = "gdis-gtk4-window";
static const guint INVALID_ATOM_INDEX = G_MAXUINT;

static GtkWidget *gdis_gtk4_window_find_model_button(GdisGtk4Window *self,
                                                     const char *path);
static GdisModel *gdis_gtk4_window_find_model(GdisGtk4Window *self,
                                              const char *path);
static GtkWidget *gdis_gtk4_window_add_loaded_model(GdisGtk4Window *self,
                                                    GdisModel *model,
                                                    gboolean select_after_add);
static void gdis_gtk4_window_remove_model(GdisGtk4Window *self,
                                          GdisModel *model);
static gboolean gdis_gtk4_window_add_model_from_path(GdisGtk4Window *self,
                                                     const char *path,
                                                     gboolean select_after_add);
static guint gdis_gtk4_window_get_image_axis_limits(const GdisModel *model,
                                                    gint negative[3],
                                                    gint positive[3]);
static void gdis_gtk4_window_set_active_model(GdisGtk4Window *self,
                                              GdisModel *model);
static void gdis_gtk4_window_update_details(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_viewer(GdisGtk4Window *self);
static void gdis_gtk4_window_sync_image_controls(GdisGtk4Window *self);
static void gdis_gtk4_window_reset_view(GdisGtk4Window *self);
static void gdis_gtk4_window_apply_axis_view(GdisGtk4Window *self,
                                             gdouble rotation_x,
                                             gdouble rotation_y);
static void gdis_gtk4_window_remember_atom_pick(GdisGtk4Window *self, guint atom_index);
static void gdis_gtk4_window_clear_atom_picks(GdisGtk4Window *self);
static void gdis_gtk4_window_clear_selected_atoms(GdisGtk4Window *self);
static void gdis_gtk4_window_apply_selection_mode(GdisGtk4Window *self, guint atom_index);
static void gdis_gtk4_window_set_click_mode(GdisGtk4Window *self,
                                            GdisClickMode click_mode);
static void gdis_gtk4_window_collect_connected_atoms(const GdisModel *model,
                                                     guint start_index,
                                                     GArray *atoms_out);
static gboolean gdis_gtk4_window_collect_fragment_path(const GdisModel *model,
                                                       guint start_index,
                                                       guint end_index,
                                                       GArray *atoms_out);
static gboolean gdis_model_build_cell_vectors(const GdisModel *model,
                                              gdouble a_vec[3],
                                              gdouble b_vec[3],
                                              gdouble c_vec[3]);
static void gdis_gtk4_window_refresh_model_buttons(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_measure_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_edit_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_diffraction_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_after_model_edit(GdisGtk4Window *self,
                                                      gboolean refresh_model_buttons);
static void gdis_gtk4_window_restore_layout(GdisGtk4Window *self);
static void gdis_gtk4_window_present_measure_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_edit_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_diffraction_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_report(GdisGtk4Window *self,
                                            const char *title,
                                            const char *body);
static gboolean gdis_gtk4_window_model_supports_diffraction(const GdisModel *model,
                                                            gchar **message_out);
static gboolean gdis_gtk4_window_parse_entry_double(GtkWidget *entry,
                                                    const char *field_name,
                                                    gdouble *value_out,
                                                    GError **error);
static gboolean gdis_gtk4_window_parse_entry_uint(GtkWidget *entry,
                                                  const char *field_name,
                                                  guint *value_out,
                                                  GError **error);
static void gdis_gtk4_window_configure_button_label(GtkWidget *button);
static gboolean gdis_gtk4_window_restore_layout_tick(GtkWidget *widget,
                                                     GdkFrameClock *frame_clock,
                                                     gpointer user_data);
static gboolean gdis_gtk4_window_get_last_two_picks(GdisGtk4Window *self,
                                                    guint *atom_index_a,
                                                    guint *atom_index_b,
                                                    GError **error);
static gboolean gdis_gtk4_window_atom_array_contains(const GArray *array,
                                                     guint atom_index);
static const char *gdis_gtk4_window_selection_mode_label(GdisSelectionMode mode);
static const char *gdis_gtk4_window_click_mode_label(GdisClickMode mode);
static const GdisAtom *gdis_gtk4_window_get_selected_atom(const GdisGtk4Window *self);
static void gdis_gtk4_window_clear_status_log(GdisGtk4Window *self);
static void on_model_button_clicked(GtkButton *button, gpointer user_data);
static void on_selection_button_clicked(GtkButton *button, gpointer user_data);
static void on_view_toggle_toggled(GtkToggleButton *button, gpointer user_data);
static void on_reset_view_clicked(GtkButton *button, gpointer user_data);
static void on_view_x_clicked(GtkButton *button, gpointer user_data);
static void on_view_y_clicked(GtkButton *button, gpointer user_data);
static void on_view_z_clicked(GtkButton *button, gpointer user_data);
static void on_apply_image_limits_clicked(GtkButton *button, gpointer user_data);
static void on_reset_image_limits_clicked(GtkButton *button, gpointer user_data);
static void on_confine_to_cell_clicked(GtkButton *button, gpointer user_data);
static void on_confine_molecules_clicked(GtkButton *button, gpointer user_data);
static void on_force_p1_clicked(GtkButton *button, gpointer user_data);
static void on_make_supercell_clicked(GtkButton *button, gpointer user_data);
static void on_viewer_drag_begin(GtkGestureDrag *gesture,
                                 gdouble start_x,
                                 gdouble start_y,
                                 gpointer user_data);
static void on_viewer_drag_update(GtkGestureDrag *gesture,
                                  gdouble offset_x,
                                  gdouble offset_y,
                                  gpointer user_data);
static gboolean on_viewer_scroll(GtkEventControllerScroll *controller,
                                 gdouble dx,
                                 gdouble dy,
                                 gpointer user_data);
static void on_viewer_click_pressed(GtkGestureClick *gesture,
                                    gint n_press,
                                    gdouble x,
                                    gdouble y,
                                    gpointer user_data);
static void on_viewer_click_released(GtkGestureClick *gesture,
                                     gint n_press,
                                     gdouble x,
                                     gdouble y,
                                     gpointer user_data);
static GtkWidget *new_toggle_button(const char *label, gboolean active);
static void present_save_dialog(GdisGtk4Window *self, GdisModel *model);
static void on_save_dialog_complete(GObject *source_object,
                                    GAsyncResult *result,
                                    gpointer user_data);

static void
gdis_gtk4_window_install_css(GtkWidget *widget)
{
  GtkCssProvider *provider;
  GdkDisplay *display;
  const char *css =
    ".gdis-root {"
    "  background: #e8edf2;"
    "}"
    ".gdis-sidebar {"
    "  background: #f5f8fb;"
    "}"
    ".gdis-section-button {"
    "  font-weight: 700;"
    "  color: #18212b;"
    "  padding: 0;"
    "}"
    ".gdis-sidebar-card {"
    "  background: #f4f7fa;"
    "  border: 1px solid #d9e1e8;"
    "  border-radius: 10px;"
    "}"
    ".gdis-sidebar-panel {"
    "  background: #fbfcfd;"
    "  border: 1px solid #dde4eb;"
    "  border-radius: 10px;"
    "}"
    ".gdis-sidebar-stack {"
    "  padding: 10px;"
    "}"
    ".gdis-model-button {"
    "  padding: 8px 10px;"
    "  border-radius: 8px;"
    "  color: #18212b;"
    "}"
    ".gdis-model-button-active {"
    "  background: #dce9f2;"
    "  color: #10202f;"
    "}"
    ".gdis-muted {"
    "  color: #5b6672;"
    "}"
    ".gdis-sidebar textview,"
    ".gdis-sidebar textview text,"
    ".gdis-sidebar button,"
    ".gdis-sidebar togglebutton,"
    ".gdis-sidebar switcher button,"
    ".gdis-toolbar button,"
    "popovermenubar,"
    "popovermenubar > box > button {"
    "  color: #18212b;"
    "}"
    ".gdis-sidebar textview,"
    ".gdis-sidebar textview text {"
    "  background: #ffffff;"
    "}"
    ".gdis-viewer-frame {"
    "  background: #05080b;"
    "}"
    ".gdis-status textview,"
    ".gdis-status textview text {"
    "  background: #fbfcfd;"
    "  color: #17212b;"
    "}";

  g_return_if_fail(GTK_IS_WIDGET(widget));

  display = gtk_widget_get_display(widget);
  provider = gtk_css_provider_new();
  gtk_css_provider_load_from_string(provider, css);
  gtk_style_context_add_provider_for_display(display,
                                             GTK_STYLE_PROVIDER(provider),
                                             GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
  g_object_unref(provider);
}

static void
gdis_gtk4_window_free(gpointer data)
{
  GdisGtk4Window *self;

  self = data;
  if (!self)
    return;

  if (self->models)
    g_ptr_array_free(self->models, TRUE);
  if (self->selected_atoms)
    g_array_free(self->selected_atoms, TRUE);
  if (self->picked_atoms)
    g_array_free(self->picked_atoms, TRUE);
  if (self->measure_tool)
    self->measure_tool->owner = NULL;
  if (self->edit_tool)
    self->edit_tool->owner = NULL;
  if (self->diffraction_tool)
    self->diffraction_tool->owner = NULL;
  g_free(self);
}

static GdisGtk4Window *
gdis_gtk4_window_lookup(GtkApplication *app)
{
  GtkWindow *window;

  window = gtk_application_get_active_window(app);
  if (!window)
    return NULL;

  return g_object_get_data(G_OBJECT(window), WINDOW_DATA_KEY);
}

static void
gdis_gtk4_window_log(GdisGtk4Window *self, const char *format, ...)
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

static void
gdis_gtk4_window_clear_status_log(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);
  g_return_if_fail(self->status_buffer != NULL);

  gtk_text_buffer_set_text(self->status_buffer, "", -1);
}

static void
gdis_text_buffer_set(GtkTextBuffer *buffer, const char *format, ...)
{
  va_list args;
  gchar *message;

  g_return_if_fail(buffer != NULL);
  g_return_if_fail(format != NULL);

  va_start(args, format);
  message = g_strdup_vprintf(format, args);
  va_end(args);

  gtk_text_buffer_set_text(buffer, message, -1);
  g_free(message);
}

static void
gdis_gtk4_window_refresh_viewer(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);

  if (self->viewer_area)
    gtk_widget_queue_draw(self->viewer_area);
}

static void
gdis_gtk4_window_refresh_model_buttons(GdisGtk4Window *self)
{
  GtkWidget *button;

  g_return_if_fail(self != NULL);

  if (!self->model_list)
    return;

  for (button = gtk_widget_get_first_child(self->model_list);
       button != NULL;
       button = gtk_widget_get_next_sibling(button))
    {
      GdisModel *model;
      g_autofree gchar *label = NULL;

      model = g_object_get_data(G_OBJECT(button), "gdis-model");
      if (!model)
        continue;

      label = g_strdup_printf("%s [%s | %u atoms]",
                              model->basename,
                              model->format_label,
                              model->atom_count);
      gtk_button_set_label(GTK_BUTTON(button), label);
      gdis_gtk4_window_configure_button_label(button);
      g_object_set_data(G_OBJECT(button), "model-path", model->path);
      gtk_widget_set_tooltip_text(button, model->path);
    }
}

static void
gdis_gtk4_window_refresh_after_model_edit(GdisGtk4Window *self,
                                          gboolean refresh_model_buttons)
{
  g_return_if_fail(self != NULL);

  if (refresh_model_buttons)
    gdis_gtk4_window_refresh_model_buttons(self);
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_refresh_diffraction_tool(self);
}

static void
gdis_gtk4_window_present_report(GdisGtk4Window *self,
                                const char *title,
                                const char *body)
{
  GtkWidget *window;
  GtkWidget *scroller;
  GtkWidget *text_view;
  GtkTextBuffer *buffer;

  g_return_if_fail(self != NULL);
  g_return_if_fail(title != NULL);
  g_return_if_fail(body != NULL);

  window = gtk_window_new();
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_modal(GTK_WINDOW(window), FALSE);
  gtk_window_set_title(GTK_WINDOW(window), title);
  gtk_window_set_default_size(GTK_WINDOW(window), 720, 520);

  scroller = gtk_scrolled_window_new();
  gtk_window_set_child(GTK_WINDOW(window), scroller);

  text_view = gtk_text_view_new();
  gtk_text_view_set_editable(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_WORD_CHAR);
  gtk_text_view_set_monospace(GTK_TEXT_VIEW(text_view), TRUE);
  gtk_text_view_set_left_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_right_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_top_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_bottom_margin(GTK_TEXT_VIEW(text_view), 12);
  buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
  gtk_text_buffer_set_text(buffer, body, -1);
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  gtk_window_present(GTK_WINDOW(window));
}

static gchar *
gdis_gtk4_window_resolve_path(const char *path)
{
  g_autofree gchar *cwd = NULL;
  g_autofree gchar *trimmed_relative = NULL;
  gchar *canonical_path;
  gchar *fallback_path;
  gchar *search_dir;
  const char *suffix;
  const char *relative_path;

  g_return_val_if_fail(path != NULL, NULL);

  if (g_file_test(path, G_FILE_TEST_EXISTS))
    return g_canonicalize_filename(path, NULL);

  cwd = g_get_current_dir();

  relative_path = path;
  while (g_str_has_prefix(relative_path, "./"))
    relative_path += 2;
  while (g_str_has_prefix(relative_path, "../"))
    relative_path += 3;

  if (!g_path_is_absolute(path) && *relative_path != '\0')
    {
      trimmed_relative = g_strdup(relative_path);
      search_dir = g_strdup(cwd);
      while (search_dir)
        {
          g_autofree gchar *candidate = NULL;

          candidate = g_build_filename(search_dir, trimmed_relative, NULL);
          if (g_file_test(candidate, G_FILE_TEST_EXISTS))
            {
              canonical_path = g_canonicalize_filename(candidate, NULL);
              g_free(search_dir);
              return canonical_path;
            }

          if (g_strcmp0(search_dir, G_DIR_SEPARATOR_S) == 0)
            {
              g_free(search_dir);
              search_dir = NULL;
            }
          else
            {
              gchar *parent_dir;

              parent_dir = g_path_get_dirname(search_dir);
              if (g_strcmp0(parent_dir, search_dir) == 0)
                {
                  g_free(parent_dir);
                  g_free(search_dir);
                  search_dir = NULL;
                }
              else
                {
                  g_free(search_dir);
                  search_dir = parent_dir;
                }
            }
        }
    }

  fallback_path = NULL;
  if (!g_path_is_absolute(path))
    {
      fallback_path = g_build_filename("..", "..", path, NULL);
      if (!g_file_test(fallback_path, G_FILE_TEST_EXISTS))
        {
          g_free(fallback_path);
          fallback_path = NULL;
        }
    }
  else
    {
      suffix = NULL;
      if (g_str_has_prefix(path, cwd) && path[strlen(cwd)] == G_DIR_SEPARATOR)
        suffix = path + strlen(cwd) + 1;
      if (suffix)
        {
          fallback_path = g_build_filename("..", "..", suffix, NULL);
          if (!g_file_test(fallback_path, G_FILE_TEST_EXISTS))
            {
              g_free(fallback_path);
              fallback_path = NULL;
            }
        }
    }

  if (!fallback_path)
    return NULL;

  canonical_path = g_canonicalize_filename(fallback_path, NULL);
  g_free(fallback_path);
  return canonical_path;
}

static const GdisAtom *
gdis_gtk4_window_get_selected_atom(const GdisGtk4Window *self)
{
  g_return_val_if_fail(self != NULL, NULL);

  if (!self->active_model)
    return NULL;

  if (self->selected_atom_index == INVALID_ATOM_INDEX)
    return NULL;

  if (self->selected_atom_index >= self->active_model->atoms->len)
    return NULL;

  return g_ptr_array_index(self->active_model->atoms, self->selected_atom_index);
}

static void
gdis_gtk4_window_clear_atom_picks(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);

  if (self->picked_atoms)
    g_array_set_size(self->picked_atoms, 0);
}

static void
gdis_gtk4_window_clear_selected_atoms(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);

  if (self->selected_atoms)
    g_array_set_size(self->selected_atoms, 0);
}

static gboolean
gdis_gtk4_window_atom_array_contains(const GArray *array, guint atom_index)
{
  guint i;

  if (!array)
    return FALSE;

  for (i = 0; i < array->len; i++)
    {
      if (g_array_index(array, guint, i) == atom_index)
        return TRUE;
    }

  return FALSE;
}

static void
gdis_gtk4_window_configure_button_label(GtkWidget *button)
{
  GtkWidget *child;

  g_return_if_fail(button != NULL);

  child = gtk_button_get_child(GTK_BUTTON(button));
  if (!GTK_IS_LABEL(child))
    return;

  gtk_label_set_xalign(GTK_LABEL(child), 0.0f);
  gtk_label_set_ellipsize(GTK_LABEL(child), PANGO_ELLIPSIZE_END);
  gtk_label_set_single_line_mode(GTK_LABEL(child), TRUE);
}

static void
gdis_gtk4_window_add_selected_atom(GdisGtk4Window *self, guint atom_index)
{
  g_return_if_fail(self != NULL);
  g_return_if_fail(self->selected_atoms != NULL);

  if (!gdis_gtk4_window_atom_array_contains(self->selected_atoms, atom_index))
    g_array_append_val(self->selected_atoms, atom_index);
}

static const char *
gdis_gtk4_window_selection_mode_label(GdisSelectionMode mode)
{
  switch (mode)
    {
    case GDIS_SELECTION_MODE_LABEL:
      return "Atom Label";
    case GDIS_SELECTION_MODE_FF_TYPE:
      return "Atom FF Type";
    case GDIS_SELECTION_MODE_ELEMENTS:
      return "Elements";
    case GDIS_SELECTION_MODE_ELEMENTS_IN_MOLECULE:
      return "Elements in Molecule";
    case GDIS_SELECTION_MODE_MOLECULES:
      return "Molecules";
    case GDIS_SELECTION_MODE_MOLECULE_FRAGMENTS:
      return "Molecule Fragments";
    case GDIS_SELECTION_MODE_REGIONS:
      return "Regions";
    case GDIS_SELECTION_MODE_ATOMS:
    default:
      return "Atoms";
    }
}

static void
gdis_gtk4_window_set_click_mode(GdisGtk4Window *self, GdisClickMode click_mode)
{
  g_return_if_fail(self != NULL);

  self->click_mode = click_mode;
}

static const char *
gdis_gtk4_window_click_mode_label(GdisClickMode mode)
{
  switch (mode)
    {
    case GDIS_CLICK_MODE_ADD_BOND:
      return "Pick Add Bond";
    case GDIS_CLICK_MODE_REMOVE_BOND:
      return "Pick Remove Bond";
    case GDIS_CLICK_MODE_SELECT:
    default:
      return "Select";
    }
}

static void
gdis_gtk4_window_collect_connected_atoms(const GdisModel *model,
                                         guint start_index,
                                         GArray *atoms_out)
{
  gboolean *visited;
  GQueue queue = G_QUEUE_INIT;

  g_return_if_fail(model != NULL);
  g_return_if_fail(model->atoms != NULL);
  g_return_if_fail(model->bonds != NULL);
  g_return_if_fail(atoms_out != NULL);

  if (start_index >= model->atoms->len)
    return;

  visited = g_new0(gboolean, model->atoms->len);
  visited[start_index] = TRUE;
  g_queue_push_tail(&queue, GUINT_TO_POINTER(start_index + 1));

  while (!g_queue_is_empty(&queue))
    {
      guint current_index;
      guint i;

      current_index = GPOINTER_TO_UINT(g_queue_pop_head(&queue)) - 1;
      g_array_append_val(atoms_out, current_index);

      for (i = 0; i < model->bonds->len; i++)
        {
          const GdisBond *bond;
          guint neighbor_index;

          bond = &g_array_index(model->bonds, GdisBond, i);
          if (bond->atom_index_a == current_index)
            neighbor_index = bond->atom_index_b;
          else if (bond->atom_index_b == current_index)
            neighbor_index = bond->atom_index_a;
          else
            continue;

          if (neighbor_index >= model->atoms->len || visited[neighbor_index])
            continue;

          visited[neighbor_index] = TRUE;
          g_queue_push_tail(&queue, GUINT_TO_POINTER(neighbor_index + 1));
        }
    }

  g_free(visited);
}

static gboolean
gdis_gtk4_window_collect_fragment_path(const GdisModel *model,
                                       guint start_index,
                                       guint end_index,
                                       GArray *atoms_out)
{
  gboolean *visited;
  gint *previous;
  GQueue queue = G_QUEUE_INIT;
  gboolean found;

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(model->atoms != NULL, FALSE);
  g_return_val_if_fail(model->bonds != NULL, FALSE);
  g_return_val_if_fail(atoms_out != NULL, FALSE);

  if (start_index >= model->atoms->len || end_index >= model->atoms->len)
    return FALSE;

  if (start_index == end_index)
    {
      g_array_append_val(atoms_out, start_index);
      return TRUE;
    }

  visited = g_new0(gboolean, model->atoms->len);
  previous = g_new(gint, model->atoms->len);
  for (guint i = 0; i < model->atoms->len; i++)
    previous[i] = -1;

  found = FALSE;
  visited[start_index] = TRUE;
  previous[start_index] = (gint) start_index;
  g_queue_push_tail(&queue, GUINT_TO_POINTER(start_index + 1));

  while (!g_queue_is_empty(&queue) && !found)
    {
      guint current_index;

      current_index = GPOINTER_TO_UINT(g_queue_pop_head(&queue)) - 1;
      for (guint i = 0; i < model->bonds->len; i++)
        {
          const GdisBond *bond;
          guint neighbor_index;

          bond = &g_array_index(model->bonds, GdisBond, i);
          if (bond->atom_index_a == current_index)
            neighbor_index = bond->atom_index_b;
          else if (bond->atom_index_b == current_index)
            neighbor_index = bond->atom_index_a;
          else
            continue;

          if (neighbor_index >= model->atoms->len || visited[neighbor_index])
            continue;

          visited[neighbor_index] = TRUE;
          previous[neighbor_index] = (gint) current_index;
          if (neighbor_index == end_index)
            {
              found = TRUE;
              break;
            }

          g_queue_push_tail(&queue, GUINT_TO_POINTER(neighbor_index + 1));
        }
    }

  if (found)
    {
      GArray *reverse_path;
      guint cursor;

      reverse_path = g_array_new(FALSE, FALSE, sizeof(guint));
      cursor = end_index;
      while (cursor != start_index)
        {
          g_array_append_val(reverse_path, cursor);
          cursor = (guint) previous[cursor];
        }
      g_array_append_val(reverse_path, start_index);

      for (gint i = (gint) reverse_path->len - 1; i >= 0; i--)
        {
          guint atom_index;

          atom_index = g_array_index(reverse_path, guint, i);
          g_array_append_val(atoms_out, atom_index);
        }

      g_array_free(reverse_path, TRUE);
    }

  g_free(visited);
  g_free(previous);
  return found;
}

static void
gdis_gtk4_window_apply_selection_mode(GdisGtk4Window *self, guint atom_index)
{
  const GdisAtom *clicked_atom;
  guint i;

  g_return_if_fail(self != NULL);

  if (!self->active_model || atom_index >= self->active_model->atoms->len)
    return;

  clicked_atom = g_ptr_array_index(self->active_model->atoms, atom_index);

  if (self->selection_mode != GDIS_SELECTION_MODE_MOLECULE_FRAGMENTS)
    self->fragment_anchor_index = INVALID_ATOM_INDEX;

  gdis_gtk4_window_clear_selected_atoms(self);

  switch (self->selection_mode)
    {
    case GDIS_SELECTION_MODE_LABEL:
      for (i = 0; i < self->active_model->atoms->len; i++)
        {
          const GdisAtom *atom;

          atom = g_ptr_array_index(self->active_model->atoms, i);
          if (g_strcmp0(atom->label, clicked_atom->label) == 0)
            gdis_gtk4_window_add_selected_atom(self, i);
        }
      break;

    case GDIS_SELECTION_MODE_FF_TYPE:
      for (i = 0; i < self->active_model->atoms->len; i++)
        {
          const GdisAtom *atom;

          atom = g_ptr_array_index(self->active_model->atoms, i);
          if (g_strcmp0(atom->ff_type, clicked_atom->ff_type) == 0)
            gdis_gtk4_window_add_selected_atom(self, i);
        }
      break;

    case GDIS_SELECTION_MODE_ELEMENTS:
      for (i = 0; i < self->active_model->atoms->len; i++)
        {
          const GdisAtom *atom;

          atom = g_ptr_array_index(self->active_model->atoms, i);
          if (g_strcmp0(atom->element, clicked_atom->element) == 0)
            gdis_gtk4_window_add_selected_atom(self, i);
        }
      break;

    case GDIS_SELECTION_MODE_ELEMENTS_IN_MOLECULE:
      {
        GArray *component;

        component = g_array_new(FALSE, FALSE, sizeof(guint));
        gdis_gtk4_window_collect_connected_atoms(self->active_model, atom_index, component);
        for (i = 0; i < component->len; i++)
          {
            guint component_index;
            const GdisAtom *atom;

            component_index = g_array_index(component, guint, i);
            atom = g_ptr_array_index(self->active_model->atoms, component_index);
            if (g_strcmp0(atom->element, clicked_atom->element) == 0)
              gdis_gtk4_window_add_selected_atom(self, component_index);
          }
        g_array_free(component, TRUE);
      }
      break;

    case GDIS_SELECTION_MODE_MOLECULES:
      gdis_gtk4_window_collect_connected_atoms(self->active_model,
                                               atom_index,
                                               self->selected_atoms);
      break;

    case GDIS_SELECTION_MODE_MOLECULE_FRAGMENTS:
      if (self->fragment_anchor_index == INVALID_ATOM_INDEX ||
          self->fragment_anchor_index >= self->active_model->atoms->len)
        {
          self->fragment_anchor_index = atom_index;
          gdis_gtk4_window_add_selected_atom(self, atom_index);
          gdis_gtk4_window_log(self,
                               "Fragment anchor set to atom #%u. Pick a second atom in the same molecule to select the bonded path.\n",
                               clicked_atom->serial);
        }
      else
        {
          gboolean found_path;

          found_path = gdis_gtk4_window_collect_fragment_path(self->active_model,
                                                              self->fragment_anchor_index,
                                                              atom_index,
                                                              self->selected_atoms);
          if (!found_path)
            {
              gdis_gtk4_window_add_selected_atom(self, self->fragment_anchor_index);
              gdis_gtk4_window_add_selected_atom(self, atom_index);
              gdis_gtk4_window_log(self,
                                   "Fragment selection could not find a bonded path between the chosen atoms, so only the endpoints were selected.\n");
            }
          self->fragment_anchor_index = INVALID_ATOM_INDEX;
        }
      break;

    case GDIS_SELECTION_MODE_REGIONS:
      gdis_gtk4_window_add_selected_atom(self, atom_index);
      break;

    case GDIS_SELECTION_MODE_ATOMS:
    default:
      gdis_gtk4_window_add_selected_atom(self, atom_index);
      break;
    }

  if (!self->selected_atoms || self->selected_atoms->len == 0)
    gdis_gtk4_window_add_selected_atom(self, atom_index);
}

static void
gdis_gtk4_window_remember_atom_pick(GdisGtk4Window *self, guint atom_index)
{
  guint i;

  g_return_if_fail(self != NULL);
  g_return_if_fail(self->picked_atoms != NULL);

  for (i = 0; i < self->picked_atoms->len; i++)
    {
      guint existing;

      existing = g_array_index(self->picked_atoms, guint, i);
      if (existing == atom_index)
        {
          g_array_remove_index(self->picked_atoms, i);
          break;
        }
    }

  g_array_append_val(self->picked_atoms, atom_index);
  while (self->picked_atoms->len > 4)
    g_array_remove_index(self->picked_atoms, 0);
}

static void
on_measure_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisMeasureTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->measure_tool = NULL;
  g_free(tool);
}

static void
on_edit_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisEditTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->edit_tool = NULL;
  g_free(tool);
}

static gchar *
gdis_gtk4_window_default_diffraction_output_name(const GdisModel *model)
{
  g_autofree gchar *name = NULL;
  gchar *dot;

  if (!model || !model->basename || !model->basename[0])
    return g_strdup("diffraction");

  name = g_strdup(model->basename);
  dot = strrchr(name, '.');
  if (dot)
    *dot = '\0';
  if (!name[0])
    return g_strdup("diffraction");

  return g_steal_pointer(&name);
}

static gboolean
gdis_gtk4_window_model_supports_diffraction(const GdisModel *model,
                                            gchar **message_out)
{
  gdouble a_vec[3];
  gdouble b_vec[3];
  gdouble c_vec[3];

  if (message_out)
    *message_out = NULL;

  if (!model)
    {
      if (message_out)
        *message_out = g_strdup("Powder diffraction needs an active model.");
      return FALSE;
    }

  if (!model->periodic || model->periodicity != 3)
    {
      if (message_out)
        *message_out = g_strdup_printf("Powder diffraction in original GDIS is for 3D periodic crystal models.\n\n"
                                       "The current model `%s` is not a 3D periodic crystal, so there is no valid unit cell to diffract.\n\n"
                                       "Try a periodic file such as `./examples/rocksalt_demo.cif` or another CIF / ARC / CAR model with full cell data.",
                                       model->basename ? model->basename : "current model");
      return FALSE;
    }

  if (!gdis_model_build_cell_vectors(model, a_vec, b_vec, c_vec))
    {
      if (message_out)
        *message_out = g_strdup_printf("The current model `%s` is marked periodic, but its unit cell is incomplete or invalid.\n\n"
                                       "Powder diffraction needs a valid 3D unit cell before it can run.",
                                       model->basename ? model->basename : "current model");
      return FALSE;
    }

  return TRUE;
}

static const char *
gdis_gtk4_window_diffraction_plot_title(GdisDiffractionRadiation radiation)
{
  switch (radiation)
    {
    case GDIS_DIFFRACTION_RADIATION_NEUTRON:
      return "Neutron diffraction";
    case GDIS_DIFFRACTION_RADIATION_ELECTRON:
      return "Electron diffraction";
    case GDIS_DIFFRACTION_RADIATION_XRAY:
    default:
      return "X-Ray diffraction";
    }
}

static void
gdis_diffraction_tool_set_report(GdisDiffractionTool *tool, const char *text)
{
  g_return_if_fail(tool != NULL);
  g_return_if_fail(tool->report_buffer != NULL);

  gtk_text_buffer_set_text(tool->report_buffer, text ? text : "", -1);
}

static gchar *
gdis_diffraction_tool_build_report(GdisGtk4Window *self,
                                   GdisDiffractionTool *tool,
                                   const GdisDiffractionSettings *settings)
{
  GString *report;
  guint peak_count;
  guint i;

  g_return_val_if_fail(self != NULL, NULL);
  g_return_val_if_fail(tool != NULL, NULL);
  g_return_val_if_fail(settings != NULL, NULL);

  report = g_string_new("");
  g_string_append_printf(report,
                         "Powder Diffraction\n"
                         "Model: %s\n"
                         "Radiation: %s\n"
                         "Broadening: %s\n"
                         "Wavelength: %.6f\n"
                         "2Theta range: %.2f to %.2f step %.2f\n"
                         "U/V/W: %.2f  %.2f  %.2f\n"
                         "Mixing parameter: %.2f\n",
                         self->active_model ? self->active_model->basename : "none",
                         gdis_diffraction_radiation_label(settings->radiation),
                         gdis_diffraction_broadening_label(settings->broadening),
                         settings->wavelength,
                         settings->theta_min,
                         tool->pattern ? tool->pattern->settings.theta_max : settings->theta_max,
                         settings->theta_step,
                         settings->u,
                         settings->v,
                         settings->w,
                         settings->asym);

  if (tool->last_export_path)
    g_string_append_printf(report, "Exported spectrum: %s\n", tool->last_export_path);

  peak_count = (tool->pattern && tool->pattern->peaks) ? tool->pattern->peaks->len : 0;
  g_string_append_printf(report, "\nGenerated peaks: %u\n", peak_count);
  if (peak_count > 0)
    g_string_append(report, "Top peaks:\n  2Theta    rel I     d(A)     h   k   l\n");

  for (i = 0; i < peak_count && i < 12; i++)
    {
      const GdisDiffractionPeak *peak;

      peak = &g_array_index(tool->pattern->peaks, GdisDiffractionPeak, i);
      g_string_append_printf(report,
                             " %7.3f  %6.2f  %7.4f   %3d %3d %3d\n",
                             peak->two_theta,
                             peak->relative_intensity,
                             peak->d_spacing,
                             peak->h,
                             peak->k,
                             peak->l);
    }

  return g_string_free(report, FALSE);
}

static gboolean
gdis_diffraction_tool_read_settings(GdisDiffractionTool *tool,
                                    GdisDiffractionSettings *settings)
{
  g_return_val_if_fail(tool != NULL, FALSE);
  g_return_val_if_fail(settings != NULL, FALSE);

  settings->radiation =
    (GdisDiffractionRadiation) gtk_drop_down_get_selected(GTK_DROP_DOWN(tool->radiation_dropdown));
  settings->broadening =
    (GdisDiffractionBroadening) gtk_drop_down_get_selected(GTK_DROP_DOWN(tool->broadening_dropdown));
  settings->wavelength = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->wavelength_spin));
  settings->asym = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->asym_spin));
  settings->theta_min = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->theta_min_spin));
  settings->theta_max = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->theta_max_spin));
  settings->theta_step = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->theta_step_spin));
  settings->u = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->u_spin));
  settings->v = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->v_spin));
  settings->w = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->w_spin));
  return TRUE;
}

static void
on_diffraction_plot_destroy(GtkWindow *window, gpointer user_data)
{
  GdisDiffractionTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  tool->plot_window = NULL;
  tool->plot_area = NULL;
}

static void
on_diffraction_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisDiffractionTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->plot_window && GTK_IS_WINDOW(tool->plot_window))
    gtk_window_destroy(GTK_WINDOW(tool->plot_window));
  if (tool->pattern)
    gdis_diffraction_pattern_free(tool->pattern);
  g_free(tool->plot_title);
  g_free(tool->last_export_path);
  if (tool->owner)
    tool->owner->diffraction_tool = NULL;
  g_free(tool);
}

static void
gdis_diffraction_plot_draw(GtkDrawingArea *area,
                           cairo_t *cr,
                           int width,
                           int height,
                           gpointer user_data)
{
  GdisDiffractionTool *tool;
  const gdouble left = 62.0;
  const gdouble top = 40.0;
  const gdouble right = 26.0;
  const gdouble bottom = 52.0;
  gdouble x_min;
  gdouble x_max;
  gdouble y_max;
  guint i;

  (void) area;

  tool = user_data;

  cairo_set_source_rgb(cr, 0.02, 0.03, 0.05);
  cairo_paint(cr);

  cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 0.92);
  cairo_set_line_width(cr, 1.2);
  cairo_rectangle(cr, left, top, width - left - right, height - top - bottom);
  cairo_stroke(cr);

  if (!tool || !tool->pattern || !tool->pattern->x_values || !tool->pattern->y_values ||
      tool->pattern->x_values->len == 0 || tool->pattern->y_values->len == 0)
    {
      cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
      cairo_set_font_size(cr, 20.0);
      cairo_move_to(cr, 28.0, 34.0);
      cairo_show_text(cr, "Powder Diffraction");
      cairo_set_font_size(cr, 15.0);
      cairo_move_to(cr, 72.0, height * 0.5);
      cairo_show_text(cr, "Execute diffraction from the setup dialog to generate a spectrum.");
      return;
    }

  x_min = tool->pattern->settings.theta_min;
  x_max = tool->pattern->settings.theta_max;
  y_max = tool->pattern->max_intensity;
  if (x_max <= x_min)
    x_max = x_min + 1.0;
  if (y_max <= 1.0e-9)
    y_max = 1.0;

  cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, 18.0);
  cairo_move_to(cr, 28.0, 30.0);
  cairo_show_text(cr, tool->plot_title ? tool->plot_title : "Powder Diffraction");

  {
    g_autofree gchar *wavelength = NULL;

    wavelength = g_strdup_printf("lambda = %.6f", tool->pattern->settings.wavelength);
    cairo_set_font_size(cr, 14.0);
    cairo_move_to(cr, width * 0.52, 30.0);
    cairo_show_text(cr, wavelength);
  }

  cairo_set_font_size(cr, 11.0);
  for (i = 0; i <= 9; i++)
    {
      const gdouble fraction = (gdouble) i / 9.0;
      const gdouble x = left + fraction * (width - left - right);
      g_autofree gchar *label = NULL;

      cairo_move_to(cr, x, height - bottom);
      cairo_line_to(cr, x, height - bottom + 6.0);
      cairo_stroke(cr);

      label = g_strdup_printf("%.2f", x_min + fraction * (x_max - x_min));
      cairo_move_to(cr, x - 12.0, height - bottom + 22.0);
      cairo_show_text(cr, label);
    }

  cairo_save(cr);
  cairo_translate(cr, 20.0, height * 0.58);
  cairo_rotate(cr, -G_PI / 2.0);
  cairo_move_to(cr, 0.0, 0.0);
  cairo_show_text(cr, "intensities (a.u.)");
  cairo_restore(cr);
  cairo_move_to(cr, width * 0.45, height - 12.0);
  cairo_show_text(cr, "2 theta (degree)");

  cairo_set_source_rgba(cr, 0.36, 0.73, 0.97, 0.32);
  cairo_set_line_width(cr, 1.0);
  for (i = 0; i < tool->pattern->peaks->len; i++)
    {
      const GdisDiffractionPeak *peak;
      gdouble peak_x;
      gdouble peak_y;

      peak = &g_array_index(tool->pattern->peaks, GdisDiffractionPeak, i);
      peak_x = left + ((peak->two_theta - x_min) / (x_max - x_min)) * (width - left - right);
      peak_y = height - bottom - (peak->intensity / y_max) * (height - top - bottom);

      cairo_move_to(cr, peak_x, height - bottom);
      cairo_line_to(cr, peak_x, peak_y);
    }
  cairo_stroke(cr);

  cairo_set_source_rgba(cr, 0.33, 0.75, 0.98, 0.95);
  cairo_set_line_width(cr, 2.0);
  for (i = 0; i < tool->pattern->x_values->len; i++)
    {
      const gdouble x_value = g_array_index(tool->pattern->x_values, gdouble, i);
      const gdouble y_value = g_array_index(tool->pattern->y_values, gdouble, i);
      const gdouble plot_x =
        left + ((x_value - x_min) / (x_max - x_min)) * (width - left - right);
      const gdouble plot_y =
        height - bottom - (y_value / y_max) * (height - top - bottom);

      if (i == 0)
        cairo_move_to(cr, plot_x, plot_y);
      else
        cairo_line_to(cr, plot_x, plot_y);
    }
  cairo_stroke(cr);
}

static void
gdis_diffraction_tool_present_plot(GdisDiffractionTool *tool)
{
  GtkWidget *window;
  GtkWidget *area;

  g_return_if_fail(tool != NULL);

  if (tool->plot_window && GTK_IS_WINDOW(tool->plot_window))
    {
      if (tool->plot_title)
        gtk_window_set_title(GTK_WINDOW(tool->plot_window), tool->plot_title);
      if (tool->plot_area)
        gtk_widget_queue_draw(tool->plot_area);
      gtk_window_present(GTK_WINDOW(tool->plot_window));
      return;
    }

  window = gtk_window_new();
  tool->plot_window = window;
  gtk_window_set_application(GTK_WINDOW(window), tool->owner->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(tool->window));
  gtk_window_set_default_size(GTK_WINDOW(window), 980, 620);
  gtk_window_set_title(GTK_WINDOW(window),
                       tool->plot_title ? tool->plot_title : "Powder Diffraction");

  area = gtk_drawing_area_new();
  tool->plot_area = area;
  gtk_widget_set_hexpand(area, TRUE);
  gtk_widget_set_vexpand(area, TRUE);
  gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(area),
                                 gdis_diffraction_plot_draw,
                                 tool,
                                 NULL);
  gtk_window_set_child(GTK_WINDOW(window), area);

  g_signal_connect(window, "destroy", G_CALLBACK(on_diffraction_plot_destroy), tool);
  gtk_window_present(GTK_WINDOW(window));
}

static void
on_diffraction_execute_clicked(GtkButton *button, gpointer user_data)
{
  GdisDiffractionTool *tool;
  GdisGtk4Window *self;
  GdisDiffractionSettings settings;
  GError *error;
  g_autofree gchar *report = NULL;

  (void) button;

  tool = user_data;
  if (!tool || !tool->owner)
    return;

  self = tool->owner;
  {
    g_autofree gchar *availability_message = NULL;

    if (!gdis_gtk4_window_model_supports_diffraction(self->active_model, &availability_message))
      {
        gdis_diffraction_tool_set_report(tool,
                                         availability_message ? availability_message :
                                         "Powder diffraction is unavailable for the current model.");
        gdis_gtk4_window_log(self, "%s\n",
                             availability_message ? availability_message :
                             "Powder diffraction is unavailable for the current model.");
        return;
      }
  }

  if (!self->active_model)
    {
      gdis_diffraction_tool_set_report(tool, "No active model loaded.");
      return;
    }

  gdis_diffraction_tool_read_settings(tool, &settings);

  if (tool->pattern)
    {
      gdis_diffraction_pattern_free(tool->pattern);
      tool->pattern = NULL;
    }
  g_clear_pointer(&tool->plot_title, g_free);
  g_clear_pointer(&tool->last_export_path, g_free);

  error = NULL;
  if (!gdis_diffraction_calculate(self->active_model, &settings, &tool->pattern, &error))
    {
      report = g_strdup_printf("Diffraction calculation failed.\n\n%s",
                               error ? error->message : "unknown error");
      gdis_diffraction_tool_set_report(tool, report);
      gdis_gtk4_window_log(self, "Diffraction failed for %s: %s\n",
                           self->active_model->basename,
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      if (tool->plot_area)
        gtk_widget_queue_draw(tool->plot_area);
      return;
    }

  tool->plot_title = g_strdup(gdis_gtk4_window_diffraction_plot_title(settings.radiation));

  {
    const char *output_name;

    output_name = gtk_editable_get_text(GTK_EDITABLE(tool->output_entry));
    if (output_name && output_name[0] != '\0')
      {
        g_autofree gchar *resolved_path = NULL;

        resolved_path = g_canonicalize_filename(output_name, g_get_current_dir());
        if (!gdis_diffraction_export(tool->pattern, resolved_path, FALSE, &error))
          {
            gdis_gtk4_window_log(self, "Diffraction export failed: %s\n",
                                 error ? error->message : "unknown error");
            g_clear_error(&error);
          }
        else
          {
            tool->last_export_path = g_strdup(resolved_path);
            gdis_gtk4_window_log(self, "Diffraction spectrum exported to %s\n", resolved_path);
          }
      }
  }

  report = gdis_diffraction_tool_build_report(self, tool, &settings);
  gdis_diffraction_tool_set_report(tool, report);
  gdis_gtk4_window_log(self, "Diffraction calculated for %s using %s.\n",
                       self->active_model->basename,
                       gdis_diffraction_radiation_label(settings.radiation));
  gdis_diffraction_tool_present_plot(tool);
}

static void
gdis_gtk4_window_refresh_measure_tool(GdisGtk4Window *self)
{
  g_autofree char *report = NULL;

  g_return_if_fail(self != NULL);

  if (!self->measure_tool || !self->measure_tool->buffer)
    return;

  report = gdis_report_measurements(self->active_model,
                                    self->picked_atoms ? (const guint *) self->picked_atoms->data : NULL,
                                    self->picked_atoms ? self->picked_atoms->len : 0,
                                    self->measure_tool->mode);
  gtk_text_buffer_set_text(self->measure_tool->buffer, report, -1);
}

static void
gdis_gtk4_window_refresh_edit_tool(GdisGtk4Window *self)
{
  GString *report;
  const GdisAtom *selected_atom;
  guint i;

  g_return_if_fail(self != NULL);

  if (!self->edit_tool || !self->edit_tool->buffer)
    return;

  report = g_string_new("");
  if (!self->active_model)
    {
      g_string_append(report, "No active model loaded.");
      gtk_text_buffer_set_text(self->edit_tool->buffer, report->str, -1);
      g_string_free(report, TRUE);
      return;
    }

  selected_atom = gdis_gtk4_window_get_selected_atom(self);

  if (selected_atom)
    {
      gchar number[64];

      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->label_entry),
                            selected_atom->label ? selected_atom->label : "");
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->element_entry),
                            selected_atom->element ? selected_atom->element : "");
      g_snprintf(number, sizeof(number), "%.6f", selected_atom->position[0]);
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->x_entry), number);
      g_snprintf(number, sizeof(number), "%.6f", selected_atom->position[1]);
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->y_entry), number);
      g_snprintf(number, sizeof(number), "%.6f", selected_atom->position[2]);
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->z_entry), number);
    }

  g_string_append_printf(report,
                         "Model Editor\nModel: %s\nFormat: %s\nAtoms: %u   Bonds: %u\nSelection Mode: %s\nClick Mode: %s\nSelected Set: %u   Pick History: %u\n\n",
                         self->active_model->basename,
                         self->active_model->format_label,
                         self->active_model->atom_count,
                         self->active_model->bond_count,
                         gdis_gtk4_window_selection_mode_label(self->selection_mode),
                         gdis_gtk4_window_click_mode_label(self->click_mode),
                         self->selected_atoms ? self->selected_atoms->len : 0,
                         self->picked_atoms ? self->picked_atoms->len : 0);
  g_string_append(report,
                  "Actions:\n"
                  "  Apply To Selected Atom: updates the selected atom from the fields above.\n"
                  "  Add Atom From Fields: creates a new atom using the current fields.\n"
                  "  Delete Selected Atom: removes the current yellow-highlighted atom.\n"
                  "  Delete Selected Set: removes the blue highlighted selection group.\n"
                  "  Delete Picked Atoms: removes the atoms in the current pick history.\n"
                  "  Add Bond Between Last 2 Picks: creates or updates an explicit bond.\n"
                  "  Remove Bond Between Last 2 Picks: removes that bond.\n"
                  "  Pick Add Bond / Pick Remove Bond: use the next 2 viewer picks as a bond edit command.\n"
                  "  Clear Picks: clears the measurement/edit pick set.\n\n");
  g_string_append(report, "Flags: S = selected atom, G = selected group, P = picked atom\n\n");
  g_string_append(report, "Flg  Serial  Elem  FFType     Label                X           Y           Z\n");

  for (i = 0; i < self->active_model->atoms->len && i < 200; i++)
    {
      GdisAtom *atom;
      gboolean selected;
      gboolean group_selected;
      gboolean picked;
      char flags[4] = {' ', ' ', ' ', '\0'};

      atom = g_ptr_array_index(self->active_model->atoms, i);
      selected = (i == self->selected_atom_index);
      group_selected = gdis_gtk4_window_atom_array_contains(self->selected_atoms, i);
      picked = gdis_gtk4_window_atom_array_contains(self->picked_atoms, i);
      if (selected)
        flags[0] = 'S';
      if (group_selected)
        flags[1] = 'P';
      if (picked)
        flags[2] = 'P';
      if (group_selected)
        flags[1] = 'G';

      g_string_append_printf(report,
                             " %3s   %5u   %-3s   %-10s %-16s %10.4f %10.4f %10.4f\n",
                             flags,
                             atom->serial,
                             atom->element,
                             atom->ff_type ? atom->ff_type : "",
                             atom->label,
                             atom->position[0],
                             atom->position[1],
                             atom->position[2]);
    }

  if (self->active_model->atoms->len > 200)
    g_string_append_printf(report,
                           "\n... %u more atoms not shown ...\n",
                           self->active_model->atom_count - 200);

  gtk_text_buffer_set_text(self->edit_tool->buffer, report->str, -1);
  g_string_free(report, TRUE);
}

static void
gdis_gtk4_window_refresh_diffraction_tool(GdisGtk4Window *self)
{
  GdisDiffractionTool *tool;
  g_autofree gchar *availability_message = NULL;
  g_autofree gchar *report = NULL;
  gboolean supported;

  g_return_if_fail(self != NULL);

  tool = self->diffraction_tool;
  if (!tool || !tool->report_buffer)
    return;

  if (tool->pattern)
    {
      gdis_diffraction_pattern_free(tool->pattern);
      tool->pattern = NULL;
    }
  g_clear_pointer(&tool->plot_title, g_free);
  g_clear_pointer(&tool->last_export_path, g_free);

  supported = gdis_gtk4_window_model_supports_diffraction(self->active_model,
                                                          &availability_message);
  if (tool->execute_button)
    gtk_widget_set_sensitive(tool->execute_button, supported);

  if (self->active_model)
    {
      const char *current_text;
      g_autofree gchar *default_name = NULL;

      current_text = gtk_editable_get_text(GTK_EDITABLE(tool->output_entry));
      if (!current_text || current_text[0] == '\0')
        {
          default_name = gdis_gtk4_window_default_diffraction_output_name(self->active_model);
          gtk_editable_set_text(GTK_EDITABLE(tool->output_entry), default_name);
        }
    }

  if (!supported)
    {
      report = g_strdup_printf("Powder Diffraction\n\n%s\n\n"
                               "The Execute button is disabled until a valid 3D periodic crystal model is active.",
                               availability_message ? availability_message : "Powder diffraction is unavailable for the current model.");
    }
  else
    {
      report = g_strdup_printf("Powder Diffraction\n"
                               "Model: %s\n\n"
                               "Choose the radiation, wavelength, broadening, and 2Theta range,\n"
                               "then click Execute to calculate a spectrum and open the plot window.",
                               self->active_model ? self->active_model->basename : "none");
    }
  gdis_diffraction_tool_set_report(tool, report);

  if (tool->plot_area)
    gtk_widget_queue_draw(tool->plot_area);
}

static gboolean
gdis_gtk4_window_parse_entry_double(GtkWidget *entry,
                                    const char *field_name,
                                    gdouble *value_out,
                                    GError **error)
{
  const char *text;
  char *endptr;
  gdouble value;

  g_return_val_if_fail(GTK_IS_EDITABLE(entry), FALSE);
  g_return_val_if_fail(field_name != NULL, FALSE);
  g_return_val_if_fail(value_out != NULL, FALSE);

  text = gtk_editable_get_text(GTK_EDITABLE(entry));
  value = g_ascii_strtod(text, &endptr);
  if (!text || text[0] == '\0' || !endptr || *endptr != '\0')
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "%s must be a valid number.",
                  field_name);
      return FALSE;
    }

  *value_out = value;
  return TRUE;
}

static gboolean
gdis_gtk4_window_parse_entry_uint(GtkWidget *entry,
                                  const char *field_name,
                                  guint *value_out,
                                  GError **error)
{
  const char *text;
  char *endptr;
  guint64 value;

  g_return_val_if_fail(GTK_IS_EDITABLE(entry), FALSE);
  g_return_val_if_fail(field_name != NULL, FALSE);
  g_return_val_if_fail(value_out != NULL, FALSE);

  text = gtk_editable_get_text(GTK_EDITABLE(entry));
  value = g_ascii_strtoull(text, &endptr, 10);
  if (!text || text[0] == '\0' || !endptr || *endptr != '\0' || value == 0 || value > 255)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "%s must be a whole number between 1 and 255.",
                  field_name);
      return FALSE;
    }

  *value_out = (guint) value;
  return TRUE;
}

static gboolean
gdis_gtk4_window_get_last_two_picks(GdisGtk4Window *self,
                                    guint *atom_index_a,
                                    guint *atom_index_b,
                                    GError **error)
{
  g_return_val_if_fail(self != NULL, FALSE);
  g_return_val_if_fail(atom_index_a != NULL, FALSE);
  g_return_val_if_fail(atom_index_b != NULL, FALSE);

  if (!self->picked_atoms || self->picked_atoms->len < 2)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Pick 2 atoms in the viewer first.");
      return FALSE;
    }

  *atom_index_a = g_array_index(self->picked_atoms, guint, self->picked_atoms->len - 2);
  *atom_index_b = g_array_index(self->picked_atoms, guint, self->picked_atoms->len - 1);
  return TRUE;
}

static void
on_measure_mode_clicked(GtkButton *button, gpointer user_data)
{
  GdisMeasureTool *tool;
  const char *mode_name;

  tool = user_data;
  mode_name = g_object_get_data(G_OBJECT(button), "measure-mode");
  if (!tool || !tool->owner || !mode_name)
    return;

  if (g_strcmp0(mode_name, "distance") == 0)
    tool->mode = GDIS_MEASURE_MODE_DISTANCE;
  else if (g_strcmp0(mode_name, "angle") == 0)
    tool->mode = GDIS_MEASURE_MODE_ANGLE;
  else if (g_strcmp0(mode_name, "torsion") == 0)
    tool->mode = GDIS_MEASURE_MODE_TORSION;
  else
    tool->mode = GDIS_MEASURE_MODE_AUTO;

  gdis_gtk4_window_refresh_measure_tool(tool->owner);
}

static void
on_clear_picks_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self)
    return;

  self->selected_atom_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_set_click_mode(self, GDIS_CLICK_MODE_SELECT);
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_log(self, "Pick history cleared.\n");
}

static void
on_delete_selected_atom_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  guint index;

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  if (self->selected_atom_index == INVALID_ATOM_INDEX ||
      self->selected_atom_index >= self->active_model->atoms->len)
    {
      gdis_gtk4_window_log(self, "Delete selected atom requested, but no atom is selected.\n");
      return;
    }

  index = self->selected_atom_index;
  error = NULL;
  if (!gdis_model_delete_atoms(self->active_model, &index, 1, &error))
    {
      gdis_gtk4_window_log(self, "Delete failed: %s\n", error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  self->selected_atom_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_refresh_model_buttons(self);
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_log(self, "Deleted the selected atom.\n");
}

static void
on_delete_picked_atoms_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  if (!self->picked_atoms || self->picked_atoms->len == 0)
    {
      gdis_gtk4_window_log(self, "Delete picked atoms requested, but the pick history is empty.\n");
      return;
    }

  error = NULL;
  if (!gdis_model_delete_atoms(self->active_model,
                               (const guint *) self->picked_atoms->data,
                               self->picked_atoms->len,
                               &error))
    {
      gdis_gtk4_window_log(self, "Delete failed: %s\n", error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  self->selected_atom_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_log(self, "Deleted atoms from the current pick history.\n");
}

static void
on_delete_selected_group_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  if (!self->selected_atoms || self->selected_atoms->len == 0)
    {
      gdis_gtk4_window_log(self, "Delete selected set requested, but nothing is selected.\n");
      return;
    }

  error = NULL;
  if (!gdis_model_delete_atoms(self->active_model,
                               (const guint *) self->selected_atoms->data,
                               self->selected_atoms->len,
                               &error))
    {
      gdis_gtk4_window_log(self, "Delete selected set failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  self->selected_atom_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_log(self, "Deleted the selected atom set.\n");
}

static void
on_edit_apply_selected_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  gdouble x;
  gdouble y;
  gdouble z;

  (void) button;

  self = user_data;
  if (!self || !self->active_model || !self->edit_tool)
    return;

  if (self->selected_atom_index == INVALID_ATOM_INDEX ||
      self->selected_atom_index >= self->active_model->atoms->len)
    {
      gdis_gtk4_window_log(self, "Edit requested, but no atom is selected.\n");
      return;
    }

  error = NULL;
  if (!gdis_gtk4_window_parse_entry_double(self->edit_tool->x_entry, "X", &x, &error) ||
      !gdis_gtk4_window_parse_entry_double(self->edit_tool->y_entry, "Y", &y, &error) ||
      !gdis_gtk4_window_parse_entry_double(self->edit_tool->z_entry, "Z", &z, &error))
    {
      gdis_gtk4_window_log(self, "Edit failed: %s\n", error ? error->message : "invalid coordinates");
      g_clear_error(&error);
      return;
    }

  if (!gdis_model_update_atom(self->active_model,
                              self->selected_atom_index,
                              gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->label_entry)),
                              gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->element_entry)),
                              x,
                              y,
                              z,
                              &error))
    {
      gdis_gtk4_window_log(self, "Edit failed: %s\n", error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_apply_selection_mode(self, self->selected_atom_index);
  gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
  gdis_gtk4_window_log(self, "Updated the selected atom.\n");
}

static void
on_edit_add_atom_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  gdouble x;
  gdouble y;
  gdouble z;

  (void) button;

  self = user_data;
  if (!self || !self->active_model || !self->edit_tool)
    return;

  error = NULL;
  if (!gdis_gtk4_window_parse_entry_double(self->edit_tool->x_entry, "X", &x, &error) ||
      !gdis_gtk4_window_parse_entry_double(self->edit_tool->y_entry, "Y", &y, &error) ||
      !gdis_gtk4_window_parse_entry_double(self->edit_tool->z_entry, "Z", &z, &error))
    {
      gdis_gtk4_window_log(self, "Add atom failed: %s\n", error ? error->message : "invalid coordinates");
      g_clear_error(&error);
      return;
    }

  if (!gdis_model_add_atom(self->active_model,
                           gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->label_entry)),
                           gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->element_entry)),
                           x,
                           y,
                           z,
                           &error))
    {
      gdis_gtk4_window_log(self, "Add atom failed: %s\n", error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  self->selected_atom_index = self->active_model->atoms->len - 1;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_apply_selection_mode(self, self->selected_atom_index);
  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_log(self, "Added a new atom to the model.\n");
}

static void
on_edit_add_bond_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  guint atom_index_a;
  guint atom_index_b;
  guint order;

  (void) button;

  self = user_data;
  if (!self || !self->active_model || !self->edit_tool)
    return;

  error = NULL;
  if (!gdis_gtk4_window_get_last_two_picks(self, &atom_index_a, &atom_index_b, &error) ||
      !gdis_gtk4_window_parse_entry_uint(self->edit_tool->bond_order_entry, "Bond order", &order, &error))
    {
      gdis_gtk4_window_log(self, "Add bond failed: %s\n", error ? error->message : "invalid input");
      g_clear_error(&error);
      return;
    }

  if (!gdis_model_add_explicit_bond(self->active_model,
                                    atom_index_a,
                                    atom_index_b,
                                    (guint8) order,
                                    &error))
    {
      gdis_gtk4_window_log(self, "Add bond failed: %s\n", error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
  gdis_gtk4_window_log(self, "Added or updated an explicit bond between the last 2 picks.\n");
}

static void
on_edit_remove_bond_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  guint atom_index_a;
  guint atom_index_b;

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  error = NULL;
  if (!gdis_gtk4_window_get_last_two_picks(self, &atom_index_a, &atom_index_b, &error))
    {
      gdis_gtk4_window_log(self, "Remove bond failed: %s\n", error ? error->message : "invalid input");
      g_clear_error(&error);
      return;
    }

  if (!gdis_model_remove_bond(self->active_model, atom_index_a, atom_index_b, &error))
    {
      gdis_gtk4_window_log(self, "Remove bond failed: %s\n", error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
  gdis_gtk4_window_log(self, "Removed the bond between the last 2 picks.\n");
}

static void
on_pick_add_bond_mode_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self)
    return;

  gdis_gtk4_window_set_click_mode(self, GDIS_CLICK_MODE_ADD_BOND);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_log(self, "Click mode set to Pick Add Bond. Pick 2 atoms in the viewer.\n");
}

static void
on_pick_remove_bond_mode_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self)
    return;

  gdis_gtk4_window_set_click_mode(self, GDIS_CLICK_MODE_REMOVE_BOND);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_log(self, "Click mode set to Pick Remove Bond. Pick 2 atoms in the viewer.\n");
}

static void
gdis_gtk4_window_present_measure_tool(GdisGtk4Window *self)
{
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *row;
  GtkWidget *button;
  GtkWidget *scroller;
  GtkWidget *text_view;
  GdisMeasureTool *tool;

  g_return_if_fail(self != NULL);

  if (self->measure_tool && GTK_IS_WINDOW(self->measure_tool->window))
    {
      gdis_gtk4_window_refresh_measure_tool(self);
      gtk_window_present(GTK_WINDOW(self->measure_tool->window));
      return;
    }

  tool = g_new0(GdisMeasureTool, 1);
  tool->owner = self;
  tool->mode = GDIS_MEASURE_MODE_AUTO;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Measurements");
  gtk_window_set_default_size(GTK_WINDOW(window), 720, 480);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_box_append(GTK_BOX(root), row);

  button = gtk_button_new_with_label("Auto");
  g_object_set_data(G_OBJECT(button), "measure-mode", "auto");
  g_signal_connect(button, "clicked", G_CALLBACK(on_measure_mode_clicked), tool);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Distance");
  g_object_set_data(G_OBJECT(button), "measure-mode", "distance");
  g_signal_connect(button, "clicked", G_CALLBACK(on_measure_mode_clicked), tool);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Angle");
  g_object_set_data(G_OBJECT(button), "measure-mode", "angle");
  g_signal_connect(button, "clicked", G_CALLBACK(on_measure_mode_clicked), tool);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Torsion");
  g_object_set_data(G_OBJECT(button), "measure-mode", "torsion");
  g_signal_connect(button, "clicked", G_CALLBACK(on_measure_mode_clicked), tool);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Clear Picks");
  g_signal_connect(button, "clicked", G_CALLBACK(on_clear_picks_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(scroller, TRUE);
  gtk_widget_set_vexpand(scroller, TRUE);
  gtk_box_append(GTK_BOX(root), scroller);

  text_view = gtk_text_view_new();
  gtk_text_view_set_editable(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_WORD_CHAR);
  gtk_text_view_set_monospace(GTK_TEXT_VIEW(text_view), TRUE);
  gtk_text_view_set_left_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_right_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_top_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_bottom_margin(GTK_TEXT_VIEW(text_view), 12);
  tool->buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  g_signal_connect(window, "destroy", G_CALLBACK(on_measure_tool_destroy), tool);

  self->measure_tool = tool;
  gdis_gtk4_window_refresh_measure_tool(self);
  gtk_window_present(GTK_WINDOW(window));
}

static void
gdis_gtk4_window_present_edit_tool(GdisGtk4Window *self)
{
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *grid;
  GtkWidget *row;
  GtkWidget *button;
  GtkWidget *label;
  GtkWidget *scroller;
  GtkWidget *text_view;
  GdisEditTool *tool;

  g_return_if_fail(self != NULL);

  if (self->edit_tool && GTK_IS_WINDOW(self->edit_tool->window))
    {
      gdis_gtk4_window_refresh_edit_tool(self);
      gtk_window_present(GTK_WINDOW(self->edit_tool->window));
      return;
    }

  tool = g_new0(GdisEditTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Model Editor");
  gtk_window_set_default_size(GTK_WINDOW(window), 920, 640);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_box_append(GTK_BOX(root), grid);

  label = gtk_label_new("Label");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  tool->label_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->label_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->label_entry, 1, 0, 1, 1);

  label = gtk_label_new("Element");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 2, 0, 1, 1);
  tool->element_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->element_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->element_entry, 3, 0, 1, 1);

  label = gtk_label_new("X");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
  tool->x_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->x_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->x_entry, 1, 1, 1, 1);

  label = gtk_label_new("Y");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 2, 1, 1, 1);
  tool->y_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->y_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->y_entry, 3, 1, 1, 1);

  label = gtk_label_new("Z");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 4, 1, 1, 1);
  tool->z_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->z_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->z_entry, 5, 1, 1, 1);

  label = gtk_label_new("Bond Order");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 4, 0, 1, 1);
  tool->bond_order_entry = gtk_entry_new();
  gtk_editable_set_text(GTK_EDITABLE(tool->bond_order_entry), "1");
  gtk_widget_set_hexpand(tool->bond_order_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->bond_order_entry, 5, 0, 1, 1);

  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_box_append(GTK_BOX(root), row);

  button = gtk_button_new_with_label("Apply To Selected Atom");
  g_signal_connect(button, "clicked", G_CALLBACK(on_edit_apply_selected_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Add Atom From Fields");
  g_signal_connect(button, "clicked", G_CALLBACK(on_edit_add_atom_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Delete Selected Atom");
  g_signal_connect(button, "clicked", G_CALLBACK(on_delete_selected_atom_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Delete Selected Set");
  g_signal_connect(button, "clicked", G_CALLBACK(on_delete_selected_group_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Delete Picked Atoms");
  g_signal_connect(button, "clicked", G_CALLBACK(on_delete_picked_atoms_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_box_append(GTK_BOX(root), row);

  button = gtk_button_new_with_label("Add Bond Between Last 2 Picks");
  g_signal_connect(button, "clicked", G_CALLBACK(on_edit_add_bond_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Remove Bond Between Last 2 Picks");
  g_signal_connect(button, "clicked", G_CALLBACK(on_edit_remove_bond_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Pick Add Bond");
  g_signal_connect(button, "clicked", G_CALLBACK(on_pick_add_bond_mode_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Pick Remove Bond");
  g_signal_connect(button, "clicked", G_CALLBACK(on_pick_remove_bond_mode_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Clear Picks");
  g_signal_connect(button, "clicked", G_CALLBACK(on_clear_picks_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(scroller, TRUE);
  gtk_widget_set_vexpand(scroller, TRUE);
  gtk_box_append(GTK_BOX(root), scroller);

  text_view = gtk_text_view_new();
  gtk_text_view_set_editable(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_NONE);
  gtk_text_view_set_monospace(GTK_TEXT_VIEW(text_view), TRUE);
  gtk_text_view_set_left_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_right_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_top_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_bottom_margin(GTK_TEXT_VIEW(text_view), 12);
  tool->buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  g_signal_connect(window, "destroy", G_CALLBACK(on_edit_tool_destroy), tool);

  self->edit_tool = tool;
  gdis_gtk4_window_refresh_edit_tool(self);
  gtk_window_present(GTK_WINDOW(window));
}

static void
gdis_gtk4_window_present_diffraction_tool(GdisGtk4Window *self)
{
  static const char *const radiation_items[] = {
    "X-Rays",
    "Neutrons",
    "Electrons",
    NULL
  };
  static const char *const broadening_items[] = {
    "Gaussian",
    "Lorentzian",
    "Pseudo-Voigt",
    NULL
  };
  GdisDiffractionTool *tool;
  GdisDiffractionSettings defaults;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *frame;
  GtkWidget *grid;
  GtkWidget *hbox;
  GtkWidget *button_row;
  GtkWidget *label;
  GtkWidget *scroller;
  GtkWidget *text_view;
  GtkStringList *radiation_model;
  GtkStringList *broadening_model;
  g_autofree gchar *default_output_name = NULL;

  g_return_if_fail(self != NULL);

  if (self->diffraction_tool && GTK_IS_WINDOW(self->diffraction_tool->window))
    {
      gdis_gtk4_window_refresh_diffraction_tool(self);
      gtk_window_present(GTK_WINDOW(self->diffraction_tool->window));
      return;
    }

  gdis_diffraction_settings_init(&defaults);

  tool = g_new0(GdisDiffractionTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Powder Diffraction");
  gtk_window_set_default_size(GTK_WINDOW(window), 760, 700);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  frame = gtk_frame_new(NULL);
  gtk_box_append(GTK_BOX(root), frame);
  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 10);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 10);
  gtk_widget_set_margin_start(grid, 10);
  gtk_widget_set_margin_end(grid, 10);
  gtk_widget_set_margin_top(grid, 10);
  gtk_widget_set_margin_bottom(grid, 10);
  gtk_frame_set_child(GTK_FRAME(frame), grid);

  label = gtk_label_new("Radiation");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  radiation_model = gtk_string_list_new(radiation_items);
  tool->radiation_dropdown = gtk_drop_down_new(G_LIST_MODEL(radiation_model), NULL);
  g_object_unref(radiation_model);
  gtk_drop_down_set_selected(GTK_DROP_DOWN(tool->radiation_dropdown), defaults.radiation);
  gtk_widget_set_hexpand(tool->radiation_dropdown, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->radiation_dropdown, 1, 0, 1, 1);

  label = gtk_label_new("Wavelength");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
  tool->wavelength_spin = gtk_spin_button_new_with_range(0.000100, 9.999999, 0.000100);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->wavelength_spin), 6);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->wavelength_spin), defaults.wavelength);
  gtk_widget_set_hexpand(tool->wavelength_spin, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->wavelength_spin, 1, 1, 1, 1);

  label = gtk_label_new("Broadening function");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 2, 1, 1);
  broadening_model = gtk_string_list_new(broadening_items);
  tool->broadening_dropdown = gtk_drop_down_new(G_LIST_MODEL(broadening_model), NULL);
  g_object_unref(broadening_model);
  gtk_drop_down_set_selected(GTK_DROP_DOWN(tool->broadening_dropdown), defaults.broadening);
  gtk_widget_set_hexpand(tool->broadening_dropdown, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->broadening_dropdown, 1, 2, 1, 1);

  label = gtk_label_new("Mixing parameter");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 3, 1, 1);
  tool->asym_spin = gtk_spin_button_new_with_range(0.0, 1.0, 0.01);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->asym_spin), 2);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->asym_spin), defaults.asym);
  gtk_widget_set_hexpand(tool->asym_spin, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->asym_spin, 1, 3, 1, 1);

  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 12);
  gtk_box_append(GTK_BOX(root), hbox);

  frame = gtk_frame_new(NULL);
  gtk_widget_set_hexpand(frame, TRUE);
  gtk_box_append(GTK_BOX(hbox), frame);
  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 10);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 10);
  gtk_widget_set_margin_start(grid, 10);
  gtk_widget_set_margin_end(grid, 10);
  gtk_widget_set_margin_top(grid, 10);
  gtk_widget_set_margin_bottom(grid, 10);
  gtk_frame_set_child(GTK_FRAME(frame), grid);

  label = gtk_label_new("2Theta min");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  tool->theta_min_spin = gtk_spin_button_new_with_range(0.0, 170.0, 0.1);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->theta_min_spin), 2);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->theta_min_spin), defaults.theta_min);
  gtk_grid_attach(GTK_GRID(grid), tool->theta_min_spin, 1, 0, 1, 1);

  label = gtk_label_new("2Theta max");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
  tool->theta_max_spin = gtk_spin_button_new_with_range(0.0, 170.0, 0.1);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->theta_max_spin), 2);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->theta_max_spin), defaults.theta_max);
  gtk_grid_attach(GTK_GRID(grid), tool->theta_max_spin, 1, 1, 1, 1);

  label = gtk_label_new("2Theta step");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 2, 1, 1);
  tool->theta_step_spin = gtk_spin_button_new_with_range(0.01, 1.0, 0.01);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->theta_step_spin), 2);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->theta_step_spin), defaults.theta_step);
  gtk_grid_attach(GTK_GRID(grid), tool->theta_step_spin, 1, 2, 1, 1);

  frame = gtk_frame_new(NULL);
  gtk_widget_set_hexpand(frame, TRUE);
  gtk_box_append(GTK_BOX(hbox), frame);
  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 10);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 10);
  gtk_widget_set_margin_start(grid, 10);
  gtk_widget_set_margin_end(grid, 10);
  gtk_widget_set_margin_top(grid, 10);
  gtk_widget_set_margin_bottom(grid, 10);
  gtk_frame_set_child(GTK_FRAME(frame), grid);

  label = gtk_label_new("U");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  tool->u_spin = gtk_spin_button_new_with_range(-9.9, 9.9, 0.05);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->u_spin), 2);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->u_spin), defaults.u);
  gtk_grid_attach(GTK_GRID(grid), tool->u_spin, 1, 0, 1, 1);

  label = gtk_label_new("V");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
  tool->v_spin = gtk_spin_button_new_with_range(-9.9, 9.9, 0.05);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->v_spin), 2);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->v_spin), defaults.v);
  gtk_grid_attach(GTK_GRID(grid), tool->v_spin, 1, 1, 1, 1);

  label = gtk_label_new("W");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 2, 1, 1);
  tool->w_spin = gtk_spin_button_new_with_range(-9.9, 9.9, 0.05);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->w_spin), 2);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->w_spin), defaults.w);
  gtk_grid_attach(GTK_GRID(grid), tool->w_spin, 1, 2, 1, 1);

  frame = gtk_frame_new(NULL);
  gtk_box_append(GTK_BOX(root), frame);
  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 10);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 10);
  gtk_widget_set_margin_start(grid, 10);
  gtk_widget_set_margin_end(grid, 10);
  gtk_widget_set_margin_top(grid, 10);
  gtk_widget_set_margin_bottom(grid, 10);
  gtk_frame_set_child(GTK_FRAME(frame), grid);

  label = gtk_label_new("Output filename");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  tool->output_entry = gtk_entry_new();
  default_output_name = gdis_gtk4_window_default_diffraction_output_name(self->active_model);
  gtk_editable_set_text(GTK_EDITABLE(tool->output_entry), default_output_name);
  gtk_widget_set_hexpand(tool->output_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->output_entry, 1, 0, 1, 1);

  scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(scroller, TRUE);
  gtk_widget_set_vexpand(scroller, TRUE);
  gtk_box_append(GTK_BOX(root), scroller);

  text_view = gtk_text_view_new();
  gtk_text_view_set_editable(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_WORD_CHAR);
  gtk_text_view_set_monospace(GTK_TEXT_VIEW(text_view), TRUE);
  gtk_text_view_set_left_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_right_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_top_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_bottom_margin(GTK_TEXT_VIEW(text_view), 12);
  tool->report_buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  button_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_widget_set_halign(button_row, GTK_ALIGN_END);
  gtk_box_append(GTK_BOX(root), button_row);

  label = gtk_button_new_with_label("Execute");
  tool->execute_button = label;
  g_signal_connect(label, "clicked", G_CALLBACK(on_diffraction_execute_clicked), tool);
  gtk_box_append(GTK_BOX(button_row), label);

  label = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(label, "clicked", G_CALLBACK(gtk_window_destroy), window);
  gtk_box_append(GTK_BOX(button_row), label);

  g_signal_connect(window, "destroy", G_CALLBACK(on_diffraction_tool_destroy), tool);

  self->diffraction_tool = tool;
  gdis_gtk4_window_refresh_diffraction_tool(self);
  gtk_window_present(GTK_WINDOW(window));
}

static void
gdis_gtk4_window_update_details(GdisGtk4Window *self)
{
  GString *summary;
  GString *content;
  GString *editing;
  GString *images;
  GString *symmetry;
  const GdisAtom *selected_atom;
  gint image_negative[3];
  gint image_positive[3];
  gint image_limits[6];
  gdouble a_vec[3];
  gdouble b_vec[3];
  gdouble c_vec[3];
  gdouble cell_area;
  gdouble cell_volume;
  gboolean have_cell_vectors;
  guint image_dims;
  guint displayed_cells;
  guint i;

  g_return_if_fail(self != NULL);

  gdis_gtk4_window_sync_image_controls(self);

  if (!self->active_model)
    {
      if (self->active_summary_buffer)
        gdis_text_buffer_set(self->active_summary_buffer,
                             "No active model\nOpen a file or choose a sample model.");
      if (self->content_buffer)
        gdis_text_buffer_set(self->content_buffer, "No active model loaded.");
      if (self->editing_buffer)
        gdis_text_buffer_set(self->editing_buffer, "No active model loaded.");
      if (self->images_buffer)
        gdis_text_buffer_set(self->images_buffer, "No active model loaded.");
      if (self->symmetry_buffer)
        gdis_text_buffer_set(self->symmetry_buffer, "No active model loaded.");
      return;
    }

  selected_atom = gdis_gtk4_window_get_selected_atom(self);
  gdis_model_get_image_limits(self->active_model, image_limits);
  image_dims = gdis_gtk4_window_get_image_axis_limits(self->active_model,
                                                      image_negative,
                                                      image_positive);
  have_cell_vectors = gdis_model_build_cell_vectors(self->active_model, a_vec, b_vec, c_vec);
  cell_area = 0.0;
  cell_volume = 0.0;
  if (have_cell_vectors)
    {
      gdouble cross_ab[3];

      cross_ab[0] = a_vec[1] * b_vec[2] - a_vec[2] * b_vec[1];
      cross_ab[1] = a_vec[2] * b_vec[0] - a_vec[0] * b_vec[2];
      cross_ab[2] = a_vec[0] * b_vec[1] - a_vec[1] * b_vec[0];
      cell_area = sqrt(cross_ab[0] * cross_ab[0] +
                       cross_ab[1] * cross_ab[1] +
                       cross_ab[2] * cross_ab[2]);
      cell_volume = fabs(cross_ab[0] * c_vec[0] +
                         cross_ab[1] * c_vec[1] +
                         cross_ab[2] * c_vec[2]);
    }
  displayed_cells = 1;
  for (i = 0; i < 3; i++)
    displayed_cells *= (guint) (image_negative[i] + image_positive[i]);

  summary = g_string_new("");
  g_string_append_printf(summary,
                         "%s\nFormat: %s\nAtoms: %u   Bonds: %u\nZoom: %.2fx   Picks: %u   Selected: %u\nSelection Mode: %s\nImage Range: %u displayed cell%s\n%s",
                         self->active_model->basename,
                         self->active_model->format_label,
                         self->active_model->atom_count,
                         self->active_model->bond_count,
                         self->zoom,
                         self->picked_atoms ? self->picked_atoms->len : 0,
                         self->selected_atoms ? self->selected_atoms->len : 0,
                         gdis_gtk4_window_selection_mode_label(self->selection_mode),
                         displayed_cells,
                         displayed_cells == 1 ? "" : "s",
                         self->active_model->path);
  if (selected_atom)
    {
      g_string_append_printf(summary,
                             "\n\nSelected: %s [%s] #%u\n(%.4f, %.4f, %.4f)",
                             selected_atom->label,
                             selected_atom->element,
                             selected_atom->serial,
                             selected_atom->position[0],
                             selected_atom->position[1],
                             selected_atom->position[2]);
    }
  else
    {
      g_string_append(summary, "\n\nSelected: none");
    }
  if (self->active_summary_buffer)
    gdis_text_buffer_set(self->active_summary_buffer, "%s", summary->str);
  g_string_free(summary, TRUE);

  content = g_string_new("");
  g_string_append_printf(content, "Title: %s\n", self->active_model->title);
  g_string_append_printf(content, "Path: %s\n", self->active_model->path);
  g_string_append_printf(content, "Format: %s\n", self->active_model->format_label);
  g_string_append_printf(content, "Atoms: %u\n", self->active_model->atom_count);
  g_string_append_printf(content, "Selection mode: %s\n",
                         gdis_gtk4_window_selection_mode_label(self->selection_mode));
  g_string_append_printf(content, "Selected set size: %u\n",
                         self->selected_atoms ? self->selected_atoms->len : 0);
  g_string_append_printf(content, "Bonds: %u", self->active_model->bond_count);
  if (self->active_model->explicit_bond_count > 0)
    g_string_append_printf(content, " (%u explicit)\n", self->active_model->explicit_bond_count);
  else
    g_string_append(content, " (inferred)\n");
  g_string_append_printf(content, "Periodic: %s", self->active_model->periodic ? "yes" : "no");
  if (self->active_model->periodic)
    g_string_append_printf(content, " (%uD)", self->active_model->periodicity);
  g_string_append_c(content, '\n');
  if (selected_atom)
    {
      g_string_append_printf(content,
                             "\nSelected atom:\n  Label: %s\n  Element: %s\n  FF Type: %s\n  Serial: %u\n  Position: %.4f %.4f %.4f\n",
                             selected_atom->label,
                             selected_atom->element,
                             selected_atom->ff_type ? selected_atom->ff_type : "",
                             selected_atom->serial,
                             selected_atom->position[0],
                             selected_atom->position[1],
                             selected_atom->position[2]);
    }
  if (self->content_buffer)
    gdis_text_buffer_set(self->content_buffer, "%s", content->str);
  g_string_free(content, TRUE);

  editing = g_string_new(" Sel  Serial  Elem  FFType     Label                X           Y           Z\n");
  for (i = 0; i < self->active_model->atoms->len && i < 32; i++)
    {
      GdisAtom *atom;
      const char *marker;

      atom = g_ptr_array_index(self->active_model->atoms, i);
      if (i == self->selected_atom_index)
        marker = "S";
      else if (gdis_gtk4_window_atom_array_contains(self->selected_atoms, i))
        marker = "G";
      else
        marker = " ";
      g_string_append_printf(editing,
                             "  %s   %5u   %-3s   %-10s %-16s %10.4f %10.4f %10.4f\n",
                             marker,
                             atom->serial,
                             atom->element,
                             atom->ff_type ? atom->ff_type : "",
                             atom->label,
                             atom->position[0],
                             atom->position[1],
                             atom->position[2]);
    }
  if (self->active_model->atoms->len > 32)
    g_string_append_printf(editing,
                           "\n... %u more atoms not shown ...\n",
                           self->active_model->atom_count - 32);
  if (selected_atom && self->selected_atom_index >= 32)
    {
      g_string_append_printf(editing,
                             "\nSelected atom outside preview table:\n  #%u  %s [%s]  %.4f %.4f %.4f\n",
                             selected_atom->serial,
                             selected_atom->label,
                             selected_atom->element,
                             selected_atom->position[0],
                             selected_atom->position[1],
                             selected_atom->position[2]);
    }
  if (self->editing_buffer)
    gdis_text_buffer_set(self->editing_buffer, "%s", editing->str);
  g_string_free(editing, TRUE);

  images = g_string_new("");
  g_string_append_printf(images, "Periodic: %s\n", self->active_model->periodic ? "yes" : "no");
  g_string_append_printf(images, "Periodic dimensions: %u\n", image_dims);
  g_string_append_printf(images,
                         "Cell lengths: a=%0.4f  b=%0.4f  c=%0.4f\n",
                         self->active_model->cell_lengths[0],
                         self->active_model->cell_lengths[1],
                         self->active_model->cell_lengths[2]);
  g_string_append_printf(images,
                         "Cell angles:  a=%0.2f  b=%0.2f  g=%0.2f\n\n",
                         self->active_model->cell_angles[0],
                         self->active_model->cell_angles[1],
                         self->active_model->cell_angles[2]);
  g_string_append_printf(images,
                         "Viewer toggles:\n  Atoms: %s\n  Bonds: %s\n  Cell: %s\n  Labels: %s\n",
                         self->show_atoms_toggle &&
                         gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_atoms_toggle)) ? "on" : "off",
                         self->show_bonds_toggle &&
                         gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_bonds_toggle)) ? "on" : "off",
                         self->show_cell_toggle &&
                         gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_cell_toggle)) ? "on" : "off",
                         self->show_labels_toggle &&
                         gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_labels_toggle)) ? "on" : "off");
  g_string_append_printf(images,
                         "\nStored image range:\n"
                         "  a: %d negative, %d positive image%s (%d total cells on a)\n"
                         "  b: %d negative, %d positive image%s (%d total cells on b)\n"
                         "  c: %d negative, %d positive image%s (%d total cells on c)\n"
                         "Displayed cell span: %u cell%s\n"
                         "Crystal tools: apply image range, confine to cell, force P1, or bake a supercell.\n",
                         image_limits[0],
                         MAX(image_limits[1] - 1, 0),
                         MAX(image_limits[1] - 1, 0) == 1 ? "" : "s",
                         image_negative[0] + image_positive[0],
                         image_limits[2],
                         MAX(image_limits[3] - 1, 0),
                         MAX(image_limits[3] - 1, 0) == 1 ? "" : "s",
                         image_negative[1] + image_positive[1],
                         image_limits[4],
                         MAX(image_limits[5] - 1, 0),
                         MAX(image_limits[5] - 1, 0) == 1 ? "" : "s",
                         image_negative[2] + image_positive[2],
                         displayed_cells,
                         displayed_cells == 1 ? "" : "s");
  if (self->images_buffer)
    gdis_text_buffer_set(self->images_buffer, "%s", images->str);
  g_string_free(images, TRUE);

  symmetry = g_string_new("");
  g_string_append_printf(symmetry, "Space group: %s\n",
                         self->active_model->space_group ? self->active_model->space_group : "unknown");
  g_string_append_printf(symmetry, "Periodic dimensions: %u\n", self->active_model->periodicity);
  g_string_append_printf(symmetry,
                         "Cell: %0.4f %0.4f %0.4f / %0.2f %0.2f %0.2f\n",
                         self->active_model->cell_lengths[0],
                         self->active_model->cell_lengths[1],
                         self->active_model->cell_lengths[2],
                         self->active_model->cell_angles[0],
                         self->active_model->cell_angles[1],
                         self->active_model->cell_angles[2]);
  if (have_cell_vectors)
    {
      if (self->active_model->periodicity == 2)
        g_string_append_printf(symmetry, "Cell area: %.4f A^2\n", cell_area);
      else if (self->active_model->periodicity >= 3)
        g_string_append_printf(symmetry, "Cell volume: %.4f A^3\n", cell_volume);
    }
  g_string_append_printf(symmetry,
                         "Displayed image span: %u cell%s\n\n",
                         displayed_cells,
                         displayed_cells == 1 ? "" : "s");
  g_string_append(symmetry,
                  "Loaded symmetry metadata:\n"
                  "  - space-group label\n"
                  "  - periodic dimensionality\n"
                  "  - unit-cell geometry\n\n"
                  "Still deferred in this GTK4 bridge:\n"
                  "  - symmetry-operator expansion\n"
                  "  - asymmetric/full-cell toggling\n"
                  "  - full legacy symmetry editor");
  if (self->symmetry_buffer)
    gdis_text_buffer_set(self->symmetry_buffer, "%s", symmetry->str);
  g_string_free(symmetry, TRUE);
}

static GdisModel *
gdis_gtk4_window_find_model(GdisGtk4Window *self, const char *path)
{
  guint i;

  g_return_val_if_fail(self != NULL, NULL);
  g_return_val_if_fail(path != NULL, NULL);

  for (i = 0; i < self->models->len; i++)
    {
      GdisModel *model;

      model = g_ptr_array_index(self->models, i);
      if (g_strcmp0(model->path, path) == 0)
        return model;
    }

  return NULL;
}

static GtkWidget *
gdis_gtk4_window_find_model_button(GdisGtk4Window *self, const char *path)
{
  GtkWidget *button;

  g_return_val_if_fail(self != NULL, NULL);
  g_return_val_if_fail(path != NULL, NULL);

  for (button = gtk_widget_get_first_child(self->model_list);
       button != NULL;
       button = gtk_widget_get_next_sibling(button))
    {
      const char *existing_path;

      existing_path = g_object_get_data(G_OBJECT(button), "model-path");
      if (g_strcmp0(existing_path, path) == 0)
        return button;
    }

  return NULL;
}

static void
gdis_gtk4_window_reset_view(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);

  self->rotation_x = -0.45;
  self->rotation_y = 0.55;
  self->zoom = 1.0;
}

static void
gdis_gtk4_window_apply_axis_view(GdisGtk4Window *self,
                                 gdouble rotation_x,
                                 gdouble rotation_y)
{
  g_return_if_fail(self != NULL);

  self->rotation_x = rotation_x;
  self->rotation_y = rotation_y;
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
}

static void
gdis_gtk4_window_set_active_model(GdisGtk4Window *self, GdisModel *model)
{
  GtkWidget *button;

  g_return_if_fail(self != NULL);

  self->active_model = model;
  self->selected_atom_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_set_click_mode(self, GDIS_CLICK_MODE_SELECT);
  gdis_gtk4_window_reset_view(self);
  gdis_gtk4_window_clear_status_log(self);

  for (button = gtk_widget_get_first_child(self->model_list);
       button != NULL;
       button = gtk_widget_get_next_sibling(button))
    gtk_widget_remove_css_class(button, "gdis-model-button-active");

  if (model)
    {
      gtk_window_set_title(GTK_WINDOW(self->window), model->basename);
      button = gdis_gtk4_window_find_model_button(self, model->path);
      if (button)
        gtk_widget_add_css_class(button, "gdis-model-button-active");
      gdis_gtk4_window_log(self,
                           "Active model: %s\n"
                           "Path: %s\n"
                           "Atoms: %u   Bonds: %u\n\n",
                           model->basename,
                           model->path,
                           model->atom_count,
                           model->bond_count);
    }
  else
    {
      gtk_window_set_title(GTK_WINDOW(self->window), "GDIS GTK4");
      gdis_gtk4_window_log(self, "Active model cleared.\n");
    }

  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_refresh_diffraction_tool(self);
}

static gboolean
gdis_gtk4_window_save_model_to_path(GdisGtk4Window *self,
                                    GdisModel *model,
                                    const char *path)
{
  GError *error;

  g_return_val_if_fail(self != NULL, FALSE);
  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(path != NULL, FALSE);

  error = NULL;
  if (!gdis_model_save(model, path, &error))
    {
      gdis_gtk4_window_log(self, "Save failed: %s\n", error ? error->message : "unknown error");
      g_clear_error(&error);
      return FALSE;
    }

  gdis_gtk4_window_refresh_model_buttons(self);
  if (self->active_model == model)
    gtk_window_set_title(GTK_WINDOW(self->window), model->basename);
  gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
  gdis_gtk4_window_log(self, "Saved model to: %s\n", model->path);
  return TRUE;
}

static gboolean
gdis_gtk4_window_add_new_empty_model(GdisGtk4Window *self)
{
  g_autofree gchar *path = NULL;
  GdisModel *model;

  g_return_val_if_fail(self != NULL, FALSE);

  self->untitled_counter++;
  path = g_strdup_printf("Untitled-%u.xyz", self->untitled_counter);
  model = gdis_model_create_empty(path, GDIS_MODEL_FORMAT_XYZ);
  if (!model)
    return FALSE;

  gdis_gtk4_window_add_loaded_model(self, model, TRUE);
  gdis_gtk4_window_log(self, "Created new model: %s\n", model->basename);
  return TRUE;
}

static gboolean
gdis_gtk4_window_add_model_from_path(GdisGtk4Window *self,
                                     const char *path,
                                     gboolean select_after_add)
{
  g_autofree gchar *resolved_path = NULL;
  GdisModel *model;

  g_return_val_if_fail(self != NULL, FALSE);
  g_return_val_if_fail(path != NULL, FALSE);

  resolved_path = gdis_gtk4_window_resolve_path(path);
  if (!resolved_path)
    {
      gdis_gtk4_window_log(self, "Path not found: %s\n", path);
      return FALSE;
    }

  if (g_file_test(resolved_path, G_FILE_TEST_IS_DIR))
    {
      gdis_gtk4_window_log(self, "Skipping directory: %s\n", resolved_path);
      return FALSE;
    }

  model = gdis_gtk4_window_find_model(self, resolved_path);
  if (!model)
    {
      GError *error;

      error = NULL;
      model = gdis_model_load(resolved_path, &error);
      if (!model)
        {
          gdis_gtk4_window_log(self, "Load failed for %s: %s\n",
                               resolved_path,
                               error ? error->message : "unknown error");
          g_clear_error(&error);
          return FALSE;
        }

      gdis_gtk4_window_add_loaded_model(self, model, select_after_add);
      gdis_gtk4_window_log(self, "Loaded model: %s (%s, %u atoms)\n",
                           model->path,
                           model->format_label,
                           model->atom_count);
    }
  else if (select_after_add)
    {
      gdis_gtk4_window_set_active_model(self, model);
    }

  return TRUE;
}

static guint
gdis_gtk4_window_get_image_axis_limits(const GdisModel *model,
                                       gint negative[3],
                                       gint positive[3])
{
  guint dims;
  guint axis;

  dims = 0;
  if (model && model->periodic)
    {
      dims = model->periodicity;
      if (dims == 0)
        dims = 3;
    }
  dims = MIN(dims, 3u);

  for (axis = 0; axis < 3; axis++)
    {
      if (negative)
        negative[axis] = 0;
      if (positive)
        positive[axis] = 1;
      if (model && axis < dims)
        {
          if (negative)
            negative[axis] = MAX(model->image_limits[2 * axis], 0);
          if (positive)
            positive[axis] = MAX(model->image_limits[2 * axis + 1], 1);
        }
    }

  return dims;
}

static void
gdis_gtk4_window_sync_image_controls(GdisGtk4Window *self)
{
  gint limits[6] = {0, 1, 0, 1, 0, 1};
  guint dims;

  g_return_if_fail(self != NULL);

  dims = gdis_gtk4_window_get_image_axis_limits(self->active_model, NULL, NULL);
  if (self->active_model)
    gdis_model_get_image_limits(self->active_model, limits);

  for (guint axis = 0; axis < 3; axis++)
    {
      gboolean axis_enabled;

      axis_enabled = (self->active_model != NULL && axis < dims);

      if (self->image_limit_spin[2 * axis])
        {
          gtk_spin_button_set_value(GTK_SPIN_BUTTON(self->image_limit_spin[2 * axis]),
                                    limits[2 * axis]);
          gtk_widget_set_sensitive(self->image_limit_spin[2 * axis], axis_enabled);
        }

      if (self->image_limit_spin[2 * axis + 1])
        {
          gtk_spin_button_set_value(GTK_SPIN_BUTTON(self->image_limit_spin[2 * axis + 1]),
                                    limits[2 * axis + 1]);
          gtk_widget_set_sensitive(self->image_limit_spin[2 * axis + 1], axis_enabled);
        }

      if (self->supercell_repeat_spin[axis])
        {
          if (!axis_enabled)
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(self->supercell_repeat_spin[axis]), 1.0);
          gtk_widget_set_sensitive(self->supercell_repeat_spin[axis], axis_enabled);
        }
    }
}

static GtkWidget *
gdis_gtk4_window_add_loaded_model(GdisGtk4Window *self,
                                  GdisModel *model,
                                  gboolean select_after_add)
{
  GtkWidget *button;
  gchar *label;

  g_return_val_if_fail(self != NULL, NULL);
  g_return_val_if_fail(model != NULL, NULL);

  g_ptr_array_add(self->models, model);
  label = g_strdup_printf("%s [%s | %u atoms]",
                          model->basename,
                          model->format_label,
                          model->atom_count);
  button = gtk_button_new_with_label(label);
  gtk_widget_add_css_class(button, "flat");
  gtk_widget_add_css_class(button, "gdis-model-button");
  gtk_widget_set_hexpand(button, TRUE);
  gtk_widget_set_halign(button, GTK_ALIGN_FILL);
  gdis_gtk4_window_configure_button_label(button);
  g_object_set_data(G_OBJECT(button), "model-path", model->path);
  g_object_set_data(G_OBJECT(button), "gdis-model", model);
  gtk_widget_set_tooltip_text(button, model->path);
  g_signal_connect(button, "clicked", G_CALLBACK(on_model_button_clicked), self);
  gtk_box_append(GTK_BOX(self->model_list), button);
  g_free(label);

  if (select_after_add)
    gdis_gtk4_window_set_active_model(self, model);

  return button;
}

static void
gdis_gtk4_window_remove_model(GdisGtk4Window *self, GdisModel *model)
{
  GtkWidget *button;
  GdisModel *replacement;
  guint index;

  g_return_if_fail(self != NULL);
  g_return_if_fail(model != NULL);

  replacement = NULL;
  index = G_MAXUINT;
  for (guint i = 0; i < self->models->len; i++)
    {
      if (g_ptr_array_index(self->models, i) == model)
        {
          index = i;
          break;
        }
    }

  if (index == G_MAXUINT)
    return;

  if (self->models->len > 1)
    {
      if (index + 1 < self->models->len)
        replacement = g_ptr_array_index(self->models, index + 1);
      else if (index > 0)
        replacement = g_ptr_array_index(self->models, index - 1);
    }

  button = gdis_gtk4_window_find_model_button(self, model->path);
  if (button)
    gtk_box_remove(GTK_BOX(self->model_list), button);

  if (self->active_model == model)
    self->active_model = NULL;

  g_ptr_array_remove_index(self->models, index);
  gdis_gtk4_window_set_active_model(self, replacement);
  gdis_gtk4_window_refresh_model_buttons(self);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_refresh_measure_tool(self);
}

static void
gdis_rotate_point(const gdouble input[3],
                  gdouble rotation_x,
                  gdouble rotation_y,
                  gdouble output[3])
{
  gdouble sin_x;
  gdouble cos_x;
  gdouble sin_y;
  gdouble cos_y;
  gdouble tmp[3];

  sin_x = sin(rotation_x);
  cos_x = cos(rotation_x);
  sin_y = sin(rotation_y);
  cos_y = cos(rotation_y);

  tmp[0] = input[0];
  tmp[1] = input[1] * cos_x - input[2] * sin_x;
  tmp[2] = input[1] * sin_x + input[2] * cos_x;

  output[0] = tmp[0] * cos_y + tmp[2] * sin_y;
  output[1] = tmp[1];
  output[2] = -tmp[0] * sin_y + tmp[2] * cos_y;
}

static gboolean
gdis_model_build_cell_vectors(const GdisModel *model,
                              gdouble a_vec[3],
                              gdouble b_vec[3],
                              gdouble c_vec[3])
{
  gdouble alpha;
  gdouble beta;
  gdouble gamma;
  gdouble sin_gamma;
  gdouble cx;
  gdouble cy;
  gdouble cz2;

  g_return_val_if_fail(model != NULL, FALSE);

  if (!model->periodic)
    return FALSE;

  if (model->cell_lengths[0] <= 1.0e-6 ||
      model->cell_lengths[1] <= 1.0e-6 ||
      model->cell_lengths[2] <= 1.0e-6)
    return FALSE;

  alpha = model->cell_angles[0] * (G_PI / 180.0);
  beta = model->cell_angles[1] * (G_PI / 180.0);
  gamma = model->cell_angles[2] * (G_PI / 180.0);
  sin_gamma = sin(gamma);
  if (fabs(sin_gamma) < 1.0e-6)
    return FALSE;

  a_vec[0] = model->cell_lengths[0];
  a_vec[1] = 0.0;
  a_vec[2] = 0.0;

  b_vec[0] = model->cell_lengths[1] * cos(gamma);
  b_vec[1] = model->cell_lengths[1] * sin_gamma;
  b_vec[2] = 0.0;

  cx = model->cell_lengths[2] * cos(beta);
  cy = model->cell_lengths[2] * (cos(alpha) - cos(beta) * cos(gamma)) / sin_gamma;
  cz2 = model->cell_lengths[2] * model->cell_lengths[2] - cx * cx - cy * cy;
  if (cz2 < 0.0)
    cz2 = 0.0;

  c_vec[0] = cx;
  c_vec[1] = cy;
  c_vec[2] = sqrt(cz2);

  return TRUE;
}

static gboolean
gdis_model_build_cell_vertices(const GdisModel *model, GdisVec3 vertices[8])
{
  gdouble a_vec[3];
  gdouble b_vec[3];
  gdouble c_vec[3];

  g_return_val_if_fail(model != NULL, FALSE);
  g_return_val_if_fail(vertices != NULL, FALSE);

  if (!gdis_model_build_cell_vectors(model, a_vec, b_vec, c_vec))
    return FALSE;

  vertices[0] = (GdisVec3) {0.0, 0.0, 0.0};
  vertices[1] = (GdisVec3) {a_vec[0], a_vec[1], a_vec[2]};
  vertices[2] = (GdisVec3) {b_vec[0], b_vec[1], b_vec[2]};
  vertices[3] = (GdisVec3) {a_vec[0] + b_vec[0], a_vec[1] + b_vec[1], a_vec[2] + b_vec[2]};
  vertices[4] = (GdisVec3) {c_vec[0], c_vec[1], c_vec[2]};
  vertices[5] = (GdisVec3) {a_vec[0] + c_vec[0], a_vec[1] + c_vec[1], a_vec[2] + c_vec[2]};
  vertices[6] = (GdisVec3) {b_vec[0] + c_vec[0], b_vec[1] + c_vec[1], b_vec[2] + c_vec[2]};
  vertices[7] = (GdisVec3) {a_vec[0] + b_vec[0] + c_vec[0],
                            a_vec[1] + b_vec[1] + c_vec[1],
                            a_vec[2] + b_vec[2] + c_vec[2]};

  return TRUE;
}

static gboolean
gdis_model_compute_bounds(const GdisModel *model,
                          gboolean include_cell,
                          gdouble min_bound[3],
                          gdouble max_bound[3])
{
  gboolean have_point;
  gint image_negative[3];
  gint image_positive[3];
  gdouble a_vec[3];
  gdouble b_vec[3];
  gdouble c_vec[3];
  gboolean have_image_vectors;
  gint ia;
  gint ib;
  gint ic;
  guint i;

  g_return_val_if_fail(model != NULL, FALSE);

  have_point = FALSE;
  have_image_vectors = FALSE;
  gdis_gtk4_window_get_image_axis_limits(model, image_negative, image_positive);
  if (model->periodic)
    have_image_vectors = gdis_model_build_cell_vectors(model, a_vec, b_vec, c_vec);
  if (!have_image_vectors)
    {
      image_negative[0] = 0;
      image_negative[1] = 0;
      image_negative[2] = 0;
      image_positive[0] = 1;
      image_positive[1] = 1;
      image_positive[2] = 1;
    }

  for (ia = -image_negative[0]; ia < image_positive[0]; ia++)
    {
      for (ib = -image_negative[1]; ib < image_positive[1]; ib++)
        {
          for (ic = -image_negative[2]; ic < image_positive[2]; ic++)
            {
              gdouble shift[3] = {0.0, 0.0, 0.0};

              if (have_image_vectors)
                {
                  shift[0] = (gdouble) ia * a_vec[0] + (gdouble) ib * b_vec[0] + (gdouble) ic * c_vec[0];
                  shift[1] = (gdouble) ia * a_vec[1] + (gdouble) ib * b_vec[1] + (gdouble) ic * c_vec[1];
                  shift[2] = (gdouble) ia * a_vec[2] + (gdouble) ib * b_vec[2] + (gdouble) ic * c_vec[2];
                }

              for (i = 0; i < model->atoms->len; i++)
                {
                  const GdisAtom *atom;
                  gdouble point[3];
                  guint axis;

                  atom = g_ptr_array_index(model->atoms, i);
                  point[0] = atom->position[0] + shift[0];
                  point[1] = atom->position[1] + shift[1];
                  point[2] = atom->position[2] + shift[2];

                  if (!have_point)
                    {
                      for (axis = 0; axis < 3; axis++)
                        {
                          min_bound[axis] = point[axis];
                          max_bound[axis] = point[axis];
                        }
                      have_point = TRUE;
                    }
                  else
                    {
                      for (axis = 0; axis < 3; axis++)
                        {
                          min_bound[axis] = MIN(min_bound[axis], point[axis]);
                          max_bound[axis] = MAX(max_bound[axis], point[axis]);
                        }
                    }
                }
            }
        }
    }

  if (include_cell && model->periodic)
    {
      GdisVec3 vertices[8];

      if (gdis_model_build_cell_vertices(model, vertices))
        {
          for (ia = -image_negative[0]; ia < image_positive[0]; ia++)
            {
              for (ib = -image_negative[1]; ib < image_positive[1]; ib++)
                {
                  for (ic = -image_negative[2]; ic < image_positive[2]; ic++)
                    {
                      gdouble shift[3] = {0.0, 0.0, 0.0};

                      if (have_image_vectors)
                        {
                          shift[0] = (gdouble) ia * a_vec[0] + (gdouble) ib * b_vec[0] + (gdouble) ic * c_vec[0];
                          shift[1] = (gdouble) ia * a_vec[1] + (gdouble) ib * b_vec[1] + (gdouble) ic * c_vec[1];
                          shift[2] = (gdouble) ia * a_vec[2] + (gdouble) ib * b_vec[2] + (gdouble) ic * c_vec[2];
                        }

                      for (i = 0; i < G_N_ELEMENTS(vertices); i++)
                        {
                          gdouble point[3] = {vertices[i].x + shift[0],
                                              vertices[i].y + shift[1],
                                              vertices[i].z + shift[2]};
                          guint axis;

                          if (!have_point)
                            {
                              for (axis = 0; axis < 3; axis++)
                                {
                                  min_bound[axis] = point[axis];
                                  max_bound[axis] = point[axis];
                                }
                              have_point = TRUE;
                            }
                          else
                            {
                              for (axis = 0; axis < 3; axis++)
                                {
                                  min_bound[axis] = MIN(min_bound[axis], point[axis]);
                                  max_bound[axis] = MAX(max_bound[axis], point[axis]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  return have_point;
}

static void
gdis_project_point(const gdouble point[3],
                   const gdouble center[3],
                   gdouble rotation_x,
                   gdouble rotation_y,
                   gdouble scale,
                   int width,
                   int height,
                   gdouble *screen_x,
                   gdouble *screen_y,
                   gdouble *depth)
{
  gdouble centered[3];
  gdouble rotated[3];

  centered[0] = point[0] - center[0];
  centered[1] = point[1] - center[1];
  centered[2] = point[2] - center[2];
  gdis_rotate_point(centered, rotation_x, rotation_y, rotated);

  if (screen_x)
    *screen_x = (width * 0.5) + rotated[0] * scale;
  if (screen_y)
    *screen_y = (height * 0.5) - rotated[1] * scale;
  if (depth)
    *depth = rotated[2];
}

static gdouble
gdis_element_draw_radius(const char *element)
{
  if (!element || !*element)
    return 0.65;

  if (g_ascii_strcasecmp(element, "H") == 0)
    return 0.32;
  if (g_ascii_strcasecmp(element, "C") == 0)
    return 0.72;
  if (g_ascii_strcasecmp(element, "N") == 0)
    return 0.68;
  if (g_ascii_strcasecmp(element, "O") == 0)
    return 0.66;
  if (g_ascii_strcasecmp(element, "S") == 0)
    return 1.02;
  if (g_ascii_strcasecmp(element, "P") == 0)
    return 1.06;
  if (g_ascii_strcasecmp(element, "Cl") == 0)
    return 0.99;
  if (g_ascii_strcasecmp(element, "Na") == 0)
    return 1.02;
  if (g_ascii_strcasecmp(element, "Ca") == 0)
    return 1.08;
  return 0.82;
}

static void
gdis_set_element_color(cairo_t *cr, const char *element, gdouble alpha)
{
  if (!element || !*element)
    {
      cairo_set_source_rgba(cr, 0.78, 0.82, 0.88, alpha);
      return;
    }

  if (g_ascii_strcasecmp(element, "H") == 0)
    cairo_set_source_rgba(cr, 0.96, 0.96, 0.96, alpha);
  else if (g_ascii_strcasecmp(element, "C") == 0)
    cairo_set_source_rgba(cr, 0.35, 0.37, 0.40, alpha);
  else if (g_ascii_strcasecmp(element, "N") == 0)
    cairo_set_source_rgba(cr, 0.20, 0.45, 0.92, alpha);
  else if (g_ascii_strcasecmp(element, "O") == 0)
    cairo_set_source_rgba(cr, 0.88, 0.18, 0.22, alpha);
  else if (g_ascii_strcasecmp(element, "S") == 0)
    cairo_set_source_rgba(cr, 0.93, 0.79, 0.24, alpha);
  else if (g_ascii_strcasecmp(element, "P") == 0)
    cairo_set_source_rgba(cr, 0.95, 0.56, 0.20, alpha);
  else if (g_ascii_strcasecmp(element, "Cl") == 0)
    cairo_set_source_rgba(cr, 0.19, 0.74, 0.30, alpha);
  else
    cairo_set_source_rgba(cr, 0.58, 0.74, 0.86, alpha);
}

static gint
gdis_projected_atom_compare(gconstpointer left, gconstpointer right)
{
  const GdisProjectedAtom *a;
  const GdisProjectedAtom *b;

  a = left;
  b = right;

  if (a->depth < b->depth)
    return -1;
  if (a->depth > b->depth)
    return 1;
  return 0;
}

static gboolean
gdis_prepare_projection(GdisGtk4Window *self,
                        int width,
                        int height,
                        GdisProjectedAtom **projected_out,
                        guint *count_out,
                        guint *cell_count_out,
                        gdouble center_out[3],
                        gdouble *scale_out)
{
  gboolean include_cell;
  gint image_negative[3];
  gint image_positive[3];
  gdouble a_vec[3];
  gdouble b_vec[3];
  gdouble c_vec[3];
  gboolean have_image_vectors;
  gdouble min_bound[3];
  gdouble max_bound[3];
  gdouble max_range;
  guint base_atom_count;
  guint cell_count;
  guint cell_index;
  gint ia;
  gint ib;
  gint ic;
  guint i;

  g_return_val_if_fail(self != NULL, FALSE);
  g_return_val_if_fail(projected_out != NULL, FALSE);
  g_return_val_if_fail(count_out != NULL, FALSE);
  g_return_val_if_fail(cell_count_out != NULL, FALSE);
  g_return_val_if_fail(center_out != NULL, FALSE);
  g_return_val_if_fail(scale_out != NULL, FALSE);

  *projected_out = NULL;
  *count_out = 0;
  *cell_count_out = 0;

  if (!self->active_model || self->active_model->atoms->len == 0)
    return FALSE;

  include_cell = self->show_cell_toggle &&
                 gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_cell_toggle));
  if (!gdis_model_compute_bounds(self->active_model, include_cell, min_bound, max_bound))
    return FALSE;

  center_out[0] = (min_bound[0] + max_bound[0]) * 0.5;
  center_out[1] = (min_bound[1] + max_bound[1]) * 0.5;
  center_out[2] = (min_bound[2] + max_bound[2]) * 0.5;

  max_range = MAX(max_bound[0] - min_bound[0], max_bound[1] - min_bound[1]);
  max_range = MAX(max_range, max_bound[2] - min_bound[2]);
  if (max_range < 1.0e-4)
    max_range = 1.0;

  base_atom_count = self->active_model->atoms->len;
  have_image_vectors = FALSE;
  gdis_gtk4_window_get_image_axis_limits(self->active_model, image_negative, image_positive);
  if (self->active_model->periodic)
    have_image_vectors = gdis_model_build_cell_vectors(self->active_model, a_vec, b_vec, c_vec);
  if (!have_image_vectors)
    {
      image_negative[0] = 0;
      image_negative[1] = 0;
      image_negative[2] = 0;
      image_positive[0] = 1;
      image_positive[1] = 1;
      image_positive[2] = 1;
    }

  cell_count = (guint) (image_negative[0] + image_positive[0]) *
               (guint) (image_negative[1] + image_positive[1]) *
               (guint) (image_negative[2] + image_positive[2]);
  if (cell_count == 0)
    cell_count = 1;

  *scale_out = (MIN(width, height) * 0.32 * self->zoom) / max_range;
  *count_out = base_atom_count * cell_count;
  *cell_count_out = cell_count;
  *projected_out = g_new0(GdisProjectedAtom, *count_out);

  cell_index = 0;
  for (ia = -image_negative[0]; ia < image_positive[0]; ia++)
    {
      for (ib = -image_negative[1]; ib < image_positive[1]; ib++)
        {
          for (ic = -image_negative[2]; ic < image_positive[2]; ic++)
            {
              gdouble shift[3] = {0.0, 0.0, 0.0};

              if (have_image_vectors)
                {
                  shift[0] = (gdouble) ia * a_vec[0] + (gdouble) ib * b_vec[0] + (gdouble) ic * c_vec[0];
                  shift[1] = (gdouble) ia * a_vec[1] + (gdouble) ib * b_vec[1] + (gdouble) ic * c_vec[1];
                  shift[2] = (gdouble) ia * a_vec[2] + (gdouble) ib * b_vec[2] + (gdouble) ic * c_vec[2];
                }

              for (i = 0; i < base_atom_count; i++)
                {
                  const GdisAtom *atom;
                  gdouble point[3];
                  guint projected_index;

                  projected_index = cell_index * base_atom_count + i;
                  atom = g_ptr_array_index(self->active_model->atoms, i);
                  point[0] = atom->position[0] + shift[0];
                  point[1] = atom->position[1] + shift[1];
                  point[2] = atom->position[2] + shift[2];
                  (*projected_out)[projected_index].atom_index = i;
                  (*projected_out)[projected_index].cell_index = cell_index;
                  (*projected_out)[projected_index].image_offset[0] = ia;
                  (*projected_out)[projected_index].image_offset[1] = ib;
                  (*projected_out)[projected_index].image_offset[2] = ic;
                  gdis_project_point(point,
                                     center_out,
                                     self->rotation_x,
                                     self->rotation_y,
                                     *scale_out,
                                     width,
                                     height,
                                     &(*projected_out)[projected_index].screen_x,
                                     &(*projected_out)[projected_index].screen_y,
                                     &(*projected_out)[projected_index].depth);
                  (*projected_out)[projected_index].radius =
                    CLAMP(*scale_out * gdis_element_draw_radius(atom->element) * 0.48,
                          2.8,
                          22.0);
                }

              cell_index++;
            }
        }
    }

  return TRUE;
}

static void
gdis_draw_placeholder(GdisGtk4Window *self, cairo_t *cr, int width, int height)
{
  const char *line_one;
  const char *line_two;

  line_one = "Open a sample model or use File > Open to load a structure.";
  line_two = "Drag to rotate, scroll to zoom, and click atoms once a model is loaded.";
  if (self && self->active_model && self->active_model->atom_count == 0)
    {
      line_one = "This model is empty. Open Edit > Edit Structure to add atoms.";
      line_two = "After adding atoms, use File > Save to write an XYZ, PDB, ARC/CAR, or CIF file.";
    }

  cairo_set_source_rgb(cr, 0.02, 0.03, 0.05);
  cairo_paint(cr);

  cairo_set_source_rgba(cr, 0.16, 0.24, 0.31, 1.0);
  cairo_rectangle(cr, 1.0, 1.0, width - 2.0, height - 2.0);
  cairo_set_line_width(cr, 1.0);
  cairo_stroke(cr);

  cairo_set_source_rgba(cr, 0.82, 0.86, 0.90, 1.0);
  cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, 22.0);
  cairo_move_to(cr, 34.0, 54.0);
  cairo_show_text(cr, "GDIS GTK4");

  cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, 14.0);
  cairo_move_to(cr, 34.0, 84.0);
  cairo_show_text(cr, line_one);
  cairo_move_to(cr, 34.0, 110.0);
  cairo_show_text(cr, line_two);
}

static void
gdis_draw_cell_edges(cairo_t *cr,
                     GdisGtk4Window *self,
                     const gdouble center[3],
                     gdouble scale,
                     int width,
                     int height)
{
  static const guint8 edges[][2] = {
    {0, 1}, {0, 2}, {0, 4},
    {1, 3}, {1, 5},
    {2, 3}, {2, 6},
    {3, 7},
    {4, 5}, {4, 6},
    {5, 7},
    {6, 7}
  };
  GdisVec3 vertices[8];
  gint image_negative[3];
  gint image_positive[3];
  gdouble a_vec[3];
  gdouble b_vec[3];
  gdouble c_vec[3];
  gboolean have_image_vectors;
  gint ia;
  gint ib;
  gint ic;
  guint i;

  if (!self->active_model || !gdis_model_build_cell_vertices(self->active_model, vertices))
    return;

  have_image_vectors = FALSE;
  gdis_gtk4_window_get_image_axis_limits(self->active_model, image_negative, image_positive);
  if (self->active_model->periodic)
    have_image_vectors = gdis_model_build_cell_vectors(self->active_model, a_vec, b_vec, c_vec);
  if (!have_image_vectors)
    {
      image_negative[0] = 0;
      image_negative[1] = 0;
      image_negative[2] = 0;
      image_positive[0] = 1;
      image_positive[1] = 1;
      image_positive[2] = 1;
    }

  cairo_save(cr);
  cairo_set_source_rgba(cr, 0.42, 0.74, 0.92, 0.70);
  cairo_set_line_width(cr, 1.4);
  for (ia = -image_negative[0]; ia < image_positive[0]; ia++)
    {
      for (ib = -image_negative[1]; ib < image_positive[1]; ib++)
        {
          for (ic = -image_negative[2]; ic < image_positive[2]; ic++)
            {
              gdouble shift[3] = {0.0, 0.0, 0.0};
              gdouble screen[8][2];

              if (have_image_vectors)
                {
                  shift[0] = (gdouble) ia * a_vec[0] + (gdouble) ib * b_vec[0] + (gdouble) ic * c_vec[0];
                  shift[1] = (gdouble) ia * a_vec[1] + (gdouble) ib * b_vec[1] + (gdouble) ic * c_vec[1];
                  shift[2] = (gdouble) ia * a_vec[2] + (gdouble) ib * b_vec[2] + (gdouble) ic * c_vec[2];
                }

              for (i = 0; i < G_N_ELEMENTS(vertices); i++)
                {
                  gdouble point[3] = {vertices[i].x + shift[0],
                                      vertices[i].y + shift[1],
                                      vertices[i].z + shift[2]};

                  gdis_project_point(point,
                                     center,
                                     self->rotation_x,
                                     self->rotation_y,
                                     scale,
                                     width,
                                     height,
                                     &screen[i][0],
                                     &screen[i][1],
                                     NULL);
                }

              for (i = 0; i < G_N_ELEMENTS(edges); i++)
                {
                  cairo_move_to(cr, screen[edges[i][0]][0], screen[edges[i][0]][1]);
                  cairo_line_to(cr, screen[edges[i][1]][0], screen[edges[i][1]][1]);
                }
            }
        }
    }
  cairo_stroke(cr);
  cairo_restore(cr);
}

static void
gdis_draw_overlay(GdisGtk4Window *self, cairo_t *cr, int width, int height)
{
  gchar *top_left;
  gchar *bottom_left;
  gint image_negative[3];
  gint image_positive[3];
  guint displayed_cells;

  g_return_if_fail(self != NULL);

  if (!self->active_model)
    return;

  gdis_gtk4_window_get_image_axis_limits(self->active_model, image_negative, image_positive);
  displayed_cells = (guint) (image_negative[0] + image_positive[0]) *
                    (guint) (image_negative[1] + image_positive[1]) *
                    (guint) (image_negative[2] + image_positive[2]);

  top_left = g_strdup_printf("%s  |  %s  |  %u atoms  |  %u bonds  |  %u cell%s",
                             self->active_model->basename,
                             self->active_model->format_label,
                             self->active_model->atom_count,
                             self->active_model->bond_count,
                             displayed_cells,
                             displayed_cells == 1 ? "" : "s");
  bottom_left = g_strdup_printf("Drag: rotate    Scroll: zoom    Click: select atom    Zoom %.2fx",
                                self->zoom);

  cairo_save(cr);
  cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, 13.0);
  cairo_set_source_rgba(cr, 0.92, 0.95, 0.97, 1.0);
  cairo_move_to(cr, 16.0, 24.0);
  cairo_show_text(cr, top_left);

  cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, 12.0);
  cairo_set_source_rgba(cr, 0.70, 0.78, 0.84, 0.95);
  cairo_move_to(cr, 16.0, height - 16.0);
  cairo_show_text(cr, bottom_left);

  if (width > 520)
    {
      cairo_move_to(cr, width - 112.0, 24.0);
      cairo_show_text(cr, "X");
      cairo_set_source_rgba(cr, 0.45, 0.70, 0.95, 0.95);
      cairo_move_to(cr, width - 94.0, 24.0);
      cairo_show_text(cr, "Y");
      cairo_set_source_rgba(cr, 0.86, 0.45, 0.40, 0.95);
      cairo_move_to(cr, width - 76.0, 24.0);
      cairo_show_text(cr, "Z");
    }
  cairo_restore(cr);

  g_free(top_left);
  g_free(bottom_left);
}

static void
viewer_draw(GtkDrawingArea *area,
            cairo_t *cr,
            int width,
            int height,
            gpointer data)
{
  GdisGtk4Window *self;
  gboolean show_atoms;
  gboolean show_bonds;
  gboolean show_cell;
  gboolean show_labels;
  gdouble center[3];
  gdouble scale;
  GdisProjectedAtom *projected;
  GdisProjectedAtom *draw_order;
  guint atom_count;
  guint base_atom_count;
  guint cell_count;
  guint i;

  (void) area;

  self = data;
  projected = NULL;
  draw_order = NULL;
  atom_count = 0;
  base_atom_count = 0;
  cell_count = 0;

  if (!self || !self->active_model)
    {
      gdis_draw_placeholder(self, cr, width, height);
      return;
    }

  show_atoms = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_atoms_toggle));
  show_bonds = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_bonds_toggle));
  show_cell = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_cell_toggle));
  show_labels = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_labels_toggle));

  cairo_set_source_rgb(cr, 0.02, 0.03, 0.05);
  cairo_paint(cr);

  cairo_set_source_rgba(cr, 0.14, 0.21, 0.28, 1.0);
  cairo_set_line_width(cr, 1.0);
  cairo_rectangle(cr, 1.0, 1.0, width - 2.0, height - 2.0);
  cairo_stroke(cr);

  if (!gdis_prepare_projection(self,
                               width,
                               height,
                               &projected,
                               &atom_count,
                               &cell_count,
                               center,
                               &scale))
    {
      gdis_draw_placeholder(self, cr, width, height);
      return;
    }
  base_atom_count = self->active_model->atoms->len;

  if (show_cell)
    gdis_draw_cell_edges(cr, self, center, scale, width, height);

  if (show_bonds)
    {
      cairo_save(cr);
      cairo_set_line_width(cr, CLAMP(scale * 0.10, 1.0, 3.2));
      cairo_set_source_rgba(cr, 0.74, 0.79, 0.84, 0.70);

      for (guint cell_index = 0; cell_index < cell_count; cell_index++)
        {
          guint cell_offset;

          cell_offset = cell_index * base_atom_count;
          for (i = 0; i < self->active_model->bonds->len; i++)
            {
              GdisBond *bond;

              bond = &g_array_index(self->active_model->bonds, GdisBond, i);
              if (bond->atom_index_a >= base_atom_count || bond->atom_index_b >= base_atom_count)
                continue;

              cairo_move_to(cr,
                            projected[cell_offset + bond->atom_index_a].screen_x,
                            projected[cell_offset + bond->atom_index_a].screen_y);
              cairo_line_to(cr,
                            projected[cell_offset + bond->atom_index_b].screen_x,
                            projected[cell_offset + bond->atom_index_b].screen_y);
            }
        }

      cairo_stroke(cr);
      cairo_restore(cr);
    }

  if (show_atoms)
    {
      draw_order = g_memdup2(projected, sizeof(*projected) * atom_count);
      qsort(draw_order, atom_count, sizeof(*draw_order), gdis_projected_atom_compare);

      for (i = 0; i < atom_count; i++)
        {
          const GdisAtom *atom;
          gboolean selected;
          gboolean group_selected;
          gdouble ring_radius;

          atom = g_ptr_array_index(self->active_model->atoms, draw_order[i].atom_index);
          selected = (draw_order[i].atom_index == self->selected_atom_index);
          group_selected = gdis_gtk4_window_atom_array_contains(self->selected_atoms,
                                                                draw_order[i].atom_index);

          cairo_save(cr);
          cairo_arc(cr,
                    draw_order[i].screen_x,
                    draw_order[i].screen_y,
                    draw_order[i].radius,
                    0.0,
                    2.0 * G_PI);
          gdis_set_element_color(cr, atom->element, 0.97);
          cairo_fill_preserve(cr);
          cairo_set_source_rgba(cr, 0.05, 0.07, 0.09, 0.95);
          cairo_set_line_width(cr, selected ? 2.2 : 1.0);
          cairo_stroke(cr);

          if (group_selected)
            {
              ring_radius = draw_order[i].radius + (selected ? 7.0 : 4.0);
              cairo_arc(cr,
                        draw_order[i].screen_x,
                        draw_order[i].screen_y,
                        ring_radius,
                        0.0,
                        2.0 * G_PI);
              cairo_set_source_rgba(cr, 0.43, 0.76, 0.97, 0.92);
              cairo_set_line_width(cr, 1.6);
              cairo_stroke(cr);
            }

          if (selected)
            {
              ring_radius = draw_order[i].radius + 9.5;
              cairo_arc(cr,
                        draw_order[i].screen_x,
                        draw_order[i].screen_y,
                        ring_radius,
                        0.0,
                        2.0 * G_PI);
              cairo_set_source_rgba(cr, 0.98, 0.85, 0.28, 0.95);
              cairo_set_line_width(cr, 2.0);
              cairo_stroke(cr);
            }

          cairo_restore(cr);
        }
    }

  if (show_labels && atom_count <= 180)
    {
      cairo_save(cr);
      cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
      cairo_set_font_size(cr, 11.0);
      cairo_set_source_rgba(cr, 0.92, 0.95, 0.98, 0.94);
      for (i = 0; i < atom_count; i++)
        {
          const GdisAtom *atom;

          if (projected[i].image_offset[0] != 0 ||
              projected[i].image_offset[1] != 0 ||
              projected[i].image_offset[2] != 0)
            continue;

          atom = g_ptr_array_index(self->active_model->atoms, projected[i].atom_index);
          cairo_move_to(cr,
                        projected[i].screen_x + projected[i].radius + 3.0,
                        projected[i].screen_y - projected[i].radius - 2.0);
          cairo_show_text(cr, atom->label && *atom->label ? atom->label : atom->element);
        }
      cairo_restore(cr);
    }

  gdis_draw_overlay(self, cr, width, height);

  g_free(draw_order);
  g_free(projected);
}

static void
gdis_gtk4_window_select_atom_at(GdisGtk4Window *self, gdouble x, gdouble y)
{
  GdisProjectedAtom *projected;
  gdouble center[3];
  gdouble scale;
  guint atom_count;
  guint cell_count;
  guint i;
  guint best_index;
  guint best_projected_index;
  gdouble best_distance2;

  g_return_if_fail(self != NULL);

  projected = NULL;
  atom_count = 0;
  cell_count = 0;
  best_index = INVALID_ATOM_INDEX;
  best_projected_index = INVALID_ATOM_INDEX;
  best_distance2 = G_MAXDOUBLE;

  if (!self->active_model || !self->viewer_area)
    return;

  if (!gdis_prepare_projection(self,
                               gtk_widget_get_width(self->viewer_area),
                               gtk_widget_get_height(self->viewer_area),
                               &projected,
                               &atom_count,
                               &cell_count,
                               center,
                               &scale))
    return;
  (void) cell_count;

  for (i = 0; i < atom_count; i++)
    {
      gdouble dx;
      gdouble dy;
      gdouble distance2;
      gdouble threshold;

      dx = projected[i].screen_x - x;
      dy = projected[i].screen_y - y;
      distance2 = dx * dx + dy * dy;
      threshold = MAX(projected[i].radius + 5.0, 12.0);

      if (distance2 <= threshold * threshold && distance2 < best_distance2)
        {
          best_distance2 = distance2;
          best_index = projected[i].atom_index;
          best_projected_index = i;
        }
    }

  if (best_index != INVALID_ATOM_INDEX)
    {
      const GdisAtom *atom;

      self->selected_atom_index = best_index;
      gdis_gtk4_window_apply_selection_mode(self, best_index);
      gdis_gtk4_window_remember_atom_pick(self, best_index);
      atom = g_ptr_array_index(self->active_model->atoms, best_index);
      gdis_gtk4_window_log(self,
                           "Selected atom in %s: %s [%s] #%u at %.4f %.4f %.4f (image offset: %d %d %d, mode: %s, selected set: %u, pick set: %u)\n",
                           self->active_model->basename,
                           atom->label,
                           atom->element,
                           atom->serial,
                           atom->position[0],
                           atom->position[1],
                           atom->position[2],
                           projected[best_projected_index].image_offset[0],
                           projected[best_projected_index].image_offset[1],
                           projected[best_projected_index].image_offset[2],
                           gdis_gtk4_window_selection_mode_label(self->selection_mode),
                           self->selected_atoms ? self->selected_atoms->len : 0,
                           self->picked_atoms ? self->picked_atoms->len : 0);

      if (self->click_mode != GDIS_CLICK_MODE_SELECT &&
          self->picked_atoms &&
          self->picked_atoms->len >= 2)
        {
          GError *error;
          guint atom_index_a;
          guint atom_index_b;

          error = NULL;
          if (gdis_gtk4_window_get_last_two_picks(self, &atom_index_a, &atom_index_b, &error))
            {
              if (self->click_mode == GDIS_CLICK_MODE_ADD_BOND)
                {
                  guint order;

                  order = 1;
                  if (self->edit_tool && self->edit_tool->bond_order_entry)
                    {
                      GError *parse_error;

                      parse_error = NULL;
                      if (gdis_gtk4_window_parse_entry_uint(self->edit_tool->bond_order_entry,
                                                            "Bond order",
                                                            &order,
                                                            &parse_error))
                        {
                          (void) parse_error;
                        }
                      else
                        {
                          g_clear_error(&parse_error);
                          order = 1;
                        }
                    }

                  if (!gdis_model_add_explicit_bond(self->active_model,
                                                    atom_index_a,
                                                    atom_index_b,
                                                    (guint8) order,
                                                    &error))
                    gdis_gtk4_window_log(self, "Pick Add Bond failed: %s\n",
                                         error ? error->message : "unknown error");
                  else
                    {
                      gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
                      gdis_gtk4_window_log(self, "Pick Add Bond completed.\n");
                    }
                }
              else if (!gdis_model_remove_bond(self->active_model,
                                               atom_index_a,
                                               atom_index_b,
                                               &error))
                {
                  gdis_gtk4_window_log(self, "Pick Remove Bond failed: %s\n",
                                       error ? error->message : "unknown error");
                }
              else
                {
                  gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
                  gdis_gtk4_window_log(self, "Pick Remove Bond completed.\n");
                }
            }

          g_clear_error(&error);
          gdis_gtk4_window_set_click_mode(self, GDIS_CLICK_MODE_SELECT);
          gdis_gtk4_window_log(self, "Click mode returned to Select.\n");
        }
    }
  else if (self->selected_atom_index != INVALID_ATOM_INDEX)
    {
      self->selected_atom_index = INVALID_ATOM_INDEX;
      gdis_gtk4_window_clear_selected_atoms(self);
      gdis_gtk4_window_log(self, "Selection cleared.\n");
    }

  g_free(projected);
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
}

static GtkWidget *
new_section_button(const char *text)
{
  GtkWidget *button;

  button = gtk_button_new_with_label(text);
  gtk_widget_add_css_class(button, "flat");
  gtk_widget_add_css_class(button, "gdis-section-button");
  gtk_widget_set_sensitive(button, FALSE);
  gtk_widget_set_halign(button, GTK_ALIGN_START);

  return button;
}

static GtkWidget *
new_readonly_text_view(GtkTextBuffer **buffer_out,
                       gboolean monospace,
                       int min_height)
{
  GtkWidget *scroller;
  GtkWidget *text_view;

  scroller = gtk_scrolled_window_new();
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                 monospace ? GTK_POLICY_AUTOMATIC : GTK_POLICY_NEVER,
                                 GTK_POLICY_AUTOMATIC);
  gtk_scrolled_window_set_min_content_height(GTK_SCROLLED_WINDOW(scroller), min_height);
  gtk_widget_set_vexpand(scroller, TRUE);

  text_view = gtk_text_view_new();
  gtk_text_view_set_editable(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view),
                              monospace ? GTK_WRAP_NONE : GTK_WRAP_WORD_CHAR);
  gtk_text_view_set_monospace(GTK_TEXT_VIEW(text_view), monospace);
  gtk_widget_set_focusable(text_view, FALSE);
  gtk_text_view_set_left_margin(GTK_TEXT_VIEW(text_view), 8);
  gtk_text_view_set_right_margin(GTK_TEXT_VIEW(text_view), 8);
  gtk_text_view_set_top_margin(GTK_TEXT_VIEW(text_view), 8);
  gtk_text_view_set_bottom_margin(GTK_TEXT_VIEW(text_view), 8);
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  if (buffer_out)
    *buffer_out = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));

  return scroller;
}

static GtkWidget *
build_content_page(GdisGtk4Window *self)
{
  GtkWidget *box;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_box_append(GTK_BOX(box), new_readonly_text_view(&self->content_buffer, FALSE, 180));
  gdis_text_buffer_set(self->content_buffer, "No active model loaded.");

  return box;
}

static GtkWidget *
build_editing_page(GdisGtk4Window *self)
{
  GtkWidget *box;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_box_append(GTK_BOX(box), new_readonly_text_view(&self->editing_buffer, TRUE, 180));
  gdis_text_buffer_set(self->editing_buffer, "No active model loaded.");

  return box;
}

static GtkWidget *
new_toggle_button(const char *label, gboolean active)
{
  GtkWidget *button;

  button = gtk_toggle_button_new_with_label(label);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), active);

  return button;
}

static GtkWidget *
build_display_page(GdisGtk4Window *self)
{
  GtkWidget *box;
  GtkWidget *grid;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);

  self->show_atoms_toggle = new_toggle_button("Atoms", TRUE);
  self->show_bonds_toggle = new_toggle_button("Bonds", TRUE);
  self->show_cell_toggle = new_toggle_button("Cell", TRUE);
  self->show_labels_toggle = new_toggle_button("Labels", FALSE);

  g_signal_connect(self->show_atoms_toggle, "toggled", G_CALLBACK(on_view_toggle_toggled), self);
  g_signal_connect(self->show_bonds_toggle, "toggled", G_CALLBACK(on_view_toggle_toggled), self);
  g_signal_connect(self->show_cell_toggle, "toggled", G_CALLBACK(on_view_toggle_toggled), self);
  g_signal_connect(self->show_labels_toggle, "toggled", G_CALLBACK(on_view_toggle_toggled), self);

  gtk_widget_set_hexpand(self->show_atoms_toggle, TRUE);
  gtk_widget_set_hexpand(self->show_bonds_toggle, TRUE);
  gtk_widget_set_hexpand(self->show_cell_toggle, TRUE);
  gtk_widget_set_hexpand(self->show_labels_toggle, TRUE);

  gtk_grid_attach(GTK_GRID(grid), self->show_atoms_toggle, 0, 0, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), self->show_bonds_toggle, 1, 0, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), self->show_cell_toggle, 0, 1, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), self->show_labels_toggle, 1, 1, 1, 1);
  gtk_box_append(GTK_BOX(box), grid);

  return box;
}

static GtkWidget *
build_images_page(GdisGtk4Window *self)
{
  GtkWidget *box;
  GtkWidget *grid;
  GtkWidget *label;
  GtkWidget *button;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_box_append(GTK_BOX(box), new_section_button("Image Range"));

  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);

  label = gtk_label_new("");
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  label = gtk_label_new("-");
  gtk_grid_attach(GTK_GRID(grid), label, 1, 0, 1, 1);
  label = gtk_label_new("+");
  gtk_grid_attach(GTK_GRID(grid), label, 2, 0, 1, 1);

  for (int axis = 0; axis < 3; axis++)
    {
      static const char *const axis_labels[] = {"a", "b", "c"};

      label = gtk_label_new(axis_labels[axis]);
      gtk_grid_attach(GTK_GRID(grid), label, 0, axis + 1, 1, 1);

      self->image_limit_spin[2 * axis] = gtk_spin_button_new_with_range(0.0, 10.0, 1.0);
      gtk_spin_button_set_digits(GTK_SPIN_BUTTON(self->image_limit_spin[2 * axis]), 0);
      gtk_grid_attach(GTK_GRID(grid), self->image_limit_spin[2 * axis], 1, axis + 1, 1, 1);

      self->image_limit_spin[2 * axis + 1] = gtk_spin_button_new_with_range(1.0, 10.0, 1.0);
      gtk_spin_button_set_digits(GTK_SPIN_BUTTON(self->image_limit_spin[2 * axis + 1]), 0);
      gtk_grid_attach(GTK_GRID(grid), self->image_limit_spin[2 * axis + 1], 2, axis + 1, 1, 1);
    }

  gtk_box_append(GTK_BOX(box), grid);

  grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  button = gtk_button_new_with_label("Apply Image Range");
  g_signal_connect(button, "clicked", G_CALLBACK(on_apply_image_limits_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 0, 0, 1, 1);
  button = gtk_button_new_with_label("Reset Image Range");
  g_signal_connect(button, "clicked", G_CALLBACK(on_reset_image_limits_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 1, 0, 1, 1);
  gtk_box_append(GTK_BOX(box), grid);

  gtk_box_append(GTK_BOX(box), new_section_button("Crystal Tools"));
  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);

  label = gtk_label_new("Supercell repeats");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  for (int axis = 0; axis < 3; axis++)
    {
      static const char *const axis_labels[] = {"a", "b", "c"};

      label = gtk_label_new(axis_labels[axis]);
      gtk_grid_attach(GTK_GRID(grid), label, axis * 2 + 1, 0, 1, 1);
      self->supercell_repeat_spin[axis] = gtk_spin_button_new_with_range(1.0, 8.0, 1.0);
      gtk_spin_button_set_digits(GTK_SPIN_BUTTON(self->supercell_repeat_spin[axis]), 0);
      gtk_grid_attach(GTK_GRID(grid), self->supercell_repeat_spin[axis], axis * 2 + 2, 0, 1, 1);
    }
  gtk_box_append(GTK_BOX(box), grid);

  grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  button = gtk_button_new_with_label("Bake Supercell");
  g_signal_connect(button, "clicked", G_CALLBACK(on_make_supercell_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 0, 0, 1, 1);
  button = gtk_button_new_with_label("Confine To Cell");
  g_signal_connect(button, "clicked", G_CALLBACK(on_confine_to_cell_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 1, 0, 1, 1);
  button = gtk_button_new_with_label("Confine Molecules");
  g_signal_connect(button, "clicked", G_CALLBACK(on_confine_molecules_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 0, 1, 1, 1);
  button = gtk_button_new_with_label("Force P1");
  g_signal_connect(button, "clicked", G_CALLBACK(on_force_p1_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 1, 1, 1, 1);
  gtk_box_append(GTK_BOX(box), grid);

  gtk_box_append(GTK_BOX(box), new_readonly_text_view(&self->images_buffer, FALSE, 180));
  gdis_text_buffer_set(self->images_buffer, "No active model loaded.");

  return box;
}

static GtkWidget *
build_symmetry_page(GdisGtk4Window *self)
{
  GtkWidget *box;
  GtkWidget *row;
  GtkWidget *button;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  button = gtk_button_new_with_label("Force P1");
  g_signal_connect(button, "clicked", G_CALLBACK(on_force_p1_clicked), self);
  gtk_box_append(GTK_BOX(row), button);
  gtk_box_append(GTK_BOX(box), row);
  gtk_box_append(GTK_BOX(box), new_readonly_text_view(&self->symmetry_buffer, FALSE, 180));
  gdis_text_buffer_set(self->symmetry_buffer, "No active model loaded.");

  return box;
}

static GtkWidget *
build_viewing_page(GdisGtk4Window *self)
{
  GtkWidget *box;
  GtkWidget *grid;
  GtkWidget *button;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_box_append(GTK_BOX(box), new_section_button("Viewing Controls"));

  grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);

  button = gtk_button_new_with_label("Reset View");
  g_signal_connect(button, "clicked", G_CALLBACK(on_reset_view_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 0, 0, 1, 1);

  button = gtk_button_new_with_label("View X");
  g_signal_connect(button, "clicked", G_CALLBACK(on_view_x_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 1, 0, 1, 1);

  button = gtk_button_new_with_label("View Y");
  g_signal_connect(button, "clicked", G_CALLBACK(on_view_y_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 0, 1, 1, 1);

  button = gtk_button_new_with_label("View Z");
  g_signal_connect(button, "clicked", G_CALLBACK(on_view_z_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 1, 1, 1, 1);

  gtk_box_append(GTK_BOX(box), grid);

  return box;
}

static GtkWidget *
build_sidebar_stack(GdisGtk4Window *self)
{
  GtkWidget *stack;
  GtkWidget *page;

  stack = gtk_stack_new();
  gtk_stack_set_transition_type(GTK_STACK(stack), GTK_STACK_TRANSITION_TYPE_SLIDE_LEFT_RIGHT);
  gtk_widget_set_vexpand(stack, TRUE);

  page = build_content_page(self);
  gtk_stack_add_titled(GTK_STACK(stack), page, "content", "Content");

  page = build_editing_page(self);
  gtk_stack_add_titled(GTK_STACK(stack), page, "editing", "Editing");

  page = build_display_page(self);
  gtk_stack_add_titled(GTK_STACK(stack), page, "display", "Display");

  page = build_images_page(self);
  gtk_stack_add_titled(GTK_STACK(stack), page, "images", "Images");

  page = build_symmetry_page(self);
  gtk_stack_add_titled(GTK_STACK(stack), page, "symmetry", "Symmetry");

  page = build_viewing_page(self);
  gtk_stack_add_titled(GTK_STACK(stack), page, "viewing", "Viewing");

  return stack;
}

static void
on_model_button_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisModel *model;

  self = user_data;
  model = g_object_get_data(G_OBJECT(button), "gdis-model");
  if (model)
    gdis_gtk4_window_set_active_model(self, model);
}

static void
on_selection_button_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  const char *name;

  self = user_data;
  name = g_object_get_data(G_OBJECT(button), "selection-name");
  if (name)
    {
      if (g_strcmp0(name, "Select : Atom Label") == 0)
        self->selection_mode = GDIS_SELECTION_MODE_LABEL;
      else if (g_strcmp0(name, "Select : Atom FF Type") == 0)
        self->selection_mode = GDIS_SELECTION_MODE_FF_TYPE;
      else if (g_strcmp0(name, "Select : Elements in Molecule") == 0)
        self->selection_mode = GDIS_SELECTION_MODE_ELEMENTS_IN_MOLECULE;
      else if (g_strcmp0(name, "Select : Molecules") == 0)
        self->selection_mode = GDIS_SELECTION_MODE_MOLECULES;
      else if (g_strcmp0(name, "Select : Molecule Fragments") == 0)
        self->selection_mode = GDIS_SELECTION_MODE_MOLECULE_FRAGMENTS;
      else if (g_strcmp0(name, "Select : Regions") == 0)
        self->selection_mode = GDIS_SELECTION_MODE_REGIONS;
      else if (g_strcmp0(name, "Select : Elements") == 0)
        self->selection_mode = GDIS_SELECTION_MODE_ELEMENTS;
      else
        self->selection_mode = GDIS_SELECTION_MODE_ATOMS;

      if (self->selection_mode != GDIS_SELECTION_MODE_MOLECULE_FRAGMENTS)
        self->fragment_anchor_index = INVALID_ATOM_INDEX;

      if (self->active_model &&
          self->selected_atom_index != INVALID_ATOM_INDEX &&
          self->selected_atom_index < self->active_model->atoms->len)
        gdis_gtk4_window_apply_selection_mode(self, self->selected_atom_index);

      gdis_gtk4_window_update_details(self);
      gdis_gtk4_window_refresh_viewer(self);
      gdis_gtk4_window_refresh_edit_tool(self);
      gdis_gtk4_window_log(self,
                           "Selection mode set to: %s\n",
                           gdis_gtk4_window_selection_mode_label(self->selection_mode));
    }
}

static void
on_view_toggle_toggled(GtkToggleButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
}

static void
on_reset_view_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  gdis_gtk4_window_reset_view(self);
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_log(self, "View reset.\n");
}

static void
on_view_x_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  gdis_gtk4_window_apply_axis_view(self, 0.0, G_PI / 2.0);
  gdis_gtk4_window_log(self, "View changed to X axis.\n");
}

static void
on_view_y_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  gdis_gtk4_window_apply_axis_view(self, -G_PI / 2.0, 0.0);
  gdis_gtk4_window_log(self, "View changed to Y axis.\n");
}

static void
on_view_z_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  gdis_gtk4_window_apply_axis_view(self, 0.0, 0.0);
  gdis_gtk4_window_log(self, "View changed to Z axis.\n");
}

static void
on_apply_image_limits_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  gint limits[6];

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  for (guint axis = 0; axis < 3; axis++)
    {
      limits[2 * axis] =
        gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->image_limit_spin[2 * axis]));
      limits[2 * axis + 1] =
        gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->image_limit_spin[2 * axis + 1]));
    }

  error = NULL;
  if (!gdis_model_set_image_limits(self->active_model, limits, &error))
    {
      gdis_gtk4_window_log(self, "Image range update failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      gdis_gtk4_window_sync_image_controls(self);
      return;
    }

  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_log(self, "Updated the model image range.\n");
}

static void
on_reset_image_limits_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  gint limits[6] = {0, 1, 0, 1, 0, 1};

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  error = NULL;
  if (!gdis_model_set_image_limits(self->active_model, limits, &error))
    {
      gdis_gtk4_window_log(self, "Image range reset failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_sync_image_controls(self);
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_log(self, "Reset the model image range to the legacy defaults.\n");
}

static void
on_confine_to_cell_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  error = NULL;
  if (!gdis_model_confine_atoms_to_cell(self->active_model, &error))
    {
      gdis_gtk4_window_log(self, "Confine to cell failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_log(self, "Confined atoms to the active unit cell.\n");
}

static void
on_confine_molecules_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  error = NULL;
  if (!gdis_model_confine_molecules_to_cell(self->active_model, &error))
    {
      gdis_gtk4_window_log(self, "Confine molecules failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_log(self, "Confined molecular fragments to the active unit cell.\n");
}

static void
on_force_p1_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  error = NULL;
  if (!gdis_model_force_p1(self->active_model, &error))
    {
      gdis_gtk4_window_log(self, "Force to P1 failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_log(self, "Forced the active model to space group P 1.\n");
}

static void
on_make_supercell_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  guint repeat_a;
  guint repeat_b;
  guint repeat_c;

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  repeat_a = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->supercell_repeat_spin[0]));
  repeat_b = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->supercell_repeat_spin[1]));
  repeat_c = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->supercell_repeat_spin[2]));

  error = NULL;
  if (!gdis_model_make_supercell(self->active_model, repeat_a, repeat_b, repeat_c, &error))
    {
      gdis_gtk4_window_log(self, "Make supercell failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  self->selected_atom_index = INVALID_ATOM_INDEX;
  self->fragment_anchor_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_clear_atom_picks(self);
  for (guint axis = 0; axis < 3; axis++)
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(self->supercell_repeat_spin[axis]), 1.0);
  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_log(self,
                       "Built a %ux%ux%u supercell from the active crystal model.\n",
                       repeat_a,
                       repeat_b,
                       repeat_c);
}

static void
on_viewer_drag_begin(GtkGestureDrag *gesture,
                     gdouble start_x,
                     gdouble start_y,
                     gpointer user_data)
{
  GdisGtk4Window *self;

  (void) gesture;
  (void) start_x;
  (void) start_y;

  self = user_data;
  self->drag_origin_x = self->rotation_x;
  self->drag_origin_y = self->rotation_y;
}

static void
on_viewer_drag_update(GtkGestureDrag *gesture,
                      gdouble offset_x,
                      gdouble offset_y,
                      gpointer user_data)
{
  GdisGtk4Window *self;

  (void) gesture;

  self = user_data;
  self->rotation_y = self->drag_origin_y + offset_x * 0.010;
  self->rotation_x = self->drag_origin_x + offset_y * 0.010;
  gdis_gtk4_window_refresh_viewer(self);
}

static gboolean
on_viewer_scroll(GtkEventControllerScroll *controller,
                 gdouble dx,
                 gdouble dy,
                 gpointer user_data)
{
  GdisGtk4Window *self;
  gdouble dominant;

  (void) controller;

  self = user_data;
  dominant = fabs(dy) > fabs(dx) ? dy : dx;
  if (dominant < 0.0)
    self->zoom *= 1.10;
  else if (dominant > 0.0)
    self->zoom /= 1.10;
  self->zoom = CLAMP(self->zoom, 0.15, 8.0);
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);

  return TRUE;
}

static void
on_viewer_click_pressed(GtkGestureClick *gesture,
                        gint n_press,
                        gdouble x,
                        gdouble y,
                        gpointer user_data)
{
  GdisGtk4Window *self;

  (void) gesture;
  (void) n_press;

  self = user_data;
  self->press_x = x;
  self->press_y = y;
}

static void
on_viewer_click_released(GtkGestureClick *gesture,
                         gint n_press,
                         gdouble x,
                         gdouble y,
                         gpointer user_data)
{
  GdisGtk4Window *self;

  (void) gesture;
  (void) n_press;

  self = user_data;
  if (fabs(x - self->press_x) <= 5.0 && fabs(y - self->press_y) <= 5.0)
    gdis_gtk4_window_select_atom_at(self, x, y);
}

static void
on_open_dialog_complete(GObject *source_object,
                        GAsyncResult *result,
                        gpointer user_data)
{
  GdisGtk4Window *self;
  GtkFileDialog *dialog;
  GFile *file;
  GError *error;
  gchar *path;

  dialog = GTK_FILE_DIALOG(source_object);
  self = user_data;

  error = NULL;
  file = gtk_file_dialog_open_finish(dialog, result, &error);
  if (!file)
    {
      if (error)
        {
          if (!g_error_matches(error, GTK_DIALOG_ERROR, GTK_DIALOG_ERROR_CANCELLED) &&
              !g_error_matches(error, G_IO_ERROR, G_IO_ERROR_CANCELLED))
            gdis_gtk4_window_log(self, "Open failed: %s\n", error->message);
          g_error_free(error);
        }
      g_object_unref(dialog);
      return;
    }

  path = g_file_get_path(file);
  if (path)
    {
      gdis_gtk4_window_add_model_from_path(self, path, TRUE);
      g_free(path);
    }

  g_object_unref(file);
  g_object_unref(dialog);
}

static void
on_save_dialog_complete(GObject *source_object,
                        GAsyncResult *result,
                        gpointer user_data)
{
  GdisSaveDialogContext *context;
  GdisGtk4Window *self;
  GtkFileDialog *dialog;
  GFile *file;
  GError *error;
  gchar *path;

  dialog = GTK_FILE_DIALOG(source_object);
  context = user_data;
  self = context ? context->owner : NULL;
  file = NULL;
  error = NULL;
  path = NULL;

  file = gtk_file_dialog_save_finish(dialog, result, &error);
  if (!file)
    {
      if (self && error)
        {
          if (!g_error_matches(error, GTK_DIALOG_ERROR, GTK_DIALOG_ERROR_CANCELLED) &&
              !g_error_matches(error, G_IO_ERROR, G_IO_ERROR_CANCELLED))
            gdis_gtk4_window_log(self, "Save failed: %s\n", error->message);
          g_error_free(error);
        }
      g_free(context);
      g_object_unref(dialog);
      return;
    }

  path = g_file_get_path(file);
  if (!self || !context->model)
    goto cleanup;

  if (!path)
    {
      gdis_gtk4_window_log(self, "Save failed: GTK did not return a local filesystem path.\n");
      goto cleanup;
    }

  gdis_gtk4_window_save_model_to_path(self, context->model, path);

cleanup:
  g_free(path);
  g_object_unref(file);
  g_free(context);
  g_object_unref(dialog);
}

static void
present_open_dialog(GdisGtk4Window *self)
{
  GtkFileDialog *dialog;

  g_return_if_fail(self != NULL);

  dialog = gtk_file_dialog_new();
  gtk_file_dialog_set_title(dialog, "Open model");
  gtk_file_dialog_open(dialog,
                       GTK_WINDOW(self->window),
                       NULL,
                       on_open_dialog_complete,
                       self);
}

static void
present_save_dialog(GdisGtk4Window *self, GdisModel *model)
{
  GtkFileDialog *dialog;
  GdisSaveDialogContext *context;
  g_autofree gchar *parent_path = NULL;
  g_autofree gchar *initial_name = NULL;

  g_return_if_fail(self != NULL);
  g_return_if_fail(model != NULL);

  dialog = gtk_file_dialog_new();
  gtk_file_dialog_set_title(dialog, "Save model");
  gtk_file_dialog_set_accept_label(dialog, "Save");

  if (model->path && g_path_is_absolute(model->path))
    {
      GFile *initial_file;

      initial_file = g_file_new_for_path(model->path);
      gtk_file_dialog_set_initial_file(dialog, initial_file);
      g_object_unref(initial_file);
    }
  else
    {
      initial_name = g_path_get_basename(model->path ? model->path : "Untitled.xyz");
      gtk_file_dialog_set_initial_name(dialog, initial_name);
    }

  if (model->path && *model->path)
    {
      parent_path = g_path_get_dirname(model->path);
      if (parent_path && g_path_is_absolute(parent_path) && g_file_test(parent_path, G_FILE_TEST_IS_DIR))
        {
          GFile *folder;

          folder = g_file_new_for_path(parent_path);
          gtk_file_dialog_set_initial_folder(dialog, folder);
          g_object_unref(folder);
        }
    }

  context = g_new0(GdisSaveDialogContext, 1);
  context->owner = self;
  context->model = model;
  gtk_file_dialog_save(dialog,
                       GTK_WINDOW(self->window),
                       NULL,
                       on_save_dialog_complete,
                       context);
}

static GtkWidget *
build_active_model_card(GdisGtk4Window *self)
{
  GtkWidget *frame;
  GtkWidget *box;

  frame = gtk_frame_new(NULL);
  gtk_widget_add_css_class(frame, "gdis-sidebar-card");

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_widget_set_margin_start(box, 10);
  gtk_widget_set_margin_end(box, 10);
  gtk_widget_set_margin_top(box, 10);
  gtk_widget_set_margin_bottom(box, 10);

  gtk_box_append(GTK_BOX(box), new_section_button("Active Model"));
  gtk_box_append(GTK_BOX(box), new_readonly_text_view(&self->active_summary_buffer, FALSE, 120));
  gdis_text_buffer_set(self->active_summary_buffer,
                       "No active model\nOpen a file or choose a sample model.");

  gtk_frame_set_child(GTK_FRAME(frame), box);

  return frame;
}

static GtkWidget *
build_model_list(GdisGtk4Window *self)
{
  GtkWidget *list;
  int i;

  list = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
  gtk_widget_set_hexpand(list, TRUE);

  self->model_list = list;

  for (i = 0; gdis_legacy_model_samples[i] != NULL; i++)
    gdis_gtk4_window_add_model_from_path(self, gdis_legacy_model_samples[i], FALSE);

  if (!self->active_model && self->models->len > 0)
    gdis_gtk4_window_set_active_model(self, g_ptr_array_index(self->models, 0));

  return list;
}

static GtkWidget *
build_sidebar(GdisGtk4Window *self)
{
  GtkWidget *box;
  GtkWidget *active_card;
  GtkWidget *footer;
  GtkWidget *mode_grid;
  GtkWidget *scroller;
  GtkWidget *stack_frame;
  GtkWidget *stack_sidebar;
  GtkWidget *stack;
  int i;

  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12);
  gtk_widget_add_css_class(box, "gdis-sidebar");
  gtk_widget_set_size_request(box, 300, -1);
  gtk_widget_set_margin_start(box, 12);
  gtk_widget_set_margin_end(box, 12);
  gtk_widget_set_margin_top(box, 12);
  gtk_widget_set_margin_bottom(box, 12);

  gtk_box_append(GTK_BOX(box), new_section_button("Models"));
  active_card = build_active_model_card(self);
  gtk_box_append(GTK_BOX(box), active_card);

  gtk_box_append(GTK_BOX(box), new_section_button("Panels"));
  stack = build_sidebar_stack(self);
  gtk_widget_set_vexpand(stack, TRUE);
  stack_sidebar = gtk_stack_sidebar_new();
  gtk_stack_sidebar_set_stack(GTK_STACK_SIDEBAR(stack_sidebar), GTK_STACK(stack));
  gtk_widget_set_hexpand(stack_sidebar, FALSE);

  build_model_list(self);
  scroller = gtk_scrolled_window_new();
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                 GTK_POLICY_NEVER,
                                 GTK_POLICY_AUTOMATIC);
  gtk_scrolled_window_set_min_content_height(GTK_SCROLLED_WINDOW(scroller), 210);
  gtk_widget_set_size_request(scroller, -1, 220);
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), self->model_list);
  gtk_box_append(GTK_BOX(box), scroller);

  gtk_box_append(GTK_BOX(box), stack_sidebar);

  stack_frame = gtk_frame_new(NULL);
  gtk_widget_add_css_class(stack_frame, "gdis-sidebar-panel");
  gtk_widget_set_vexpand(stack_frame, TRUE);
  gtk_widget_add_css_class(stack, "gdis-sidebar-stack");
  gtk_frame_set_child(GTK_FRAME(stack_frame), stack);
  gtk_box_append(GTK_BOX(box), stack_frame);

  gtk_box_append(GTK_BOX(box), new_section_button("Selection Modes"));
  mode_grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(mode_grid), 6);
  gtk_grid_set_row_spacing(GTK_GRID(mode_grid), 6);
  for (i = 0; gdis_legacy_selection_modes[i] != NULL && i < 8; i++)
    {
      GtkWidget *button;

      button = gtk_button_new_with_label(gdis_legacy_selection_modes[i] + 9);
      gtk_widget_add_css_class(button, "flat");
      g_object_set_data_full(G_OBJECT(button), "selection-name",
                             g_strdup(gdis_legacy_selection_modes[i]), g_free);
      g_signal_connect(button, "clicked",
                       G_CALLBACK(on_selection_button_clicked), self);
      gtk_widget_set_tooltip_text(button, gdis_legacy_selection_modes[i]);
      gtk_widget_set_hexpand(button, TRUE);
      gdis_gtk4_window_configure_button_label(button);
      if (g_strcmp0(gdis_legacy_selection_modes[i], "Select : Regions") == 0)
        {
          gtk_widget_set_sensitive(button, FALSE);
          gtk_widget_set_tooltip_text(button,
                                      "Region selection needs region-labelled model data, which is not loaded by the current GTK4 bridge yet.");
        }
      gtk_grid_attach(GTK_GRID(mode_grid), button, i % 2, i / 2, 1, 1);
    }
  gtk_box_append(GTK_BOX(box), mode_grid);

  footer = gtk_button_new_with_label("Legacy reference: ../legacy_snapshot");
  gtk_widget_add_css_class(footer, "flat");
  gtk_widget_add_css_class(footer, "gdis-muted");
  gtk_widget_set_sensitive(footer, FALSE);
  gtk_box_append(GTK_BOX(box), footer);

  return box;
}

static GtkWidget *
build_viewer_area(GdisGtk4Window *self)
{
  GtkWidget *frame;
  GtkWidget *drawing_area;
  GtkGesture *drag;
  GtkGesture *click;
  GtkEventController *scroll;

  frame = gtk_frame_new(NULL);
  gtk_widget_add_css_class(frame, "gdis-viewer-frame");
  gtk_widget_set_hexpand(frame, TRUE);
  gtk_widget_set_vexpand(frame, TRUE);
  gtk_widget_set_size_request(frame, 820, 540);

  drawing_area = gtk_drawing_area_new();
  gtk_widget_set_hexpand(drawing_area, TRUE);
  gtk_widget_set_vexpand(drawing_area, TRUE);
  gtk_drawing_area_set_content_width(GTK_DRAWING_AREA(drawing_area), 820);
  gtk_drawing_area_set_content_height(GTK_DRAWING_AREA(drawing_area), 540);
  gtk_widget_set_cursor_from_name(drawing_area, "crosshair");
  self->viewer_area = drawing_area;
  gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawing_area), viewer_draw, self, NULL);

  drag = gtk_gesture_drag_new();
  gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag), GDK_BUTTON_PRIMARY);
  g_signal_connect(drag, "drag-begin", G_CALLBACK(on_viewer_drag_begin), self);
  g_signal_connect(drag, "drag-update", G_CALLBACK(on_viewer_drag_update), self);
  gtk_widget_add_controller(drawing_area, GTK_EVENT_CONTROLLER(drag));

  click = gtk_gesture_click_new();
  gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click), GDK_BUTTON_PRIMARY);
  g_signal_connect(click, "pressed", G_CALLBACK(on_viewer_click_pressed), self);
  g_signal_connect(click, "released", G_CALLBACK(on_viewer_click_released), self);
  gtk_widget_add_controller(drawing_area, GTK_EVENT_CONTROLLER(click));

  scroll = gtk_event_controller_scroll_new(GTK_EVENT_CONTROLLER_SCROLL_BOTH_AXES |
                                           GTK_EVENT_CONTROLLER_SCROLL_DISCRETE);
  g_signal_connect(scroll, "scroll", G_CALLBACK(on_viewer_scroll), self);
  gtk_widget_add_controller(drawing_area, scroll);

  gtk_frame_set_child(GTK_FRAME(frame), drawing_area);

  return frame;
}

static GtkWidget *
build_status_view(GdisGtk4Window *self)
{
  GtkWidget *scroller;
  GtkWidget *text_view;

  scroller = gtk_scrolled_window_new();
  gtk_widget_add_css_class(scroller, "gdis-status");
  gtk_widget_set_vexpand(scroller, TRUE);
  gtk_widget_set_hexpand(scroller, TRUE);
  gtk_widget_set_size_request(scroller, 760, 150);
  gtk_scrolled_window_set_min_content_height(GTK_SCROLLED_WINDOW(scroller), 150);

  text_view = gtk_text_view_new();
  gtk_text_view_set_editable(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_WORD_CHAR);
  gtk_text_view_set_monospace(GTK_TEXT_VIEW(text_view), TRUE);
  gtk_text_view_set_left_margin(GTK_TEXT_VIEW(text_view), 8);
  gtk_text_view_set_right_margin(GTK_TEXT_VIEW(text_view), 8);
  gtk_text_view_set_top_margin(GTK_TEXT_VIEW(text_view), 8);
  gtk_text_view_set_bottom_margin(GTK_TEXT_VIEW(text_view), 8);
  self->status_buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));

  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  return scroller;
}

static GtkWidget *
build_main_content(GdisGtk4Window *self)
{
  GtkWidget *main_paned;
  GtkWidget *right_paned;
  GtkWidget *viewer;
  GtkWidget *status;

  main_paned = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
  gtk_widget_set_hexpand(main_paned, TRUE);
  gtk_widget_set_vexpand(main_paned, TRUE);
  gtk_widget_set_size_request(main_paned, 1180, 760);
  gtk_paned_set_wide_handle(GTK_PANED(main_paned), TRUE);
  self->main_paned = main_paned;

  right_paned = gtk_paned_new(GTK_ORIENTATION_VERTICAL);
  gtk_widget_set_hexpand(right_paned, TRUE);
  gtk_widget_set_vexpand(right_paned, TRUE);
  gtk_widget_set_size_request(right_paned, 820, 720);
  gtk_paned_set_wide_handle(GTK_PANED(right_paned), TRUE);
  self->right_paned = right_paned;

  viewer = build_viewer_area(self);
  status = build_status_view(self);

  gtk_paned_set_start_child(GTK_PANED(main_paned), build_sidebar(self));
  gtk_paned_set_end_child(GTK_PANED(main_paned), right_paned);
  gtk_paned_set_resize_start_child(GTK_PANED(main_paned), FALSE);
  gtk_paned_set_shrink_start_child(GTK_PANED(main_paned), FALSE);
  gtk_paned_set_resize_end_child(GTK_PANED(main_paned), TRUE);
  gtk_paned_set_shrink_end_child(GTK_PANED(main_paned), FALSE);
  gtk_paned_set_position(GTK_PANED(main_paned), 320);

  gtk_paned_set_start_child(GTK_PANED(right_paned), viewer);
  gtk_paned_set_end_child(GTK_PANED(right_paned), status);
  gtk_paned_set_resize_start_child(GTK_PANED(right_paned), TRUE);
  gtk_paned_set_shrink_start_child(GTK_PANED(right_paned), FALSE);
  gtk_paned_set_resize_end_child(GTK_PANED(right_paned), FALSE);
  gtk_paned_set_shrink_end_child(GTK_PANED(right_paned), FALSE);
  gtk_paned_set_position(GTK_PANED(right_paned), 690);

  return main_paned;
}

static void
gdis_gtk4_window_restore_layout(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);

  if (self->main_paned)
    gtk_paned_set_position(GTK_PANED(self->main_paned), 320);

  if (self->right_paned)
    gtk_paned_set_position(GTK_PANED(self->right_paned), 690);
}

static gboolean
gdis_gtk4_window_restore_layout_tick(GtkWidget *widget,
                                     GdkFrameClock *frame_clock,
                                     gpointer user_data)
{
  GdisGtk4Window *self;
  gint width;
  gint height;

  (void) widget;
  (void) frame_clock;

  self = user_data;
  if (!self)
    return G_SOURCE_REMOVE;

  self->layout_restore_tick_id = 0;
  if (!GTK_IS_WINDOW(self->window))
    return G_SOURCE_REMOVE;

  width = gtk_widget_get_width(self->window);
  height = gtk_widget_get_height(self->window);
  if (width < 1000 || height < 700)
    gtk_window_set_default_size(GTK_WINDOW(self->window), 1480, 920);

  gdis_gtk4_window_restore_layout(self);

  return G_SOURCE_REMOVE;
}

static void
action_activate(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  const char *action_name;

  self = user_data;
  action_name = g_object_get_data(G_OBJECT(button), "action-name");
  if (!action_name)
    return;

  g_action_group_activate_action(G_ACTION_GROUP(self->app), action_name, NULL);
}

static GtkWidget *
build_toolbar(GdisGtk4Window *self)
{
  GtkWidget *toolbar;
  int i;

  toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_widget_add_css_class(toolbar, "gdis-toolbar");
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

static GMenuModel *
build_menu_bar_model(void)
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
  g_menu_append(file_menu, "Save As", "app.save-as");
  g_menu_append(file_menu, "Close", "app.close-model");
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

static void
action_dispatch(GSimpleAction *action, GVariant *parameter, gpointer user_data)
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
      gdis_gtk4_window_log(self, "Quit requested from GTK4 rebuild.\n");
      g_application_quit(G_APPLICATION(app));
      return;
    }

  if (g_strcmp0(name, "open") == 0)
    {
      present_open_dialog(self);
      return;
    }

  if (g_strcmp0(name, "save") == 0)
    {
      if (!self->active_model)
        {
          gdis_gtk4_window_log(self, "Save requested, but no active model is loaded.\n");
          return;
        }

      if (self->active_model->path &&
          g_file_test(self->active_model->path, G_FILE_TEST_EXISTS))
        {
          gdis_gtk4_window_save_model_to_path(self,
                                              self->active_model,
                                              self->active_model->path);
        }
      else
        {
          present_save_dialog(self, self->active_model);
        }
      return;
    }

  if (g_strcmp0(name, "save-as") == 0)
    {
      if (!self->active_model)
        {
          gdis_gtk4_window_log(self, "Save As requested, but no active model is loaded.\n");
          return;
        }

      present_save_dialog(self, self->active_model);
      return;
    }

  if (g_strcmp0(name, "close-model") == 0)
    {
      if (!self->active_model)
        {
          gdis_gtk4_window_log(self, "Close requested, but no active model is loaded.\n");
          return;
        }

      gdis_gtk4_window_log(self, "Closing model: %s\n", self->active_model->path);
      gdis_gtk4_window_remove_model(self, self->active_model);
      return;
    }

  if (g_strcmp0(name, "reset-view") == 0)
    {
      gdis_gtk4_window_reset_view(self);
      gdis_gtk4_window_update_details(self);
      gdis_gtk4_window_refresh_viewer(self);
      gdis_gtk4_window_log(self, "View reset.\n");
      return;
    }

  if (g_strcmp0(name, "about") == 0)
    {
      gdis_gtk4_window_log(self,
                           "GTK4 rebuild for GDIS.\n"
                           "Legacy-derived loaders now support XYZ, PDB, ARC/CAR, and CIF metadata.\n");
      return;
    }

  if (g_strcmp0(name, "measure") == 0)
    {
      gdis_gtk4_window_present_measure_tool(self);
      gdis_gtk4_window_log(self, "Opened measurement tool.\n");
      return;
    }

  if (g_strcmp0(name, "surface") == 0)
    {
      g_autofree char *report = NULL;

      report = gdis_report_surface(self->active_model);
      gdis_gtk4_window_present_report(self, "Surface Explorer", report);
      gdis_gtk4_window_log(self, "Opened surface explorer.\n");
      return;
    }

  if (g_strcmp0(name, "diffraction") == 0)
    {
      g_autofree gchar *availability_message = NULL;

      if (!gdis_gtk4_window_model_supports_diffraction(self->active_model, &availability_message))
        {
          gdis_gtk4_window_log(self, "%s\n",
                               availability_message ? availability_message :
                               "Powder diffraction is unavailable for the current model.");
          gdis_gtk4_window_present_report(self,
                                          "Powder Diffraction",
                                          availability_message ? availability_message :
                                          "Powder diffraction is unavailable for the current model.");
          return;
        }

      gdis_gtk4_window_present_diffraction_tool(self);
      gdis_gtk4_window_log(self, "Opened powder diffraction tool.\n");
      return;
    }

  if (g_strcmp0(name, "isosurface") == 0)
    {
      g_autofree char *report = NULL;

      report = gdis_report_isosurface(self->active_model);
      gdis_gtk4_window_present_report(self, "Iso-surface Status", report);
      gdis_gtk4_window_log(self, "Opened iso-surface status report.\n");
      return;
    }

  if (g_strcmp0(name, "render") == 0)
    {
      g_autofree char *report = NULL;

      report = g_strdup_printf("Render Controls\n\n"
                               "Viewer toggles:\n"
                               "  Atoms: %s\n"
                               "  Bonds: %s\n"
                               "  Cell: %s\n"
                               "  Labels: %s\n\n"
                               "Controls:\n"
                               "  Drag: rotate\n"
                               "  Scroll: zoom\n"
                               "  Click: pick atom\n"
                               "  View X / Y / Z: preset orientations\n\n"
                               "This GTK4 rebuild currently uses a lightweight native renderer.\n"
                               "The full legacy display/render property dialog is a separate port.",
                               gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_atoms_toggle)) ? "on" : "off",
                               gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_bonds_toggle)) ? "on" : "off",
                               gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_cell_toggle)) ? "on" : "off",
                               gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_labels_toggle)) ? "on" : "off");
      gdis_gtk4_window_present_report(self, "Render Controls", report);
      gdis_gtk4_window_log(self, "Opened render control summary.\n");
      return;
    }

  if (g_strcmp0(name, "new-model") == 0 ||
      g_strcmp0(name, "edit") == 0)
    {
      if (g_strcmp0(name, "edit") == 0)
        {
          gdis_gtk4_window_present_edit_tool(self);
          gdis_gtk4_window_log(self, "Opened model editor.\n");
        }
      else
        {
          if (!gdis_gtk4_window_add_new_empty_model(self))
            {
              gdis_gtk4_window_log(self, "Could not create a new empty model.\n");
              return;
            }

          gdis_gtk4_window_present_edit_tool(self);
          gdis_gtk4_window_log(self, "Opened model editor for the new empty model.\n");
        }
      return;
    }

  gdis_gtk4_window_log(self, "Activated action: %s\n", name);
}

static void
install_actions(GtkApplication *app)
{
  const GActionEntry entries[] = {
    {.name = "open", .activate = action_dispatch},
    {.name = "save", .activate = action_dispatch},
    {.name = "save-as", .activate = action_dispatch},
    {.name = "close-model", .activate = action_dispatch},
    {.name = "quit", .activate = action_dispatch},
    {.name = "new-model", .activate = action_dispatch},
    {.name = "edit", .activate = action_dispatch},
    {.name = "render", .activate = action_dispatch},
    {.name = "measure", .activate = action_dispatch},
    {.name = "isosurface", .activate = action_dispatch},
    {.name = "surface", .activate = action_dispatch},
    {.name = "diffraction", .activate = action_dispatch},
    {.name = "reset-view", .activate = action_dispatch},
    {.name = "about", .activate = action_dispatch}
  };

  g_action_map_add_action_entries(G_ACTION_MAP(app),
                                  entries,
                                  G_N_ELEMENTS(entries),
                                  app);
}

GdisGtk4Window *
gdis_gtk4_window_new(GtkApplication *app)
{
  GdisGtk4Window *self;
  GtkWidget *root;
  GtkWidget *menu_bar;
  GMenuModel *menu_model;

  g_return_val_if_fail(GTK_IS_APPLICATION(app), NULL);

  self = g_new0(GdisGtk4Window, 1);
  self->app = app;
  self->models = g_ptr_array_new_with_free_func((GDestroyNotify) gdis_model_free);
  self->selected_atoms = g_array_new(FALSE, FALSE, sizeof(guint));
  self->picked_atoms = g_array_new(FALSE, FALSE, sizeof(guint));
  self->selected_atom_index = INVALID_ATOM_INDEX;
  self->selection_mode = GDIS_SELECTION_MODE_ATOMS;
  self->click_mode = GDIS_CLICK_MODE_SELECT;
  self->fragment_anchor_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_reset_view(self);

  install_actions(app);

  self->window = gtk_application_window_new(app);
  gtk_window_set_title(GTK_WINDOW(self->window), "GDIS GTK4");
  gtk_window_set_default_size(GTK_WINDOW(self->window), 1480, 920);
  gdis_gtk4_window_install_css(self->window);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_widget_add_css_class(root, "gdis-root");
  gtk_widget_set_size_request(root, 1180, 760);
  gtk_window_set_child(GTK_WINDOW(self->window), root);

  menu_model = build_menu_bar_model();
  menu_bar = gtk_popover_menu_bar_new_from_model(menu_model);
  g_object_unref(menu_model);

  gtk_box_append(GTK_BOX(root), menu_bar);
  gtk_box_append(GTK_BOX(root), build_toolbar(self));
  gtk_box_append(GTK_BOX(root), gtk_separator_new(GTK_ORIENTATION_HORIZONTAL));
  gtk_box_append(GTK_BOX(root), build_main_content(self));
  gdis_gtk4_window_restore_layout(self);

  g_object_set_data_full(G_OBJECT(self->window), WINDOW_DATA_KEY, self, gdis_gtk4_window_free);

  gdis_gtk4_window_log(self, "GTK4 rebuild initialized.\n");
  gdis_gtk4_window_log(self, "Legacy-derived model loaders are active for XYZ, PDB, ARC/CAR, and CIF.\n");
  gdis_gtk4_window_log(self, "Renderer supports rotate, zoom, model switching, and atom picking.\n");

  return self;
}

void
gdis_gtk4_window_present(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);
  g_return_if_fail(GTK_IS_WINDOW(self->window));

  gtk_window_set_default_size(GTK_WINDOW(self->window), 1480, 920);
  gdis_gtk4_window_restore_layout(self);
  gtk_window_present(GTK_WINDOW(self->window));
  if (self->layout_restore_tick_id == 0)
    {
      self->layout_restore_tick_id =
        gtk_widget_add_tick_callback(self->window,
                                     gdis_gtk4_window_restore_layout_tick,
                                     self,
                                     NULL);
    }
}

GtkWindow *
gdis_gtk4_window_get_window(GdisGtk4Window *self)
{
  g_return_val_if_fail(self != NULL, NULL);
  g_return_val_if_fail(GTK_IS_WINDOW(self->window), NULL);

  return GTK_WINDOW(self->window);
}

void
gdis_gtk4_window_open_startup_paths(GdisGtk4Window *self,
                                    const char *const *paths)
{
  int i;

  g_return_if_fail(self != NULL);

  if (!paths)
    return;

  for (i = 0; paths[i] != NULL; i++)
    gdis_gtk4_window_add_model_from_path(self, paths[i], TRUE);
}
