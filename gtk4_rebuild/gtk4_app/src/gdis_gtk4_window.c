#include "gdis_gtk4_window.h"

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <glib/gstdio.h>

#include "gdis_elements.h"
#include "gdis_legacy_map.h"
#include "gdis_isosurface.h"
#include "gdis_macos_menu.h"
#include "gdis_model.h"
#include "gdis_restoration.h"
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
  gdouble screen[3][2];
  gdouble depth;
  gdouble shade;
} GdisProjectedTriangle;

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
  GtkWidget *capture_toggle;
  GdisMeasureMode mode;
} GdisMeasureTool;

typedef struct _GdisEditTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkTextBuffer *buffer;
  GtkWidget *label_entry;
  GtkWidget *element_entry;
  GtkWidget *ff_type_entry;
  GtkWidget *region_entry;
  GtkWidget *x_entry;
  GtkWidget *y_entry;
  GtkWidget *z_entry;
  GtkWidget *bond_order_entry;
} GdisEditTool;

typedef struct _GdisDisplayTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkWidget *atoms_toggle;
  GtkWidget *bonds_toggle;
  GtkWidget *cell_toggle;
  GtkWidget *labels_toggle;
} GdisDisplayTool;

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

typedef struct _GdisSurfaceTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkWidget *h_spin;
  GtkWidget *k_spin;
  GtkWidget *l_spin;
  GtkWidget *shift_spin;
  GtkWidget *region_a_spin;
  GtkWidget *region_b_spin;
  GtkWidget *repeat_a_spin;
  GtkWidget *repeat_b_spin;
  GtkWidget *vacuum_spin;
  GtkWidget *execute_button;
  GtkTextBuffer *report_buffer;
} GdisSurfaceTool;

typedef struct _GdisIsosurfaceTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkWidget *mode_dropdown;
  GtkWidget *grid_spin;
  GtkWidget *blur_spin;
  GtkWidget *value_spin;
  GtkWidget *execute_button;
  GtkWidget *clear_button;
  GtkTextBuffer *report_buffer;
} GdisIsosurfaceTool;

typedef struct _GdisAnimationTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkWidget *play_button;
  GtkWidget *frame_scale;
  GtkWidget *delay_spin;
  GtkWidget *step_spin;
  GtkWidget *loop_toggle;
  GtkWidget *preserve_connectivity_toggle;
  GtkWidget *preserve_scale_toggle;
  GtkWidget *confine_none_toggle;
  GtkWidget *confine_atoms_toggle;
  GtkWidget *confine_molecules_toggle;
  GtkLabel *summary_label;
  GtkLabel *processing_label;
  GtkLabel *rendering_label;
  guint timer_id;
  gboolean suppress_scale_signal;
} GdisAnimationTool;

typedef struct _GdisRecordingTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkWidget *output_entry;
  GtkWidget *prefix_entry;
  GtkWidget *width_spin;
  GtkWidget *height_spin;
  GtkWidget *movie_enable_toggle;
  GtkWidget *movie_mp4_toggle;
  GtkWidget *movie_gif_toggle;
  GtkWidget *fps_spin;
  GtkWidget *keep_frames_toggle;
  GtkWidget *movie_button;
  GtkLabel *summary_label;
} GdisRecordingTool;

typedef struct _GdisZmatrixTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkWidget *scope_dropdown;
  GtkWidget *row_spin;
  GtkWidget *distance_entry;
  GtkWidget *angle_entry;
  GtkWidget *torsion_entry;
  GtkLabel *row_label;
  GtkTextBuffer *buffer;
  GArray *scope;
  GArray *rows;
  gboolean suppress_row_signal;
} GdisZmatrixTool;

typedef struct _GdisDislocationTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkWidget *type_dropdown;
  GtkWidget *line_entries[3];
  GtkWidget *burgers_entries[3];
  GtkWidget *origin_entries[3];
  GtkWidget *radius_spin;
  GtkWidget *poisson_spin;
  GtkWidget *selection_toggle;
  GtkTextBuffer *buffer;
} GdisDislocationTool;

typedef struct _GdisDockingTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkWidget *output_entry;
  GtkWidget *name_entry;
  GtkWidget *grid_spins[3];
  GtkWidget *rotation_spins[3];
  GtkWidget *span_spins[3];
  GtkWidget *preview_spin;
  GtkTextBuffer *buffer;
} GdisDockingTool;

typedef struct _GdisQboxSpeciesRow
{
  gchar *element;
  GtkWidget *alias_entry;
  GtkWidget *pseudo_entry;
} GdisQboxSpeciesRow;

typedef struct _GdisQboxTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkWidget *job_entry;
  GtkWidget *workdir_entry;
  GtkWidget *input_entry;
  GtkWidget *output_entry;
  GtkWidget *save_entry;
  GtkWidget *restart_entry;
  GtkWidget *exec_entry;
  GtkWidget *launcher_entry;
  GtkWidget *pseudo_dir_entry;
  GtkWidget *xc_entry;
  GtkWidget *wf_dyn_entry;
  GtkWidget *atoms_dyn_entry;
  GtkWidget *scf_tol_entry;
  GtkWidget *ecut_spin;
  GtkWidget *charge_spin;
  GtkWidget *padding_spin;
  GtkWidget *kmesh_spins[3];
  GtkWidget *ionic_steps_spin;
  GtkWidget *scf_steps_spin;
  GtkWidget *density_update_spin;
  GtkWidget *use_cell_toggle;
  GtkWidget *atomic_density_toggle;
  GtkWidget *randomize_toggle;
  GtkWidget *species_scroller;
  GtkLabel *summary_label;
  GtkTextBuffer *deck_buffer;
  GtkTextBuffer *report_buffer;
  GPtrArray *species_rows;
  GdisModel *source_model;
  gchar *last_generated_input;
  gboolean suppress_input_signal;
  gboolean editor_dirty;
} GdisQboxTool;

typedef struct _GdisPeriodicTableTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkLabel *title_label;
  GtkLabel *info_label;
  GtkWidget *copy_button;
  GtkWidget *apply_button;
  guint selected_atomic_number;
} GdisPeriodicTableTool;

typedef struct _GdisTaskManagerTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkTextBuffer *buffer;
} GdisTaskManagerTool;

typedef struct _GdisExecutablePathsTool
{
  struct _GdisGtk4Window *owner;
  GtkWidget *window;
  GtkTextBuffer *preview_buffer;
  GtkWidget **entries;
  guint entry_count;
  gchar *focus_backend;
} GdisExecutablePathsTool;

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

typedef enum
{
  GDIS_DRAG_MODE_NONE = 0,
  GDIS_DRAG_MODE_ROTATE,
  GDIS_DRAG_MODE_BOX_SELECT
} GdisDragMode;

typedef enum
{
  GDIS_ANIMATION_SOURCE_NONE = 0,
  GDIS_ANIMATION_SOURCE_MODEL_FRAMES,
  GDIS_ANIMATION_SOURCE_SESSION_MODELS
} GdisAnimationSourceType;

typedef enum
{
  GDIS_ANIMATION_CONFINE_NONE = 0,
  GDIS_ANIMATION_CONFINE_ATOMS,
  GDIS_ANIMATION_CONFINE_MOLECULES
} GdisAnimationConfineMode;

typedef struct
{
  GdisMeasureMode mode;
  guint atom_count;
  guint atom_indices[4];
} GdisMeasurementRecord;

typedef struct
{
  GdisModel *model_snapshot;
  GArray *selected_atoms;
  GArray *picked_atoms;
  GPtrArray *measurement_records;
  guint selected_atom_index;
  guint fragment_anchor_index;
} GdisUndoEntry;

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
  GHashTable *measurement_records;
  GHashTable *undo_stacks;
  GHashTable *iso_surfaces;
  GdisDisplayTool *display_tool;
  GdisMeasureTool *measure_tool;
  GdisEditTool *edit_tool;
  GdisDiffractionTool *diffraction_tool;
  GdisSurfaceTool *surface_tool;
  GdisIsosurfaceTool *isosurface_tool;
  GdisAnimationTool *animation_tool;
  GdisRecordingTool *recording_tool;
  GdisZmatrixTool *zmatrix_tool;
  GdisDislocationTool *dislocation_tool;
  GdisDockingTool *docking_tool;
  GdisQboxTool *qbox_tool;
  GdisPeriodicTableTool *periodic_table_tool;
  GdisTaskManagerTool *task_manager_tool;
  GdisExecutablePathsTool *exec_paths_tool;
  GdisSelectionMode selection_mode;
  GdisClickMode click_mode;
  GHashTable *executable_paths;
  guint fragment_anchor_index;
  guint untitled_counter;
  guint layout_restore_tick_id;
  GdkModifierType press_modifiers;
  GdisDragMode drag_mode;
  gdouble rotation_x;
  gdouble rotation_y;
  gdouble drag_origin_x;
  gdouble drag_origin_y;
  gdouble press_x;
  gdouble press_y;
  gdouble drag_current_x;
  gdouble drag_current_y;
  gdouble zoom;
  gchar *qbox_last_summary;
  gchar *qbox_last_workdir;
  gchar *qbox_last_input_path;
  gchar *qbox_last_output_path;
  gchar *qbox_last_stderr_path;
  gchar *qbox_last_save_path;
};

static const char *const WINDOW_DATA_KEY = "gdis-gtk4-window";
static const guint INVALID_ATOM_INDEX = G_MAXUINT;
static const gdouble GDIS_ANGSTROM_TO_BOHR = 1.8897261254578281;
static const char *const GDIS_QBOX_DEFAULT_PSEUDO_SOURCE_URL = "http://quantum-simulation.org/potentials/sg15_oncv/xml";
static const gint GDIS_MAIN_WINDOW_WIDTH = 1360;
static const gint GDIS_MAIN_WINDOW_HEIGHT = 820;
static const gint GDIS_MAIN_ROOT_WIDTH = 1120;
static const gint GDIS_MAIN_ROOT_HEIGHT = 680;
static const gint GDIS_VIEWER_MIN_WIDTH = 800;
static const gint GDIS_VIEWER_MIN_HEIGHT = 500;
static const gint GDIS_STATUS_MIN_HEIGHT = 130;
static const gint GDIS_MAIN_SIDEBAR_WIDTH = 320;
static const gint GDIS_RIGHT_PANED_POSITION = 560;
static gboolean qbox_startup_write_consumed = FALSE;
static gboolean qbox_startup_run_consumed = FALSE;

typedef struct _GdisExecutableSpec GdisExecutableSpec;

static GtkWidget *gdis_gtk4_window_find_model_button(GdisGtk4Window *self,
                                                     const char *path);
static GdisModel *gdis_gtk4_window_find_model(GdisGtk4Window *self,
                                              const char *path);
static gint gdis_gtk4_window_compare_models(const GdisModel *left,
                                            const GdisModel *right);
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
static void gdis_gtk4_window_select_all_atoms(GdisGtk4Window *self);
static void gdis_gtk4_window_invert_selected_atoms(GdisGtk4Window *self);
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
static void gdis_gtk4_window_refresh_surface_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_isosurface_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_animation_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_recording_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_zmatrix_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_qbox_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_periodic_table_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_task_manager_tool(GdisGtk4Window *self);
static void gdis_qbox_clear_last_run_state(GdisGtk4Window *self);
static void gdis_qbox_set_last_run_state(GdisGtk4Window *self,
                                         const char *workdir,
                                         const char *input_path,
                                         const char *output_path,
                                         const char *stderr_path,
                                         const char *save_path);
static gboolean gdis_qbox_import_result_into_active_model(GdisGtk4Window *self,
                                                          const char *result_path,
                                                          GError **error);
static gchar *gdis_qbox_build_results_report(GdisGtk4Window *self);
static gchar *gdis_qbox_make_unique_continue_name(const char *workdir,
                                                  const char *base_stem,
                                                  const char *extension);
static gchar *gdis_qbox_extract_last_etotal(const char *text);
static gchar *gdis_qbox_tail_text(const char *text, gsize max_chars);
static void on_qbox_use_last_xml_clicked(GtkButton *button, gpointer user_data);
static void on_qbox_continue_clicked(GtkButton *button, gpointer user_data);
static void on_qbox_import_result_clicked(GtkButton *button, gpointer user_data);
static void on_qbox_results_clicked(GtkButton *button, gpointer user_data);
static void on_qbox_apply_mp_mesh_clicked(GtkButton *button, gpointer user_data);
static void on_qbox_append_hw5_phonons_clicked(GtkButton *button, gpointer user_data);
static void on_qbox_append_hw5_homo_lumo_clicked(GtkButton *button, gpointer user_data);
static void gdis_gtk4_window_refresh_executable_paths_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_sync_display_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_refresh_after_model_edit(GdisGtk4Window *self,
                                                      gboolean refresh_model_buttons);
static void gdis_gtk4_window_restore_layout(GdisGtk4Window *self);
static void gdis_gtk4_window_present_display_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_measure_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_edit_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_diffraction_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_surface_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_isosurface_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_animation_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_recording_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_zmatrix_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_dislocation_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_docking_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_qbox_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_periodic_table_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_task_manager_tool(GdisGtk4Window *self);
static void gdis_gtk4_window_present_executable_paths_tool(GdisGtk4Window *self,
                                                           const char *focus_backend);
static const GdisExecutableSpec *gdis_executable_spec_for_id(const char *id);
static gchar *gdis_executable_detect_path(const GdisExecutableSpec *spec);
static void gdis_gtk4_window_present_report(GdisGtk4Window *self,
                                            const char *title,
                                            const char *body);
static void gdis_gtk4_window_present_placeholder_feature(GdisGtk4Window *self,
                                                         const char *title,
                                                         const char *body);
static GdisIsoSurface *gdis_gtk4_window_get_isosurface(GdisGtk4Window *self,
                                                       GdisModel *model,
                                                       gboolean create_entry);
static void gdis_gtk4_window_set_isosurface(GdisGtk4Window *self,
                                            GdisModel *model,
                                            GdisIsoSurface *surface);
static void gdis_gtk4_window_clear_isosurface(GdisGtk4Window *self,
                                              GdisModel *model,
                                              const char *reason);
static void gdis_draw_isosurface(cairo_t *cr,
                                 GdisGtk4Window *self,
                                 const GdisIsoSurface *surface,
                                 const gdouble center[3],
                                 gdouble scale,
                                 int width,
                                 int height);
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
static gboolean gdis_gtk4_window_parse_region_entry(GtkWidget *entry,
                                                    gint *value_out,
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
static gint gdis_gtk4_window_find_atom_index_in_array(const GArray *array,
                                                      guint atom_index);
static GArray *gdis_gtk4_window_copy_atom_array(const GArray *array);
static void gdis_gtk4_window_toggle_atom_in_array(GArray *array, guint atom_index);
static GPtrArray *gdis_measurement_record_array_clone(const GPtrArray *records);
static GPtrArray *gdis_gtk4_window_get_undo_stack(GdisGtk4Window *self,
                                                  GdisModel *model,
                                                  gboolean create);
static void gdis_gtk4_window_update_undo_action(GdisGtk4Window *self);
static gboolean gdis_gtk4_window_push_undo_snapshot(GdisGtk4Window *self,
                                                    const char *reason);
static void gdis_gtk4_window_discard_undo_snapshot(GdisGtk4Window *self);
static gboolean gdis_gtk4_window_perform_undo(GdisGtk4Window *self);
static gint gdis_gtk4_window_get_model_index(GdisGtk4Window *self,
                                             GdisModel *model);
static guint gdis_gtk4_window_hit_test_atom_at(GdisGtk4Window *self,
                                               gdouble x,
                                               gdouble y,
                                               GdisProjectedAtom *projected_out);
static gboolean gdis_gtk4_window_apply_box_selection(GdisGtk4Window *self,
                                                     gdouble x0,
                                                     gdouble y0,
                                                     gdouble x1,
                                                     gdouble y1,
                                                     gboolean toggle_existing);
static GdisAnimationSourceType gdis_gtk4_window_get_animation_source(GdisGtk4Window *self,
                                                                     guint *count_out,
                                                                     gint *active_index_out);
static gboolean gdis_gtk4_window_should_capture_measure_picks(const GdisGtk4Window *self);
static gboolean gdis_gtk4_window_should_record_viewer_picks(const GdisGtk4Window *self);
static gboolean gdis_gtk4_window_set_picks_from_array(GdisGtk4Window *self,
                                                      const GArray *atoms,
                                                      guint limit,
                                                      const char *log_message);
static GPtrArray *gdis_gtk4_window_get_measurement_records(GdisGtk4Window *self,
                                                           GdisModel *model,
                                                           gboolean create);
static void gdis_gtk4_window_clear_saved_measurements(GdisGtk4Window *self,
                                                      GdisModel *model,
                                                      const char *reason);
static char *gdis_gtk4_window_build_measurement_tool_report(GdisGtk4Window *self);
static const char *gdis_gtk4_window_selection_mode_label(GdisSelectionMode mode);
static const char *gdis_gtk4_window_click_mode_label(GdisClickMode mode);
static const GdisAtom *gdis_gtk4_window_get_selected_atom(const GdisGtk4Window *self);
static void gdis_gtk4_window_clear_status_log(GdisGtk4Window *self);
static void on_model_button_clicked(GtkButton *button, gpointer user_data);
static void on_selection_button_clicked(GtkButton *button, gpointer user_data);
static void on_view_toggle_toggled(GtkToggleButton *button, gpointer user_data);
static void on_display_tool_toggle_toggled(GtkToggleButton *button, gpointer user_data);
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
static void on_viewer_secondary_drag_begin(GtkGestureDrag *gesture,
                                           gdouble start_x,
                                           gdouble start_y,
                                           gpointer user_data);
static void on_viewer_drag_update(GtkGestureDrag *gesture,
                                  gdouble offset_x,
                                  gdouble offset_y,
                                  gpointer user_data);
static void on_viewer_secondary_drag_update(GtkGestureDrag *gesture,
                                            gdouble offset_x,
                                            gdouble offset_y,
                                            gpointer user_data);
static void on_viewer_drag_end(GtkGestureDrag *gesture,
                               gdouble offset_x,
                               gdouble offset_y,
                               gpointer user_data);
static void on_viewer_secondary_drag_end(GtkGestureDrag *gesture,
                                         gdouble offset_x,
                                         gdouble offset_y,
                                         gpointer user_data);
static gboolean on_viewer_scroll(GtkEventControllerScroll *controller,
                                 gdouble dx,
                                 gdouble dy,
                                 gpointer user_data);
static gboolean gdis_gtk4_window_export_view_png(GdisGtk4Window *self,
                                                 const char *path,
                                                 gint width,
                                                 gint height,
                                                 GError **error);
static gboolean gdis_gtk4_window_animation_preserve_connectivity(const GdisGtk4Window *self);
static gboolean gdis_gtk4_window_animation_preserve_scale(const GdisGtk4Window *self);
static GdisAnimationConfineMode gdis_gtk4_window_get_animation_confine_mode(const GdisGtk4Window *self);
static void gdis_gtk4_window_apply_animation_processing(GdisGtk4Window *self);
static gboolean gdis_gtk4_window_export_animation_sequence(GdisGtk4Window *self,
                                                           const char *output_dir,
                                                           const char *prefix,
                                                           gint width,
                                                           gint height,
                                                           guint *frame_count_out,
                                                           GError **error);
static gboolean gdis_gtk4_window_export_movie(GdisGtk4Window *self,
                                              const char *output_dir,
                                              const char *prefix,
                                              gint width,
                                              gint height,
                                              guint fps,
                                              guint format_index,
                                              gboolean keep_frames,
                                              gchar **movie_path_out,
                                              guint *frame_count_out,
                                              GError **error);
static void gdis_zmatrix_tool_clear_working_set(GdisZmatrixTool *tool);
static void gdis_gtk4_window_refresh_zmatrix_row_controls(GdisGtk4Window *self);
static gboolean gdis_gtk4_window_commit_zmatrix_row_edits(GdisGtk4Window *self,
                                                          GError **error);
static void on_animation_processing_changed(GtkWidget *button, gpointer user_data);
static void on_recording_settings_changed(GtkWidget *widget, gpointer user_data);
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
static GtkWidget *new_section_button(const char *text);
static GtkWidget *new_toggle_button(const char *label, gboolean active);
static void present_save_dialog(GdisGtk4Window *self, GdisModel *model);
static void on_save_dialog_complete(GObject *source_object,
                                    GAsyncResult *result,
                                    gpointer user_data);
static void on_measure_capture_toggled(GtkToggleButton *button, gpointer user_data);
static void on_measure_use_selection_clicked(GtkButton *button, gpointer user_data);
static void on_measure_add_current_clicked(GtkButton *button, gpointer user_data);
static void on_measure_delete_last_clicked(GtkButton *button, gpointer user_data);
static void on_measure_clear_saved_clicked(GtkButton *button, gpointer user_data);
static void on_edit_apply_region_clicked(GtkButton *button, gpointer user_data);
static void on_use_selection_as_picks_clicked(GtkButton *button, gpointer user_data);
static void on_clear_selection_clicked(GtkButton *button, gpointer user_data);

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
  if (self->measurement_records)
    g_hash_table_unref(self->measurement_records);
  if (self->undo_stacks)
    g_hash_table_unref(self->undo_stacks);
  if (self->iso_surfaces)
    g_hash_table_unref(self->iso_surfaces);
  if (self->executable_paths)
    g_hash_table_unref(self->executable_paths);
  if (self->measure_tool)
    self->measure_tool->owner = NULL;
  if (self->edit_tool)
    self->edit_tool->owner = NULL;
  if (self->display_tool)
    self->display_tool->owner = NULL;
  if (self->diffraction_tool)
    self->diffraction_tool->owner = NULL;
  if (self->surface_tool)
    self->surface_tool->owner = NULL;
  if (self->isosurface_tool)
    self->isosurface_tool->owner = NULL;
  if (self->animation_tool)
    self->animation_tool->owner = NULL;
  if (self->recording_tool)
    self->recording_tool->owner = NULL;
  if (self->zmatrix_tool)
    self->zmatrix_tool->owner = NULL;
  if (self->dislocation_tool)
    self->dislocation_tool->owner = NULL;
  if (self->docking_tool)
    self->docking_tool->owner = NULL;
  if (self->qbox_tool)
    self->qbox_tool->owner = NULL;
  if (self->periodic_table_tool)
    self->periodic_table_tool->owner = NULL;
  if (self->task_manager_tool)
    self->task_manager_tool->owner = NULL;
  if (self->exec_paths_tool)
    self->exec_paths_tool->owner = NULL;
  g_clear_pointer(&self->qbox_last_summary, g_free);
  g_clear_pointer(&self->qbox_last_workdir, g_free);
  g_clear_pointer(&self->qbox_last_input_path, g_free);
  g_clear_pointer(&self->qbox_last_output_path, g_free);
  g_clear_pointer(&self->qbox_last_stderr_path, g_free);
  g_clear_pointer(&self->qbox_last_save_path, g_free);
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
  g_return_if_fail(format != NULL);

  va_start(args, format);
  message = g_strdup_vprintf(format, args);
  va_end(args);

  if (self->status_buffer != NULL)
    {
      gtk_text_buffer_get_end_iter(self->status_buffer, &end);
      gtk_text_buffer_insert(self->status_buffer, &end, message, -1);
    }
  else
    {
      g_printerr("%s", message);
    }

  g_free(message);
}

static void
gdis_gtk4_window_clear_status_log(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);

  if (self->status_buffer == NULL)
    return;

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
  if (self->active_model)
    gdis_gtk4_window_clear_isosurface(self,
                                      self->active_model,
                                      "Saved iso-surface preview was cleared because the model geometry changed.\n");
  gdis_gtk4_window_sync_display_tool(self);
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_refresh_diffraction_tool(self);
  gdis_gtk4_window_refresh_surface_tool(self);
  gdis_gtk4_window_refresh_isosurface_tool(self);
  gdis_gtk4_window_refresh_animation_tool(self);
  gdis_gtk4_window_refresh_recording_tool(self);
  gdis_gtk4_window_refresh_zmatrix_tool(self);
  gdis_gtk4_window_refresh_qbox_tool(self);
  gdis_gtk4_window_refresh_periodic_table_tool(self);
  gdis_gtk4_window_refresh_task_manager_tool(self);
  gdis_gtk4_window_update_undo_action(self);
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

static void
gdis_gtk4_window_present_placeholder_feature(GdisGtk4Window *self,
                                             const char *title,
                                             const char *body)
{
  g_autofree gchar *report = NULL;

  g_return_if_fail(self != NULL);
  g_return_if_fail(title != NULL);
  g_return_if_fail(body != NULL);

  report = g_strdup_printf(
    "%s\n\n"
    "%s\n\n"
    "The menu position now matches the legacy GDIS layout, but this module itself is not restored yet in the GTK4 rebuild.",
    title,
    body);
  gdis_gtk4_window_present_report(self, title, report);
}

static GdisIsoSurface *
gdis_gtk4_window_get_isosurface(GdisGtk4Window *self,
                                GdisModel *model,
                                gboolean create_entry)
{
  (void) create_entry;

  g_return_val_if_fail(self != NULL, NULL);

  if (!self->iso_surfaces || !model)
    return NULL;

  return g_hash_table_lookup(self->iso_surfaces, model);
}

static void
gdis_gtk4_window_set_isosurface(GdisGtk4Window *self,
                                GdisModel *model,
                                GdisIsoSurface *surface)
{
  g_return_if_fail(self != NULL);

  if (!self->iso_surfaces || !model)
    {
      g_clear_pointer(&surface, gdis_isosurface_free);
      return;
    }

  if (!surface)
    g_hash_table_remove(self->iso_surfaces, model);
  else
    g_hash_table_replace(self->iso_surfaces, model, surface);

  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_isosurface_tool(self);
}

static void
gdis_gtk4_window_clear_isosurface(GdisGtk4Window *self,
                                  GdisModel *model,
                                  const char *reason)
{
  g_return_if_fail(self != NULL);

  if (!self->iso_surfaces || !model)
    return;

  if (g_hash_table_remove(self->iso_surfaces, model) && reason)
    gdis_gtk4_window_log(self, "%s", reason);

  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_isosurface_tool(self);
}

static void
gdis_gtk4_window_sync_display_tool(GdisGtk4Window *self)
{
  GdisDisplayTool *tool;

  g_return_if_fail(self != NULL);

  tool = self->display_tool;
  if (!tool)
    return;

  if (tool->atoms_toggle && self->show_atoms_toggle)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->atoms_toggle),
                                 gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_atoms_toggle)));
  if (tool->bonds_toggle && self->show_bonds_toggle)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->bonds_toggle),
                                 gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_bonds_toggle)));
  if (tool->cell_toggle && self->show_cell_toggle)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->cell_toggle),
                                 gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_cell_toggle)));
  if (tool->labels_toggle && self->show_labels_toggle)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->labels_toggle),
                                 gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->show_labels_toggle)));
}

static void
on_display_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisDisplayTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->display_tool = NULL;

  g_free(tool);
}

static void
on_display_tool_toggle_toggled(GtkToggleButton *button, gpointer user_data)
{
  GdisDisplayTool *tool;
  GdisGtk4Window *self;
  const char *toggle_name;
  gboolean active;

  tool = user_data;
  if (!tool || !tool->owner)
    return;

  self = tool->owner;
  toggle_name = g_object_get_data(G_OBJECT(button), "display-toggle");
  active = gtk_toggle_button_get_active(button);
  if (!toggle_name)
    return;

  if (g_strcmp0(toggle_name, "atoms") == 0 && self->show_atoms_toggle)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(self->show_atoms_toggle), active);
  else if (g_strcmp0(toggle_name, "bonds") == 0 && self->show_bonds_toggle)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(self->show_bonds_toggle), active);
  else if (g_strcmp0(toggle_name, "cell") == 0 && self->show_cell_toggle)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(self->show_cell_toggle), active);
  else if (g_strcmp0(toggle_name, "labels") == 0 && self->show_labels_toggle)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(self->show_labels_toggle), active);
}

static void
gdis_gtk4_window_present_display_tool(GdisGtk4Window *self)
{
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *grid;
  GtkWidget *button;
  GtkWidget *label;
  GdisDisplayTool *tool;

  g_return_if_fail(self != NULL);

  if (self->display_tool && GTK_IS_WINDOW(self->display_tool->window))
    {
      gdis_gtk4_window_sync_display_tool(self);
      gtk_window_present(GTK_WINDOW(self->display_tool->window));
      return;
    }

  tool = g_new0(GdisDisplayTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Display Properties");
  gtk_window_set_default_size(GTK_WINDOW(window), 560, 360);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  gtk_box_append(GTK_BOX(root), new_section_button("Viewer Toggles"));

  grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);

  tool->atoms_toggle = new_toggle_button("Atoms", TRUE);
  tool->bonds_toggle = new_toggle_button("Bonds", TRUE);
  tool->cell_toggle = new_toggle_button("Cell", TRUE);
  tool->labels_toggle = new_toggle_button("Labels", FALSE);

  g_object_set_data(G_OBJECT(tool->atoms_toggle), "display-toggle", "atoms");
  g_object_set_data(G_OBJECT(tool->bonds_toggle), "display-toggle", "bonds");
  g_object_set_data(G_OBJECT(tool->cell_toggle), "display-toggle", "cell");
  g_object_set_data(G_OBJECT(tool->labels_toggle), "display-toggle", "labels");

  g_signal_connect(tool->atoms_toggle, "toggled", G_CALLBACK(on_display_tool_toggle_toggled), tool);
  g_signal_connect(tool->bonds_toggle, "toggled", G_CALLBACK(on_display_tool_toggle_toggled), tool);
  g_signal_connect(tool->cell_toggle, "toggled", G_CALLBACK(on_display_tool_toggle_toggled), tool);
  g_signal_connect(tool->labels_toggle, "toggled", G_CALLBACK(on_display_tool_toggle_toggled), tool);

  gtk_widget_set_hexpand(tool->atoms_toggle, TRUE);
  gtk_widget_set_hexpand(tool->bonds_toggle, TRUE);
  gtk_widget_set_hexpand(tool->cell_toggle, TRUE);
  gtk_widget_set_hexpand(tool->labels_toggle, TRUE);

  gtk_grid_attach(GTK_GRID(grid), tool->atoms_toggle, 0, 0, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), tool->bonds_toggle, 1, 0, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), tool->cell_toggle, 0, 1, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), tool->labels_toggle, 1, 1, 1, 1);
  gtk_box_append(GTK_BOX(root), grid);

  gtk_box_append(GTK_BOX(root), new_section_button("View Presets"));

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

  button = gtk_button_new_with_label("Reset Model Images");
  g_signal_connect(button, "clicked", G_CALLBACK(on_reset_image_limits_clicked), self);
  gtk_widget_set_hexpand(button, TRUE);
  gtk_grid_attach(GTK_GRID(grid), button, 0, 2, 2, 1);
  gtk_box_append(GTK_BOX(root), grid);

  label = gtk_label_new("This GTK4 window restores the legacy everyday display subset: viewer toggles, view presets, and model-image reset. The deeper camera, lights, colours, OpenGL, and POVRay pages are still pending.");
  gtk_label_set_wrap(GTK_LABEL(label), TRUE);
  gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
  gtk_widget_add_css_class(label, "gdis-muted");
  gtk_box_append(GTK_BOX(root), label);

  g_signal_connect(window, "destroy", G_CALLBACK(on_display_tool_destroy), tool);

  self->display_tool = tool;
  gdis_gtk4_window_sync_display_tool(self);
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

static gint
gdis_gtk4_window_find_atom_index_in_array(const GArray *array, guint atom_index)
{
  g_return_val_if_fail(array != NULL, -1);

  for (guint i = 0; i < array->len; i++)
    {
      if (g_array_index(array, guint, i) == atom_index)
        return (gint) i;
    }

  return -1;
}

static GArray *
gdis_gtk4_window_copy_atom_array(const GArray *array)
{
  GArray *copy;

  copy = g_array_new(FALSE, FALSE, sizeof(guint));
  if (array && array->len > 0)
    g_array_append_vals(copy, array->data, array->len);

  return copy;
}

static void
gdis_gtk4_window_toggle_atom_in_array(GArray *array, guint atom_index)
{
  g_return_if_fail(array != NULL);

  for (guint i = 0; i < array->len; i++)
    {
      if (g_array_index(array, guint, i) == atom_index)
        {
          g_array_remove_index(array, i);
          return;
        }
    }

  g_array_append_val(array, atom_index);
}

static void
gdis_measurement_record_free(gpointer data)
{
  g_free(data);
}

static void
gdis_uint_array_free(gpointer data)
{
  GArray *array;

  array = data;
  if (array)
    g_array_free(array, TRUE);
}

static void
gdis_measurement_record_array_free(gpointer data)
{
  GPtrArray *records;

  records = data;
  if (!records)
    return;

  g_ptr_array_free(records, TRUE);
}

static GPtrArray *
gdis_measurement_record_array_clone(const GPtrArray *records)
{
  GPtrArray *copy;

  copy = g_ptr_array_new_with_free_func(gdis_measurement_record_free);
  if (!records)
    return copy;

  for (guint i = 0; i < records->len; i++)
    {
      const GdisMeasurementRecord *record;
      GdisMeasurementRecord *record_copy;

      record = g_ptr_array_index((GPtrArray *) records, i);
      record_copy = g_new(GdisMeasurementRecord, 1);
      *record_copy = *record;
      g_ptr_array_add(copy, record_copy);
    }

  return copy;
}

static void
gdis_undo_entry_free(gpointer data)
{
  GdisUndoEntry *entry;

  entry = data;
  if (!entry)
    return;

  g_clear_pointer(&entry->model_snapshot, gdis_model_free);
  g_clear_pointer(&entry->selected_atoms, gdis_uint_array_free);
  g_clear_pointer(&entry->picked_atoms, gdis_uint_array_free);
  g_clear_pointer(&entry->measurement_records, g_ptr_array_unref);
  g_free(entry);
}

static GPtrArray *
gdis_gtk4_window_get_undo_stack(GdisGtk4Window *self,
                                GdisModel *model,
                                gboolean create)
{
  GPtrArray *stack;

  g_return_val_if_fail(self != NULL, NULL);

  if (!self->undo_stacks || !model)
    return NULL;

  stack = g_hash_table_lookup(self->undo_stacks, model);
  if (!stack && create)
    {
      stack = g_ptr_array_new_with_free_func(gdis_undo_entry_free);
      g_hash_table_insert(self->undo_stacks, model, stack);
    }

  return stack;
}

static void
gdis_gtk4_window_update_undo_action(GdisGtk4Window *self)
{
  GAction *action;
  GPtrArray *stack;
  gboolean enabled;

  g_return_if_fail(self != NULL);

  action = g_action_map_lookup_action(G_ACTION_MAP(self->app), "undo");
  if (!action || !G_IS_SIMPLE_ACTION(action))
    return;

  stack = gdis_gtk4_window_get_undo_stack(self, self->active_model, FALSE);
  enabled = (stack != NULL && stack->len > 0);
  g_simple_action_set_enabled(G_SIMPLE_ACTION(action), enabled);
}

static gboolean
gdis_gtk4_window_push_undo_snapshot(GdisGtk4Window *self, const char *reason)
{
  GPtrArray *stack;
  GPtrArray *records;
  GdisUndoEntry *entry;

  g_return_val_if_fail(self != NULL, FALSE);

  if (!self->active_model)
    return FALSE;

  stack = gdis_gtk4_window_get_undo_stack(self, self->active_model, TRUE);
  if (!stack)
    return FALSE;

  entry = g_new0(GdisUndoEntry, 1);
  entry->model_snapshot = gdis_model_clone(self->active_model);
  entry->selected_atoms = gdis_gtk4_window_copy_atom_array(self->selected_atoms);
  entry->picked_atoms = gdis_gtk4_window_copy_atom_array(self->picked_atoms);
  records = gdis_gtk4_window_get_measurement_records(self, self->active_model, FALSE);
  entry->measurement_records = gdis_measurement_record_array_clone(records);
  entry->selected_atom_index = self->selected_atom_index;
  entry->fragment_anchor_index = self->fragment_anchor_index;
  g_ptr_array_add(stack, entry);
  while (stack->len > 32)
    g_ptr_array_remove_index(stack, 0);

  if (reason)
    gdis_gtk4_window_log(self, "%s", reason);

  gdis_gtk4_window_update_undo_action(self);
  return TRUE;
}

static void
gdis_gtk4_window_discard_undo_snapshot(GdisGtk4Window *self)
{
  GPtrArray *stack;

  g_return_if_fail(self != NULL);

  stack = gdis_gtk4_window_get_undo_stack(self, self->active_model, FALSE);
  if (stack && stack->len > 0)
    g_ptr_array_remove_index(stack, stack->len - 1);

  gdis_gtk4_window_update_undo_action(self);
}

static gboolean
gdis_gtk4_window_perform_undo(GdisGtk4Window *self)
{
  GPtrArray *stack;
  GdisUndoEntry *entry;
  GError *error;

  g_return_val_if_fail(self != NULL, FALSE);

  if (!self->active_model)
    return FALSE;

  stack = gdis_gtk4_window_get_undo_stack(self, self->active_model, FALSE);
  if (!stack || stack->len == 0)
    return FALSE;

  entry = g_ptr_array_index(stack, stack->len - 1);
  error = NULL;
  if (!gdis_model_copy_from(self->active_model, entry->model_snapshot, &error))
    {
      gdis_gtk4_window_log(self, "Undo failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return FALSE;
    }

  if (self->measurement_records)
    {
      GPtrArray *records_copy;

      records_copy = gdis_measurement_record_array_clone(entry->measurement_records);
      g_hash_table_insert(self->measurement_records, self->active_model, records_copy);
    }

  gdis_gtk4_window_clear_selected_atoms(self);
  if (entry->selected_atoms && entry->selected_atoms->len > 0)
    g_array_append_vals(self->selected_atoms,
                        entry->selected_atoms->data,
                        entry->selected_atoms->len);

  gdis_gtk4_window_clear_atom_picks(self);
  if (entry->picked_atoms && entry->picked_atoms->len > 0)
    g_array_append_vals(self->picked_atoms,
                        entry->picked_atoms->data,
                        entry->picked_atoms->len);

  self->selected_atom_index = entry->selected_atom_index;
  self->fragment_anchor_index = entry->fragment_anchor_index;
  if (self->active_model->atoms &&
      self->selected_atom_index >= self->active_model->atoms->len)
    self->selected_atom_index = INVALID_ATOM_INDEX;
  if (self->fragment_anchor_index != INVALID_ATOM_INDEX &&
      (!self->active_model->atoms ||
       self->fragment_anchor_index >= self->active_model->atoms->len))
    self->fragment_anchor_index = INVALID_ATOM_INDEX;

  g_ptr_array_remove_index(stack, stack->len - 1);
  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_log(self, "Undid the most recent model edit.\n");
  gdis_gtk4_window_update_undo_action(self);
  return TRUE;
}

static GPtrArray *
gdis_gtk4_window_get_measurement_records(GdisGtk4Window *self,
                                         GdisModel *model,
                                         gboolean create)
{
  GPtrArray *records;

  g_return_val_if_fail(self != NULL, NULL);

  if (!self->measurement_records || !model)
    return NULL;

  records = g_hash_table_lookup(self->measurement_records, model);
  if (!records && create)
    {
      records = g_ptr_array_new_with_free_func(gdis_measurement_record_free);
      g_hash_table_insert(self->measurement_records, model, records);
    }

  return records;
}

static void
gdis_gtk4_window_clear_saved_measurements(GdisGtk4Window *self,
                                          GdisModel *model,
                                          const char *reason)
{
  GPtrArray *records;

  g_return_if_fail(self != NULL);

  records = gdis_gtk4_window_get_measurement_records(self, model, FALSE);
  if (records)
    g_ptr_array_set_size(records, 0);

  if (reason)
    gdis_gtk4_window_log(self, "%s", reason);
}

static gboolean
gdis_gtk4_window_should_capture_measure_picks(const GdisGtk4Window *self)
{
  g_return_val_if_fail(self != NULL, FALSE);

  return self->measure_tool &&
         self->measure_tool->capture_toggle &&
         gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->measure_tool->capture_toggle));
}

static gboolean
gdis_gtk4_window_should_record_viewer_picks(const GdisGtk4Window *self)
{
  g_return_val_if_fail(self != NULL, FALSE);

  if (self->click_mode != GDIS_CLICK_MODE_SELECT)
    return TRUE;

  return gdis_gtk4_window_should_capture_measure_picks(self);
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

static void
gdis_gtk4_window_finish_selection_change(GdisGtk4Window *self, const char *log_message)
{
  g_return_if_fail(self != NULL);

  self->fragment_anchor_index = INVALID_ATOM_INDEX;

  if (!self->active_model || !self->selected_atoms || self->selected_atoms->len == 0)
    {
      self->selected_atom_index = INVALID_ATOM_INDEX;
    }
  else if (self->selected_atom_index == INVALID_ATOM_INDEX ||
           self->selected_atom_index >= self->active_model->atoms->len ||
           !gdis_gtk4_window_atom_array_contains(self->selected_atoms, self->selected_atom_index))
    {
      self->selected_atom_index = g_array_index(self->selected_atoms, guint, 0);
    }

  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);

  if (log_message)
    gdis_gtk4_window_log(self, "%s", log_message);
}

static void
gdis_gtk4_window_select_all_atoms(GdisGtk4Window *self)
{
  guint i;

  g_return_if_fail(self != NULL);

  if (!self->active_model)
    {
      gdis_gtk4_window_log(self, "Select all requested, but no active model is loaded.\n");
      return;
    }

  if (!self->active_model->atoms || self->active_model->atoms->len == 0)
    {
      gdis_gtk4_window_log(self, "Select all requested, but the active model has no atoms.\n");
      return;
    }

  gdis_gtk4_window_clear_selected_atoms(self);
  for (i = 0; i < self->active_model->atoms->len; i++)
    gdis_gtk4_window_add_selected_atom(self, i);

  gdis_gtk4_window_finish_selection_change(self, "Selected all atoms in the active model.\n");
}

static void
gdis_gtk4_window_invert_selected_atoms(GdisGtk4Window *self)
{
  gboolean *selected_mask;
  guint atom_count;
  guint i;

  g_return_if_fail(self != NULL);

  if (!self->active_model)
    {
      gdis_gtk4_window_log(self, "Invert selection requested, but no active model is loaded.\n");
      return;
    }

  if (!self->active_model->atoms || self->active_model->atoms->len == 0)
    {
      gdis_gtk4_window_log(self, "Invert selection requested, but the active model has no atoms.\n");
      return;
    }

  atom_count = self->active_model->atoms->len;
  selected_mask = g_new0(gboolean, atom_count);

  for (i = 0; self->selected_atoms && i < self->selected_atoms->len; i++)
    {
      guint atom_index;

      atom_index = g_array_index(self->selected_atoms, guint, i);
      if (atom_index < atom_count)
        selected_mask[atom_index] = TRUE;
    }

  gdis_gtk4_window_clear_selected_atoms(self);
  for (i = 0; i < atom_count; i++)
    {
      if (!selected_mask[i])
        gdis_gtk4_window_add_selected_atom(self, i);
    }
  g_free(selected_mask);

  gdis_gtk4_window_finish_selection_change(self, "Inverted the current atom selection.\n");
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
      if (clicked_atom->region < 0)
        {
          gdis_gtk4_window_add_selected_atom(self, atom_index);
          gdis_gtk4_window_log(self,
                               "The clicked atom has no legacy region label, so region selection fell back to that atom only.\n");
        }
      else
        {
          for (i = 0; i < self->active_model->atoms->len; i++)
            {
              const GdisAtom *atom;

              atom = g_ptr_array_index(self->active_model->atoms, i);
              if (atom->region == clicked_atom->region)
                gdis_gtk4_window_add_selected_atom(self, i);
            }
        }
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
gdis_gtk4_window_append_measurement_line(GString *report,
                                         const GdisModel *model,
                                         GdisMeasureMode mode,
                                         const guint *atom_indices,
                                         guint count,
                                         const char *prefix)
{
  guint used_indices[4] = {0, 0, 0, 0};
  guint used_count;
  GdisMeasureMode resolved_mode;
  gdouble value;
  GError *error;

  g_return_if_fail(report != NULL);

  error = NULL;
  if (!gdis_measure_calculate(model,
                              mode,
                              atom_indices,
                              count,
                              used_indices,
                              &used_count,
                              &resolved_mode,
                              &value,
                              &error))
    {
      g_string_append_printf(report,
                             "%s%s\n",
                             prefix ? prefix : "",
                             error ? error->message : "Measurement unavailable.");
      g_clear_error(&error);
      return;
    }

  g_string_append_printf(report, "%s", prefix ? prefix : "");
  switch (resolved_mode)
    {
    case GDIS_MEASURE_MODE_DISTANCE:
      g_string_append_printf(report,
                             "Distance %s-%s = %.5f A\n",
                             ((GdisAtom *) g_ptr_array_index(model->atoms, used_indices[0]))->label,
                             ((GdisAtom *) g_ptr_array_index(model->atoms, used_indices[1]))->label,
                             value);
      break;
    case GDIS_MEASURE_MODE_ANGLE:
      g_string_append_printf(report,
                             "Angle %s-%s-%s = %.3f deg\n",
                             ((GdisAtom *) g_ptr_array_index(model->atoms, used_indices[0]))->label,
                             ((GdisAtom *) g_ptr_array_index(model->atoms, used_indices[1]))->label,
                             ((GdisAtom *) g_ptr_array_index(model->atoms, used_indices[2]))->label,
                             value);
      break;
    case GDIS_MEASURE_MODE_TORSION:
      g_string_append_printf(report,
                             "Torsion %s-%s-%s-%s = %.3f deg\n",
                             ((GdisAtom *) g_ptr_array_index(model->atoms, used_indices[0]))->label,
                             ((GdisAtom *) g_ptr_array_index(model->atoms, used_indices[1]))->label,
                             ((GdisAtom *) g_ptr_array_index(model->atoms, used_indices[2]))->label,
                             ((GdisAtom *) g_ptr_array_index(model->atoms, used_indices[3]))->label,
                             value);
      break;
    case GDIS_MEASURE_MODE_AUTO:
    default:
      break;
    }
}

static char *
gdis_gtk4_window_build_measurement_tool_report(GdisGtk4Window *self)
{
  GString *report;
  const GArray *source_atoms;
  const char *source_name;
  GPtrArray *records;

  g_return_val_if_fail(self != NULL, g_strdup("Measurement tool unavailable."));

  report = g_string_new("");
  if (!self->active_model)
    {
      g_string_append(report, "No active model loaded.");
      return g_string_free(report, FALSE);
    }

  source_atoms = NULL;
  source_name = "none";
  if (self->picked_atoms && self->picked_atoms->len > 0)
    {
      source_atoms = self->picked_atoms;
      source_name = "viewer picks";
    }
  else if (self->selected_atoms && self->selected_atoms->len > 0)
    {
      source_atoms = self->selected_atoms;
      source_name = "current selection";
    }

  g_string_append_printf(report,
                         "Measurements\nModel: %s\nMode: %s\nCapture Viewer Picks: %s\nCurrent source: %s\n\n",
                         self->active_model->basename,
                         self->measure_tool->mode == GDIS_MEASURE_MODE_DISTANCE ? "Distance" :
                         self->measure_tool->mode == GDIS_MEASURE_MODE_ANGLE ? "Angle" :
                         self->measure_tool->mode == GDIS_MEASURE_MODE_TORSION ? "Torsion" : "Auto",
                         gdis_gtk4_window_should_capture_measure_picks(self) ? "on" : "off",
                         source_name);

  if (source_atoms && source_atoms->len > 0)
    {
      guint start_index;

      start_index = source_atoms->len > 4 ? source_atoms->len - 4 : 0;
      for (guint i = start_index; i < source_atoms->len; i++)
        {
          guint atom_index;
          const GdisAtom *atom;

          atom_index = g_array_index(source_atoms, guint, i);
          if (atom_index >= self->active_model->atoms->len)
            continue;
          atom = g_ptr_array_index(self->active_model->atoms, atom_index);
          g_string_append_printf(report,
                                 "Atom %u: %s [%s] #%u  %.5f %.5f %.5f\n",
                                 i - start_index + 1,
                                 atom->label,
                                 atom->element,
                                 atom->serial,
                                 atom->position[0],
                                 atom->position[1],
                                 atom->position[2]);
        }
      g_string_append(report, "\nCurrent result:\n");
      gdis_gtk4_window_append_measurement_line(report,
                                               self->active_model,
                                               self->measure_tool->mode,
                                               (const guint *) source_atoms->data,
                                               source_atoms->len,
                                               "  ");
    }
  else
    {
      g_string_append(report,
                      "No measurement atoms are active yet.\n"
                      "Use Shift+click or box selection, or enable Capture Viewer Picks.\n\n");
    }

  g_string_append(report, "\nSaved measurements:\n");
  records = gdis_gtk4_window_get_measurement_records(self, self->active_model, FALSE);
  if (!records || records->len == 0)
    {
      g_string_append(report, "  none\n");
    }
  else
    {
      for (guint i = 0; i < records->len; i++)
        {
          const GdisMeasurementRecord *record;
          g_autofree gchar *prefix = NULL;

          record = g_ptr_array_index(records, i);
          prefix = g_strdup_printf("  %u. ", i + 1);
          gdis_gtk4_window_append_measurement_line(report,
                                                   self->active_model,
                                                   record->mode,
                                                   record->atom_indices,
                                                   record->atom_count,
                                                   prefix);
        }
    }

  g_string_append(report,
                  "\nActions:\n"
                  "  - Capture Viewer Picks: click atoms without disturbing the selection set.\n"
                  "  - Use Selection: load the current selection into the pick list.\n"
                  "  - Add Current: store the current distance/angle/torsion.\n"
                  "  - Delete Last / Clear Saved: manage the stored list.\n");

  return g_string_free(report, FALSE);
}

static void
gdis_gtk4_window_refresh_measure_tool(GdisGtk4Window *self)
{
  g_autofree char *report = NULL;

  g_return_if_fail(self != NULL);

  if (!self->measure_tool || !self->measure_tool->buffer)
    return;

  report = gdis_gtk4_window_build_measurement_tool_report(self);
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
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->ff_type_entry),
                            selected_atom->ff_type ? selected_atom->ff_type : "");
      g_snprintf(number, sizeof(number), "%d", selected_atom->region);
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->region_entry),
                            selected_atom->region >= 0 ? number : "");
      g_snprintf(number, sizeof(number), "%.6f", selected_atom->position[0]);
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->x_entry), number);
      g_snprintf(number, sizeof(number), "%.6f", selected_atom->position[1]);
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->y_entry), number);
      g_snprintf(number, sizeof(number), "%.6f", selected_atom->position[2]);
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->z_entry), number);
    }
  else
    {
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->label_entry), "");
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->element_entry), "");
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->ff_type_entry), "");
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->region_entry), "");
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->x_entry), "0.0");
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->y_entry), "0.0");
      gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->z_entry), "0.0");
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
                  "  Apply Region To Selected Set: writes the Region field onto every selected atom.\n"
                  "  Delete Selected Atom: removes the current yellow-highlighted atom.\n"
                  "  Delete Selected Set: removes the blue highlighted selection group.\n"
                  "  Delete Picked Atoms: removes the atoms in the current pick history.\n"
                  "  Add Bond Between Last 2 Picks: creates or updates an explicit bond.\n"
                  "  Remove Bond Between Last 2 Picks: removes that bond.\n"
                  "  Pick Add Bond / Pick Remove Bond: use the next 2 viewer picks as a bond edit command.\n"
                  "  Use Selection As Picks: copies up to the last 4 selected atoms into the pick list.\n"
                  "  Clear Picks: clears the measurement/edit pick set.\n\n");
  g_string_append(report, "Flags: S = selected atom, G = selected group, P = picked atom\n\n");
  g_string_append(report, "Flg  Serial  Elem  FFType     Region  Label                X           Y           Z\n");

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
                             " %3s   %5u   %-3s   %-10s %6d  %-16s %10.4f %10.4f %10.4f\n",
                             flags,
                             atom->serial,
                             atom->element,
                             atom->ff_type ? atom->ff_type : "",
                             atom->region,
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
gdis_gtk4_window_parse_region_entry(GtkWidget *entry,
                                    gint *value_out,
                                    GError **error)
{
  const char *text;
  gchar *endptr;
  glong value;

  g_return_val_if_fail(GTK_IS_EDITABLE(entry), FALSE);
  g_return_val_if_fail(value_out != NULL, FALSE);

  text = gtk_editable_get_text(GTK_EDITABLE(entry));
  if (!text || text[0] == '\0')
    {
      *value_out = -1;
      return TRUE;
    }

  value = g_ascii_strtoll(text, &endptr, 10);
  if (!endptr || *endptr != '\0' || value < -1 || value > 3)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Region must be blank or an integer between -1 and 3.");
      return FALSE;
    }

  *value_out = (gint) value;
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
      if (self->selected_atoms && self->selected_atoms->len == 2)
        {
          *atom_index_a = g_array_index(self->selected_atoms, guint, 0);
          *atom_index_b = g_array_index(self->selected_atoms, guint, 1);
          return TRUE;
        }

      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Pick 2 atoms in the viewer first, or select exactly 2 atoms.");
      return FALSE;
    }

  *atom_index_a = g_array_index(self->picked_atoms, guint, self->picked_atoms->len - 2);
  *atom_index_b = g_array_index(self->picked_atoms, guint, self->picked_atoms->len - 1);
  return TRUE;
}

static gboolean
gdis_gtk4_window_set_picks_from_array(GdisGtk4Window *self,
                                      const GArray *atoms,
                                      guint limit,
                                      const char *log_message)
{
  guint start_index;

  g_return_val_if_fail(self != NULL, FALSE);

  if (!atoms || atoms->len == 0)
    {
      gdis_gtk4_window_log(self, "No atoms are available for picking.\n");
      return FALSE;
    }

  gdis_gtk4_window_clear_atom_picks(self);
  start_index = atoms->len > limit ? atoms->len - limit : 0;
  for (guint i = start_index; i < atoms->len; i++)
    {
      guint atom_index;

      atom_index = g_array_index(atoms, guint, i);
      if (self->active_model && atom_index < self->active_model->atoms->len)
        g_array_append_val(self->picked_atoms, atom_index);
    }

  if (log_message)
    gdis_gtk4_window_log(self, "%s", log_message);

  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
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
on_measure_capture_toggled(GtkToggleButton *button, gpointer user_data)
{
  GdisMeasureTool *tool;

  tool = user_data;
  if (!tool || !tool->owner)
    return;

  gdis_gtk4_window_log(tool->owner,
                       "Measurement viewer-pick capture %s.\n",
                       gtk_toggle_button_get_active(button) ? "enabled" : "disabled");
  gdis_gtk4_window_refresh_measure_tool(tool->owner);
}

static void
on_measure_use_selection_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self)
    return;

  if (!self->selected_atoms || self->selected_atoms->len == 0)
    {
      gdis_gtk4_window_log(self, "Use Selection failed: no atoms are selected.\n");
      return;
    }

  gdis_gtk4_window_set_picks_from_array(self,
                                        self->selected_atoms,
                                        4,
                                        "Loaded the current selection into the measurement pick list.\n");
}

static void
on_measure_add_current_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  const GArray *source_atoms;
  GdisMeasurementRecord *record;
  guint used_indices[4] = {0, 0, 0, 0};
  guint used_count;
  GdisMeasureMode resolved_mode;
  gdouble value;
  GError *error;
  GPtrArray *records;

  (void) button;

  self = user_data;
  if (!self || !self->active_model || !self->measure_tool)
    return;

  source_atoms = NULL;
  if (self->picked_atoms && self->picked_atoms->len > 0)
    source_atoms = self->picked_atoms;
  else if (self->selected_atoms && self->selected_atoms->len > 0)
    source_atoms = self->selected_atoms;

  if (!source_atoms)
    {
      gdis_gtk4_window_log(self, "Add Current failed: no picked or selected atoms are available.\n");
      return;
    }

  error = NULL;
  if (!gdis_measure_calculate(self->active_model,
                              self->measure_tool->mode,
                              (const guint *) source_atoms->data,
                              source_atoms->len,
                              used_indices,
                              &used_count,
                              &resolved_mode,
                              &value,
                              &error))
    {
      gdis_gtk4_window_log(self, "Add Current failed: %s\n",
                           error ? error->message : "measurement unavailable");
      g_clear_error(&error);
      return;
    }
  (void) value;

  records = gdis_gtk4_window_get_measurement_records(self, self->active_model, TRUE);
  record = g_new0(GdisMeasurementRecord, 1);
  record->mode = resolved_mode;
  record->atom_count = used_count;
  memcpy(record->atom_indices, used_indices, sizeof(used_indices));
  g_ptr_array_add(records, record);

  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_log(self, "Saved the current measurement.\n");
}

static void
on_measure_delete_last_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GPtrArray *records;

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  records = gdis_gtk4_window_get_measurement_records(self, self->active_model, FALSE);
  if (!records || records->len == 0)
    {
      gdis_gtk4_window_log(self, "Delete Last failed: there are no saved measurements.\n");
      return;
    }

  g_ptr_array_remove_index(records, records->len - 1);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_log(self, "Removed the last saved measurement.\n");
}

static void
on_measure_clear_saved_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self || !self->active_model)
    return;

  gdis_gtk4_window_clear_saved_measurements(self,
                                            self->active_model,
                                            "Cleared the saved measurement list for the active model.\n");
  gdis_gtk4_window_refresh_measure_tool(self);
}

static void
on_clear_picks_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self)
    return;

  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_set_click_mode(self, GDIS_CLICK_MODE_SELECT);
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_log(self, "Pick history cleared.\n");
}

static void
on_use_selection_as_picks_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self)
    return;

  if (!self->selected_atoms || self->selected_atoms->len == 0)
    {
      gdis_gtk4_window_log(self, "Use Selection As Picks failed: no atoms are selected.\n");
      return;
    }

  gdis_gtk4_window_set_picks_from_array(self,
                                        self->selected_atoms,
                                        4,
                                        "Loaded the current selection into the shared pick list.\n");
}

static void
on_clear_selection_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self)
    return;

  self->selected_atom_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_update_details(self);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_log(self, "Selection cleared.\n");
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
  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_delete_atoms(self->active_model, &index, 1, &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
      gdis_gtk4_window_log(self, "Delete failed: %s\n", error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  self->selected_atom_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_clear_saved_measurements(self,
                                            self->active_model,
                                            "Saved measurements were cleared because atom indices changed.\n");
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
  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_delete_atoms(self->active_model,
                               (const guint *) self->picked_atoms->data,
                               self->picked_atoms->len,
                               &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
      gdis_gtk4_window_log(self, "Delete failed: %s\n", error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  self->selected_atom_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_clear_saved_measurements(self,
                                            self->active_model,
                                            "Saved measurements were cleared because atom indices changed.\n");
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
  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_delete_atoms(self->active_model,
                               (const guint *) self->selected_atoms->data,
                               self->selected_atoms->len,
                               &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
      gdis_gtk4_window_log(self, "Delete selected set failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  self->selected_atom_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_clear_saved_measurements(self,
                                            self->active_model,
                                            "Saved measurements were cleared because atom indices changed.\n");
  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_log(self, "Deleted the selected atom set.\n");
}

static void
on_edit_apply_selected_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  gint region;
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
      !gdis_gtk4_window_parse_entry_double(self->edit_tool->z_entry, "Z", &z, &error) ||
      !gdis_gtk4_window_parse_region_entry(self->edit_tool->region_entry, &region, &error))
    {
      gdis_gtk4_window_log(self, "Edit failed: %s\n", error ? error->message : "invalid coordinates");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_update_atom(self->active_model,
                              self->selected_atom_index,
                              gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->label_entry)),
                              gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->element_entry)),
                              gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->ff_type_entry)),
                              region,
                              x,
                              y,
                              z,
                              &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
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
  gint region;
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
      !gdis_gtk4_window_parse_entry_double(self->edit_tool->z_entry, "Z", &z, &error) ||
      !gdis_gtk4_window_parse_region_entry(self->edit_tool->region_entry, &region, &error))
    {
      gdis_gtk4_window_log(self, "Add atom failed: %s\n", error ? error->message : "invalid coordinates");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_add_atom(self->active_model,
                           gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->label_entry)),
                           gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->element_entry)),
                           gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->ff_type_entry)),
                           region,
                           x,
                           y,
                           z,
                           &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
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
on_edit_apply_region_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  gint region;

  (void) button;

  self = user_data;
  if (!self || !self->active_model || !self->edit_tool)
    return;

  if (!self->selected_atoms || self->selected_atoms->len == 0)
    {
      gdis_gtk4_window_log(self, "Apply Region failed: no atoms are selected.\n");
      return;
    }

  error = NULL;
  if (!gdis_gtk4_window_parse_region_entry(self->edit_tool->region_entry, &region, &error))
    {
      gdis_gtk4_window_log(self, "Apply Region failed: %s\n",
                           error ? error->message : "invalid region");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  for (guint i = 0; i < self->selected_atoms->len; i++)
    {
      guint atom_index;
      GdisAtom *atom;

      atom_index = g_array_index(self->selected_atoms, guint, i);
      if (atom_index >= self->active_model->atoms->len)
        continue;

      atom = g_ptr_array_index(self->active_model->atoms, atom_index);
      atom->region = region;
    }

  gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
  gdis_gtk4_window_log(self, "Updated the region for the selected atom set.\n");
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

  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_add_explicit_bond(self->active_model,
                                    atom_index_a,
                                    atom_index_b,
                                    (guint8) order,
                                    &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
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

  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_remove_bond(self->active_model, atom_index_a, atom_index_b, &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
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

  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_set_click_mode(self, GDIS_CLICK_MODE_ADD_BOND);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
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

  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_set_click_mode(self, GDIS_CLICK_MODE_REMOVE_BOND);
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
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

  tool->capture_toggle = gtk_toggle_button_new_with_label("Capture Viewer Picks");
  g_signal_connect(tool->capture_toggle, "toggled", G_CALLBACK(on_measure_capture_toggled), tool);
  gtk_box_append(GTK_BOX(row), tool->capture_toggle);

  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_box_append(GTK_BOX(root), row);

  button = gtk_button_new_with_label("Use Selection");
  g_signal_connect(button, "clicked", G_CALLBACK(on_measure_use_selection_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Add Current");
  g_signal_connect(button, "clicked", G_CALLBACK(on_measure_add_current_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Delete Last");
  g_signal_connect(button, "clicked", G_CALLBACK(on_measure_delete_last_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Clear Saved");
  g_signal_connect(button, "clicked", G_CALLBACK(on_measure_clear_saved_clicked), self);
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

  label = gtk_label_new("FF Type");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 4, 0, 1, 1);
  tool->ff_type_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->ff_type_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->ff_type_entry, 5, 0, 1, 1);

  label = gtk_label_new("Region");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 6, 0, 1, 1);
  tool->region_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->region_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->region_entry, 7, 0, 1, 1);

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
  gtk_grid_attach(GTK_GRID(grid), label, 6, 1, 1, 1);
  tool->bond_order_entry = gtk_entry_new();
  gtk_editable_set_text(GTK_EDITABLE(tool->bond_order_entry), "1");
  gtk_widget_set_hexpand(tool->bond_order_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->bond_order_entry, 7, 1, 1, 1);

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

  button = gtk_button_new_with_label("Apply Region To Selected Set");
  g_signal_connect(button, "clicked", G_CALLBACK(on_edit_apply_region_clicked), self);
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

  button = gtk_button_new_with_label("Use Selection As Picks");
  g_signal_connect(button, "clicked", G_CALLBACK(on_use_selection_as_picks_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Clear Picks");
  g_signal_connect(button, "clicked", G_CALLBACK(on_clear_picks_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Clear Selection");
  g_signal_connect(button, "clicked", G_CALLBACK(on_clear_selection_clicked), self);
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
on_surface_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisSurfaceTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->surface_tool = NULL;
  g_free(tool);
}

static void
gdis_gtk4_window_refresh_surface_tool(GdisGtk4Window *self)
{
  GdisSurfaceTool *tool;
  gboolean enabled;
  g_autofree gchar *report = NULL;

  g_return_if_fail(self != NULL);

  tool = self->surface_tool;
  if (!tool)
    return;

  enabled = (self->active_model != NULL &&
             self->active_model->periodic &&
             self->active_model->periodicity >= 3);
  if (tool->execute_button)
    gtk_widget_set_sensitive(tool->execute_button, enabled);

  report = gdis_report_surface(self->active_model);
  if (tool->report_buffer)
    gtk_text_buffer_set_text(tool->report_buffer, report ? report : "Surface tool unavailable.", -1);
}

static void
on_surface_execute_clicked(GtkButton *button, gpointer user_data)
{
  GdisSurfaceTool *tool;
  GdisGtk4Window *self;
  GdisModel *surface_model;
  GError *error;
  gchar *summary;
  gint h;
  gint k;
  gint l;
  gdouble shift;
  guint region_a;
  guint region_b;
  guint repeat_a;
  guint repeat_b;
  gdouble vacuum;

  (void) button;

  tool = user_data;
  self = tool ? tool->owner : NULL;
  if (!self || !self->active_model || !tool)
    return;

  h = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->h_spin));
  k = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->k_spin));
  l = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->l_spin));
  shift = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->shift_spin));
  region_a = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->region_a_spin));
  region_b = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->region_b_spin));
  repeat_a = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->repeat_a_spin));
  repeat_b = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->repeat_b_spin));
  vacuum = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->vacuum_spin));

  surface_model = NULL;
  summary = NULL;
  error = NULL;
  if (!gdis_model_build_surface_slab(self->active_model,
                                     h,
                                     k,
                                     l,
                                     shift,
                                     region_a,
                                     region_b,
                                     repeat_a,
                                     repeat_b,
                                     vacuum,
                                     &surface_model,
                                     &summary,
                                     &error))
    {
      const char *message;

      message = error ? error->message : "unknown error";
      gdis_gtk4_window_log(self, "Surface build failed: %s\n", message);
      if (tool->report_buffer)
        gtk_text_buffer_set_text(tool->report_buffer, message, -1);
      g_clear_error(&error);
      g_free(summary);
      return;
    }

  gdis_gtk4_window_add_loaded_model(self, surface_model, TRUE);
  if (summary && tool->report_buffer)
    gtk_text_buffer_set_text(tool->report_buffer, summary, -1);
  gdis_gtk4_window_log(self,
                       "Built a surface slab: (%d %d %d), layers %u + %u, repeats %u x %u, vacuum %.2f A.\n",
                       h,
                       k,
                       l,
                       region_a,
                       region_b,
                       repeat_a,
                       repeat_b,
                       vacuum);
  g_free(summary);
}

static void
gdis_gtk4_window_present_surface_tool(GdisGtk4Window *self)
{
  GdisSurfaceTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *grid;
  GtkWidget *row;
  GtkWidget *label;
  GtkWidget *button;
  GtkWidget *scroller;
  GtkWidget *text_view;

  g_return_if_fail(self != NULL);

  if (self->surface_tool && GTK_IS_WINDOW(self->surface_tool->window))
    {
      gdis_gtk4_window_refresh_surface_tool(self);
      gtk_window_present(GTK_WINDOW(self->surface_tool->window));
      return;
    }

  tool = g_new0(GdisSurfaceTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Surface Builder");
  gtk_window_set_default_size(GTK_WINDOW(window), 760, 560);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(grid), 10);
  gtk_grid_set_row_spacing(GTK_GRID(grid), 10);
  gtk_box_append(GTK_BOX(root), grid);

  label = gtk_label_new("Miller h k l");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_grid_attach(GTK_GRID(grid), row, 1, 0, 1, 1);
  tool->h_spin = gtk_spin_button_new_with_range(-9.0, 9.0, 1.0);
  tool->k_spin = gtk_spin_button_new_with_range(-9.0, 9.0, 1.0);
  tool->l_spin = gtk_spin_button_new_with_range(-9.0, 9.0, 1.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->h_spin), 1.0);
  gtk_box_append(GTK_BOX(row), tool->h_spin);
  gtk_box_append(GTK_BOX(row), tool->k_spin);
  gtk_box_append(GTK_BOX(row), tool->l_spin);

  label = gtk_label_new("Shift (Dhkl units)");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
  tool->shift_spin = gtk_spin_button_new_with_range(-4.0, 4.0, 0.05);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->shift_spin), 2);
  gtk_grid_attach(GTK_GRID(grid), tool->shift_spin, 1, 1, 1, 1);

  label = gtk_label_new("Region layers A / B");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 2, 1, 1);
  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_grid_attach(GTK_GRID(grid), row, 1, 2, 1, 1);
  tool->region_a_spin = gtk_spin_button_new_with_range(0.0, 48.0, 1.0);
  tool->region_b_spin = gtk_spin_button_new_with_range(0.0, 48.0, 1.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->region_a_spin), 4.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->region_b_spin), 4.0);
  gtk_box_append(GTK_BOX(row), tool->region_a_spin);
  gtk_box_append(GTK_BOX(row), tool->region_b_spin);

  label = gtk_label_new("In-plane repeats A / B");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 3, 1, 1);
  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_grid_attach(GTK_GRID(grid), row, 1, 3, 1, 1);
  tool->repeat_a_spin = gtk_spin_button_new_with_range(1.0, 24.0, 1.0);
  tool->repeat_b_spin = gtk_spin_button_new_with_range(1.0, 24.0, 1.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->repeat_a_spin), 1.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->repeat_b_spin), 1.0);
  gtk_box_append(GTK_BOX(row), tool->repeat_a_spin);
  gtk_box_append(GTK_BOX(row), tool->repeat_b_spin);

  label = gtk_label_new("Vacuum thickness");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 4, 1, 1);
  tool->vacuum_spin = gtk_spin_button_new_with_range(0.0, 80.0, 0.5);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->vacuum_spin), 1);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->vacuum_spin), 12.0);
  gtk_grid_attach(GTK_GRID(grid), tool->vacuum_spin, 1, 4, 1, 1);

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

  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_widget_set_halign(row, GTK_ALIGN_END);
  gtk_box_append(GTK_BOX(root), row);

  button = gtk_button_new_with_label("Build Surface");
  tool->execute_button = button;
  g_signal_connect(button, "clicked", G_CALLBACK(on_surface_execute_clicked), tool);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_destroy), window);
  gtk_box_append(GTK_BOX(row), button);

  g_signal_connect(window, "destroy", G_CALLBACK(on_surface_tool_destroy), tool);

  self->surface_tool = tool;
  gdis_gtk4_window_refresh_surface_tool(self);
  gtk_window_present(GTK_WINDOW(window));
}

static void
on_isosurface_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisIsosurfaceTool *tool;
  GdisGtk4Window *owner;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  owner = tool->owner;
  if (owner)
    owner->isosurface_tool = NULL;

  /* The window's child widgets can still emit teardown-time signals that
   * reference this tool, so invalidate the pointers now and let the window's
   * object data own the final g_free() once destruction fully completes.
   */
  tool->owner = NULL;
  tool->window = NULL;
  tool->mode_dropdown = NULL;
  tool->grid_spin = NULL;
  tool->blur_spin = NULL;
  tool->value_spin = NULL;
  tool->execute_button = NULL;
  tool->clear_button = NULL;
  tool->report_buffer = NULL;
}

static void
gdis_gtk4_window_refresh_isosurface_tool(GdisGtk4Window *self)
{
  GdisIsosurfaceTool *tool;
  GdisIsoSurface *surface;
  gboolean has_model;
  gboolean supported_model;
  guint selected_mode;
  gboolean needs_selection;
  gboolean has_selection;
  g_autofree gchar *report = NULL;

  g_return_if_fail(self != NULL);

  tool = self->isosurface_tool;
  if (!tool)
    return;

  has_model = (self->active_model != NULL);
  selected_mode = tool->mode_dropdown ?
    gtk_drop_down_get_selected(GTK_DROP_DOWN(tool->mode_dropdown)) :
    GDIS_ISOSURFACE_MODE_MOLECULAR;
  needs_selection = (selected_mode == GDIS_ISOSURFACE_MODE_PROMOLECULE ||
                     selected_mode == GDIS_ISOSURFACE_MODE_HIRSHFELD);
  has_selection = (self->selected_atoms && self->selected_atoms->len > 0u);
  supported_model = has_model && (!needs_selection || has_selection);
  if (tool->execute_button)
    gtk_widget_set_sensitive(tool->execute_button, supported_model);

  if (tool->clear_button)
    gtk_widget_set_sensitive(tool->clear_button,
                             has_model && gdis_gtk4_window_get_isosurface(self, self->active_model, FALSE) != NULL);

  if (!tool->report_buffer)
    return;

  if (!has_model)
    {
      gtk_text_buffer_set_text(tool->report_buffer,
                               "No active model loaded.\nOpen a structure first, then generate an iso-surface.",
                               -1);
      return;
    }

  if (needs_selection && !has_selection)
    {
      gtk_text_buffer_set_text(tool->report_buffer,
                               "This iso-surface mode follows the legacy selected-atom workflow.\n\n"
                               "Select one or more atoms in the main viewer first, then generate the surface.",
                               -1);
      return;
    }

  surface = gdis_gtk4_window_get_isosurface(self, self->active_model, FALSE);
  if (surface && surface->summary)
    {
      gtk_text_buffer_set_text(tool->report_buffer, surface->summary, -1);
      return;
    }

  report = g_strdup(
    "GTK4 Iso-surfaces\n\n"
    "Restored now:\n"
    "  - Molecular surface\n"
    "  - Promolecule isosurface\n"
    "  - Hirshfeld-style selected-fragment surface\n"
    "  - Analytic electron-density approximation for the whole model\n\n"
    "Still pending:\n"
    "  - True volumetric-grid electron density import/render\n\n"
    "Generate a surface to preview it directly in the main viewer.");
  gtk_text_buffer_set_text(tool->report_buffer, report, -1);
}

static void
on_isosurface_clear_clicked(GtkButton *button, gpointer user_data)
{
  GdisIsosurfaceTool *tool;
  GdisGtk4Window *self;

  (void) button;

  tool = user_data;
  self = tool ? tool->owner : NULL;
  if (!self || !self->active_model)
    return;

  gdis_gtk4_window_clear_isosurface(self,
                                    self->active_model,
                                    "Cleared the active iso-surface preview.\n");
}

static void
on_isosurface_mode_changed(GObject *object, GParamSpec *pspec, gpointer user_data)
{
  GdisIsosurfaceTool *tool;
  GdisGtk4Window *self;
  guint selected_mode;
  gdouble default_value;

  (void) object;
  (void) pspec;

  tool = user_data;
  self = tool ? tool->owner : NULL;
  if (!tool || !self || !tool->value_spin)
    return;

  selected_mode = gtk_drop_down_get_selected(GTK_DROP_DOWN(tool->mode_dropdown));
  default_value = 0.080;
  if (selected_mode == GDIS_ISOSURFACE_MODE_HIRSHFELD)
    default_value = 0.500;
  else if (selected_mode == GDIS_ISOSURFACE_MODE_ELECTRON_DENSITY)
    default_value = 0.100;

  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->value_spin), default_value);
  gdis_gtk4_window_refresh_isosurface_tool(self);
}

static void
on_isosurface_execute_clicked(GtkButton *button, gpointer user_data)
{
  GdisIsosurfaceTool *tool;
  GdisGtk4Window *self;
  GdisIsosurfaceSettings settings;
  GdisIsoSurface *surface;
  GError *error;
  guint selected_mode;

  (void) button;

  tool = user_data;
  self = tool ? tool->owner : NULL;
  if (!self || !self->active_model || !tool)
    return;

  gdis_isosurface_settings_init(&settings);
  selected_mode = gtk_drop_down_get_selected(GTK_DROP_DOWN(tool->mode_dropdown));
  settings.mode = (GdisIsosurfaceMode) selected_mode;
  settings.grid_size = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->grid_spin));
  settings.blur = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->blur_spin));
  settings.isovalue = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->value_spin));
  settings.selected_atoms = self->selected_atoms;

  surface = NULL;
  error = NULL;
  if (!gdis_isosurface_generate(self->active_model, &settings, &surface, &error))
    {
      const char *message;

      message = error ? error->message : "unknown error";
      gdis_gtk4_window_log(self, "Iso-surface generation failed: %s\n", message);
      if (tool->report_buffer)
        gtk_text_buffer_set_text(tool->report_buffer, message, -1);
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_set_isosurface(self, self->active_model, surface);
  if (surface->summary && tool->report_buffer)
    gtk_text_buffer_set_text(tool->report_buffer, surface->summary, -1);
  gdis_gtk4_window_log(self,
                       "Generated %s with %u triangles.\n",
                       gdis_isosurface_mode_label(settings.mode),
                       surface->triangles ? surface->triangles->len : 0u);
}

static void
gdis_gtk4_window_present_isosurface_tool(GdisGtk4Window *self)
{
  static const char *const mode_items[] = {
    "Molecular surface",
    "Promolecule isosurface",
    "Hirshfeld surface",
    "Electron density",
    NULL
  };
  GdisIsosurfaceTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *grid;
  GtkWidget *label;
  GtkWidget *button_row;
  GtkWidget *button;
  GtkWidget *scroller;
  GtkWidget *text_view;
  GtkStringList *mode_model;

  g_return_if_fail(self != NULL);

  if (self->isosurface_tool && GTK_IS_WINDOW(self->isosurface_tool->window))
    {
      gdis_gtk4_window_refresh_isosurface_tool(self);
      gtk_window_present(GTK_WINDOW(self->isosurface_tool->window));
      return;
    }

  tool = g_new0(GdisIsosurfaceTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  g_object_set_data_full(G_OBJECT(window), "gdis-isosurface-tool", tool, g_free);
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_hide_on_close(GTK_WINDOW(window), TRUE);
  gtk_window_set_title(GTK_WINDOW(window), "Iso-surfaces");
  gtk_window_set_default_size(GTK_WINDOW(window), 760, 520);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(grid), 10);
  gtk_grid_set_row_spacing(GTK_GRID(grid), 10);
  gtk_box_append(GTK_BOX(root), grid);

  label = gtk_label_new("Iso-surface type");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  mode_model = gtk_string_list_new(mode_items);
  tool->mode_dropdown = gtk_drop_down_new(G_LIST_MODEL(mode_model), NULL);
  g_object_unref(mode_model);
  gtk_drop_down_set_selected(GTK_DROP_DOWN(tool->mode_dropdown), GDIS_ISOSURFACE_MODE_MOLECULAR);
  gtk_widget_set_hexpand(tool->mode_dropdown, TRUE);
  g_signal_connect(tool->mode_dropdown,
                   "notify::selected",
                   G_CALLBACK(on_isosurface_mode_changed),
                   tool);
  gtk_grid_attach(GTK_GRID(grid), tool->mode_dropdown, 1, 0, 1, 1);

  label = gtk_label_new("Grid size");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
  tool->grid_spin = gtk_spin_button_new_with_range(0.12, 2.50, 0.02);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->grid_spin), 2);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->grid_spin), 0.40);
  gtk_grid_attach(GTK_GRID(grid), tool->grid_spin, 1, 1, 1, 1);

  label = gtk_label_new("Blur / smoothing");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 2, 1, 1);
  tool->blur_spin = gtk_spin_button_new_with_range(0.00, 2.00, 0.05);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->blur_spin), 2);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->blur_spin), 0.35);
  gtk_grid_attach(GTK_GRID(grid), tool->blur_spin, 1, 2, 1, 1);

  label = gtk_label_new("Isovalue");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 3, 1, 1);
  tool->value_spin = gtk_spin_button_new_with_range(0.001, 5.000, 0.005);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(tool->value_spin), 3);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->value_spin), 0.080);
  gtk_grid_attach(GTK_GRID(grid), tool->value_spin, 1, 3, 1, 1);

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

  button = gtk_button_new_with_label("Generate");
  tool->execute_button = button;
  g_signal_connect(button, "clicked", G_CALLBACK(on_isosurface_execute_clicked), tool);
  gtk_box_append(GTK_BOX(button_row), button);

  button = gtk_button_new_with_label("Clear");
  tool->clear_button = button;
  g_signal_connect(button, "clicked", G_CALLBACK(on_isosurface_clear_clicked), tool);
  gtk_box_append(GTK_BOX(button_row), button);

  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_close), window);
  gtk_box_append(GTK_BOX(button_row), button);

  g_signal_connect(window, "destroy", G_CALLBACK(on_isosurface_tool_destroy), tool);

  self->isosurface_tool = tool;
  gdis_gtk4_window_refresh_isosurface_tool(self);
  gtk_window_present(GTK_WINDOW(window));
}

static void
on_animation_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisAnimationTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->timer_id != 0u)
    g_source_remove(tool->timer_id);
  if (tool->owner)
    tool->owner->animation_tool = NULL;
  g_free(tool);
}

static void
gdis_gtk4_window_get_animation_source_label(GdisAnimationSourceType source,
                                            gchar *buffer,
                                            gsize buffer_len)
{
  g_return_if_fail(buffer != NULL);

  switch (source)
    {
    case GDIS_ANIMATION_SOURCE_MODEL_FRAMES:
      g_strlcpy(buffer, "active model frames", buffer_len);
      break;
    case GDIS_ANIMATION_SOURCE_SESSION_MODELS:
      g_strlcpy(buffer, "loaded model list", buffer_len);
      break;
    case GDIS_ANIMATION_SOURCE_NONE:
    default:
      g_strlcpy(buffer, "none", buffer_len);
      break;
    }
}

static gboolean
gdis_gtk4_window_animation_preserve_connectivity(const GdisGtk4Window *self)
{
  g_return_val_if_fail(self != NULL, FALSE);

  return self->animation_tool &&
         self->animation_tool->preserve_connectivity_toggle &&
         gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->animation_tool->preserve_connectivity_toggle));
}

static gboolean
gdis_gtk4_window_animation_preserve_scale(const GdisGtk4Window *self)
{
  g_return_val_if_fail(self != NULL, TRUE);

  return !self->animation_tool ||
         !self->animation_tool->preserve_scale_toggle ||
         gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->animation_tool->preserve_scale_toggle));
}

static GdisAnimationConfineMode
gdis_gtk4_window_get_animation_confine_mode(const GdisGtk4Window *self)
{
  g_return_val_if_fail(self != NULL, GDIS_ANIMATION_CONFINE_NONE);

  if (!self->animation_tool)
    return GDIS_ANIMATION_CONFINE_NONE;
  if (self->animation_tool->confine_atoms_toggle &&
      gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->animation_tool->confine_atoms_toggle)))
    return GDIS_ANIMATION_CONFINE_ATOMS;
  if (self->animation_tool->confine_molecules_toggle &&
      gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->animation_tool->confine_molecules_toggle)))
    return GDIS_ANIMATION_CONFINE_MOLECULES;
  return GDIS_ANIMATION_CONFINE_NONE;
}

static void
gdis_gtk4_window_apply_animation_processing(GdisGtk4Window *self)
{
  GError *error;

  g_return_if_fail(self != NULL);

  if (!self->active_model)
    return;

  error = NULL;
  switch (gdis_gtk4_window_get_animation_confine_mode(self))
    {
    case GDIS_ANIMATION_CONFINE_ATOMS:
      gdis_model_confine_atoms_to_cell(self->active_model, &error);
      break;
    case GDIS_ANIMATION_CONFINE_MOLECULES:
      gdis_model_confine_molecules_to_cell(self->active_model, &error);
      break;
    case GDIS_ANIMATION_CONFINE_NONE:
    default:
      break;
    }

  if (error)
    {
      gdis_gtk4_window_log(self, "Animation processing warning: %s\n", error->message);
      g_clear_error(&error);
    }
}

static GdisAnimationSourceType
gdis_gtk4_window_get_animation_source(GdisGtk4Window *self,
                                      guint *count_out,
                                      gint *active_index_out)
{
  guint count;
  gint active_index;
  GdisAnimationSourceType source;

  g_return_val_if_fail(self != NULL, GDIS_ANIMATION_SOURCE_NONE);

  source = GDIS_ANIMATION_SOURCE_NONE;
  count = 0u;
  active_index = -1;

  if (self->active_model)
    {
      count = gdis_model_get_frame_count(self->active_model);
      if (count > 1u)
        {
          source = GDIS_ANIMATION_SOURCE_MODEL_FRAMES;
          active_index = (gint) gdis_model_get_current_frame_index(self->active_model);
        }
    }

  if (source == GDIS_ANIMATION_SOURCE_NONE && self->models->len > 1u)
    {
      source = GDIS_ANIMATION_SOURCE_SESSION_MODELS;
      count = self->models->len;
      active_index = gdis_gtk4_window_get_model_index(self, self->active_model);
    }

  if (source == GDIS_ANIMATION_SOURCE_NONE && self->active_model)
    {
      count = 1u;
      active_index = 0;
    }

  if (count_out)
    *count_out = count;
  if (active_index_out)
    *active_index_out = active_index;
  return source;
}

static void
gdis_gtk4_window_animation_step_to_index(GdisGtk4Window *self,
                                         gint index,
                                         gboolean log_change)
{
  GdisAnimationSourceType source;
  GdisModel *model;
  gdouble rotation_x;
  gdouble rotation_y;
  gdouble zoom;
  guint source_count;
  GtkTextIter start;
  GtkTextIter end;
  gchar source_label[32];
  GError *error;
  GArray *saved_bonds;
  guint saved_explicit_bond_count;
  guint saved_atom_count;
  g_autofree gchar *status_text = NULL;

  g_return_if_fail(self != NULL);

  source = gdis_gtk4_window_get_animation_source(self, &source_count, NULL);
  if (source == GDIS_ANIMATION_SOURCE_NONE ||
      index < 0 ||
      (guint) index >= source_count)
    return;

  rotation_x = self->rotation_x;
  rotation_y = self->rotation_y;
  zoom = self->zoom;
  saved_bonds = NULL;
  saved_explicit_bond_count = 0u;
  saved_atom_count = self->active_model ? self->active_model->atom_count : 0u;
  if (self->status_buffer)
    {
      gtk_text_buffer_get_start_iter(self->status_buffer, &start);
      gtk_text_buffer_get_end_iter(self->status_buffer, &end);
      status_text = gtk_text_buffer_get_text(self->status_buffer, &start, &end, FALSE);
    }

  error = NULL;
  if (source == GDIS_ANIMATION_SOURCE_MODEL_FRAMES)
    {
      if (gdis_gtk4_window_animation_preserve_connectivity(self) &&
          self->active_model &&
          self->active_model->bonds)
        {
          saved_bonds = g_array_sized_new(FALSE,
                                          FALSE,
                                          sizeof(GdisBond),
                                          self->active_model->bonds->len);
          if (self->active_model->bonds->len > 0u)
            g_array_append_vals(saved_bonds,
                                self->active_model->bonds->data,
                                self->active_model->bonds->len);
          saved_explicit_bond_count = self->active_model->explicit_bond_count;
        }

      if (self->active_model &&
          (gint) gdis_model_get_current_frame_index(self->active_model) == index)
        {
          if (saved_bonds)
            g_array_free(saved_bonds, TRUE);
          gdis_gtk4_window_refresh_animation_tool(self);
          return;
        }

      if (!self->active_model || !gdis_model_set_frame_index(self->active_model, (guint) index, &error))
        {
          if (saved_bonds)
            g_array_free(saved_bonds, TRUE);
          gdis_gtk4_window_log(self, "Animation frame change failed: %s\n",
                               error ? error->message : "unknown error");
          g_clear_error(&error);
          return;
        }
      model = self->active_model;
      gdis_gtk4_window_apply_animation_processing(self);
      if (saved_bonds)
        {
          if (model->atom_count == saved_atom_count)
            {
              if (model->bonds)
                g_array_free(model->bonds, TRUE);
              model->bonds = saved_bonds;
              saved_bonds = NULL;
              model->explicit_bond_count = MIN(saved_explicit_bond_count, model->bonds->len);
              model->bond_count = model->bonds->len;
            }
          else
            {
              g_array_free(saved_bonds, TRUE);
              saved_bonds = NULL;
            }
        }
      gdis_gtk4_window_update_details(self);
    }
  else
    {
      model = g_ptr_array_index(self->models, index);
      if (self->active_model == model)
        {
          if (saved_bonds)
            g_array_free(saved_bonds, TRUE);
          gdis_gtk4_window_refresh_animation_tool(self);
          return;
        }
      gdis_gtk4_window_set_active_model(self, model);
      gdis_gtk4_window_apply_animation_processing(self);
    }

  if (gdis_gtk4_window_animation_preserve_scale(self))
    {
      self->rotation_x = rotation_x;
      self->rotation_y = rotation_y;
      self->zoom = zoom;
    }
  gdis_gtk4_window_refresh_viewer(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_refresh_zmatrix_tool(self);

  if (self->status_buffer && status_text)
    {
      gtk_text_buffer_set_text(self->status_buffer, status_text, -1);
      if (log_change)
        {
          gdis_gtk4_window_get_animation_source_label(source, source_label, sizeof(source_label));
          gdis_gtk4_window_log(self,
                               "Animation %s: %d / %u (%s)\n",
                               source == GDIS_ANIMATION_SOURCE_MODEL_FRAMES ? "frame" : "step",
                               index + 1,
                               source_count,
                               source_label);
        }
    }

  gdis_gtk4_window_refresh_animation_tool(self);
  gdis_gtk4_window_refresh_recording_tool(self);
}

static gboolean
on_animation_tick(gpointer user_data)
{
  GdisAnimationTool *tool;
  GdisGtk4Window *self;
  GdisAnimationSourceType source;
  gint current_index;
  gint step;
  gint next_index;
  guint source_count;
  gboolean loop;

  tool = user_data;
  self = tool ? tool->owner : NULL;
  if (!tool || !self || !tool->window)
    return G_SOURCE_REMOVE;

  source = gdis_gtk4_window_get_animation_source(self, &source_count, &current_index);
  if (source == GDIS_ANIMATION_SOURCE_NONE || source_count < 2u)
    {
      tool->timer_id = 0u;
      gdis_gtk4_window_refresh_animation_tool(self);
      return G_SOURCE_REMOVE;
    }

  step = (gint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->step_spin));
  loop = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(tool->loop_toggle));
  next_index = current_index + MAX(step, 1);
  if (next_index >= (gint) source_count)
    {
      if (!loop)
        {
          tool->timer_id = 0u;
          gdis_gtk4_window_refresh_animation_tool(self);
          return G_SOURCE_REMOVE;
        }
      next_index %= (gint) source_count;
    }

  gdis_gtk4_window_animation_step_to_index(self, next_index, FALSE);
  return G_SOURCE_CONTINUE;
}

static void
gdis_gtk4_window_refresh_animation_tool(GdisGtk4Window *self)
{
  GdisAnimationTool *tool;
  GdisAnimationSourceType source;
  guint source_count;
  gint active_index;
  gchar source_label[32];
  const char *confine_label;
  gboolean can_play;
  g_autofree gchar *summary = NULL;
  g_autofree gchar *processing = NULL;
  g_autofree gchar *rendering = NULL;
  g_autofree gchar *ffmpeg_path = NULL;

  g_return_if_fail(self != NULL);

  tool = self->animation_tool;
  if (!tool)
    return;

  source = gdis_gtk4_window_get_animation_source(self, &source_count, &active_index);
  gdis_gtk4_window_get_animation_source_label(source, source_label, sizeof(source_label));
  can_play = (source_count > 1u);
  switch (gdis_gtk4_window_get_animation_confine_mode(self))
    {
    case GDIS_ANIMATION_CONFINE_ATOMS:
      confine_label = "Confine atoms to PBC";
      break;
    case GDIS_ANIMATION_CONFINE_MOLECULES:
      confine_label = "Confine molecules to PBC";
      break;
    case GDIS_ANIMATION_CONFINE_NONE:
    default:
      confine_label = "Cell confinement off";
      break;
    }

  if (tool->frame_scale)
    {
      tool->suppress_scale_signal = TRUE;
      gtk_range_set_range(GTK_RANGE(tool->frame_scale),
                          1.0,
                          MAX((gdouble) source_count, 1.0));
      gtk_range_set_value(GTK_RANGE(tool->frame_scale),
                          active_index >= 0 ? (gdouble) active_index + 1.0 : 1.0);
      gtk_widget_set_sensitive(tool->frame_scale, source_count > 0u);
      tool->suppress_scale_signal = FALSE;
    }

  if (tool->play_button)
    {
      gtk_widget_set_sensitive(tool->play_button, can_play);
      gtk_button_set_label(GTK_BUTTON(tool->play_button),
                           tool->timer_id != 0u ? "Pause" : "Play");
    }

  summary = g_strdup_printf(
    "Animation source: %s\n"
    "Loaded models: %u\n"
    "Available steps: %u\n"
    "Current frame: %d / %u\n"
    "Active model: %s\n\n"
    "Delay: %d ms   Step size: %d   Loop: %s\n\n"
    "%s",
    source_label,
    self->models->len,
    source_count,
    active_index >= 0 ? active_index + 1 : 0,
    source_count,
    self->active_model ? self->active_model->basename : "none",
    tool->delay_spin ? gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->delay_spin)) : 0,
    tool->step_spin ? gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->step_spin)) : 1,
    (tool->loop_toggle && gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(tool->loop_toggle))) ? "on" : "off",
    source == GDIS_ANIMATION_SOURCE_MODEL_FRAMES
      ? "GTK4 animation is now driving the active model's internal frame set."
      : (source == GDIS_ANIMATION_SOURCE_SESSION_MODELS
         ? "GTK4 animation is using the loaded model list because the active model is single-frame."
         : "Load a trajectory-capable model or multiple models to enable playback."));
  gtk_label_set_text(tool->summary_label, summary);

  processing = g_strdup_printf(
    "Per-cycle processing\n"
    "  Connectivity: %s\n"
    "  View scale/rotation: %s\n"
    "  Cell handling: %s\n\n"
    "These controls mirror the old Linux dialog. Connectivity preservation is applied when stepping active-model frames.",
    gdis_gtk4_window_animation_preserve_connectivity(self) ? "preserve current bonding" : "recompute as frame changes",
    gdis_gtk4_window_animation_preserve_scale(self) ? "preserve current zoom/rotation" : "allow per-step reset",
    confine_label);
  if (tool->processing_label)
    gtk_label_set_text(tool->processing_label, processing);

  ffmpeg_path = g_find_program_in_path("ffmpeg");
  rendering = g_strdup_printf(
    "Rendering / recording\n"
    "  Snapshot and PNG sequence export are available now.\n"
    "  Movie backend: %s\n"
    "  Use Record... to open the GTK4 rendering/export panel with MP4/GIF controls.",
    ffmpeg_path ? ffmpeg_path : "ffmpeg not found");
  if (tool->rendering_label)
    gtk_label_set_text(tool->rendering_label, rendering);
}

static void
on_animation_play_clicked(GtkButton *button, gpointer user_data)
{
  GdisAnimationTool *tool;
  guint interval_ms;

  (void) button;

  tool = user_data;
  if (!tool || !tool->owner)
    return;

  if (tool->timer_id != 0u)
    {
      g_source_remove(tool->timer_id);
      tool->timer_id = 0u;
      gdis_gtk4_window_refresh_animation_tool(tool->owner);
      return;
    }

  interval_ms = (guint) MAX(80, gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->delay_spin)));
  tool->timer_id = g_timeout_add(interval_ms, on_animation_tick, tool);
  gdis_gtk4_window_refresh_animation_tool(tool->owner);
}

static void
on_animation_prev_clicked(GtkButton *button, gpointer user_data)
{
  GdisAnimationTool *tool;
  GdisAnimationSourceType source;
  guint source_count;
  gint active_index;
  gint step;
  gint next_index;

  (void) button;

  tool = user_data;
  if (!tool || !tool->owner)
    return;

  source = gdis_gtk4_window_get_animation_source(tool->owner, &source_count, &active_index);
  if (source == GDIS_ANIMATION_SOURCE_NONE || source_count == 0u)
    return;
  if (active_index < 0)
    active_index = 0;
  step = (gint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->step_spin));
  next_index = active_index - MAX(step, 1);
  while (next_index < 0)
    next_index += (gint) source_count;
  gdis_gtk4_window_animation_step_to_index(tool->owner, next_index, TRUE);
}

static void
on_animation_next_clicked(GtkButton *button, gpointer user_data)
{
  GdisAnimationTool *tool;

  (void) button;

  tool = user_data;
  if (!tool || !tool->owner)
    return;

  on_animation_tick(tool);
}

static void
on_animation_scale_value_changed(GtkRange *range, gpointer user_data)
{
  GdisAnimationTool *tool;
  gint index;

  tool = user_data;
  if (!tool || !tool->owner || tool->suppress_scale_signal)
    return;

  index = (gint) llround(gtk_range_get_value(range)) - 1;
  gdis_gtk4_window_animation_step_to_index(tool->owner, index, TRUE);
}

static void
on_animation_processing_changed(GtkWidget *button, gpointer user_data)
{
  GdisAnimationTool *tool;

  (void) button;

  tool = user_data;
  if (!tool || !tool->owner)
    return;

  gdis_gtk4_window_refresh_animation_tool(tool->owner);
}

static void
gdis_gtk4_window_present_animation_tool(GdisGtk4Window *self)
{
  GdisAnimationTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *controls;
  GtkWidget *button;
  GtkWidget *label;
  GtkWidget *frame;
  GtkWidget *top_frame;
  GtkWidget *switcher;
  GtkWidget *stack;
  GtkWidget *page;
  GtkWidget *section;
  GtkWidget *grid;
  GtkWidget *transport_frame;
  GtkWidget *actions;
  GtkWidget *record_row;
  GtkWidget *toggle;
  GtkCheckButton *group_head;

  g_return_if_fail(self != NULL);

  if (self->animation_tool && GTK_IS_WINDOW(self->animation_tool->window))
    {
      gdis_gtk4_window_refresh_animation_tool(self);
      gtk_window_present(GTK_WINDOW(self->animation_tool->window));
      return;
    }

  tool = g_new0(GdisAnimationTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Animation");
  gtk_window_set_default_size(GTK_WINDOW(window), 720, 560);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  top_frame = gtk_frame_new(NULL);
  gtk_box_append(GTK_BOX(root), top_frame);

  page = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_widget_set_margin_start(page, 10);
  gtk_widget_set_margin_end(page, 10);
  gtk_widget_set_margin_top(page, 10);
  gtk_widget_set_margin_bottom(page, 10);
  gtk_frame_set_child(GTK_FRAME(top_frame), page);

  switcher = gtk_stack_switcher_new();
  gtk_widget_set_halign(switcher, GTK_ALIGN_START);
  gtk_box_append(GTK_BOX(page), switcher);

  stack = gtk_stack_new();
  gtk_stack_set_transition_type(GTK_STACK(stack), GTK_STACK_TRANSITION_TYPE_SLIDE_LEFT_RIGHT);
  gtk_stack_set_hhomogeneous(GTK_STACK(stack), FALSE);
  gtk_stack_set_vhomogeneous(GTK_STACK(stack), FALSE);
  gtk_widget_set_vexpand(stack, TRUE);
  gtk_stack_switcher_set_stack(GTK_STACK_SWITCHER(switcher), GTK_STACK(stack));
  gtk_box_append(GTK_BOX(page), stack);

  page = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
  section = gtk_frame_new("Playback");
  gtk_box_append(GTK_BOX(page), section);
  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 10);
  gtk_widget_set_margin_start(grid, 10);
  gtk_widget_set_margin_end(grid, 10);
  gtk_widget_set_margin_top(grid, 10);
  gtk_widget_set_margin_bottom(grid, 10);
  gtk_frame_set_child(GTK_FRAME(section), grid);

  label = gtk_label_new("Delay (ms)");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  tool->delay_spin = gtk_spin_button_new_with_range(80.0, 5000.0, 20.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->delay_spin), 350.0);
  gtk_grid_attach(GTK_GRID(grid), tool->delay_spin, 1, 0, 1, 1);

  label = gtk_label_new("Step size");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 2, 0, 1, 1);
  tool->step_spin = gtk_spin_button_new_with_range(1.0, 24.0, 1.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->step_spin), 1.0);
  gtk_grid_attach(GTK_GRID(grid), tool->step_spin, 3, 0, 1, 1);

  tool->loop_toggle = gtk_check_button_new_with_label("Loop");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->loop_toggle), TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->loop_toggle, 0, 1, 2, 1);

  tool->preserve_connectivity_toggle = gtk_check_button_new_with_label("Don't recalculate connectivity");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->preserve_connectivity_toggle), TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->preserve_connectivity_toggle, 2, 1, 2, 1);

  tool->preserve_scale_toggle = gtk_check_button_new_with_label("Don't recalculate scale");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->preserve_scale_toggle), TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->preserve_scale_toggle, 0, 2, 2, 1);

  frame = gtk_frame_new("Control summary");
  gtk_box_append(GTK_BOX(page), frame);
  tool->summary_label = GTK_LABEL(gtk_label_new(""));
  gtk_label_set_wrap(tool->summary_label, TRUE);
  gtk_label_set_xalign(tool->summary_label, 0.0f);
  gtk_widget_set_margin_start(GTK_WIDGET(tool->summary_label), 10);
  gtk_widget_set_margin_end(GTK_WIDGET(tool->summary_label), 10);
  gtk_widget_set_margin_top(GTK_WIDGET(tool->summary_label), 10);
  gtk_widget_set_margin_bottom(GTK_WIDGET(tool->summary_label), 10);
  gtk_frame_set_child(GTK_FRAME(frame), GTK_WIDGET(tool->summary_label));
  gtk_stack_add_titled(GTK_STACK(stack), page, "control", "Control");

  page = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
  section = gtk_frame_new("Per-cycle processing");
  gtk_box_append(GTK_BOX(page), section);
  controls = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_widget_set_margin_start(controls, 10);
  gtk_widget_set_margin_end(controls, 10);
  gtk_widget_set_margin_top(controls, 10);
  gtk_widget_set_margin_bottom(controls, 10);
  gtk_frame_set_child(GTK_FRAME(section), controls);

  group_head = GTK_CHECK_BUTTON(gtk_check_button_new_with_label("Cell confinement off"));
  tool->confine_none_toggle = GTK_WIDGET(group_head);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->confine_none_toggle), TRUE);
  gtk_box_append(GTK_BOX(controls), tool->confine_none_toggle);

  toggle = gtk_check_button_new_with_label("Confine atoms to PBC");
  gtk_check_button_set_group(GTK_CHECK_BUTTON(toggle), group_head);
  tool->confine_atoms_toggle = toggle;
  gtk_box_append(GTK_BOX(controls), tool->confine_atoms_toggle);

  toggle = gtk_check_button_new_with_label("Confine molecules to PBC");
  gtk_check_button_set_group(GTK_CHECK_BUTTON(toggle), group_head);
  tool->confine_molecules_toggle = toggle;
  gtk_box_append(GTK_BOX(controls), tool->confine_molecules_toggle);

  frame = gtk_frame_new("Processing summary");
  gtk_box_append(GTK_BOX(page), frame);
  tool->processing_label = GTK_LABEL(gtk_label_new(""));
  gtk_label_set_wrap(tool->processing_label, TRUE);
  gtk_label_set_xalign(tool->processing_label, 0.0f);
  gtk_widget_set_margin_start(GTK_WIDGET(tool->processing_label), 10);
  gtk_widget_set_margin_end(GTK_WIDGET(tool->processing_label), 10);
  gtk_widget_set_margin_top(GTK_WIDGET(tool->processing_label), 10);
  gtk_widget_set_margin_bottom(GTK_WIDGET(tool->processing_label), 10);
  gtk_frame_set_child(GTK_FRAME(frame), GTK_WIDGET(tool->processing_label));
  gtk_stack_add_titled(GTK_STACK(stack), page, "processing", "Processing");

  page = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
  frame = gtk_frame_new("Rendering");
  gtk_box_append(GTK_BOX(page), frame);
  controls = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
  gtk_widget_set_margin_start(controls, 10);
  gtk_widget_set_margin_end(controls, 10);
  gtk_widget_set_margin_top(controls, 10);
  gtk_widget_set_margin_bottom(controls, 10);
  gtk_frame_set_child(GTK_FRAME(frame), controls);

  tool->rendering_label = GTK_LABEL(gtk_label_new(""));
  gtk_label_set_wrap(tool->rendering_label, TRUE);
  gtk_label_set_xalign(tool->rendering_label, 0.0f);
  gtk_box_append(GTK_BOX(controls), GTK_WIDGET(tool->rendering_label));

  record_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_box_append(GTK_BOX(controls), record_row);

  button = gtk_button_new_with_label("Open Export Panel...");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gdis_gtk4_window_present_recording_tool), self);
  gtk_box_append(GTK_BOX(record_row), button);
  gtk_stack_add_titled(GTK_STACK(stack), page, "rendering", "Rendering");

  transport_frame = gtk_frame_new(NULL);
  gtk_box_append(GTK_BOX(root), transport_frame);
  page = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
  gtk_widget_set_margin_start(page, 10);
  gtk_widget_set_margin_end(page, 10);
  gtk_widget_set_margin_top(page, 10);
  gtk_widget_set_margin_bottom(page, 10);
  gtk_frame_set_child(GTK_FRAME(transport_frame), page);

  label = gtk_label_new("Frame");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_box_append(GTK_BOX(page), label);

  tool->frame_scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, 1.0, 1.0, 1.0);
  gtk_scale_set_digits(GTK_SCALE(tool->frame_scale), 0);
  gtk_scale_set_draw_value(GTK_SCALE(tool->frame_scale), TRUE);
  g_signal_connect(tool->frame_scale,
                   "value-changed",
                   G_CALLBACK(on_animation_scale_value_changed),
                   tool);
  gtk_box_append(GTK_BOX(page), tool->frame_scale);

  controls = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_box_append(GTK_BOX(page), controls);

  button = gtk_button_new_with_label("Previous");
  g_signal_connect(button, "clicked", G_CALLBACK(on_animation_prev_clicked), tool);
  gtk_box_append(GTK_BOX(controls), button);

  button = gtk_button_new_with_label("Play");
  tool->play_button = button;
  g_signal_connect(button, "clicked", G_CALLBACK(on_animation_play_clicked), tool);
  gtk_box_append(GTK_BOX(controls), button);

  button = gtk_button_new_with_label("Next");
  g_signal_connect(button, "clicked", G_CALLBACK(on_animation_next_clicked), tool);
  gtk_box_append(GTK_BOX(controls), button);

  actions = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_widget_set_halign(actions, GTK_ALIGN_END);
  gtk_box_append(GTK_BOX(root), actions);

  button = gtk_button_new_with_label("Record...");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gdis_gtk4_window_present_recording_tool), self);
  gtk_box_append(GTK_BOX(actions), button);

  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_destroy), window);
  gtk_box_append(GTK_BOX(actions), button);

  g_signal_connect(tool->loop_toggle, "toggled", G_CALLBACK(on_animation_processing_changed), tool);
  g_signal_connect(tool->delay_spin, "value-changed", G_CALLBACK(on_animation_processing_changed), tool);
  g_signal_connect(tool->step_spin, "value-changed", G_CALLBACK(on_animation_processing_changed), tool);
  g_signal_connect(tool->preserve_connectivity_toggle, "toggled", G_CALLBACK(on_animation_processing_changed), tool);
  g_signal_connect(tool->preserve_scale_toggle, "toggled", G_CALLBACK(on_animation_processing_changed), tool);
  g_signal_connect(tool->confine_none_toggle, "toggled", G_CALLBACK(on_animation_processing_changed), tool);
  g_signal_connect(tool->confine_atoms_toggle, "toggled", G_CALLBACK(on_animation_processing_changed), tool);
  g_signal_connect(tool->confine_molecules_toggle, "toggled", G_CALLBACK(on_animation_processing_changed), tool);

  g_signal_connect(window, "destroy", G_CALLBACK(on_animation_tool_destroy), tool);

  self->animation_tool = tool;
  gdis_gtk4_window_refresh_animation_tool(self);
  gtk_window_present(GTK_WINDOW(window));
}

static void
on_recording_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisRecordingTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->recording_tool = NULL;
  g_free(tool);
}

static void
gdis_gtk4_window_refresh_recording_tool(GdisGtk4Window *self)
{
  GdisRecordingTool *tool;
  GdisAnimationSourceType source;
  guint source_count;
  gint active_index;
  gchar source_label[32];
  g_autofree gchar *summary = NULL;
  g_autofree gchar *default_output = NULL;
  g_autofree gchar *default_prefix = NULL;
  g_autofree gchar *ffmpeg_path = NULL;
  const char *current_output;
  const char *current_prefix;
  const char *format_label;
  gboolean movie_enabled;
  gboolean ffmpeg_available;
  guint fps;
  gint viewer_width;
  gint viewer_height;

  g_return_if_fail(self != NULL);

  tool = self->recording_tool;
  if (!tool)
    return;

  source = gdis_gtk4_window_get_animation_source(self, &source_count, &active_index);
  gdis_gtk4_window_get_animation_source_label(source, source_label, sizeof(source_label));
  viewer_width = self->viewer_area ? gtk_widget_get_width(self->viewer_area) : 1280;
  viewer_height = self->viewer_area ? gtk_widget_get_height(self->viewer_area) : 900;

  current_output = tool->output_entry ? gtk_editable_get_text(GTK_EDITABLE(tool->output_entry)) : NULL;
  if (!current_output || current_output[0] == '\0')
    {
      default_output = g_build_filename(g_get_current_dir(), "gdis_exports", NULL);
      gtk_editable_set_text(GTK_EDITABLE(tool->output_entry), default_output);
    }

  current_prefix = tool->prefix_entry ? gtk_editable_get_text(GTK_EDITABLE(tool->prefix_entry)) : NULL;
  if ((!current_prefix || current_prefix[0] == '\0') && self->active_model)
    {
      default_prefix = g_strdup_printf("%s_capture",
                                       self->active_model->basename ? self->active_model->basename : "gdis");
      gtk_editable_set_text(GTK_EDITABLE(tool->prefix_entry), default_prefix);
    }

  if (tool->width_spin && gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->width_spin)) <= 1)
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->width_spin), MAX(viewer_width, 640));
  if (tool->height_spin && gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->height_spin)) <= 1)
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->height_spin), MAX(viewer_height, 480));
  if (tool->fps_spin && gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->fps_spin)) <= 0)
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->fps_spin), 12.0);

  ffmpeg_path = g_find_program_in_path("ffmpeg");
  ffmpeg_available = (ffmpeg_path != NULL);
  movie_enabled = tool->movie_enable_toggle &&
                  gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(tool->movie_enable_toggle));
  format_label = (tool->movie_gif_toggle &&
                  gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(tool->movie_gif_toggle)))
                   ? "Animated GIF"
                   : "MP4 (H.264)";
  fps = tool->fps_spin ? (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->fps_spin)) : 12u;

  if (tool->movie_mp4_toggle)
    gtk_widget_set_sensitive(tool->movie_mp4_toggle, ffmpeg_available && movie_enabled);
  if (tool->movie_gif_toggle)
    gtk_widget_set_sensitive(tool->movie_gif_toggle, ffmpeg_available && movie_enabled);
  if (tool->fps_spin)
    gtk_widget_set_sensitive(tool->fps_spin, ffmpeg_available && movie_enabled);
  if (tool->keep_frames_toggle)
    gtk_widget_set_sensitive(tool->keep_frames_toggle, ffmpeg_available && movie_enabled);
  if (tool->movie_button)
    gtk_widget_set_sensitive(tool->movie_button, ffmpeg_available && movie_enabled);

  summary = g_strdup_printf(
    "Recording source: %s\n"
    "Available steps: %u\n"
    "Current step: %d / %u\n"
    "Active model: %s\n"
    "Viewer size: %d x %d\n"
    "Movie backend: %s\n"
    "Create movie: %s\n"
    "Movie format: %s at %u fps\n\n"
    "Snapshot exports the current view.\n"
    "Sequence exports one PNG per animation step, using active model frames when available.\n"
    "Movie export assembles those frames with ffmpeg%s.",
    source_label,
    source_count,
    active_index >= 0 ? active_index + 1 : 0,
    source_count,
    self->active_model ? self->active_model->basename : "none",
    viewer_width,
    viewer_height,
    ffmpeg_path ? ffmpeg_path : "not available",
    movie_enabled ? "on" : "off",
    format_label,
    MAX(fps, 1u),
    ffmpeg_path ? "" : " once ffmpeg is installed");
  gtk_label_set_text(tool->summary_label, summary);
}

static void
on_recording_settings_changed(GtkWidget *widget, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) widget;

  self = user_data;
  if (!self)
    return;

  gdis_gtk4_window_refresh_recording_tool(self);
}

static void
on_recording_snapshot_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  const char *output_dir;
  const char *prefix;
  gint width;
  gint height;
  g_autofree gchar *basename = NULL;
  gchar *path;
  GError *error;

  (void) button;

  self = user_data;
  if (!self || !self->recording_tool)
    return;

  output_dir = gtk_editable_get_text(GTK_EDITABLE(self->recording_tool->output_entry));
  prefix = gtk_editable_get_text(GTK_EDITABLE(self->recording_tool->prefix_entry));
  width = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->recording_tool->width_spin));
  height = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->recording_tool->height_spin));

  if (!output_dir || output_dir[0] == '\0')
    {
      gdis_gtk4_window_log(self, "Recording snapshot failed: choose an output directory first.\n");
      return;
    }

  if (g_mkdir_with_parents(output_dir, 0755) != 0)
    {
      gdis_gtk4_window_log(self, "Recording snapshot failed: could not create %s.\n", output_dir);
      return;
    }

  basename = (prefix && prefix[0] != '\0') ? g_strdup_printf("%s_snapshot.png", prefix)
                                           : g_strdup("gdis_snapshot.png");
  path = g_build_filename(output_dir, basename, NULL);
  error = NULL;
  if (!gdis_gtk4_window_export_view_png(self, path, width, height, &error))
    {
      gdis_gtk4_window_log(self, "Recording snapshot failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      g_free(path);
      return;
    }

  gdis_gtk4_window_log(self, "Saved snapshot: %s\n", path);
  g_free(path);
}

static void
on_recording_sequence_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  gint width;
  gint height;
  const char *output_dir;
  const char *prefix;
  GError *error;
  guint frame_count;

  (void) button;

  self = user_data;
  if (!self || !self->recording_tool)
    return;

  output_dir = gtk_editable_get_text(GTK_EDITABLE(self->recording_tool->output_entry));
  prefix = gtk_editable_get_text(GTK_EDITABLE(self->recording_tool->prefix_entry));
  width = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->recording_tool->width_spin));
  height = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->recording_tool->height_spin));
  error = NULL;
  if (!gdis_gtk4_window_export_animation_sequence(self,
                                                  output_dir,
                                                  prefix,
                                                  width,
                                                  height,
                                                  &frame_count,
                                                  &error))
    {
      gdis_gtk4_window_log(self, "Recording sequence failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_log(self,
                       "Saved PNG sequence to %s using %u animation step%s.\n",
                       output_dir,
                       frame_count,
                       frame_count == 1u ? "" : "s");
  gdis_gtk4_window_refresh_recording_tool(self);
}

static void
on_recording_movie_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  const char *output_dir;
  const char *prefix;
  gint width;
  gint height;
  guint fps;
  guint format_index;
  gboolean keep_frames;
  guint frame_count;
  g_autofree gchar *movie_path = NULL;
  GError *error;

  (void) button;

  self = user_data;
  if (!self || !self->recording_tool)
    return;

  output_dir = gtk_editable_get_text(GTK_EDITABLE(self->recording_tool->output_entry));
  prefix = gtk_editable_get_text(GTK_EDITABLE(self->recording_tool->prefix_entry));
  width = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->recording_tool->width_spin));
  height = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->recording_tool->height_spin));
  fps = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(self->recording_tool->fps_spin));
  format_index = (self->recording_tool->movie_gif_toggle &&
                  gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->recording_tool->movie_gif_toggle)))
                   ? 1u
                   : 0u;
  keep_frames = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->recording_tool->keep_frames_toggle));

  if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(self->recording_tool->movie_enable_toggle)))
    {
      gdis_gtk4_window_log(self, "Recording movie failed: enable Create movie first.\n");
      return;
    }

  error = NULL;
  if (!gdis_gtk4_window_export_movie(self,
                                     output_dir,
                                     prefix,
                                     width,
                                     height,
                                     MAX(fps, 1u),
                                     format_index,
                                     keep_frames,
                                     &movie_path,
                                     &frame_count,
                                     &error))
    {
      gdis_gtk4_window_log(self, "Recording movie failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_log(self,
                       "Created movie %s from %u animation step%s.\n",
                       movie_path,
                       frame_count,
                       frame_count == 1u ? "" : "s");
  gdis_gtk4_window_refresh_recording_tool(self);
}

static void
gdis_gtk4_window_present_recording_tool(GdisGtk4Window *self)
{
  GdisRecordingTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *grid;
  GtkWidget *label;
  GtkWidget *button_row;
  GtkWidget *button;
  GtkWidget *frame;
  GtkWidget *section;
  GtkWidget *controls;
  GtkCheckButton *group_head;

  g_return_if_fail(self != NULL);

  if (self->recording_tool && GTK_IS_WINDOW(self->recording_tool->window))
    {
      gdis_gtk4_window_refresh_recording_tool(self);
      gtk_window_present(GTK_WINDOW(self->recording_tool->window));
      return;
    }

  tool = g_new0(GdisRecordingTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Recording");
  gtk_window_set_default_size(GTK_WINDOW(window), 700, 460);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  frame = gtk_frame_new("Output");
  gtk_box_append(GTK_BOX(root), frame);
  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_widget_set_margin_start(grid, 10);
  gtk_widget_set_margin_end(grid, 10);
  gtk_widget_set_margin_top(grid, 10);
  gtk_widget_set_margin_bottom(grid, 10);
  gtk_frame_set_child(GTK_FRAME(frame), grid);

  label = gtk_label_new("Output directory");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  tool->output_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->output_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->output_entry, 1, 0, 3, 1);

  label = gtk_label_new("Prefix");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
  tool->prefix_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->prefix_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->prefix_entry, 1, 1, 3, 1);

  label = gtk_label_new("Width");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 2, 1, 1);
  tool->width_spin = gtk_spin_button_new_with_range(320.0, 4096.0, 16.0);
  gtk_grid_attach(GTK_GRID(grid), tool->width_spin, 1, 2, 1, 1);

  label = gtk_label_new("Height");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 2, 2, 1, 1);
  tool->height_spin = gtk_spin_button_new_with_range(240.0, 4096.0, 16.0);
  gtk_grid_attach(GTK_GRID(grid), tool->height_spin, 3, 2, 1, 1);

  frame = gtk_frame_new("Rendering");
  gtk_box_append(GTK_BOX(root), frame);
  section = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_widget_set_margin_start(section, 10);
  gtk_widget_set_margin_end(section, 10);
  gtk_widget_set_margin_top(section, 10);
  gtk_widget_set_margin_bottom(section, 10);
  gtk_frame_set_child(GTK_FRAME(frame), section);

  tool->movie_enable_toggle = gtk_check_button_new_with_label("Create movie");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->movie_enable_toggle), TRUE);
  gtk_box_append(GTK_BOX(section), tool->movie_enable_toggle);

  controls = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
  gtk_widget_set_margin_start(controls, 16);
  gtk_box_append(GTK_BOX(section), controls);

  group_head = GTK_CHECK_BUTTON(gtk_check_button_new_with_label("MP4 (H.264)"));
  tool->movie_mp4_toggle = GTK_WIDGET(group_head);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->movie_mp4_toggle), TRUE);
  gtk_box_append(GTK_BOX(controls), tool->movie_mp4_toggle);

  tool->movie_gif_toggle = gtk_check_button_new_with_label("Animated GIF");
  gtk_check_button_set_group(GTK_CHECK_BUTTON(tool->movie_gif_toggle), group_head);
  gtk_box_append(GTK_BOX(controls), tool->movie_gif_toggle);

  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_box_append(GTK_BOX(controls), grid);

  label = gtk_label_new("Delay / fps");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  tool->fps_spin = gtk_spin_button_new_with_range(1.0, 60.0, 1.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->fps_spin), 12.0);
  gtk_grid_attach(GTK_GRID(grid), tool->fps_spin, 1, 0, 1, 1);

  tool->keep_frames_toggle = gtk_check_button_new_with_label("Keep PNG frames after movie export");
  gtk_grid_attach(GTK_GRID(grid), tool->keep_frames_toggle, 0, 1, 2, 1);

  frame = gtk_frame_new("Export summary");
  gtk_box_append(GTK_BOX(root), frame);
  tool->summary_label = GTK_LABEL(gtk_label_new(""));
  gtk_label_set_wrap(tool->summary_label, TRUE);
  gtk_label_set_xalign(tool->summary_label, 0.0f);
  gtk_widget_set_margin_start(GTK_WIDGET(tool->summary_label), 10);
  gtk_widget_set_margin_end(GTK_WIDGET(tool->summary_label), 10);
  gtk_widget_set_margin_top(GTK_WIDGET(tool->summary_label), 10);
  gtk_widget_set_margin_bottom(GTK_WIDGET(tool->summary_label), 10);
  gtk_frame_set_child(GTK_FRAME(frame), GTK_WIDGET(tool->summary_label));

  button_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_widget_set_halign(button_row, GTK_ALIGN_END);
  gtk_box_append(GTK_BOX(root), button_row);

  button = gtk_button_new_with_label("Snapshot");
  g_signal_connect(button, "clicked", G_CALLBACK(on_recording_snapshot_clicked), self);
  gtk_box_append(GTK_BOX(button_row), button);

  button = gtk_button_new_with_label("Sequence");
  g_signal_connect(button, "clicked", G_CALLBACK(on_recording_sequence_clicked), self);
  gtk_box_append(GTK_BOX(button_row), button);

  button = gtk_button_new_with_label("Movie");
  tool->movie_button = button;
  g_signal_connect(button, "clicked", G_CALLBACK(on_recording_movie_clicked), self);
  gtk_box_append(GTK_BOX(button_row), button);

  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_destroy), window);
  gtk_box_append(GTK_BOX(button_row), button);

  g_signal_connect(tool->movie_enable_toggle, "toggled", G_CALLBACK(on_recording_settings_changed), self);
  g_signal_connect(tool->movie_mp4_toggle, "toggled", G_CALLBACK(on_recording_settings_changed), self);
  g_signal_connect(tool->movie_gif_toggle, "toggled", G_CALLBACK(on_recording_settings_changed), self);
  g_signal_connect(tool->fps_spin, "value-changed", G_CALLBACK(on_recording_settings_changed), self);
  g_signal_connect(tool->keep_frames_toggle, "toggled", G_CALLBACK(on_recording_settings_changed), self);

  g_signal_connect(window, "destroy", G_CALLBACK(on_recording_tool_destroy), tool);

  self->recording_tool = tool;
  gdis_gtk4_window_refresh_recording_tool(self);
  gtk_window_present(GTK_WINDOW(window));
}

static void
on_zmatrix_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisZmatrixTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  gdis_zmatrix_tool_clear_working_set(tool);
  if (tool->owner)
    tool->owner->zmatrix_tool = NULL;
  g_free(tool);
}

static void
gdis_zmatrix_tool_clear_working_set(GdisZmatrixTool *tool)
{
  g_return_if_fail(tool != NULL);

  g_clear_pointer(&tool->scope, gdis_uint_array_free);
  if (tool->rows)
    {
      g_array_free(tool->rows, TRUE);
      tool->rows = NULL;
    }
}

static void
gdis_gtk4_window_refresh_zmatrix_row_controls(GdisGtk4Window *self)
{
  GdisZmatrixTool *tool;
  gint row_index;

  g_return_if_fail(self != NULL);

  tool = self->zmatrix_tool;
  if (!tool || !tool->row_spin || !tool->distance_entry || !tool->angle_entry ||
      !tool->torsion_entry || !tool->row_label)
    return;

  tool->suppress_row_signal = TRUE;
  if (!tool->rows || tool->rows->len == 0u || !self->active_model)
    {
      gtk_label_set_text(tool->row_label, "No editable Z-matrix rows are available.");
      gtk_editable_set_text(GTK_EDITABLE(tool->distance_entry), "");
      gtk_editable_set_text(GTK_EDITABLE(tool->angle_entry), "");
      gtk_editable_set_text(GTK_EDITABLE(tool->torsion_entry), "");
      gtk_widget_set_sensitive(tool->row_spin, FALSE);
      gtk_widget_set_sensitive(tool->distance_entry, FALSE);
      gtk_widget_set_sensitive(tool->angle_entry, FALSE);
      gtk_widget_set_sensitive(tool->torsion_entry, FALSE);
      tool->suppress_row_signal = FALSE;
      return;
    }

  gtk_widget_set_sensitive(tool->row_spin, TRUE);
  row_index = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->row_spin));
  row_index = CLAMP(row_index, 1, (gint) tool->rows->len) - 1;
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->row_spin), row_index + 1);

  {
    const GdisZmatrixRow *row;
    const GdisAtom *atom;
    gint bond_ref;
    gint angle_ref;
    gint torsion_ref;
    gchar number[64];
    g_autofree gchar *label = NULL;
    g_autofree gchar *bond_text = NULL;
    g_autofree gchar *angle_text = NULL;
    g_autofree gchar *torsion_text = NULL;

    row = &g_array_index(tool->rows, GdisZmatrixRow, row_index);
    atom = g_ptr_array_index(self->active_model->atoms, row->atom_index);
    bond_ref = (row->distance_ref != G_MAXUINT)
               ? gdis_gtk4_window_find_atom_index_in_array(tool->scope, row->distance_ref) + 1
               : 0;
    angle_ref = (row->angle_ref != G_MAXUINT)
                ? gdis_gtk4_window_find_atom_index_in_array(tool->scope, row->angle_ref) + 1
                : 0;
    torsion_ref = (row->torsion_ref != G_MAXUINT)
                  ? gdis_gtk4_window_find_atom_index_in_array(tool->scope, row->torsion_ref) + 1
                  : 0;
    bond_text = (bond_ref > 0) ? g_strdup_printf("bond %d", bond_ref) : g_strdup("bond -");
    angle_text = (angle_ref > 0) ? g_strdup_printf("angle %d", angle_ref) : g_strdup("angle -");
    torsion_text = (torsion_ref > 0) ? g_strdup_printf("torsion %d", torsion_ref) : g_strdup("torsion -");
    label = g_strdup_printf("Row %d: %s  |  refs: %s, %s, %s",
                            row_index + 1,
                            atom->label && atom->label[0] != '\0' ? atom->label : atom->element,
                            bond_text,
                            angle_text,
                            torsion_text);
    gtk_label_set_text(tool->row_label, label);

    if (row->distance_ref != G_MAXUINT)
      {
        g_snprintf(number, sizeof(number), "%.6f", row->distance);
        gtk_editable_set_text(GTK_EDITABLE(tool->distance_entry), number);
      }
    else
      {
        gtk_editable_set_text(GTK_EDITABLE(tool->distance_entry), "");
      }

    if (row->angle_ref != G_MAXUINT)
      {
        g_snprintf(number, sizeof(number), "%.6f", row->angle);
        gtk_editable_set_text(GTK_EDITABLE(tool->angle_entry), number);
      }
    else
      {
        gtk_editable_set_text(GTK_EDITABLE(tool->angle_entry), "");
      }

    if (row->torsion_ref != G_MAXUINT)
      {
        g_snprintf(number, sizeof(number), "%.6f", row->torsion);
        gtk_editable_set_text(GTK_EDITABLE(tool->torsion_entry), number);
      }
    else
      {
        gtk_editable_set_text(GTK_EDITABLE(tool->torsion_entry), "");
      }

    gtk_widget_set_sensitive(tool->distance_entry, row->distance_ref != G_MAXUINT);
    gtk_widget_set_sensitive(tool->angle_entry, row->angle_ref != G_MAXUINT);
    gtk_widget_set_sensitive(tool->torsion_entry, row->torsion_ref != G_MAXUINT);
  }

  tool->suppress_row_signal = FALSE;
}

static gboolean
gdis_gtk4_window_commit_zmatrix_row_edits(GdisGtk4Window *self, GError **error)
{
  GdisZmatrixTool *tool;
  GdisZmatrixRow *row;
  gint row_index;
  g_autofree gchar *report = NULL;

  g_return_val_if_fail(self != NULL, FALSE);

  tool = self->zmatrix_tool;
  if (!tool || !tool->rows || tool->rows->len == 0u || !self->active_model)
    return TRUE;

  row_index = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->row_spin));
  row_index = CLAMP(row_index, 1, (gint) tool->rows->len) - 1;
  row = &g_array_index(tool->rows, GdisZmatrixRow, row_index);

  if (row->distance_ref != G_MAXUINT)
    {
      if (!gdis_gtk4_window_parse_entry_double(tool->distance_entry, "Bond length", &row->distance, error))
        return FALSE;
      if (row->distance <= 0.0)
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_FAILED,
                      "Bond length must be greater than zero.");
          return FALSE;
        }
    }

  if (row->angle_ref != G_MAXUINT &&
      !gdis_gtk4_window_parse_entry_double(tool->angle_entry, "Angle", &row->angle, error))
    return FALSE;

  if (row->torsion_ref != G_MAXUINT &&
      !gdis_gtk4_window_parse_entry_double(tool->torsion_entry, "Torsion", &row->torsion, error))
    return FALSE;

  report = gdis_format_zmatrix_rows(self->active_model,
                                    tool->scope,
                                    tool->rows,
                                    gtk_drop_down_get_selected(GTK_DROP_DOWN(tool->scope_dropdown)) == 1u);
  gtk_text_buffer_set_text(tool->buffer, report ? report : "Z-matrix report unavailable.", -1);
  return TRUE;
}

static void
gdis_gtk4_window_refresh_zmatrix_tool(GdisGtk4Window *self)
{
  GdisZmatrixTool *tool;
  gboolean use_selection;
  GArray *scope;
  GArray *rows;
  gint wanted_row;
  GError *error;
  g_autofree gchar *report = NULL;

  g_return_if_fail(self != NULL);

  tool = self->zmatrix_tool;
  if (!tool || !tool->buffer)
    return;

  if (!self->active_model)
    {
      gdis_zmatrix_tool_clear_working_set(tool);
      gtk_text_buffer_set_text(tool->buffer, "No active model loaded.", -1);
      gdis_gtk4_window_refresh_zmatrix_row_controls(self);
      return;
    }

  use_selection = gtk_drop_down_get_selected(GTK_DROP_DOWN(tool->scope_dropdown)) == 1u;
  wanted_row = tool->row_spin ? gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->row_spin)) : 1;
  error = NULL;
  if (!gdis_build_zmatrix_rows(self->active_model,
                               self->selected_atoms,
                               use_selection,
                               &scope,
                               &rows,
                               &error))
    {
      gdis_zmatrix_tool_clear_working_set(tool);
      gtk_text_buffer_set_text(tool->buffer,
                               error ? error->message : "Z-matrix report is unavailable.",
                               -1);
      g_clear_error(&error);
      gdis_gtk4_window_refresh_zmatrix_row_controls(self);
      return;
    }

  gdis_zmatrix_tool_clear_working_set(tool);
  tool->scope = scope;
  tool->rows = rows;
  report = gdis_format_zmatrix_rows(self->active_model, tool->scope, tool->rows, use_selection);
  gtk_text_buffer_set_text(tool->buffer, report ? report : "Z-matrix report is unavailable.", -1);

  tool->suppress_row_signal = TRUE;
  gtk_spin_button_set_range(GTK_SPIN_BUTTON(tool->row_spin), 1.0, (gdouble) MAX(tool->rows->len, 1u));
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->row_spin),
                            (gdouble) CLAMP(wanted_row, 1, (gint) tool->rows->len));
  tool->suppress_row_signal = FALSE;
  gdis_gtk4_window_refresh_zmatrix_row_controls(self);
}

static void
on_zmatrix_generate_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self)
    return;

  gdis_gtk4_window_refresh_zmatrix_tool(self);
  gdis_gtk4_window_log(self, "Rebuilt the working Z-matrix from the current model.\n");
}

static void
on_zmatrix_scope_changed(GObject *object, GParamSpec *pspec, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) object;
  (void) pspec;

  self = user_data;
  if (!self)
    return;

  gdis_gtk4_window_refresh_zmatrix_tool(self);
}

static void
on_zmatrix_row_changed(GtkSpinButton *spin, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) spin;

  self = user_data;
  if (!self || !self->zmatrix_tool || self->zmatrix_tool->suppress_row_signal)
    return;

  gdis_gtk4_window_refresh_zmatrix_row_controls(self);
}

static void
on_zmatrix_entry_changed(GtkEditable *editable, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;

  (void) editable;

  self = user_data;
  if (!self || !self->zmatrix_tool || self->zmatrix_tool->suppress_row_signal)
    return;

  error = NULL;
  if (!gdis_gtk4_window_commit_zmatrix_row_edits(self, &error))
    g_clear_error(&error);
}

static void
on_zmatrix_apply_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GError *error;
  g_autofree gchar *summary = NULL;

  (void) button;

  self = user_data;
  if (!self || !self->active_model || !self->zmatrix_tool || !self->zmatrix_tool->rows)
    return;

  error = NULL;
  if (!gdis_gtk4_window_commit_zmatrix_row_edits(self, &error))
    {
      gdis_gtk4_window_log(self, "Z-matrix apply failed: %s\n",
                           error ? error->message : "invalid internal coordinates");
      g_clear_error(&error);
      return;
    }

  if (!gdis_gtk4_window_push_undo_snapshot(self, NULL))
    return;

  if (!gdis_apply_zmatrix_rows(self->active_model,
                               self->zmatrix_tool->scope,
                               self->zmatrix_tool->rows,
                               &summary,
                               &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
      gdis_gtk4_window_log(self, "Z-matrix apply failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
  gdis_gtk4_window_log(self, "%s\n", summary ? summary : "Recomputed geometry from the current Z-matrix rows.");
  gdis_gtk4_window_refresh_zmatrix_tool(self);
}

static void
on_zmatrix_build_selection_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self || !self->zmatrix_tool)
    return;

  if (!self->selected_atoms || self->selected_atoms->len == 0u)
    {
      gdis_gtk4_window_log(self, "Build Z-matrix from selection failed: select one molecule or atom set first.\n");
      return;
    }

  gtk_drop_down_set_selected(GTK_DROP_DOWN(self->zmatrix_tool->scope_dropdown), 1u);
  gdis_gtk4_window_refresh_zmatrix_tool(self);
  gdis_gtk4_window_log(self, "Rebuilt the working Z-matrix from the current selection order.\n");
}

static void
gdis_gtk4_window_present_zmatrix_tool(GdisGtk4Window *self)
{
  GdisZmatrixTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *row;
  GtkWidget *editor_frame;
  GtkWidget *editor_grid;
  GtkWidget *units_label;
  GtkWidget *button;
  GtkWidget *scroller;
  GtkWidget *text_view;
  GtkStringList *scope_model;
  const char *const scope_items[] = {"Whole model", "Current selection", NULL};

  g_return_if_fail(self != NULL);

  if (self->zmatrix_tool && GTK_IS_WINDOW(self->zmatrix_tool->window))
    {
      gdis_gtk4_window_refresh_zmatrix_tool(self);
      gtk_window_present(GTK_WINDOW(self->zmatrix_tool->window));
      return;
    }

  tool = g_new0(GdisZmatrixTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Z-matrix");
  gtk_window_set_default_size(GTK_WINDOW(window), 920, 640);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_box_append(GTK_BOX(root), row);

  scope_model = gtk_string_list_new(scope_items);
  tool->scope_dropdown = gtk_drop_down_new(G_LIST_MODEL(scope_model), NULL);
  g_object_unref(scope_model);
  gtk_drop_down_set_selected(GTK_DROP_DOWN(tool->scope_dropdown),
                             (self->selected_atoms && self->selected_atoms->len > 0u) ? 1u : 0u);
  g_signal_connect(tool->scope_dropdown,
                   "notify::selected",
                   G_CALLBACK(on_zmatrix_scope_changed),
                   self);
  gtk_box_append(GTK_BOX(row), tool->scope_dropdown);

  button = gtk_button_new_with_label("Generate");
  g_signal_connect(button, "clicked", G_CALLBACK(on_zmatrix_generate_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Build From Selection");
  g_signal_connect(button, "clicked", G_CALLBACK(on_zmatrix_build_selection_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Recompute Geometry");
  g_signal_connect(button, "clicked", G_CALLBACK(on_zmatrix_apply_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_destroy), window);
  gtk_box_append(GTK_BOX(row), button);

  editor_frame = gtk_frame_new("Row Editor");
  gtk_box_append(GTK_BOX(root), editor_frame);

  editor_grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(editor_grid), 8);
  gtk_grid_set_column_spacing(GTK_GRID(editor_grid), 8);
  gtk_widget_set_margin_start(editor_grid, 10);
  gtk_widget_set_margin_end(editor_grid, 10);
  gtk_widget_set_margin_top(editor_grid, 10);
  gtk_widget_set_margin_bottom(editor_grid, 10);
  gtk_frame_set_child(GTK_FRAME(editor_frame), editor_grid);

  units_label = gtk_label_new("Distance units: Angstrom   Angle units: degrees");
  gtk_widget_set_halign(units_label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(editor_grid), units_label, 0, 0, 4, 1);

  button = gtk_label_new("Row");
  gtk_widget_set_halign(button, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(editor_grid), button, 0, 1, 1, 1);

  tool->row_spin = gtk_spin_button_new_with_range(1.0, 1.0, 1.0);
  gtk_widget_set_hexpand(tool->row_spin, FALSE);
  gtk_grid_attach(GTK_GRID(editor_grid), tool->row_spin, 1, 1, 1, 1);

  tool->row_label = GTK_LABEL(gtk_label_new("No editable Z-matrix rows are available."));
  gtk_label_set_wrap(tool->row_label, TRUE);
  gtk_label_set_xalign(tool->row_label, 0.0f);
  gtk_grid_attach(GTK_GRID(editor_grid), GTK_WIDGET(tool->row_label), 2, 1, 2, 1);

  button = gtk_label_new("Bond");
  gtk_widget_set_halign(button, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(editor_grid), button, 0, 2, 1, 1);
  tool->distance_entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(editor_grid), tool->distance_entry, 1, 2, 1, 1);

  button = gtk_label_new("Angle");
  gtk_widget_set_halign(button, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(editor_grid), button, 2, 2, 1, 1);
  tool->angle_entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(editor_grid), tool->angle_entry, 3, 2, 1, 1);

  button = gtk_label_new("Torsion");
  gtk_widget_set_halign(button, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(editor_grid), button, 0, 3, 1, 1);
  tool->torsion_entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(editor_grid), tool->torsion_entry, 1, 3, 1, 1);

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

  g_signal_connect(tool->row_spin, "value-changed", G_CALLBACK(on_zmatrix_row_changed), self);
  g_signal_connect(tool->distance_entry, "changed", G_CALLBACK(on_zmatrix_entry_changed), self);
  g_signal_connect(tool->angle_entry, "changed", G_CALLBACK(on_zmatrix_entry_changed), self);
  g_signal_connect(tool->torsion_entry, "changed", G_CALLBACK(on_zmatrix_entry_changed), self);
  g_signal_connect(window, "destroy", G_CALLBACK(on_zmatrix_tool_destroy), tool);

  self->zmatrix_tool = tool;
  gdis_gtk4_window_refresh_zmatrix_tool(self);
  gtk_window_present(GTK_WINDOW(window));
}

static void
on_dislocation_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisDislocationTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->dislocation_tool = NULL;
  g_free(tool);
}

static void
on_dislocation_execute_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisDislocationSettings settings;
  GdisDislocationTool *tool;
  gchar *summary;
  GError *error;

  (void) button;

  self = user_data;
  if (!self || !self->active_model || !self->dislocation_tool)
    return;

  tool = self->dislocation_tool;
  gdis_dislocation_settings_init(&settings);
  settings.type = (GdisDislocationType) gtk_drop_down_get_selected(GTK_DROP_DOWN(tool->type_dropdown));
  settings.radius = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->radius_spin));
  settings.poisson_ratio = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->poisson_spin));
  settings.selected_atoms = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(tool->selection_toggle))
                              ? self->selected_atoms
                              : NULL;

  error = NULL;
  for (guint axis = 0; axis < 3; axis++)
    {
      if (!gdis_gtk4_window_parse_entry_double(tool->line_entries[axis], "Line direction",
                                               &settings.line_direction[axis], &error) ||
          !gdis_gtk4_window_parse_entry_double(tool->burgers_entries[axis], "Burgers vector",
                                               &settings.burgers_vector[axis], &error) ||
          !gdis_gtk4_window_parse_entry_double(tool->origin_entries[axis], "Origin",
                                               &settings.origin[axis], &error))
        {
          gtk_text_buffer_set_text(tool->buffer,
                                   error ? error->message : "Invalid dislocation input.",
                                   -1);
          g_clear_error(&error);
          return;
        }
    }

  if (!gdis_gtk4_window_push_undo_snapshot(self, NULL))
    return;

  summary = NULL;
  if (!gdis_model_apply_dislocation(self->active_model, &settings, &summary, &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
      gtk_text_buffer_set_text(tool->buffer,
                               error ? error->message : "Dislocation build failed.",
                               -1);
      gdis_gtk4_window_log(self, "Dislocation build failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gtk_text_buffer_set_text(tool->buffer, summary ? summary : "Dislocation applied.", -1);
  g_free(summary);
  gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
  gdis_gtk4_window_log(self, "Applied the dislocation transform to the active model.\n");
}

static void
gdis_gtk4_window_present_dislocation_tool(GdisGtk4Window *self)
{
  GdisDislocationTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *grid;
  GtkWidget *label;
  GtkWidget *row;
  GtkWidget *button;
  GtkWidget *scroller;
  GtkWidget *text_view;
  GtkStringList *type_model;
  const char *const type_items[] = {"Screw", "Edge (experimental)", NULL};
  GdisDislocationSettings defaults;

  g_return_if_fail(self != NULL);

  if (self->dislocation_tool && GTK_IS_WINDOW(self->dislocation_tool->window))
    {
      gtk_window_present(GTK_WINDOW(self->dislocation_tool->window));
      return;
    }

  gdis_dislocation_settings_init(&defaults);
  tool = g_new0(GdisDislocationTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Dislocations");
  gtk_window_set_default_size(GTK_WINDOW(window), 840, 520);

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

  label = gtk_label_new("Type");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  type_model = gtk_string_list_new(type_items);
  tool->type_dropdown = gtk_drop_down_new(G_LIST_MODEL(type_model), NULL);
  g_object_unref(type_model);
  gtk_grid_attach(GTK_GRID(grid), tool->type_dropdown, 1, 0, 2, 1);

  for (guint axis = 0; axis < 3; axis++)
    {
      static const char *const axis_labels[] = {"X", "Y", "Z"};

      label = gtk_label_new(axis_labels[axis]);
      gtk_widget_set_halign(label, GTK_ALIGN_START);
      gtk_grid_attach(GTK_GRID(grid), label, axis + 1, 1, 1, 1);
    }

  label = gtk_label_new("Line direction");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 2, 1, 1);
  label = gtk_label_new("Burgers vector");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 3, 1, 1);
  label = gtk_label_new("Origin");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 4, 1, 1);

  for (guint axis = 0; axis < 3; axis++)
    {
      tool->line_entries[axis] = gtk_entry_new();
      tool->burgers_entries[axis] = gtk_entry_new();
      tool->origin_entries[axis] = gtk_entry_new();
      gtk_editable_set_text(GTK_EDITABLE(tool->line_entries[axis]),
                            axis == 2 ? "1.0" : "0.0");
      gtk_editable_set_text(GTK_EDITABLE(tool->burgers_entries[axis]),
                            axis == 2 ? "1.0" : "0.0");
      gtk_editable_set_text(GTK_EDITABLE(tool->origin_entries[axis]), "0.0");
      gtk_grid_attach(GTK_GRID(grid), tool->line_entries[axis], axis + 1, 2, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), tool->burgers_entries[axis], axis + 1, 3, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), tool->origin_entries[axis], axis + 1, 4, 1, 1);
    }

  label = gtk_label_new("Radius");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 5, 1, 1);
  tool->radius_spin = gtk_spin_button_new_with_range(0.0, 200.0, 0.5);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->radius_spin), 0.0);
  gtk_grid_attach(GTK_GRID(grid), tool->radius_spin, 1, 5, 1, 1);

  label = gtk_label_new("Poisson");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 2, 5, 1, 1);
  tool->poisson_spin = gtk_spin_button_new_with_range(0.01, 0.49, 0.01);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->poisson_spin), defaults.poisson_ratio);
  gtk_grid_attach(GTK_GRID(grid), tool->poisson_spin, 3, 5, 1, 1);

  tool->selection_toggle = gtk_check_button_new_with_label("Apply to current selection only");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tool->selection_toggle),
                               self->selected_atoms && self->selected_atoms->len > 0u);
  gtk_box_append(GTK_BOX(root), tool->selection_toggle);

  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_widget_set_halign(row, GTK_ALIGN_END);
  gtk_box_append(GTK_BOX(root), row);
  button = gtk_button_new_with_label("Execute");
  g_signal_connect(button, "clicked", G_CALLBACK(on_dislocation_execute_clicked), self);
  gtk_box_append(GTK_BOX(row), button);
  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_destroy), window);
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
  gtk_text_buffer_set_text(tool->buffer,
                           "Configure the line direction, Burgers vector, and optional radius filter,\nthen Execute to apply a native GTK4 dislocation transform to the active model.",
                           -1);
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  g_signal_connect(window, "destroy", G_CALLBACK(on_dislocation_tool_destroy), tool);

  self->dislocation_tool = tool;
  gtk_window_present(GTK_WINDOW(window));
}

static void
on_docking_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisDockingTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->docking_tool = NULL;
  g_free(tool);
}

static void
on_docking_generate_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisDockingTool *tool;
  GdisDockingSettings settings;
  GError *error;
  gchar *summary;

  (void) button;

  self = user_data;
  if (!self || !self->active_model || !self->docking_tool)
    return;

  tool = self->docking_tool;
  gdis_docking_settings_init(&settings);
  settings.output_dir = gtk_editable_get_text(GTK_EDITABLE(tool->output_entry));
  settings.project_name = gtk_editable_get_text(GTK_EDITABLE(tool->name_entry));
  settings.preview_limit = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->preview_spin));
  settings.selected_atoms = self->selected_atoms;
  for (guint axis = 0; axis < 3; axis++)
    {
      settings.grid_counts[axis] =
        (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->grid_spins[axis]));
      settings.rotation_counts[axis] =
        (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->rotation_spins[axis]));
      settings.translation_span[axis] =
        gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->span_spins[axis]));
    }

  error = NULL;
  summary = NULL;
  if (!gdis_generate_docking_project(self->active_model, &settings, &summary, &error))
    {
      gtk_text_buffer_set_text(tool->buffer,
                               error ? error->message : "Docking project generation failed.",
                               -1);
      gdis_gtk4_window_log(self, "Docking project generation failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gtk_text_buffer_set_text(tool->buffer, summary ? summary : "Docking project generated.", -1);
  g_free(summary);
  gdis_gtk4_window_log(self, "Generated the docking project from the current selection.\n");
}

static void
gdis_gtk4_window_present_docking_tool(GdisGtk4Window *self)
{
  GdisDockingTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *grid;
  GtkWidget *label;
  GtkWidget *row;
  GtkWidget *button;
  GtkWidget *scroller;
  GtkWidget *text_view;
  g_autofree gchar *default_output = NULL;
  g_autofree gchar *default_name = NULL;

  g_return_if_fail(self != NULL);

  if (self->docking_tool && GTK_IS_WINDOW(self->docking_tool->window))
    {
      gtk_window_present(GTK_WINDOW(self->docking_tool->window));
      return;
    }

  tool = g_new0(GdisDockingTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Docking");
  gtk_window_set_default_size(GTK_WINDOW(window), 860, 540);

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

  default_output = g_build_filename(g_get_current_dir(), "gdis_docking", NULL);
  default_name = g_strdup_printf("%s_dock",
                                 self->active_model && self->active_model->basename
                                   ? self->active_model->basename
                                   : "project");

  label = gtk_label_new("Output root");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  tool->output_entry = gtk_entry_new();
  gtk_editable_set_text(GTK_EDITABLE(tool->output_entry), default_output);
  gtk_widget_set_hexpand(tool->output_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->output_entry, 1, 0, 5, 1);

  label = gtk_label_new("Project name");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
  tool->name_entry = gtk_entry_new();
  gtk_editable_set_text(GTK_EDITABLE(tool->name_entry), default_name);
  gtk_grid_attach(GTK_GRID(grid), tool->name_entry, 1, 1, 5, 1);

  for (guint axis = 0; axis < 3; axis++)
    {
      static const char *const axis_labels[] = {"X", "Y", "Z"};

      label = gtk_label_new(axis_labels[axis]);
      gtk_widget_set_halign(label, GTK_ALIGN_START);
      gtk_grid_attach(GTK_GRID(grid), label, axis + 1, 2, 1, 1);
    }

  label = gtk_label_new("Grid counts");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 3, 1, 1);
  label = gtk_label_new("Rotation counts");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 4, 1, 1);
  label = gtk_label_new("Translation span");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 5, 1, 1);

  for (guint axis = 0; axis < 3; axis++)
    {
      tool->grid_spins[axis] = gtk_spin_button_new_with_range(1.0, 25.0, 1.0);
      tool->rotation_spins[axis] = gtk_spin_button_new_with_range(1.0, 36.0, 1.0);
      tool->span_spins[axis] = gtk_spin_button_new_with_range(0.0, 100.0, 0.5);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->grid_spins[axis]), axis < 2 ? 3.0 : 1.0);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->rotation_spins[axis]), axis == 2 ? 4.0 : 1.0);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->span_spins[axis]), axis < 2 ? 6.0 : 2.0);
      gtk_grid_attach(GTK_GRID(grid), tool->grid_spins[axis], axis + 1, 3, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), tool->rotation_spins[axis], axis + 1, 4, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), tool->span_spins[axis], axis + 1, 5, 1, 1);
    }

  label = gtk_label_new("Preview limit");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 6, 1, 1);
  tool->preview_spin = gtk_spin_button_new_with_range(1.0, 500.0, 1.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->preview_spin), 24.0);
  gtk_grid_attach(GTK_GRID(grid), tool->preview_spin, 1, 6, 1, 1);

  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_widget_set_halign(row, GTK_ALIGN_END);
  gtk_box_append(GTK_BOX(root), row);
  button = gtk_button_new_with_label("Generate");
  g_signal_connect(button, "clicked", G_CALLBACK(on_docking_generate_clicked), self);
  gtk_box_append(GTK_BOX(row), button);
  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_destroy), window);
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
  gtk_text_buffer_set_text(tool->buffer,
                           "Select a fragment in the viewer first, then Generate to create a native GTK4 docking pose ensemble project.",
                           -1);
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  g_signal_connect(window, "destroy", G_CALLBACK(on_docking_tool_destroy), tool);

  self->docking_tool = tool;
  gtk_window_present(GTK_WINDOW(window));
}

static gchar *
gdis_qbox_normalize_symbol(const char *text)
{
  g_autofree gchar *trimmed = NULL;
  gchar *normalized;

  if (!text || text[0] == '\0')
    return g_strdup("X");

  trimmed = g_strstrip(g_strdup(text));
  if (trimmed[0] == '\0')
    return g_strdup("X");

  normalized = g_ascii_strdown(trimmed, -1);
  normalized[0] = (gchar) g_ascii_toupper(normalized[0]);
  for (guint i = 1; normalized[i] != '\0'; i++)
    normalized[i] = (gchar) g_ascii_tolower(normalized[i]);

  return normalized;
}

static gchar *
gdis_qbox_entry_text(GtkWidget *entry)
{
  g_autofree gchar *text = NULL;

  g_return_val_if_fail(GTK_IS_EDITABLE(entry), g_strdup(""));

  text = g_strdup(gtk_editable_get_text(GTK_EDITABLE(entry)));
  g_strstrip(text);
  return g_steal_pointer(&text);
}

static gchar *
gdis_qbox_entry_text_or_default(GtkWidget *entry, const char *fallback)
{
  g_autofree gchar *text = NULL;

  text = gdis_qbox_entry_text(entry);
  if (!text[0])
    return g_strdup(fallback ? fallback : "");
  return g_steal_pointer(&text);
}

static gboolean
gdis_qbox_is_remote_uri(const char *text)
{
  return (text &&
          (g_str_has_prefix(text, "http://") ||
           g_str_has_prefix(text, "https://") ||
           g_str_has_prefix(text, "ftp://") ||
           g_str_has_prefix(text, "file://")));
}

static gchar *
gdis_qbox_symbol_alias(const char *symbol)
{
  g_autofree gchar *lower = NULL;
  GString *alias;

  lower = g_ascii_strdown(symbol ? symbol : "x", -1);
  alias = g_string_new("");
  for (guint i = 0; lower[i] != '\0'; i++)
    {
      if (g_ascii_isalnum(lower[i]))
        g_string_append_c(alias, lower[i]);
    }
  if (alias->len == 0)
    g_string_append(alias, "x");
  return g_string_free(alias, FALSE);
}

static gchar *
gdis_qbox_sanitize_path_component(const char *text)
{
  GString *value;

  value = g_string_new("");
  if (text)
    {
      for (guint i = 0; text[i] != '\0'; i++)
        {
          const guchar ch = (guchar) text[i];

          if (g_ascii_isalnum(ch) || ch == '.' || ch == '-' || ch == '_')
            g_string_append_c(value, (gchar) ch);
          else
            g_string_append_c(value, '_');
        }
    }
  if (value->len == 0)
    g_string_append(value, "qbox_job");
  return g_string_free(value, FALSE);
}

static gchar *
gdis_qbox_job_slug(const char *text)
{
  return gdis_qbox_sanitize_path_component(text);
}

static gchar *
gdis_qbox_model_base_name(const GdisModel *model)
{
  g_autofree gchar *basename = NULL;
  gchar *dot;

  if (!model || !model->basename)
    return g_strdup("qbox_job");

  basename = g_strdup(model->basename);
  dot = strrchr(basename, '.');
  if (dot)
    *dot = '\0';
  return gdis_qbox_job_slug(basename);
}

static gchar *
gdis_qbox_default_pseudo_cache_dir(void)
{
  return g_build_filename(g_get_user_data_dir(),
                          "gdis",
                          "qbox-pseudos",
                          "sg15_oncv",
                          "xml",
                          NULL);
}

static gchar *
gdis_qbox_default_pseudo_source(void)
{
  g_autofree gchar *cache_dir = NULL;

  cache_dir = gdis_qbox_default_pseudo_cache_dir();
  if (g_file_test(cache_dir, G_FILE_TEST_IS_DIR))
    return g_steal_pointer(&cache_dir);

  return g_strdup(GDIS_QBOX_DEFAULT_PSEUDO_SOURCE_URL);
}

static gboolean
gdis_qbox_download_remote_file(const char *uri,
                               const char *dest_path,
                               GError **error)
{
  g_autofree gchar *curl_path = NULL;
  gint wait_status = 0;
  gchar *stdout_text = NULL;
  gchar *stderr_text = NULL;
  gchar *argv[9];

  g_return_val_if_fail(uri != NULL, FALSE);
  g_return_val_if_fail(dest_path != NULL, FALSE);

  curl_path = g_find_program_in_path("curl");
  if (!curl_path)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not find 'curl', which is needed to download Qbox pseudopotentials.");
      return FALSE;
    }

  argv[0] = curl_path;
  argv[1] = "-L";
  argv[2] = "--fail";
  argv[3] = "--silent";
  argv[4] = "--show-error";
  argv[5] = "-o";
  argv[6] = (gchar *) dest_path;
  argv[7] = (gchar *) uri;
  argv[8] = NULL;

  if (!g_spawn_sync(NULL,
                    argv,
                    NULL,
                    G_SPAWN_SEARCH_PATH,
                    NULL,
                    NULL,
                    &stdout_text,
                    &stderr_text,
                    &wait_status,
                    error))
    {
      g_free(stdout_text);
      g_free(stderr_text);
      return FALSE;
    }

  if (!g_spawn_check_wait_status(wait_status, error))
    {
      g_prefix_error(error, "Could not download '%s': ", uri);
      if (stderr_text && stderr_text[0])
        {
          g_autofree gchar *previous = g_strdup((*error)->message);
          g_clear_error(error);
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_IO,
                      "%s%s",
                      previous ? previous : "Download failed. ",
                      stderr_text);
        }
      g_free(stdout_text);
      g_free(stderr_text);
      return FALSE;
    }

  g_free(stdout_text);
  g_free(stderr_text);
  return TRUE;
}

static gchar *
gdis_qbox_guess_pseudo_path(const char *pseudo_dir,
                            const char *symbol,
                            const char *xc_name)
{
  g_autofree gchar *trimmed = NULL;
  g_autofree gchar *normalized_symbol = NULL;
  g_autofree gchar *upper_xc = NULL;
  static const char *const patterns[] = {
    "%s_HSCV_%s-1.0.xml",
    "%s_HSCV_%s-1.1.xml",
    "%s_HSCV_%s.xml",
    "%s_ONCV_%s-1.0.xml",
    "%s_ONCV_%s-1.2.xml",
    "%s_ONCV_%s.xml",
    "%s_%s.xml",
    "%s.xml",
    NULL
  };

  if (!pseudo_dir || pseudo_dir[0] == '\0')
    return NULL;

  trimmed = g_strdup(pseudo_dir);
  g_strstrip(trimmed);
  while (g_str_has_suffix(trimmed, "/"))
    trimmed[strlen(trimmed) - 1] = '\0';
  if (trimmed[0] == '\0')
    return NULL;

  normalized_symbol = gdis_qbox_normalize_symbol(symbol);
  upper_xc = g_ascii_strup((xc_name && xc_name[0]) ? xc_name : "PBE", -1);

  if (gdis_qbox_is_remote_uri(trimmed))
    {
      g_autofree gchar *basename = NULL;

      basename = g_strdup_printf("%s_ONCV_%s-1.2.xml", normalized_symbol, upper_xc);
      if (g_strrstr(trimmed, "sg15_oncv") || g_str_has_suffix(trimmed, "/xml"))
        return g_strdup_printf("%s/%s", trimmed, basename);

      return g_strdup_printf("%s/%s/%s_HSCV_%s-1.0.xml",
                             trimmed,
                             normalized_symbol,
                             normalized_symbol,
                             upper_xc);
    }

  if (!g_file_test(trimmed, G_FILE_TEST_IS_DIR))
    return NULL;

  for (guint i = 0; patterns[i] != NULL; i++)
    {
      g_autofree gchar *basename = NULL;
      g_autofree gchar *path = NULL;
      g_autofree gchar *nested = NULL;

      basename = g_strdup_printf(patterns[i], normalized_symbol, upper_xc);
      path = g_build_filename(trimmed, basename, NULL);
      if (g_file_test(path, G_FILE_TEST_EXISTS))
        return g_steal_pointer(&path);

      nested = g_build_filename(trimmed, normalized_symbol, basename, NULL);
      if (g_file_test(nested, G_FILE_TEST_EXISTS))
        return g_steal_pointer(&nested);
    }

  return NULL;
}

static gchar *
gdis_qbox_lookup_executable_path(GdisGtk4Window *self)
{
  const GdisExecutableSpec *spec;
  const char *stored_path;

  g_return_val_if_fail(self != NULL, NULL);

  stored_path = self->executable_paths ?
    g_hash_table_lookup(self->executable_paths, "qbox") :
    NULL;
  if (stored_path && stored_path[0] != '\0')
    return g_strdup(stored_path);

  spec = gdis_executable_spec_for_id("qbox");
  return spec ? gdis_executable_detect_path(spec) : NULL;
}

static gchar *
gdis_qbox_default_workdir(GdisGtk4Window *self, const char *job_slug)
{
  g_autofree gchar *cwd = NULL;

  (void) self;

  cwd = g_get_current_dir();
  return g_build_filename(cwd, "qbox_jobs", job_slug && job_slug[0] ? job_slug : "qbox_job", NULL);
}

static GPtrArray *
gdis_qbox_collect_unique_species(const GdisModel *model)
{
  GHashTable *seen;
  GPtrArray *species;

  g_return_val_if_fail(model != NULL, NULL);

  seen = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
  species = g_ptr_array_new_with_free_func(g_free);
  for (guint i = 0; i < model->atoms->len; i++)
    {
      const GdisAtom *atom;
      g_autofree gchar *symbol = NULL;

      atom = g_ptr_array_index(model->atoms, i);
      symbol = gdis_qbox_normalize_symbol(atom->element && atom->element[0] ? atom->element : atom->label);
      if (g_hash_table_contains(seen, symbol))
        continue;
      g_hash_table_add(seen, g_strdup(symbol));
      g_ptr_array_add(species, g_strdup(symbol));
    }
  g_hash_table_unref(seen);
  return species;
}

static guint
gdis_qbox_count_unique_species(const GdisModel *model)
{
  g_autoptr(GPtrArray) species = NULL;

  species = gdis_qbox_collect_unique_species(model);
  return species ? species->len : 0u;
}

static gboolean
gdis_qbox_build_cell_and_shift(const GdisModel *model,
                               gboolean prefer_model_cell,
                               gdouble padding_angstrom,
                               gdouble cell_bohr[9],
                               gdouble shift_angstrom[3],
                               gchar **cell_mode_out,
                               GError **error)
{
  gboolean used_model_cell;

  g_return_val_if_fail(model != NULL, FALSE);

  memset(cell_bohr, 0, sizeof(gdouble) * 9);
  shift_angstrom[0] = 0.0;
  shift_angstrom[1] = 0.0;
  shift_angstrom[2] = 0.0;
  used_model_cell = FALSE;

  if (prefer_model_cell && model->periodic)
    {
      gdouble a_vec[3];
      gdouble b_vec[3];
      gdouble c_vec[3];

      if (gdis_model_build_cell_vectors(model, a_vec, b_vec, c_vec))
        {
          cell_bohr[0] = a_vec[0] * GDIS_ANGSTROM_TO_BOHR;
          cell_bohr[1] = a_vec[1] * GDIS_ANGSTROM_TO_BOHR;
          cell_bohr[2] = a_vec[2] * GDIS_ANGSTROM_TO_BOHR;
          cell_bohr[3] = b_vec[0] * GDIS_ANGSTROM_TO_BOHR;
          cell_bohr[4] = b_vec[1] * GDIS_ANGSTROM_TO_BOHR;
          cell_bohr[5] = b_vec[2] * GDIS_ANGSTROM_TO_BOHR;
          cell_bohr[6] = c_vec[0] * GDIS_ANGSTROM_TO_BOHR;
          cell_bohr[7] = c_vec[1] * GDIS_ANGSTROM_TO_BOHR;
          cell_bohr[8] = c_vec[2] * GDIS_ANGSTROM_TO_BOHR;
          used_model_cell = TRUE;
        }
    }

  if (!used_model_cell)
    {
      gdouble min_pos[3];
      gdouble max_pos[3];
      gboolean have_atoms;

      have_atoms = FALSE;
      for (guint i = 0; i < model->atoms->len; i++)
        {
          const GdisAtom *atom;

          atom = g_ptr_array_index(model->atoms, i);
          if (!have_atoms)
            {
              memcpy(min_pos, atom->position, sizeof(min_pos));
              memcpy(max_pos, atom->position, sizeof(max_pos));
              have_atoms = TRUE;
            }
          else
            {
              for (guint axis = 0; axis < 3; axis++)
                {
                  min_pos[axis] = MIN(min_pos[axis], atom->position[axis]);
                  max_pos[axis] = MAX(max_pos[axis], atom->position[axis]);
                }
            }
        }

      if (!have_atoms)
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_FAILED,
                      "The active model does not contain any atoms to export to Qbox.");
          return FALSE;
        }

      for (guint axis = 0; axis < 3; axis++)
        {
          gdouble span_angstrom;

          span_angstrom = MAX(max_pos[axis] - min_pos[axis], 1.0e-3) + 2.0 * padding_angstrom;
          cell_bohr[axis * 3 + axis] = span_angstrom * GDIS_ANGSTROM_TO_BOHR;
          shift_angstrom[axis] = padding_angstrom - min_pos[axis];
        }
    }

  if (cell_mode_out)
    *cell_mode_out = g_strdup(used_model_cell ? "Model periodic cell" : "Bounding box cell");
  return TRUE;
}

static void
gdis_qbox_species_row_free(gpointer data)
{
  GdisQboxSpeciesRow *row;

  row = data;
  if (!row)
    return;

  g_clear_pointer(&row->element, g_free);
  g_free(row);
}

static void
gdis_qbox_tool_set_generated_text(GdisQboxTool *tool, const char *text)
{
  g_return_if_fail(tool != NULL);

  tool->suppress_input_signal = TRUE;
  gtk_text_buffer_set_text(tool->deck_buffer, text ? text : "", -1);
  tool->suppress_input_signal = FALSE;

  g_free(tool->last_generated_input);
  tool->last_generated_input = g_strdup(text ? text : "");
  tool->editor_dirty = FALSE;
}

static GdisQboxSpeciesRow *
gdis_qbox_find_species_row(GdisQboxTool *tool, const char *symbol)
{
  g_return_val_if_fail(tool != NULL, NULL);

  if (!tool->species_rows)
    return NULL;

  for (guint i = 0; i < tool->species_rows->len; i++)
    {
      GdisQboxSpeciesRow *row;

      row = g_ptr_array_index(tool->species_rows, i);
      if (g_strcmp0(row->element, symbol) == 0)
        return row;
    }

  return NULL;
}

static gchar *
gdis_qbox_restart_reference(GdisQboxTool *tool)
{
  g_autofree gchar *restart_text = NULL;
  g_autofree gchar *resolved = NULL;

  restart_text = gdis_qbox_entry_text(tool->restart_entry);
  if (!restart_text[0])
    return NULL;
  if (gdis_qbox_is_remote_uri(restart_text))
    return g_steal_pointer(&restart_text);

  resolved = gdis_gtk4_window_resolve_path(restart_text);
  if (resolved)
    return g_strdup("restart.xml");

  return g_strdup("<missing-restart-xml>");
}

static gchar *
gdis_qbox_species_reference(GdisQboxSpeciesRow *row)
{
  g_autofree gchar *pseudo_text = NULL;
  g_autofree gchar *resolved = NULL;
  g_autofree gchar *basename = NULL;

  g_return_val_if_fail(row != NULL, g_strdup("<missing-species>"));

  pseudo_text = gdis_qbox_entry_text(row->pseudo_entry);
  if (!pseudo_text[0])
    return g_strdup_printf("<set-pseudo-for-%s>", row->element ? row->element : "X");
  if (gdis_qbox_is_remote_uri(pseudo_text))
    return g_steal_pointer(&pseudo_text);

  resolved = gdis_gtk4_window_resolve_path(pseudo_text);
  if (!resolved)
    return g_strdup_printf("<missing-pseudo-for-%s>", row->element ? row->element : "X");

  basename = g_path_get_basename(resolved);
  return gdis_qbox_sanitize_path_component(basename);
}

static void
gdis_qbox_append_candidate_uri(GPtrArray *candidates, const char *uri)
{
  g_return_if_fail(candidates != NULL);

  if (!uri || !uri[0])
    return;
  for (guint i = 0; i < candidates->len; i++)
    {
      const char *existing = g_ptr_array_index(candidates, i);

      if (g_strcmp0(existing, uri) == 0)
        return;
    }
  g_ptr_array_add(candidates, g_strdup(uri));
}

static gchar *
gdis_qbox_remote_pseudo_fallback(const char *symbol, const char *xc_name)
{
  return gdis_qbox_guess_pseudo_path(GDIS_QBOX_DEFAULT_PSEUDO_SOURCE_URL,
                                     symbol,
                                     xc_name);
}

static gboolean
gdis_qbox_prepare_local_pseudo_library(GdisGtk4Window *self,
                                       GdisQboxTool *tool,
                                       guint *downloaded_out,
                                       guint *relinked_out,
                                       GString *warnings,
                                       GError **error)
{
  g_autofree gchar *cache_dir = NULL;
  g_autofree gchar *pseudo_dir = NULL;
  g_autofree gchar *xc_name = NULL;
  guint downloaded = 0u;
  guint relinked = 0u;

  g_return_val_if_fail(self != NULL, FALSE);
  g_return_val_if_fail(tool != NULL, FALSE);

  cache_dir = gdis_qbox_default_pseudo_cache_dir();
  if (g_mkdir_with_parents(cache_dir, 0755) != 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not create the Qbox pseudo cache '%s': %s",
                  cache_dir,
                  g_strerror(errno));
      return FALSE;
    }

  pseudo_dir = gdis_qbox_entry_text(tool->pseudo_dir_entry);
  xc_name = gdis_qbox_entry_text_or_default(tool->xc_entry, "PBE");

  for (guint i = 0; tool->species_rows && i < tool->species_rows->len; i++)
    {
      GdisQboxSpeciesRow *row;
      g_autofree gchar *current = NULL;
      g_autofree gchar *resolved = NULL;
      g_autofree gchar *guessed = NULL;
      g_autoptr(GPtrArray) candidates = NULL;
      gboolean linked = FALSE;

      row = g_ptr_array_index(tool->species_rows, i);
      current = gdis_qbox_entry_text(row->pseudo_entry);
      resolved = current[0] && !gdis_qbox_is_remote_uri(current) ?
        gdis_gtk4_window_resolve_path(current) : NULL;
      if (resolved)
        continue;

      candidates = g_ptr_array_new_with_free_func(g_free);
      if (current[0] && gdis_qbox_is_remote_uri(current))
        gdis_qbox_append_candidate_uri(candidates, current);

      if (pseudo_dir[0] && gdis_qbox_is_remote_uri(pseudo_dir))
        {
          guessed = gdis_qbox_guess_pseudo_path(pseudo_dir, row->element, xc_name);
          gdis_qbox_append_candidate_uri(candidates, guessed);
        }

      g_clear_pointer(&guessed, g_free);
      guessed = gdis_qbox_remote_pseudo_fallback(row->element, xc_name);
      gdis_qbox_append_candidate_uri(candidates, guessed);

      for (guint candidate_index = 0; candidate_index < candidates->len; candidate_index++)
        {
          const char *candidate_uri;
          g_autofree gchar *basename = NULL;
          g_autofree gchar *safe_name = NULL;
          g_autofree gchar *dest_path = NULL;
          GError *download_error = NULL;

          candidate_uri = g_ptr_array_index(candidates, candidate_index);
          basename = g_path_get_basename(candidate_uri);
          safe_name = gdis_qbox_sanitize_path_component(basename);
          dest_path = g_build_filename(cache_dir, safe_name, NULL);

          if (!g_file_test(dest_path, G_FILE_TEST_EXISTS))
            {
              if (!gdis_qbox_download_remote_file(candidate_uri, dest_path, &download_error))
                {
                  if (warnings)
                    g_string_append_printf(warnings,
                                           "Could not download pseudo for %s from %s: %s\n",
                                           row->element,
                                           candidate_uri,
                                           download_error ? download_error->message : "unknown error");
                  g_clear_error(&download_error);
                  continue;
                }
              downloaded++;
            }

          gtk_editable_set_text(GTK_EDITABLE(row->pseudo_entry), dest_path);
          linked = TRUE;
          relinked++;
          break;
        }

      if (!linked)
        {
          if (warnings)
            g_string_append_printf(warnings,
                                   "No usable pseudo source was found for %s.\n",
                                   row->element);
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_IO,
                      "Could not prepare a local Qbox pseudo for %s.",
                      row->element);
          return FALSE;
        }
    }

  gtk_editable_set_text(GTK_EDITABLE(tool->pseudo_dir_entry), cache_dir);
  if (downloaded_out)
    *downloaded_out = downloaded;
  if (relinked_out)
    *relinked_out = relinked;
  return TRUE;
}

static void
gdis_qbox_tool_rebuild_species_rows(GdisGtk4Window *self,
                                    GdisQboxTool *tool)
{
  GtkWidget *grid;

  g_return_if_fail(self != NULL);
  g_return_if_fail(tool != NULL);

  g_clear_pointer(&tool->species_rows, g_ptr_array_unref);
  tool->species_rows = g_ptr_array_new_with_free_func(gdis_qbox_species_row_free);

  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_widget_set_margin_start(grid, 10);
  gtk_widget_set_margin_end(grid, 10);
  gtk_widget_set_margin_top(grid, 10);
  gtk_widget_set_margin_bottom(grid, 10);

  if (!self->active_model)
    {
      GtkWidget *label;

      label = gtk_label_new("Load a model to populate the Qbox species table.");
      gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
      gtk_label_set_wrap(GTK_LABEL(label), TRUE);
      gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
      gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(tool->species_scroller), grid);
      return;
    }

  gtk_grid_attach(GTK_GRID(grid), gtk_label_new("Element"), 0, 0, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), gtk_label_new("Species"), 1, 0, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), gtk_label_new("Pseudo URI / File"), 2, 0, 1, 1);

  g_autoptr(GPtrArray) species = gdis_qbox_collect_unique_species(self->active_model);
  for (guint i = 0; species && i < species->len; i++)
    {
      GdisQboxSpeciesRow *row;
      const char *element;
      GtkWidget *label;
      g_autofree gchar *alias = NULL;
      g_autofree gchar *guess = NULL;

      row = g_new0(GdisQboxSpeciesRow, 1);
      element = g_ptr_array_index(species, i);
      row->element = g_strdup(element);

      label = gtk_label_new(element);
      gtk_widget_set_halign(label, GTK_ALIGN_START);
      gtk_grid_attach(GTK_GRID(grid), label, 0, (gint) i + 1, 1, 1);

      row->alias_entry = gtk_entry_new();
      alias = gdis_qbox_symbol_alias(element);
      gtk_editable_set_text(GTK_EDITABLE(row->alias_entry), alias);
      gtk_grid_attach(GTK_GRID(grid), row->alias_entry, 1, (gint) i + 1, 1, 1);

      row->pseudo_entry = gtk_entry_new();
      gtk_widget_set_hexpand(row->pseudo_entry, TRUE);
      guess = gdis_qbox_guess_pseudo_path(gtk_editable_get_text(GTK_EDITABLE(tool->pseudo_dir_entry)),
                                          element,
                                          gtk_editable_get_text(GTK_EDITABLE(tool->xc_entry)));
      gtk_editable_set_text(GTK_EDITABLE(row->pseudo_entry), guess ? guess : "");
      gtk_grid_attach(GTK_GRID(grid), row->pseudo_entry, 2, (gint) i + 1, 1, 1);

      g_ptr_array_add(tool->species_rows, row);
    }

  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(tool->species_scroller), grid);
}

static gboolean
gdis_qbox_build_input_deck(GdisGtk4Window *self,
                           GdisQboxTool *tool,
                           gchar **deck_out,
                           gchar **cell_mode_out,
                           GError **error)
{
  g_autofree gchar *xc_name = NULL;
  g_autofree gchar *wf_dyn = NULL;
  g_autofree gchar *atoms_dyn = NULL;
  g_autofree gchar *scf_tol = NULL;
  g_autofree gchar *restart_reference = NULL;
  gdouble cell_bohr[9];
  gdouble shift_angstrom[3];
  gdouble padding_angstrom;
  guint scf_steps;
  guint ionic_steps;
  guint density_update;
  gint charge;
  GString *deck;
  GHashTable *element_counts;
  gboolean use_atomic_density;
  gboolean use_model_cell;
  gboolean randomize_wf;
  gboolean locked_atoms;
  gboolean have_placeholders;

  g_return_val_if_fail(self != NULL, FALSE);
  g_return_val_if_fail(tool != NULL, FALSE);
  g_return_val_if_fail(deck_out != NULL, FALSE);

  if (!self->active_model)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Load a model before generating a Qbox input deck.");
      return FALSE;
    }

  xc_name = gdis_qbox_entry_text_or_default(tool->xc_entry, "PBE");
  wf_dyn = gdis_qbox_entry_text_or_default(tool->wf_dyn_entry, "PSDA");
  atoms_dyn = gdis_qbox_entry_text_or_default(tool->atoms_dyn_entry, "LOCKED");
  scf_tol = gdis_qbox_entry_text_or_default(tool->scf_tol_entry, "1.0e-8");
  restart_reference = gdis_qbox_restart_reference(tool);
  padding_angstrom = gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->padding_spin));
  scf_steps = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->scf_steps_spin));
  ionic_steps = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->ionic_steps_spin));
  density_update = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->density_update_spin));
  charge = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->charge_spin));
  use_atomic_density = gtk_check_button_get_active(GTK_CHECK_BUTTON(tool->atomic_density_toggle));
  use_model_cell = gtk_check_button_get_active(GTK_CHECK_BUTTON(tool->use_cell_toggle));
  randomize_wf = gtk_check_button_get_active(GTK_CHECK_BUTTON(tool->randomize_toggle));
  locked_atoms = (g_ascii_strcasecmp(atoms_dyn, "LOCKED") == 0 || ionic_steps == 0u);
  have_placeholders = FALSE;

  scf_steps = MAX(scf_steps, 1u);
  density_update = MAX(density_update, 1u);

  if (!gdis_qbox_build_cell_and_shift(self->active_model,
                                      use_model_cell,
                                      padding_angstrom,
                                      cell_bohr,
                                      shift_angstrom,
                                      cell_mode_out,
                                      error))
    return FALSE;

  deck = g_string_new(
    "# Qbox deck generated by the GDIS GTK4 Qbox tool.\n"
    "# Coordinates and cell vectors are written in bohr.\n"
    "# Local pseudo files and local restart XML files are staged into the job directory on Write Input / Run.\n");

  if (restart_reference && restart_reference[0])
    {
      if (restart_reference[0] == '<')
        have_placeholders = TRUE;
      g_string_append_printf(deck,
                             "# Restart source from the controls window.\n"
                             "load %s\n\n",
                             restart_reference);
    }
  else
    {
      GHashTable *written_species;

      g_string_append_printf(deck,
                             "set cell %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
                             cell_bohr[0], cell_bohr[1], cell_bohr[2],
                             cell_bohr[3], cell_bohr[4], cell_bohr[5],
                             cell_bohr[6], cell_bohr[7], cell_bohr[8]);

      written_species = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
      for (guint i = 0; tool->species_rows && i < tool->species_rows->len; i++)
        {
          GdisQboxSpeciesRow *row;
          g_autofree gchar *alias = NULL;
          g_autofree gchar *reference = NULL;

          row = g_ptr_array_index(tool->species_rows, i);
          alias = gdis_qbox_entry_text(row->alias_entry);
          if (!alias[0])
            {
              g_free(alias);
              alias = gdis_qbox_symbol_alias(row->element);
              gtk_editable_set_text(GTK_EDITABLE(row->alias_entry), alias);
            }
          reference = gdis_qbox_species_reference(row);
          if (reference[0] == '<')
            have_placeholders = TRUE;
          if (!g_hash_table_contains(written_species, row->element))
            {
              g_hash_table_add(written_species, g_strdup(row->element));
              g_string_append_printf(deck, "species %s %s\n", alias, reference);
            }
        }
      g_hash_table_unref(written_species);

      element_counts = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
      for (guint i = 0; i < self->active_model->atoms->len; i++)
        {
          const GdisAtom *atom;
          g_autofree gchar *symbol = NULL;
          GdisQboxSpeciesRow *row;
          g_autofree gchar *alias = NULL;
          guint count;
          gdouble x_bohr;
          gdouble y_bohr;
          gdouble z_bohr;

          atom = g_ptr_array_index(self->active_model->atoms, i);
          symbol = gdis_qbox_normalize_symbol(atom->element && atom->element[0] ? atom->element : atom->label);
          row = gdis_qbox_find_species_row(tool, symbol);
          alias = row ? gdis_qbox_entry_text(row->alias_entry) : gdis_qbox_symbol_alias(symbol);
          if (!alias[0])
            {
              g_free(alias);
              alias = gdis_qbox_symbol_alias(symbol);
            }

          count = GPOINTER_TO_UINT(g_hash_table_lookup(element_counts, symbol)) + 1u;
          g_hash_table_replace(element_counts, g_strdup(symbol), GUINT_TO_POINTER(count));

          x_bohr = (atom->position[0] + shift_angstrom[0]) * GDIS_ANGSTROM_TO_BOHR;
          y_bohr = (atom->position[1] + shift_angstrom[1]) * GDIS_ANGSTROM_TO_BOHR;
          z_bohr = (atom->position[2] + shift_angstrom[2]) * GDIS_ANGSTROM_TO_BOHR;

          g_string_append_printf(deck,
                                 "atom %s%u %s %.10g %.10g %.10g\n",
                                 symbol,
                                 count,
                                 alias,
                                 x_bohr,
                                 y_bohr,
                                 z_bohr);
        }
      g_hash_table_unref(element_counts);
      g_string_append_c(deck, '\n');
    }

  if (have_placeholders)
    g_string_append(deck, "# Replace any <...> placeholders before running Qbox.\n");
  g_string_append_printf(deck, "set xc %s\n", xc_name);
  g_string_append_printf(deck, "set ecut %.1f\n", gtk_spin_button_get_value(GTK_SPIN_BUTTON(tool->ecut_spin)));
  g_string_append_printf(deck, "set scf_tol %s\n", scf_tol);
  g_string_append_printf(deck, "set wf_dyn %s\n", wf_dyn);
  g_string_append_printf(deck, "set atoms_dyn %s\n", atoms_dyn);
  g_string_append(deck, "set cell_dyn LOCKED\n");
  if (charge != 0)
    g_string_append_printf(deck, "set net_charge %d\n", charge);
  if (randomize_wf && !(restart_reference && restart_reference[0]))
    g_string_append(deck, "randomize_wf\n");
  g_string_append_printf(deck,
                         "run %s%u %u %u\n",
                         use_atomic_density ? "-atomic_density " : "",
                         locked_atoms ? 0u : ionic_steps,
                         scf_steps,
                         density_update);
  if (gtk_editable_get_text(GTK_EDITABLE(tool->save_entry))[0])
    g_string_append_printf(deck, "save %s\n", gtk_editable_get_text(GTK_EDITABLE(tool->save_entry)));

  *deck_out = g_string_free(deck, FALSE);
  return TRUE;
}

static gchar *
gdis_qbox_build_monkhorst_pack_block(guint nx, guint ny, guint nz)
{
  GString *block;
  gdouble weight;

  nx = MAX(nx, 1u);
  ny = MAX(ny, 1u);
  nz = MAX(nz, 1u);

  weight = 1.0 / ((gdouble) nx * (gdouble) ny * (gdouble) nz);
  block = g_string_new(
    "# Monkhorst-Pack mesh generated by the GTK4 Qbox helper.\n"
    "kpoint delete 0 0 0\n");

  for (guint ix = 1; ix <= nx; ix++)
    {
      gdouble kx = (2.0 * (gdouble) ix - (gdouble) nx - 1.0) / (2.0 * (gdouble) nx);

      for (guint iy = 1; iy <= ny; iy++)
        {
          gdouble ky = (2.0 * (gdouble) iy - (gdouble) ny - 1.0) / (2.0 * (gdouble) ny);

          for (guint iz = 1; iz <= nz; iz++)
            {
              gdouble kz = (2.0 * (gdouble) iz - (gdouble) nz - 1.0) / (2.0 * (gdouble) nz);

              g_string_append_printf(block,
                                     "kpoint add %.10f %.10f %.10f %.16f\n",
                                     kx,
                                     ky,
                                     kz,
                                     weight);
            }
        }
    }

  return g_string_free(block, FALSE);
}

static gchar *
gdis_qbox_apply_monkhorst_pack_to_text(const char *deck_text,
                                       guint nx,
                                       guint ny,
                                       guint nz)
{
  gchar **lines;
  g_autofree gchar *block = NULL;
  GString *result;
  gboolean inserted;

  block = gdis_qbox_build_monkhorst_pack_block(nx, ny, nz);
  lines = g_strsplit(deck_text ? deck_text : "", "\n", -1);
  result = g_string_new("");
  inserted = FALSE;

  for (guint i = 0; lines[i] != NULL; i++)
    {
      g_autofree gchar *trimmed = NULL;
      const gchar *line = lines[i];

      trimmed = g_strdup(line);
      g_strstrip(trimmed);

      if (trimmed[0] != '\0' && g_str_has_prefix(trimmed, "kpoint "))
        continue;

      if (!inserted &&
          (g_str_has_prefix(trimmed, "randomize_wf") ||
           g_str_has_prefix(trimmed, "run ") ||
           g_str_has_prefix(trimmed, "save ")))
        {
          if (result->len > 0 && result->str[result->len - 1] != '\n')
            g_string_append_c(result, '\n');
          g_string_append(result, block);
          inserted = TRUE;
        }

      g_string_append(result, line);
      g_string_append_c(result, '\n');
    }

  if (!inserted)
    {
      if (result->len > 0 && result->str[result->len - 1] != '\n')
        g_string_append_c(result, '\n');
      if (result->len > 0)
        g_string_append_c(result, '\n');
      g_string_append(result, block);
    }

  g_strfreev(lines);
  return g_string_free(result, FALSE);
}

static gchar *
gdis_qbox_text_buffer_contents(GtkTextBuffer *buffer)
{
  GtkTextIter start;
  GtkTextIter end;

  g_return_val_if_fail(GTK_IS_TEXT_BUFFER(buffer), NULL);

  gtk_text_buffer_get_bounds(buffer, &start, &end);
  return gtk_text_buffer_get_text(buffer, &start, &end, FALSE);
}

static void
gdis_qbox_set_report(GtkTextBuffer *buffer, const char *text)
{
  g_return_if_fail(GTK_IS_TEXT_BUFFER(buffer));

  gtk_text_buffer_set_text(buffer,
                           text && text[0] != '\0' ? text :
                           "No Qbox report yet.",
                           -1);
}

static void
on_qbox_deck_buffer_changed(GtkTextBuffer *buffer, gpointer user_data)
{
  GdisQboxTool *tool;
  g_autofree gchar *current_text = NULL;

  (void) buffer;

  tool = user_data;
  if (!tool || tool->suppress_input_signal)
    return;

  current_text = gdis_qbox_text_buffer_contents(tool->deck_buffer);
  tool->editor_dirty = (g_strcmp0(current_text, tool->last_generated_input) != 0);
}

static void
gdis_qbox_append_text_block(GdisQboxTool *tool, const char *text)
{
  GtkTextIter end;

  g_return_if_fail(tool != NULL);
  g_return_if_fail(GTK_IS_TEXT_BUFFER(tool->deck_buffer));
  g_return_if_fail(text != NULL);

  tool->suppress_input_signal = TRUE;
  gtk_text_buffer_get_end_iter(tool->deck_buffer, &end);
  if (!gtk_text_iter_is_start(&end))
    gtk_text_buffer_insert(tool->deck_buffer, &end, "\n\n", -1);
  gtk_text_buffer_insert(tool->deck_buffer, &end, text, -1);
  tool->suppress_input_signal = FALSE;
  tool->editor_dirty = TRUE;
}

static void
on_qbox_apply_mp_mesh_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;
  guint nx;
  guint ny;
  guint nz;
  g_autofree gchar *deck_text = NULL;
  g_autofree gchar *updated_text = NULL;
  g_autofree gchar *report = NULL;

  (void) button;

  self = user_data;
  if (!self || !self->qbox_tool)
    return;

  tool = self->qbox_tool;
  nx = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->kmesh_spins[0]));
  ny = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->kmesh_spins[1]));
  nz = (guint) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(tool->kmesh_spins[2]));
  deck_text = gdis_qbox_text_buffer_contents(tool->deck_buffer);
  updated_text = gdis_qbox_apply_monkhorst_pack_to_text(deck_text, nx, ny, nz);

  tool->suppress_input_signal = TRUE;
  gtk_text_buffer_set_text(tool->deck_buffer, updated_text, -1);
  tool->suppress_input_signal = FALSE;
  tool->editor_dirty = TRUE;

  report = g_strdup_printf("Applied a %ux%ux%u Monkhorst-Pack mesh to the current deck.\n"
                           "Existing kpoint lines were replaced and the generated block was inserted before the first run/save command.",
                           nx,
                           ny,
                           nz);
  gdis_qbox_set_report(tool->report_buffer, report);
  gdis_gtk4_window_log(self,
                       "Inserted a %ux%ux%u Monkhorst-Pack mesh into the Qbox deck.\n",
                       nx,
                       ny,
                       nz);
}

static void
on_qbox_append_hw5_phonons_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;
  GString *block;

  (void) button;

  self = user_data;
  tool = self ? self->qbox_tool : NULL;
  if (!tool || !self || !self->active_model || !self->active_model->atoms)
    return;

  block = g_string_new(
    "# Homework 5 frozen-phonon finite-difference block\n"
    "# Displacement step: 0.05 bohr\n");

  for (guint i = 0; i < self->active_model->atoms->len; i++)
    {
      const GdisAtom *atom;
      const char *atom_name;

      atom = g_ptr_array_index(self->active_model->atoms, i);
      atom_name = (atom->label && atom->label[0]) ? atom->label : atom->element;
      if (!atom_name || !atom_name[0])
        atom_name = "atom";

      g_string_append_printf(block,
                             "\n# %s along x\n"
                             "move %s by 0.05 0 0\n"
                             "run 0 40\n"
                             "move %s by -0.05 0 0\n"
                             "move %s by -0.05 0 0\n"
                             "run 0 40\n"
                             "move %s by 0.05 0 0\n",
                             atom_name,
                             atom_name,
                             atom_name,
                             atom_name,
                             atom_name);
      g_string_append_printf(block,
                             "\n# %s along y\n"
                             "move %s by 0 0.05 0\n"
                             "run 0 40\n"
                             "move %s by 0 -0.05 0\n"
                             "move %s by 0 -0.05 0\n"
                             "run 0 40\n"
                             "move %s by 0 0.05 0\n",
                             atom_name,
                             atom_name,
                             atom_name,
                             atom_name,
                             atom_name);
      g_string_append_printf(block,
                             "\n# %s along z\n"
                             "move %s by 0 0 0.05\n"
                             "run 0 40\n"
                             "move %s by 0 0 -0.05\n"
                             "move %s by 0 0 -0.05\n"
                             "run 0 40\n"
                             "move %s by 0 0 0.05\n",
                             atom_name,
                             atom_name,
                             atom_name,
                             atom_name,
                             atom_name);
    }

  gdis_qbox_append_text_block(tool, block->str);
  gdis_qbox_set_report(tool->report_buffer,
                       "Appended the Homework 5 frozen-phonon displacement block.\n"
                       "This follows the original +/- 0.05 bohr pattern for every atom and axis, ready to be run from the current Qbox deck.");
  gdis_gtk4_window_log(self,
                       "Appended the Homework 5 frozen-phonon block to the Qbox deck.\n");
  g_string_free(block, TRUE);
}

static void
on_qbox_append_hw5_homo_lumo_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;
  g_autofree gchar *restart_reference = NULL;
  g_autofree gchar *save_text = NULL;
  g_autofree gchar *last_save_basename = NULL;
  g_autofree gchar *load_ref = NULL;
  GString *block;

  (void) button;

  self = user_data;
  tool = self ? self->qbox_tool : NULL;
  if (!tool || !self)
    return;

  restart_reference = gdis_qbox_restart_reference(tool);
  save_text = gdis_qbox_entry_text(tool->save_entry);
  if (self->qbox_last_save_path && self->qbox_last_save_path[0])
    last_save_basename = g_path_get_basename(self->qbox_last_save_path);

  if (restart_reference && restart_reference[0] &&
      !g_str_has_prefix(restart_reference, "<"))
    load_ref = g_strdup(restart_reference);
  else if (save_text && save_text[0])
    load_ref = g_strdup(save_text);
  else if (last_save_basename && last_save_basename[0])
    load_ref = g_strdup(last_save_basename);
  else
    load_ref = g_strdup("<set-restart-xml>");

  block = g_string_new(
    "# Homework 5 HOMO/LUMO continuation block\n");
  g_string_append_printf(block,
                         "load %s\n"
                         "set nempty 1\n"
                         "set wf_dyn JD\n"
                         "run 0 120\n"
                         "plot -wf 4 HOMO.cube\n"
                         "plot -wf 5 LUMO.cube\n"
                         "save ch4_homo_lumo.xml\n",
                         load_ref);

  gdis_qbox_append_text_block(tool, block->str);
  gdis_qbox_set_report(tool->report_buffer,
                       "Appended the Homework 5 HOMO/LUMO continuation block.\n"
                       "The block loads the current restart XML, adds one empty state, switches to JD, and writes HOMO.cube and LUMO.cube.");
  gdis_gtk4_window_log(self,
                       "Appended the Homework 5 HOMO/LUMO block to the Qbox deck.\n");
  g_string_free(block, TRUE);
}

static void
gdis_qbox_set_last_summary(GdisGtk4Window *self, const char *summary)
{
  g_return_if_fail(self != NULL);

  g_free(self->qbox_last_summary);
  self->qbox_last_summary = g_strdup(summary ? summary : "");
}

static void
gdis_qbox_clear_last_run_state(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);

  g_clear_pointer(&self->qbox_last_workdir, g_free);
  g_clear_pointer(&self->qbox_last_input_path, g_free);
  g_clear_pointer(&self->qbox_last_output_path, g_free);
  g_clear_pointer(&self->qbox_last_stderr_path, g_free);
  g_clear_pointer(&self->qbox_last_save_path, g_free);
}

static void
gdis_qbox_set_last_run_state(GdisGtk4Window *self,
                             const char *workdir,
                             const char *input_path,
                             const char *output_path,
                             const char *stderr_path,
                             const char *save_path)
{
  g_return_if_fail(self != NULL);

  gdis_qbox_clear_last_run_state(self);
  self->qbox_last_workdir = g_strdup(workdir);
  self->qbox_last_input_path = g_strdup(input_path);
  self->qbox_last_output_path = g_strdup(output_path);
  self->qbox_last_stderr_path = g_strdup(stderr_path);
  self->qbox_last_save_path = g_strdup(save_path);
}

static gboolean
gdis_qbox_import_result_into_active_model(GdisGtk4Window *self,
                                          const char *result_path,
                                          GError **error)
{
  g_autofree gchar *resolved_path = NULL;
  g_autoptr(GdisModel) loaded = NULL;

  g_return_val_if_fail(self != NULL, FALSE);
  g_return_val_if_fail(result_path != NULL, FALSE);

  if (!self->active_model)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "There is no active model to replace with the Qbox result.");
      return FALSE;
    }

  resolved_path = gdis_gtk4_window_resolve_path(result_path);
  if (!resolved_path)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not locate the Qbox result XML '%s'.",
                  result_path);
      return FALSE;
    }

  loaded = gdis_model_load(resolved_path, error);
  if (!loaded)
    return FALSE;

  if (!gdis_model_copy_from(self->active_model, loaded, error))
    return FALSE;

  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_set_active_model(self, self->active_model);
  gdis_gtk4_window_log(self, "Imported Qbox result into the active model: %s\n", resolved_path);
  return TRUE;
}

static gchar *
gdis_qbox_build_results_report(GdisGtk4Window *self)
{
  GString *report;
  g_autofree gchar *stdout_text = NULL;
  g_autofree gchar *stderr_text = NULL;
  g_autofree gchar *stdout_tail = NULL;
  g_autofree gchar *stderr_tail = NULL;
  g_autofree gchar *last_etotal = NULL;

  g_return_val_if_fail(self != NULL, g_strdup("No Qbox results available."));

  report = g_string_new("");
  if (!self->qbox_last_output_path && !self->qbox_last_save_path)
    {
      g_string_append(report,
                      "No completed Qbox run is recorded in this session yet.\n"
                      "Run a Qbox job first, then reopen Results.");
      return g_string_free(report, FALSE);
    }

  g_string_append(report, "Qbox Session Results\n");
  if (self->qbox_last_summary && self->qbox_last_summary[0])
    g_string_append_printf(report, "Last activity: %s\n", self->qbox_last_summary);
  if (self->qbox_last_workdir)
    g_string_append_printf(report, "Working directory: %s\n", self->qbox_last_workdir);
  if (self->qbox_last_input_path)
    g_string_append_printf(report, "Input deck: %s\n", self->qbox_last_input_path);
  if (self->qbox_last_output_path)
    g_string_append_printf(report, "stdout report: %s\n", self->qbox_last_output_path);
  if (self->qbox_last_stderr_path)
    g_string_append_printf(report, "stderr report: %s\n", self->qbox_last_stderr_path);
  if (self->qbox_last_save_path)
    g_string_append_printf(report, "Saved XML: %s%s\n",
                           self->qbox_last_save_path,
                           g_file_test(self->qbox_last_save_path, G_FILE_TEST_EXISTS) ? "" : " (missing)");

  if (self->qbox_last_output_path)
    g_file_get_contents(self->qbox_last_output_path, &stdout_text, NULL, NULL);
  if (self->qbox_last_stderr_path)
    g_file_get_contents(self->qbox_last_stderr_path, &stderr_text, NULL, NULL);

  last_etotal = gdis_qbox_extract_last_etotal(stdout_text);
  if (last_etotal && last_etotal[0])
    g_string_append_printf(report, "Last <etotal>: %s\n", last_etotal);

  stdout_tail = gdis_qbox_tail_text(stdout_text, 12000u);
  stderr_tail = gdis_qbox_tail_text(stderr_text, 6000u);

  if (stdout_tail && stdout_tail[0])
    g_string_append_printf(report, "\nRecent stdout:\n%s", stdout_tail);
  if (stderr_tail && stderr_tail[0])
    g_string_append_printf(report, "\n\nRecent stderr:\n%s", stderr_tail);

  return g_string_free(report, FALSE);
}

static gchar *
gdis_qbox_make_unique_continue_name(const char *workdir,
                                    const char *base_stem,
                                    const char *extension)
{
  guint index;

  g_return_val_if_fail(base_stem != NULL, NULL);
  g_return_val_if_fail(extension != NULL, NULL);

  for (index = 1u; index < 1000u; index++)
    {
      g_autofree gchar *candidate = NULL;
      g_autofree gchar *candidate_path = NULL;

      if (index == 1u)
        candidate = g_strdup_printf("%s-continue.%s", base_stem, extension);
      else
        candidate = g_strdup_printf("%s-continue-%u.%s", base_stem, index, extension);

      if (!workdir || !workdir[0])
        return g_steal_pointer(&candidate);

      candidate_path = g_build_filename(workdir, candidate, NULL);
      if (!g_file_test(candidate_path, G_FILE_TEST_EXISTS))
        return g_steal_pointer(&candidate);
    }

  return g_strdup_printf("%s-continue-final.%s", base_stem, extension);
}

static void
on_qbox_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisQboxTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->qbox_tool = NULL;
  g_clear_pointer(&tool->species_rows, g_ptr_array_unref);
  g_clear_pointer(&tool->last_generated_input, g_free);
  g_free(tool);
}

static void
gdis_gtk4_window_refresh_qbox_tool(GdisGtk4Window *self)
{
  GdisQboxTool *tool;
  GString *summary;
  guint unique_species;
  gboolean model_changed;

  g_return_if_fail(self != NULL);

  tool = self->qbox_tool;
  if (!tool)
    return;

  model_changed = (tool->source_model != self->active_model);
  if (model_changed)
    {
      g_autofree gchar *default_job = gdis_qbox_model_base_name(self->active_model);
      g_autofree gchar *job_slug = NULL;
      g_autofree gchar *default_workdir = NULL;
      g_autofree gchar *default_exec = NULL;

      gtk_editable_set_text(GTK_EDITABLE(tool->job_entry), default_job);
      job_slug = gdis_qbox_job_slug(default_job);
      default_workdir = gdis_qbox_default_workdir(self, job_slug);
      gtk_editable_set_text(GTK_EDITABLE(tool->workdir_entry), default_workdir);

      {
        g_autofree gchar *input_name = g_strdup_printf("%s.i", job_slug);
        g_autofree gchar *output_name = g_strdup_printf("%s.r", job_slug);
        g_autofree gchar *save_name = g_strdup_printf("%s.xml", job_slug);
        gtk_editable_set_text(GTK_EDITABLE(tool->input_entry), input_name);
        gtk_editable_set_text(GTK_EDITABLE(tool->output_entry), output_name);
        gtk_editable_set_text(GTK_EDITABLE(tool->save_entry), save_name);
      }

      gtk_editable_set_text(GTK_EDITABLE(tool->restart_entry), "");
      gtk_editable_set_text(GTK_EDITABLE(tool->launcher_entry), "");
      {
        g_autofree gchar *default_pseudo_source = gdis_qbox_default_pseudo_source();
        gtk_editable_set_text(GTK_EDITABLE(tool->pseudo_dir_entry), default_pseudo_source);
      }
      gtk_editable_set_text(GTK_EDITABLE(tool->xc_entry), "PBE");
      gtk_editable_set_text(GTK_EDITABLE(tool->wf_dyn_entry), "PSDA");
      gtk_editable_set_text(GTK_EDITABLE(tool->atoms_dyn_entry), "LOCKED");
      gtk_editable_set_text(GTK_EDITABLE(tool->scf_tol_entry), "1.0e-8");

      gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->ecut_spin), 35.0);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->charge_spin), 0.0);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->padding_spin), 8.0);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->ionic_steps_spin), 0.0);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->scf_steps_spin), 40.0);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->density_update_spin), 10.0);
      gtk_check_button_set_active(GTK_CHECK_BUTTON(tool->use_cell_toggle),
                                  self->active_model && self->active_model->periodic);
      gtk_check_button_set_active(GTK_CHECK_BUTTON(tool->atomic_density_toggle),
                                  !(self->active_model && self->active_model->periodic));
      gtk_check_button_set_active(GTK_CHECK_BUTTON(tool->randomize_toggle), TRUE);

      default_exec = gdis_qbox_lookup_executable_path(self);
      gtk_editable_set_text(GTK_EDITABLE(tool->exec_entry), default_exec ? default_exec : "");

      gdis_qbox_tool_rebuild_species_rows(self, tool);
      tool->source_model = self->active_model;
    }
  else if (!gtk_editable_get_text(GTK_EDITABLE(tool->exec_entry))[0])
    {
      g_autofree gchar *detected_exec = gdis_qbox_lookup_executable_path(self);
      gtk_editable_set_text(GTK_EDITABLE(tool->exec_entry), detected_exec ? detected_exec : "");
    }

  gtk_widget_set_sensitive(tool->use_cell_toggle,
                           self->active_model && self->active_model->periodic);
  if (!(self->active_model && self->active_model->periodic))
    gtk_check_button_set_active(GTK_CHECK_BUTTON(tool->use_cell_toggle), FALSE);

  if (self->active_model && (!tool->editor_dirty || model_changed))
    {
      g_autofree gchar *deck = NULL;
      g_autofree gchar *cell_mode = NULL;
      GError *error = NULL;

      if (gdis_qbox_build_input_deck(self, tool, &deck, &cell_mode, &error))
        gdis_qbox_tool_set_generated_text(tool, deck);
      g_clear_error(&error);
    }

  summary = g_string_new("");
  if (!self->active_model)
    {
      g_string_append(summary,
                      "No active model.\n"
                      "The current Qbox deck remains editable, but Regenerate needs a loaded structure.");
    }
  else
    {
      unique_species = gdis_qbox_count_unique_species(self->active_model);
      g_string_append_printf(summary,
                             "Active model: %s\nAtoms: %u   Unique species: %u\nCell source: %s\nExecutable: %s\nDeck status: %s",
                             self->active_model->basename,
                             self->active_model->atom_count,
                             unique_species,
                             gtk_check_button_get_active(GTK_CHECK_BUTTON(tool->use_cell_toggle)) &&
                             self->active_model->periodic ? "Model periodic cell" : "Generated bounding box",
                             gtk_editable_get_text(GTK_EDITABLE(tool->exec_entry))[0] ?
                               gtk_editable_get_text(GTK_EDITABLE(tool->exec_entry)) :
                               "(not set)",
                             tool->editor_dirty ? "manual edits preserved" : "synchronized");
      if (self->qbox_last_summary && self->qbox_last_summary[0])
        g_string_append_printf(summary, "\nLast activity: %s", self->qbox_last_summary);
      if (self->qbox_last_save_path && self->qbox_last_save_path[0])
        g_string_append_printf(summary, "\nLast XML: %s%s",
                               self->qbox_last_save_path,
                               g_file_test(self->qbox_last_save_path, G_FILE_TEST_EXISTS) ? "" : " (missing)");
      if (self->qbox_last_output_path && self->qbox_last_output_path[0])
        g_string_append_printf(summary, "\nLast report: %s", self->qbox_last_output_path);
    }
  gtk_label_set_text(tool->summary_label, summary->str);
  g_string_free(summary, TRUE);
}

static void
on_qbox_regenerate_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;
  g_autofree gchar *deck = NULL;
  g_autofree gchar *cell_mode = NULL;
  GError *error;
  GString *report;

  (void) button;

  self = user_data;
  tool = self ? self->qbox_tool : NULL;
  if (!tool)
    return;

  error = NULL;
  if (!gdis_qbox_build_input_deck(self, tool, &deck, &cell_mode, &error))
    {
      gdis_qbox_set_report(tool->report_buffer, error ? error->message : "Qbox deck generation failed.");
      gdis_gtk4_window_log(self, "Qbox deck generation failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_qbox_tool_set_generated_text(tool, deck);

  report = g_string_new("");
  g_string_append_printf(report,
                         "Generated a Qbox starter deck for %s.\n"
                         "Cell source: %s\n"
                         "Input file: %s\n"
                         "Output capture: %s\n"
                         "Restart save: %s\n\n"
                         "The editor below is live: you can keep this generated deck or customize it before writing or running.",
                         self->active_model ? self->active_model->basename : "(none)",
                         cell_mode ? cell_mode : "unknown",
                         gtk_editable_get_text(GTK_EDITABLE(tool->input_entry)),
                         gtk_editable_get_text(GTK_EDITABLE(tool->output_entry)),
                         gtk_editable_get_text(GTK_EDITABLE(tool->save_entry)));
  gdis_qbox_set_report(tool->report_buffer, report->str);
  g_string_free(report, TRUE);
  gdis_gtk4_window_refresh_qbox_tool(self);
  gdis_gtk4_window_log(self, "Generated the Qbox deck for %s.\n",
                       self->active_model ? self->active_model->basename : "the active model");
}

static void
on_qbox_detect_exec_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  g_autofree gchar *path = NULL;

  (void) button;

  self = user_data;
  if (!self || !self->qbox_tool)
    return;

  path = gdis_qbox_lookup_executable_path(self);
  gtk_editable_set_text(GTK_EDITABLE(self->qbox_tool->exec_entry), path ? path : "");
  gdis_gtk4_window_refresh_qbox_tool(self);
}

static void
on_qbox_exec_paths_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;

  (void) button;

  self = user_data;
  if (!self)
    return;

  gdis_gtk4_window_present_executable_paths_tool(self, "qbox");
}

static gboolean
gdis_qbox_copy_local_file(const char *source_path,
                          const char *dest_path,
                          GError **error)
{
  GFile *source;
  GFile *dest;
  gboolean ok;

  if (g_strcmp0(source_path, dest_path) == 0)
    return TRUE;

  source = g_file_new_for_path(source_path);
  dest = g_file_new_for_path(dest_path);
  ok = g_file_copy(source, dest, G_FILE_COPY_OVERWRITE, NULL, NULL, NULL, error);
  g_object_unref(source);
  g_object_unref(dest);
  return ok;
}

static void
gdis_qbox_build_paths(GdisQboxTool *tool,
                      gchar **workdir_out,
                      gchar **input_path_out,
                      gchar **input_name_out,
                      gchar **output_path_out,
                      gchar **stderr_path_out,
                      GError **error)
{
  g_autofree gchar *workdir_text = NULL;
  g_autofree gchar *input_name = NULL;
  g_autofree gchar *output_name = NULL;
  g_autofree gchar *workdir = NULL;

  g_return_if_fail(tool != NULL);

  workdir_text = gdis_qbox_entry_text(tool->workdir_entry);
  input_name = gdis_qbox_entry_text(tool->input_entry);
  output_name = gdis_qbox_entry_text(tool->output_entry);
  if (!workdir_text[0] || !input_name[0] || !output_name[0])
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Qbox working directory, input filename, and output filename must not be empty.");
      return;
    }

  workdir = g_canonicalize_filename(workdir_text, g_get_current_dir());
  if (g_mkdir_with_parents(workdir, 0755) != 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not create '%s': %s",
                  workdir,
                  g_strerror(errno));
      return;
    }

  if (workdir_out)
    *workdir_out = g_steal_pointer(&workdir);
  if (input_path_out)
    *input_path_out = g_build_filename(workdir_out ? *workdir_out : workdir, input_name, NULL);
  if (input_name_out)
    *input_name_out = g_steal_pointer(&input_name);
  if (output_path_out)
    *output_path_out = g_build_filename(workdir_out ? *workdir_out : workdir, output_name, NULL);
  if (stderr_path_out)
    {
      g_autofree gchar *stderr_name = g_strdup_printf("%s.stderr", output_name);
      *stderr_path_out = g_build_filename(workdir_out ? *workdir_out : workdir, stderr_name, NULL);
    }
}

static void
gdis_qbox_validate_run_buffer(const char *deck_text, GError **error)
{
  if (!deck_text || !deck_text[0])
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "The Qbox deck editor is empty.");
      return;
    }

  if (strstr(deck_text, "<set-pseudo-for-") ||
      strstr(deck_text, "<missing-pseudo-for-") ||
      strstr(deck_text, "<missing-restart-xml>"))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "The Qbox deck still contains unresolved placeholders. Fill the species table or edit the deck before running.");
    }
}

static gboolean
gdis_qbox_stage_assets(GdisGtk4Window *self,
                       GdisQboxTool *tool,
                       const char *workdir,
                       guint *copied_files_out,
                       GString *warnings,
                       GError **error)
{
  guint copied_files;
  g_autofree gchar *restart_text = NULL;

  g_return_val_if_fail(self != NULL, FALSE);
  g_return_val_if_fail(tool != NULL, FALSE);
  g_return_val_if_fail(workdir != NULL, FALSE);

  copied_files = 0u;

  restart_text = gdis_qbox_entry_text(tool->restart_entry);
  if (restart_text[0] && !gdis_qbox_is_remote_uri(restart_text))
    {
      g_autofree gchar *resolved = NULL;
      g_autofree gchar *dest = NULL;

      resolved = gdis_gtk4_window_resolve_path(restart_text);
      if (!resolved)
        {
          g_set_error(error,
                      GDIS_MODEL_ERROR,
                      GDIS_MODEL_ERROR_IO,
                      "Could not locate the restart XML '%s'.",
                      restart_text);
          return FALSE;
        }
      dest = g_build_filename(workdir, "restart.xml", NULL);
      if (!gdis_qbox_copy_local_file(resolved, dest, error))
        return FALSE;
      copied_files++;
    }

  for (guint i = 0; tool->species_rows && i < tool->species_rows->len; i++)
    {
      GdisQboxSpeciesRow *row;
      g_autofree gchar *pseudo_text = NULL;
      g_autofree gchar *resolved = NULL;
      g_autofree gchar *dest_name = NULL;
      g_autofree gchar *dest = NULL;

      row = g_ptr_array_index(tool->species_rows, i);
      pseudo_text = gdis_qbox_entry_text(row->pseudo_entry);
      if (!pseudo_text[0] || gdis_qbox_is_remote_uri(pseudo_text))
        continue;

      resolved = gdis_gtk4_window_resolve_path(pseudo_text);
      if (!resolved)
        {
          if (warnings)
            g_string_append_printf(warnings,
                                   "Could not stage pseudo file for %s from '%s'.\n",
                                   row->element,
                                   pseudo_text);
          continue;
        }

      dest_name = gdis_qbox_species_reference(row);
      if (dest_name[0] == '<')
        continue;
      dest = g_build_filename(workdir, dest_name, NULL);
      if (!gdis_qbox_copy_local_file(resolved, dest, error))
        return FALSE;
      copied_files++;
    }

  if (copied_files_out)
    *copied_files_out = copied_files;
  return TRUE;
}

static gboolean
gdis_qbox_write_prepared_input(GdisGtk4Window *self,
                               GdisQboxTool *tool,
                               gboolean validate_for_run,
                               gchar **workdir_out,
                               gchar **input_path_out,
                               gchar **input_name_out,
                               gchar **output_path_out,
                               gchar **stderr_path_out,
                               guint *copied_files_out,
                               guint *localized_pseudos_out,
                               GString *warnings,
                               GError **error)
{
  g_autofree gchar *deck_text = NULL;
  GError *local_error = NULL;
  guint localized_pseudos = 0u;

  g_return_val_if_fail(self != NULL, FALSE);
  g_return_val_if_fail(tool != NULL, FALSE);

  if (!gdis_qbox_prepare_local_pseudo_library(self,
                                              tool,
                                              &localized_pseudos,
                                              NULL,
                                              warnings,
                                              &local_error))
    {
      g_propagate_error(error, local_error);
      return FALSE;
    }

  if (!tool->editor_dirty && self->active_model)
    {
      g_autofree gchar *deck = NULL;
      g_autofree gchar *cell_mode = NULL;

      if (gdis_qbox_build_input_deck(self, tool, &deck, &cell_mode, &local_error))
        {
          gdis_qbox_tool_set_generated_text(tool, deck);
          g_clear_error(&local_error);
        }
    }

  deck_text = gdis_qbox_text_buffer_contents(tool->deck_buffer);
  if ((!deck_text || !deck_text[0]) && self->active_model)
    {
      g_autofree gchar *deck = NULL;
      g_autofree gchar *cell_mode = NULL;

      if (gdis_qbox_build_input_deck(self, tool, &deck, &cell_mode, &local_error))
        {
          gdis_qbox_tool_set_generated_text(tool, deck);
          g_clear_error(&local_error);
          g_free(deck_text);
          deck_text = gdis_qbox_text_buffer_contents(tool->deck_buffer);
        }
    }

  if (validate_for_run)
    {
      gdis_qbox_validate_run_buffer(deck_text, &local_error);
      if (local_error)
        {
          g_propagate_error(error, local_error);
          return FALSE;
        }
    }

  if (!deck_text || !deck_text[0])
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "The Qbox deck editor is empty. Regenerate a deck or type one before writing.");
      return FALSE;
    }

  gdis_qbox_build_paths(tool,
                        workdir_out,
                        input_path_out,
                        input_name_out,
                        output_path_out,
                        stderr_path_out,
                        &local_error);
  if (local_error)
    {
      g_propagate_error(error, local_error);
      return FALSE;
    }

  if (!gdis_qbox_stage_assets(self,
                              tool,
                              *workdir_out,
                              copied_files_out,
                              warnings,
                              &local_error))
    {
      g_propagate_error(error, local_error);
      return FALSE;
    }

  if (!g_file_set_contents(*input_path_out, deck_text, -1, &local_error))
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not write '%s': %s",
                  *input_path_out,
                  local_error ? local_error->message : "unknown error");
      g_clear_error(&local_error);
      return FALSE;
    }

  if (localized_pseudos_out)
    *localized_pseudos_out = localized_pseudos;
  return TRUE;
}

static gboolean
gdis_qbox_run_shell_capture(const char *working_directory,
                            const char *command,
                            const char *stdout_path,
                            const char *stderr_path,
                            gchar **stdout_text_out,
                            gchar **stderr_text_out,
                            GError **error)
{
  g_autofree gchar *quoted_stdout = NULL;
  g_autofree gchar *quoted_stderr = NULL;
  g_autofree gchar *redirected_command = NULL;
  const gchar *argv[] = {"/bin/sh", "-lc", NULL, NULL};
  gchar *shell_stdout = NULL;
  gchar *shell_stderr = NULL;
  gchar *stdout_text = NULL;
  gchar *stderr_text = NULL;
  gint wait_status = 0;

  g_return_val_if_fail(command != NULL, FALSE);
  g_return_val_if_fail(stdout_path != NULL, FALSE);
  g_return_val_if_fail(stderr_path != NULL, FALSE);

  quoted_stdout = g_shell_quote(stdout_path);
  quoted_stderr = g_shell_quote(stderr_path);
  redirected_command = g_strdup_printf("exec %s > %s 2> %s",
                                       command,
                                       quoted_stdout,
                                       quoted_stderr);
  argv[2] = redirected_command;

  g_file_set_contents(stdout_path, "", -1, NULL);
  g_file_set_contents(stderr_path, "", -1, NULL);

  if (!g_spawn_sync(working_directory,
                    (gchar **) argv,
                    NULL,
                    G_SPAWN_SEARCH_PATH,
                    NULL,
                    NULL,
                    &shell_stdout,
                    &shell_stderr,
                    &wait_status,
                    error))
    {
      g_free(shell_stdout);
      g_free(shell_stderr);
      return FALSE;
    }

  if (!g_file_get_contents(stdout_path, &stdout_text, NULL, NULL))
    stdout_text = g_strdup("");
  if (!g_file_get_contents(stderr_path, &stderr_text, NULL, NULL))
    stderr_text = g_strdup("");

  if (shell_stdout && shell_stdout[0])
    {
      g_autofree gchar *merged = g_strconcat(stdout_text ? stdout_text : "",
                                             shell_stdout,
                                             NULL);
      g_free(stdout_text);
      stdout_text = g_steal_pointer(&merged);
      g_file_set_contents(stdout_path, stdout_text, -1, NULL);
    }
  if (shell_stderr && shell_stderr[0])
    {
      g_autofree gchar *merged = g_strconcat(stderr_text ? stderr_text : "",
                                             shell_stderr,
                                             NULL);
      g_free(stderr_text);
      stderr_text = g_steal_pointer(&merged);
      g_file_set_contents(stderr_path, stderr_text, -1, NULL);
    }

  if (!g_spawn_check_wait_status(wait_status, error))
    {
      if (stdout_text_out)
        *stdout_text_out = stdout_text;
      else
        g_free(stdout_text);
      if (stderr_text_out)
        *stderr_text_out = stderr_text;
      else
        g_free(stderr_text);
      g_free(shell_stdout);
      g_free(shell_stderr);
      return FALSE;
    }

  if (stdout_text_out)
    *stdout_text_out = stdout_text;
  else
    g_free(stdout_text);
  if (stderr_text_out)
    *stderr_text_out = stderr_text;
  else
    g_free(stderr_text);
  g_free(shell_stdout);
  g_free(shell_stderr);
  return TRUE;
}

static gchar *
gdis_qbox_extract_last_etotal(const char *text)
{
  const gchar *start;
  const gchar *end;

  if (!text)
    return NULL;

  start = g_strrstr(text, "<etotal>");
  end = g_strrstr(text, "</etotal>");
  if (!start || !end || end <= start)
    return NULL;

  start += strlen("<etotal>");
  return g_strndup(start, (gsize) (end - start));
}

static gchar *
gdis_qbox_tail_text(const char *text, gsize max_chars)
{
  gsize length;

  if (!text)
    return g_strdup("");

  length = strlen(text);
  if (length <= max_chars)
    return g_strdup(text);
  return g_strdup(text + (length - max_chars));
}

static void
on_qbox_guess_pseudos_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;
  guint filled = 0u;

  (void) button;

  self = user_data;
  tool = self ? self->qbox_tool : NULL;
  if (!tool || !tool->species_rows)
    return;

  for (guint i = 0; i < tool->species_rows->len; i++)
    {
      GdisQboxSpeciesRow *row;
      g_autofree gchar *current = NULL;
      g_autofree gchar *guess = NULL;

      row = g_ptr_array_index(tool->species_rows, i);
      current = gdis_qbox_entry_text(row->pseudo_entry);
      if (current[0] && !g_str_has_prefix(current, "<"))
        continue;

      guess = gdis_qbox_guess_pseudo_path(gtk_editable_get_text(GTK_EDITABLE(tool->pseudo_dir_entry)),
                                          row->element,
                                          gtk_editable_get_text(GTK_EDITABLE(tool->xc_entry)));
      if (guess && guess[0])
        {
          gtk_editable_set_text(GTK_EDITABLE(row->pseudo_entry), guess);
          filled++;
        }
    }

  if (filled > 0u)
    {
      gdis_qbox_set_report(tool->report_buffer,
                           "Filled missing pseudo entries from the current pseudo dir/URL setting.\n"
                           "Regenerate the deck if you want the updated URIs written into the editor.");
      gdis_gtk4_window_log(self, "Updated %u Qbox pseudo mapping%s from the current pseudo dir/URL.\n",
                           filled,
                           filled == 1u ? "" : "s");
    }
  else
    {
      gdis_qbox_set_report(tool->report_buffer,
                           "No additional pseudo files were auto-detected.\n"
                           "Edit the species table directly if your library uses different names.");
    }
}

static void
on_qbox_setup_local_pseudos_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;
  GString *warnings;
  GError *error = NULL;
  guint downloaded = 0u;
  guint relinked = 0u;
  GString *report;

  (void) button;

  self = user_data;
  tool = self ? self->qbox_tool : NULL;
  if (!tool)
    return;

  warnings = g_string_new("");
  if (!gdis_qbox_prepare_local_pseudo_library(self,
                                              tool,
                                              &downloaded,
                                              &relinked,
                                              warnings,
                                              &error))
    {
      gdis_qbox_set_report(tool->report_buffer,
                           error ? error->message : "Could not prepare the local Qbox pseudo library.");
      gdis_gtk4_window_log(self, "Qbox local pseudo setup failed: %s\n",
                           error ? error->message : "unknown error");
      g_string_free(warnings, TRUE);
      g_clear_error(&error);
      return;
    }

  if (!tool->editor_dirty && self->active_model)
    on_qbox_regenerate_clicked(NULL, self);

  report = g_string_new("");
  g_string_append_printf(report,
                         "Prepared the local Qbox pseudo library.\n"
                         "Cache directory: %s\n"
                         "Species linked: %u\n"
                         "Fresh downloads: %u\n",
                         gtk_editable_get_text(GTK_EDITABLE(tool->pseudo_dir_entry)),
                         relinked,
                         downloaded);
  if (warnings->len > 0)
    g_string_append_printf(report, "\nWarnings:\n%s", warnings->str);
  gdis_qbox_set_report(tool->report_buffer, report->str);
  g_string_free(report, TRUE);
  g_string_free(warnings, TRUE);
  gdis_qbox_set_last_summary(self, downloaded > 0u ? "local pseudos downloaded" : "local pseudos ready");
  gdis_gtk4_window_refresh_qbox_tool(self);
  gdis_gtk4_window_log(self,
                       "Prepared the local Qbox pseudo library at %s (%u download%s).\n",
                       gtk_editable_get_text(GTK_EDITABLE(tool->pseudo_dir_entry)),
                       downloaded,
                       downloaded == 1u ? "" : "s");
}

static void
on_qbox_write_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;
  g_autofree gchar *workdir = NULL;
  g_autofree gchar *input_path = NULL;
  g_autofree gchar *input_name = NULL;
  g_autofree gchar *output_path = NULL;
  g_autofree gchar *stderr_path = NULL;
  g_autofree gchar *save_path = NULL;
  GString *warnings;
  GError *error = NULL;
  guint copied_files = 0u;
  guint localized_pseudos = 0u;
  GString *report;

  (void) button;

  self = user_data;
  tool = self ? self->qbox_tool : NULL;
  if (!tool)
    return;

  if (gtk_editable_get_text(GTK_EDITABLE(tool->exec_entry))[0])
    g_hash_table_replace(self->executable_paths,
                         g_strdup("qbox"),
                         g_strdup(gtk_editable_get_text(GTK_EDITABLE(tool->exec_entry))));

  warnings = g_string_new("");
  if (!gdis_qbox_write_prepared_input(self,
                                      tool,
                                      FALSE,
                                      &workdir,
                                      &input_path,
                                      &input_name,
                                      &output_path,
                                      &stderr_path,
                                      &copied_files,
                                      &localized_pseudos,
                                      warnings,
                                      &error))
    {
      gdis_qbox_set_report(tool->report_buffer,
                           error ? error->message : "Could not write the Qbox input deck.");
      gdis_gtk4_window_log(self, "Qbox input write failed: %s\n",
                           error ? error->message : "unknown error");
      g_string_free(warnings, TRUE);
      g_clear_error(&error);
      return;
    }

  report = g_string_new("");
  if (gtk_editable_get_text(GTK_EDITABLE(tool->save_entry))[0])
    save_path = g_build_filename(workdir,
                                 gtk_editable_get_text(GTK_EDITABLE(tool->save_entry)),
                                 NULL);
  g_string_append_printf(report,
                         "Wrote the Qbox input deck.\n"
                         "Working directory: %s\n"
                         "Input file: %s\n"
                         "Planned output capture: %s\n"
                         "Planned save XML: %s\n"
                         "Staged local assets: %u\n"
                         "Localized pseudos: %u\n",
                         workdir,
                         input_path,
                         output_path,
                         save_path ? save_path : "(not set)",
                         copied_files,
                         localized_pseudos);
  if (warnings->len > 0)
    g_string_append_printf(report, "\nWarnings:\n%s", warnings->str);
  gdis_qbox_set_report(tool->report_buffer, report->str);
  g_string_free(report, TRUE);
  g_string_free(warnings, TRUE);

  gdis_qbox_set_last_summary(self, "input deck written");
  gdis_qbox_set_last_run_state(self,
                               workdir,
                               input_path,
                               output_path,
                               stderr_path,
                               save_path);
  gdis_gtk4_window_refresh_qbox_tool(self);
  gdis_gtk4_window_refresh_executable_paths_tool(self);
  gdis_gtk4_window_refresh_task_manager_tool(self);
  gdis_gtk4_window_log(self, "Wrote the Qbox input deck: %s\n", input_path);
}

static void
on_qbox_use_last_xml_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;

  (void) button;

  self = user_data;
  tool = self ? self->qbox_tool : NULL;
  if (!tool)
    return;

  if (!self->qbox_last_save_path || !g_file_test(self->qbox_last_save_path, G_FILE_TEST_EXISTS))
    {
      gdis_qbox_set_report(tool->report_buffer,
                           "There is no saved Qbox XML from this session yet. Run a job that writes a save XML first.");
      return;
    }

  gtk_editable_set_text(GTK_EDITABLE(tool->restart_entry), self->qbox_last_save_path);
  gtk_check_button_set_active(GTK_CHECK_BUTTON(tool->randomize_toggle), FALSE);
  gdis_qbox_set_report(tool->report_buffer,
                       "Loaded the most recent Qbox save XML into Restart XML.\n"
                       "Randomize_wf was turned off for a continuation-friendly deck.");
  gdis_qbox_set_last_summary(self, "last XML selected for restart");
  gdis_gtk4_window_refresh_qbox_tool(self);
}

static void
on_qbox_continue_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;
  const char *workdir_text;
  g_autofree gchar *base_stem = NULL;
  g_autofree gchar *next_input = NULL;
  g_autofree gchar *next_output = NULL;
  g_autofree gchar *next_save = NULL;

  (void) button;

  self = user_data;
  tool = self ? self->qbox_tool : NULL;
  if (!tool)
    return;

  if (!self->qbox_last_save_path || !g_file_test(self->qbox_last_save_path, G_FILE_TEST_EXISTS))
    {
      gdis_qbox_set_report(tool->report_buffer,
                           "There is no saved Qbox XML to continue from yet. Run a job that writes a save XML first.");
      return;
    }

  workdir_text = gtk_editable_get_text(GTK_EDITABLE(tool->workdir_entry));
  if (self->qbox_last_workdir && self->qbox_last_workdir[0])
    gtk_editable_set_text(GTK_EDITABLE(tool->workdir_entry), self->qbox_last_workdir);

  if (gtk_editable_get_text(GTK_EDITABLE(tool->job_entry))[0])
    base_stem = gdis_qbox_job_slug(gtk_editable_get_text(GTK_EDITABLE(tool->job_entry)));
  else if (self->qbox_last_save_path)
    {
      g_autofree gchar *basename = g_path_get_basename(self->qbox_last_save_path);
      gchar *dot = strrchr(basename, '.');
      if (dot)
        *dot = '\0';
      base_stem = gdis_qbox_job_slug(basename);
    }
  if (!base_stem || !base_stem[0])
    base_stem = g_strdup("qbox_job");

  next_input = gdis_qbox_make_unique_continue_name(self->qbox_last_workdir ? self->qbox_last_workdir : workdir_text,
                                                   base_stem,
                                                   "i");
  next_output = gdis_qbox_make_unique_continue_name(self->qbox_last_workdir ? self->qbox_last_workdir : workdir_text,
                                                    base_stem,
                                                    "r");
  next_save = gdis_qbox_make_unique_continue_name(self->qbox_last_workdir ? self->qbox_last_workdir : workdir_text,
                                                  base_stem,
                                                  "xml");

  gtk_editable_set_text(GTK_EDITABLE(tool->restart_entry), self->qbox_last_save_path);
  gtk_editable_set_text(GTK_EDITABLE(tool->input_entry), next_input);
  gtk_editable_set_text(GTK_EDITABLE(tool->output_entry), next_output);
  gtk_editable_set_text(GTK_EDITABLE(tool->save_entry), next_save);
  {
    g_autofree gchar *job_name = g_strdup(next_save);
    gchar *dot = strrchr(job_name, '.');
    if (dot)
      *dot = '\0';
    gtk_editable_set_text(GTK_EDITABLE(tool->job_entry), job_name);
  }
  gtk_check_button_set_active(GTK_CHECK_BUTTON(tool->randomize_toggle), FALSE);

  on_qbox_regenerate_clicked(NULL, self);
  gdis_qbox_set_last_summary(self, "continue deck prepared");
  gdis_gtk4_window_log(self, "Prepared a Qbox continuation deck from: %s\n", self->qbox_last_save_path);
}

static void
on_qbox_import_result_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;
  GError *error = NULL;

  (void) button;

  self = user_data;
  tool = self ? self->qbox_tool : NULL;
  if (!tool)
    return;

  if (!self->qbox_last_save_path || !g_file_test(self->qbox_last_save_path, G_FILE_TEST_EXISTS))
    {
      gdis_qbox_set_report(tool->report_buffer,
                           "There is no saved Qbox XML to import yet. Run a job that writes a save XML first.");
      return;
    }

  if (!gdis_qbox_import_result_into_active_model(self, self->qbox_last_save_path, &error))
    {
      gdis_qbox_set_report(tool->report_buffer,
                           error ? error->message : "Could not import the saved Qbox XML.");
      gdis_gtk4_window_log(self, "Qbox result import failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_qbox_set_report(tool->report_buffer,
                       "Imported the most recent Qbox result XML into the active model.\n"
                       "The viewer and model summary were refreshed from the saved Qbox state.");
  gdis_qbox_set_last_summary(self, "result imported");
  gdis_gtk4_window_refresh_qbox_tool(self);
}

static void
on_qbox_results_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  g_autofree gchar *report = NULL;

  (void) button;

  self = user_data;
  if (!self)
    return;

  report = gdis_qbox_build_results_report(self);
  gdis_gtk4_window_present_report(self, "Qbox Results", report);
}

static void
on_qbox_run_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  GdisQboxTool *tool;
  g_autofree gchar *workdir = NULL;
  g_autofree gchar *input_path = NULL;
  g_autofree gchar *input_name = NULL;
  g_autofree gchar *output_path = NULL;
  g_autofree gchar *stderr_path = NULL;
  g_autofree gchar *save_path = NULL;
  g_autofree gchar *stdout_text = NULL;
  g_autofree gchar *stderr_text = NULL;
  g_autofree gchar *tail = NULL;
  g_autofree gchar *last_etotal = NULL;
  g_autofree gchar *launcher = NULL;
  g_autofree gchar *quoted_exec = NULL;
  g_autofree gchar *quoted_input = NULL;
  g_autofree gchar *command = NULL;
  GString *warnings;
  GError *error = NULL;
  guint copied_files = 0u;
  guint localized_pseudos = 0u;
  gboolean success;
  GString *report;

  (void) button;

  self = user_data;
  tool = self ? self->qbox_tool : NULL;
  if (!tool)
    return;

  if (!gtk_editable_get_text(GTK_EDITABLE(tool->exec_entry))[0])
    {
      gdis_qbox_set_report(tool->report_buffer,
                           "Set the Qbox executable path first. You can type it here, Detect it, or open Executable Paths.");
      return;
    }

  g_hash_table_replace(self->executable_paths,
                       g_strdup("qbox"),
                       g_strdup(gtk_editable_get_text(GTK_EDITABLE(tool->exec_entry))));

  warnings = g_string_new("");
  if (!gdis_qbox_write_prepared_input(self,
                                      tool,
                                      TRUE,
                                      &workdir,
                                      &input_path,
                                      &input_name,
                                      &output_path,
                                      &stderr_path,
                                      &copied_files,
                                      &localized_pseudos,
                                      warnings,
                                      &error))
    {
      gdis_qbox_set_report(tool->report_buffer,
                           error ? error->message : "Could not prepare the Qbox run.");
      gdis_gtk4_window_log(self, "Qbox run setup failed: %s\n",
                           error ? error->message : "unknown error");
      g_string_free(warnings, TRUE);
      g_clear_error(&error);
      return;
    }

  launcher = gdis_qbox_entry_text(tool->launcher_entry);
  quoted_exec = g_shell_quote(gtk_editable_get_text(GTK_EDITABLE(tool->exec_entry)));
  quoted_input = g_shell_quote(input_name);
  if (launcher[0])
    command = g_strdup_printf("%s %s %s", launcher, quoted_exec, quoted_input);
  else
    command = g_strdup_printf("%s %s", quoted_exec, quoted_input);

  success = gdis_qbox_run_shell_capture(workdir,
                                        command,
                                        output_path,
                                        stderr_path,
                                        &stdout_text,
                                        &stderr_text,
                                        &error);

  if (gtk_editable_get_text(GTK_EDITABLE(tool->save_entry))[0])
    save_path = g_build_filename(workdir,
                                 gtk_editable_get_text(GTK_EDITABLE(tool->save_entry)),
                                 NULL);

  tail = gdis_qbox_tail_text(stdout_text ? stdout_text : stderr_text, 12000u);
  last_etotal = gdis_qbox_extract_last_etotal(stdout_text);

  report = g_string_new("");
  g_string_append_printf(report,
                         "%s Qbox run.\n"
                         "Command: %s\n"
                         "Working directory: %s\n"
                         "Input file: %s\n"
                         "Output file: %s\n"
                         "Saved XML: %s\n"
                         "Staged local assets: %u\n"
                         "Localized pseudos: %u\n",
                         success ? "Completed" : "Failed",
                         command,
                         workdir,
                         input_path,
                         output_path,
                         save_path ? save_path : "(not set)",
                         copied_files,
                         localized_pseudos);
  if (stderr_text && stderr_text[0])
    g_string_append_printf(report, "stderr capture: %s\n", stderr_path);
  if (last_etotal && last_etotal[0])
    g_string_append_printf(report, "Last <etotal>: %s\n", last_etotal);
  if (warnings->len > 0)
    g_string_append_printf(report, "\nWarnings:\n%s", warnings->str);
  if (!success)
    g_string_append_printf(report,
                           "\nFailure reason: %s\n",
                           error ? error->message : "unknown error");
  if (tail && tail[0])
    g_string_append_printf(report, "\nRecent output:\n%s", tail);

  gdis_qbox_set_last_run_state(self,
                               workdir,
                               input_path,
                               output_path,
                               stderr_path,
                               save_path);

  if (success)
    {
      GError *import_error = NULL;

      if (save_path && g_file_test(save_path, G_FILE_TEST_EXISTS))
        {
          if (gdis_qbox_import_result_into_active_model(self, save_path, &import_error))
            {
              gdis_qbox_set_last_summary(self, "run completed and result imported");
              g_string_append_printf(report,
                                     "\nResult import: loaded %s back into the active model.\n",
                                     save_path);
            }
          else
            {
              gdis_qbox_set_last_summary(self, "run completed but result import failed");
              g_string_append_printf(report,
                                     "\nResult import failed: %s\n",
                                     import_error ? import_error->message : "unknown error");
              g_clear_error(&import_error);
            }
        }
      else
        {
          gdis_qbox_set_last_summary(self, "run completed");
          if (save_path)
            g_string_append_printf(report,
                                   "\nResult import skipped because %s was not created.\n",
                                   save_path);
        }
      gdis_gtk4_window_log(self, "Qbox run completed: %s\n", output_path);
    }
  else
    {
      gdis_qbox_set_last_summary(self, "run failed");
      gdis_gtk4_window_log(self, "Qbox run failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
    }

  gdis_qbox_set_report(tool->report_buffer, report->str);
  g_string_free(report, TRUE);
  g_string_free(warnings, TRUE);

  gdis_gtk4_window_refresh_qbox_tool(self);
  gdis_gtk4_window_refresh_executable_paths_tool(self);
  gdis_gtk4_window_refresh_task_manager_tool(self);
}

static void
gdis_qbox_maybe_run_startup_write(GdisGtk4Window *self)
{
  const gchar *flag;

  g_return_if_fail(self != NULL);

  if (qbox_startup_write_consumed || !self->qbox_tool)
    return;

  flag = g_getenv("GDIS_GTK4_QBOX_AUTOWRITE");
  if (!flag || !flag[0] || g_strcmp0(flag, "0") == 0 ||
      g_ascii_strcasecmp(flag, "false") == 0)
    return;

  qbox_startup_write_consumed = TRUE;
  on_qbox_write_clicked(NULL, self);
}

static void
gdis_qbox_maybe_run_startup_run(GdisGtk4Window *self)
{
  const gchar *flag;

  g_return_if_fail(self != NULL);

  if (qbox_startup_run_consumed || !self->qbox_tool)
    return;

  flag = g_getenv("GDIS_GTK4_QBOX_AUTORUN");
  if (!flag || !flag[0] || g_strcmp0(flag, "0") == 0 ||
      g_ascii_strcasecmp(flag, "false") == 0)
    return;

  qbox_startup_run_consumed = TRUE;
  on_qbox_run_clicked(NULL, self);
}

static void
gdis_gtk4_window_present_qbox_tool(GdisGtk4Window *self)
{
  GdisQboxTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *top_box;
  GtkWidget *top_scroller;
  GtkWidget *frame;
  GtkWidget *grid;
  GtkWidget *label;
  GtkWidget *row;
  GtkWidget *button;
  GtkWidget *scroller;
  GtkWidget *text_view;
  GtkWidget *paned;

  g_return_if_fail(self != NULL);

  if (self->qbox_tool && GTK_IS_WINDOW(self->qbox_tool->window))
    {
      gdis_gtk4_window_refresh_qbox_tool(self);
      gtk_window_present(GTK_WINDOW(self->qbox_tool->window));
      return;
    }

  tool = g_new0(GdisQboxTool, 1);
  tool->owner = self;

  window = GTK_WIDGET(gtk_application_window_new(self->app));
  tool->window = window;
#ifdef __APPLE__
  gtk_window_set_modal(GTK_WINDOW(window), FALSE);
#else
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
#endif
  gtk_window_set_title(GTK_WINDOW(window), "Qbox");
  gtk_window_set_default_size(GTK_WINDOW(window), 1260, 980);
  gtk_window_set_resizable(GTK_WINDOW(window), TRUE);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  top_scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(top_scroller, TRUE);
  gtk_widget_set_vexpand(top_scroller, FALSE);
  gtk_widget_set_size_request(top_scroller, -1, 360);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(top_scroller),
                                 GTK_POLICY_NEVER,
                                 GTK_POLICY_AUTOMATIC);
  gtk_box_append(GTK_BOX(root), top_scroller);

  top_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(top_scroller), top_box);

  tool->summary_label = GTK_LABEL(gtk_label_new(""));
  gtk_label_set_wrap(tool->summary_label, TRUE);
  gtk_label_set_xalign(tool->summary_label, 0.0f);
  gtk_box_append(GTK_BOX(top_box), GTK_WIDGET(tool->summary_label));

  frame = gtk_frame_new("Qbox Setup");
  gtk_box_append(GTK_BOX(top_box), frame);

  grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_widget_set_margin_start(grid, 10);
  gtk_widget_set_margin_end(grid, 10);
  gtk_widget_set_margin_top(grid, 10);
  gtk_widget_set_margin_bottom(grid, 10);
  gtk_frame_set_child(GTK_FRAME(frame), grid);

  label = gtk_label_new("Job name");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
  tool->job_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->job_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->job_entry, 1, 0, 1, 1);

  label = gtk_label_new("Working directory");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 2, 0, 1, 1);
  tool->workdir_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->workdir_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->workdir_entry, 3, 0, 3, 1);

  label = gtk_label_new("Input file");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
  tool->input_entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(grid), tool->input_entry, 1, 1, 1, 1);

  label = gtk_label_new("Output file");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 2, 1, 1, 1);
  tool->output_entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(grid), tool->output_entry, 3, 1, 1, 1);

  label = gtk_label_new("Save XML");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 4, 1, 1, 1);
  tool->save_entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(grid), tool->save_entry, 5, 1, 1, 1);

  label = gtk_label_new("Executable");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 2, 1, 1);
  tool->exec_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->exec_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->exec_entry, 1, 2, 3, 1);
  button = gtk_button_new_with_label("Detect");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_detect_exec_clicked), self);
  gtk_grid_attach(GTK_GRID(grid), button, 4, 2, 1, 1);
  button = gtk_button_new_with_label("Executable Paths");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_exec_paths_clicked), self);
  gtk_grid_attach(GTK_GRID(grid), button, 5, 2, 1, 1);

  label = gtk_label_new("Launcher prefix");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 3, 1, 1);
  tool->launcher_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->launcher_entry, TRUE);
  gtk_widget_set_tooltip_text(tool->launcher_entry, "Optional, for example: mpirun -np 4");
  gtk_grid_attach(GTK_GRID(grid), tool->launcher_entry, 1, 3, 5, 1);

  label = gtk_label_new("Pseudo dir / URL");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 4, 1, 1);
  tool->pseudo_dir_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->pseudo_dir_entry, TRUE);
  gtk_grid_attach(GTK_GRID(grid), tool->pseudo_dir_entry, 1, 4, 3, 1);
  button = gtk_button_new_with_label("Setup Local Pseudos");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_setup_local_pseudos_clicked), self);
  gtk_grid_attach(GTK_GRID(grid), button, 4, 4, 2, 1);

  label = gtk_label_new("Restart XML");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 5, 1, 1);
  tool->restart_entry = gtk_entry_new();
  gtk_widget_set_hexpand(tool->restart_entry, TRUE);
  gtk_widget_set_tooltip_text(tool->restart_entry, "Optional local restart XML or URI");
  gtk_grid_attach(GTK_GRID(grid), tool->restart_entry, 1, 5, 3, 1);
  button = gtk_button_new_with_label("Use Last XML");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_use_last_xml_clicked), self);
  gtk_grid_attach(GTK_GRID(grid), button, 4, 5, 1, 1);
  button = gtk_button_new_with_label("Continue Last");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_continue_clicked), self);
  gtk_grid_attach(GTK_GRID(grid), button, 5, 5, 1, 1);

  label = gtk_label_new("XC");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 6, 1, 1);
  tool->xc_entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(grid), tool->xc_entry, 1, 6, 1, 1);

  label = gtk_label_new("WF dyn");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 2, 6, 1, 1);
  tool->wf_dyn_entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(grid), tool->wf_dyn_entry, 3, 6, 1, 1);

  label = gtk_label_new("Atoms dyn");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 4, 6, 1, 1);
  tool->atoms_dyn_entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(grid), tool->atoms_dyn_entry, 5, 6, 1, 1);

  label = gtk_label_new("SCF tol");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 7, 1, 1);
  tool->scf_tol_entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(grid), tool->scf_tol_entry, 1, 7, 1, 1);

  label = gtk_label_new("Ecut");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 2, 7, 1, 1);
  tool->ecut_spin = gtk_spin_button_new_with_range(5.0, 500.0, 1.0);
  gtk_grid_attach(GTK_GRID(grid), tool->ecut_spin, 3, 7, 1, 1);

  label = gtk_label_new("Net charge");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 4, 7, 1, 1);
  tool->charge_spin = gtk_spin_button_new_with_range(-20.0, 20.0, 1.0);
  gtk_grid_attach(GTK_GRID(grid), tool->charge_spin, 5, 7, 1, 1);

  label = gtk_label_new("SCF steps");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 8, 1, 1);
  tool->scf_steps_spin = gtk_spin_button_new_with_range(1.0, 5000.0, 1.0);
  gtk_grid_attach(GTK_GRID(grid), tool->scf_steps_spin, 1, 8, 1, 1);

  label = gtk_label_new("Ionic steps");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 2, 8, 1, 1);
  tool->ionic_steps_spin = gtk_spin_button_new_with_range(0.0, 2000.0, 1.0);
  gtk_grid_attach(GTK_GRID(grid), tool->ionic_steps_spin, 3, 8, 1, 1);

  label = gtk_label_new("Density update");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 4, 8, 1, 1);
  tool->density_update_spin = gtk_spin_button_new_with_range(1.0, 500.0, 1.0);
  gtk_grid_attach(GTK_GRID(grid), tool->density_update_spin, 5, 8, 1, 1);

  label = gtk_label_new("Padding (A)");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 9, 1, 1);
  tool->padding_spin = gtk_spin_button_new_with_range(1.0, 40.0, 0.5);
  gtk_grid_attach(GTK_GRID(grid), tool->padding_spin, 1, 9, 1, 1);

  tool->use_cell_toggle = gtk_check_button_new_with_label("Use model periodic cell when available");
  gtk_grid_attach(GTK_GRID(grid), tool->use_cell_toggle, 2, 9, 2, 1);
  tool->atomic_density_toggle = gtk_check_button_new_with_label("Use -atomic_density startup");
  gtk_grid_attach(GTK_GRID(grid), tool->atomic_density_toggle, 4, 9, 2, 1);

  label = gtk_label_new("MP mesh");
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), label, 0, 10, 1, 1);
  for (guint i = 0; i < 3; i++)
    {
      tool->kmesh_spins[i] = gtk_spin_button_new_with_range(1.0, 24.0, 1.0);
      gtk_widget_set_hexpand(tool->kmesh_spins[i], FALSE);
      gtk_grid_attach(GTK_GRID(grid), tool->kmesh_spins[i], 1 + (gint) i, 10, 1, 1);
    }
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->kmesh_spins[0]), 4.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->kmesh_spins[1]), 4.0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(tool->kmesh_spins[2]), 1.0);
  gtk_widget_set_tooltip_text(tool->kmesh_spins[0], "Monkhorst-Pack divisions along reciprocal a*.");
  gtk_widget_set_tooltip_text(tool->kmesh_spins[1], "Monkhorst-Pack divisions along reciprocal b*.");
  gtk_widget_set_tooltip_text(tool->kmesh_spins[2], "Monkhorst-Pack divisions along reciprocal c*.");
  button = gtk_button_new_with_label("Apply MP Mesh");
  gtk_widget_set_tooltip_text(button,
                              "Replace any existing kpoint lines in the current deck with a generated Monkhorst-Pack mesh.");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_apply_mp_mesh_clicked), self);
  gtk_grid_attach(GTK_GRID(grid), button, 4, 10, 2, 1);

  tool->randomize_toggle = gtk_check_button_new_with_label("Add randomize_wf for fresh starts");
  gtk_grid_attach(GTK_GRID(grid), tool->randomize_toggle, 2, 11, 4, 1);

  frame = gtk_frame_new("Species Mapping");
  gtk_box_append(GTK_BOX(top_box), frame);
  tool->species_scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(tool->species_scroller, TRUE);
  gtk_widget_set_size_request(tool->species_scroller, -1, 220);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(tool->species_scroller),
                                 GTK_POLICY_AUTOMATIC,
                                 GTK_POLICY_AUTOMATIC);
  gtk_frame_set_child(GTK_FRAME(frame), tool->species_scroller);

  paned = gtk_paned_new(GTK_ORIENTATION_VERTICAL);
  gtk_widget_set_vexpand(paned, TRUE);
  gtk_paned_set_wide_handle(GTK_PANED(paned), TRUE);
  gtk_paned_set_resize_start_child(GTK_PANED(paned), TRUE);
  gtk_paned_set_shrink_start_child(GTK_PANED(paned), FALSE);
  gtk_paned_set_resize_end_child(GTK_PANED(paned), TRUE);
  gtk_paned_set_shrink_end_child(GTK_PANED(paned), FALSE);
  gtk_paned_set_position(GTK_PANED(paned), 520);
  gtk_box_append(GTK_BOX(root), paned);

  scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(scroller, TRUE);
  gtk_widget_set_vexpand(scroller, TRUE);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                 GTK_POLICY_AUTOMATIC,
                                 GTK_POLICY_AUTOMATIC);
  gtk_paned_set_start_child(GTK_PANED(paned), scroller);

  text_view = gtk_text_view_new();
  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_NONE);
  gtk_text_view_set_monospace(GTK_TEXT_VIEW(text_view), TRUE);
  gtk_text_view_set_left_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_right_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_top_margin(GTK_TEXT_VIEW(text_view), 12);
  gtk_text_view_set_bottom_margin(GTK_TEXT_VIEW(text_view), 12);
  tool->deck_buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(scroller, TRUE);
  gtk_widget_set_vexpand(scroller, TRUE);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                 GTK_POLICY_AUTOMATIC,
                                 GTK_POLICY_AUTOMATIC);
  gtk_paned_set_end_child(GTK_PANED(paned), scroller);

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

  row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_widget_set_halign(row, GTK_ALIGN_END);
  gtk_box_append(GTK_BOX(root), row);

  button = gtk_button_new_with_label("Guess Pseudos");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_guess_pseudos_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Regenerate");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_regenerate_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Write Input");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_write_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Run");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_run_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Import Result");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_import_result_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Results");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_results_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("HW5 Frozen Phonons");
  gtk_widget_set_tooltip_text(button,
                              "Append the original Homework 5 +/- 0.05 bohr displacement block for every atom and axis.");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_append_hw5_phonons_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("HW5 HOMO/LUMO");
  gtk_widget_set_tooltip_text(button,
                              "Append the Homework 5 continuation block for nempty/JD plus HOMO.cube and LUMO.cube export.");
  g_signal_connect(button, "clicked", G_CALLBACK(on_qbox_append_hw5_homo_lumo_clicked), self);
  gtk_box_append(GTK_BOX(row), button);

  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_destroy), window);
  gtk_box_append(GTK_BOX(row), button);

  g_signal_connect(tool->deck_buffer, "changed", G_CALLBACK(on_qbox_deck_buffer_changed), tool);
  g_signal_connect(window, "destroy", G_CALLBACK(on_qbox_tool_destroy), tool);

  self->qbox_tool = tool;
  gdis_qbox_set_report(tool->report_buffer,
                       "Generate a Qbox starter deck for the active model, then edit, write, or run it here.\n"
                       "This GTK4 tool stages local pseudo files into the job directory and preserves manual deck edits.");
  gdis_gtk4_window_refresh_qbox_tool(self);
  if (self->active_model)
    on_qbox_regenerate_clicked(NULL, self);
  gdis_qbox_maybe_run_startup_write(self);
  gdis_qbox_maybe_run_startup_run(self);
  gtk_window_present(GTK_WINDOW(window));
}

static void
on_periodic_table_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisPeriodicTableTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->periodic_table_tool = NULL;
  g_free(tool);
}

static void
gdis_gtk4_window_select_periodic_element(GdisGtk4Window *self, guint atomic_number)
{
  GdisPeriodicTableTool *tool;
  const GdisElementInfo *element;
  gboolean can_apply;
  g_autofree gchar *title = NULL;
  g_autofree gchar *info = NULL;

  g_return_if_fail(self != NULL);

  tool = self->periodic_table_tool;
  if (!tool)
    return;

  element = gdis_element_lookup_atomic_number(atomic_number);
  tool->selected_atomic_number = element ? element->atomic_number : 6u;
  element = gdis_element_lookup_atomic_number(tool->selected_atomic_number);

  title = g_strdup_printf("%s (%s)", element->name, element->symbol);
  info = g_strdup_printf(
    "Atomic number: %u\n"
    "Family: %s\n"
    "Covalent radius: %.2f A\n"
    "Van der Waals radius: %.2f A\n"
    "Viewer radius: %.2f A\n"
    "Table position: row %u, column %u\n\n"
    "Use 'Copy to Editor' to populate the GTK4 editing tool, or 'Apply to Selected Atom' to retag the currently selected atom.",
    element->atomic_number,
    gdis_element_family_label(element->family),
    element->covalent_radius,
    element->vdw_radius,
    element->draw_radius,
    element->row,
    element->column);
  gtk_label_set_text(tool->title_label, title);
  gtk_label_set_text(tool->info_label, info);

  can_apply = (self->active_model != NULL &&
               self->selected_atom_index != INVALID_ATOM_INDEX &&
               self->selected_atom_index < self->active_model->atoms->len);
  gtk_widget_set_sensitive(tool->apply_button, can_apply);
  gtk_widget_set_sensitive(tool->copy_button, TRUE);
}

static void
on_periodic_element_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  guint atomic_number;

  self = user_data;
  atomic_number = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(button), "atomic-number"));
  if (!self || atomic_number == 0u)
    return;

  gdis_gtk4_window_select_periodic_element(self, atomic_number);
}

static void
on_periodic_copy_clicked(GtkButton *button, gpointer user_data)
{
  GdisPeriodicTableTool *tool;
  GdisGtk4Window *self;
  const GdisElementInfo *element;

  (void) button;

  tool = user_data;
  self = tool ? tool->owner : NULL;
  if (!self)
    return;

  element = gdis_element_lookup_atomic_number(tool->selected_atomic_number ? tool->selected_atomic_number : 6u);
  gdis_gtk4_window_present_edit_tool(self);
  if (!self->edit_tool)
    return;

  gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->element_entry), element->symbol);
  if (gtk_editable_get_text(GTK_EDITABLE(self->edit_tool->label_entry))[0] == '\0')
    gtk_editable_set_text(GTK_EDITABLE(self->edit_tool->label_entry), element->symbol);
  gdis_gtk4_window_log(self, "Copied %s into the GTK4 editing tool.\n", element->symbol);
  gdis_gtk4_window_refresh_periodic_table_tool(self);
}

static void
on_periodic_apply_clicked(GtkButton *button, gpointer user_data)
{
  GdisPeriodicTableTool *tool;
  GdisGtk4Window *self;
  const GdisElementInfo *element;
  const GdisAtom *atom;
  GError *error;

  (void) button;

  tool = user_data;
  self = tool ? tool->owner : NULL;
  if (!self || !self->active_model)
    return;

  if (self->selected_atom_index == INVALID_ATOM_INDEX ||
      self->selected_atom_index >= self->active_model->atoms->len)
    return;

  element = gdis_element_lookup_atomic_number(tool->selected_atomic_number ? tool->selected_atomic_number : 6u);
  atom = g_ptr_array_index(self->active_model->atoms, self->selected_atom_index);
  if (!gdis_gtk4_window_push_undo_snapshot(self, "Apply element from periodic table"))
    return;

  error = NULL;
  if (!gdis_model_update_atom(self->active_model,
                              self->selected_atom_index,
                              atom->label && atom->label[0] ? atom->label : element->symbol,
                              element->symbol,
                              atom->ff_type,
                              atom->region,
                              atom->position[0],
                              atom->position[1],
                              atom->position[2],
                              &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
      gdis_gtk4_window_log(self, "Periodic table apply failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  gdis_gtk4_window_refresh_after_model_edit(self, TRUE);
  gdis_gtk4_window_log(self,
                       "Applied element %s to atom %u.\n",
                       element->symbol,
                       atom->serial);
}

static void
gdis_gtk4_window_refresh_periodic_table_tool(GdisGtk4Window *self)
{
  GdisPeriodicTableTool *tool;

  g_return_if_fail(self != NULL);

  tool = self->periodic_table_tool;
  if (!tool)
    return;

  if (tool->selected_atomic_number == 0u)
    tool->selected_atomic_number = 6u;
  gdis_gtk4_window_select_periodic_element(self, tool->selected_atomic_number);
}

static void
gdis_gtk4_window_present_periodic_table_tool(GdisGtk4Window *self)
{
  GdisPeriodicTableTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *content;
  GtkWidget *scroller;
  GtkWidget *grid;
  GtkWidget *right;
  GtkWidget *button;
  GtkWidget *frame;

  g_return_if_fail(self != NULL);

  if (self->periodic_table_tool && GTK_IS_WINDOW(self->periodic_table_tool->window))
    {
      gdis_gtk4_window_refresh_periodic_table_tool(self);
      gtk_window_present(GTK_WINDOW(self->periodic_table_tool->window));
      return;
    }

  tool = g_new0(GdisPeriodicTableTool, 1);
  tool->owner = self;
  tool->selected_atomic_number = 6u;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Periodic Table");
  gtk_window_set_default_size(GTK_WINDOW(window), 980, 520);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  content = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 12);
  gtk_box_append(GTK_BOX(root), content);

  scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(scroller, TRUE);
  gtk_widget_set_vexpand(scroller, TRUE);
  gtk_box_append(GTK_BOX(content), scroller);

  grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(grid), 6);
  gtk_grid_set_row_spacing(GTK_GRID(grid), 6);
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), grid);

  for (guint atomic_number = 1; atomic_number <= gdis_element_count(); atomic_number++)
    {
      const GdisElementInfo *element;

      element = gdis_element_lookup_atomic_number(atomic_number);
      button = gtk_button_new_with_label(element->symbol);
      gtk_widget_set_size_request(button, 46, 34);
      gtk_widget_set_tooltip_text(button, element->name);
      g_object_set_data(G_OBJECT(button), "atomic-number", GUINT_TO_POINTER(atomic_number));
      g_signal_connect(button, "clicked", G_CALLBACK(on_periodic_element_clicked), self);
      gtk_grid_attach(GTK_GRID(grid),
                      button,
                      (gint) element->column - 1,
                      (gint) element->row - 1,
                      1,
                      1);
    }

  right = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
  gtk_widget_set_size_request(right, 280, -1);
  gtk_box_append(GTK_BOX(content), right);

  tool->title_label = GTK_LABEL(gtk_label_new(""));
  gtk_label_set_xalign(tool->title_label, 0.0f);
  gtk_box_append(GTK_BOX(right), GTK_WIDGET(tool->title_label));

  frame = gtk_frame_new("Element details");
  gtk_box_append(GTK_BOX(right), frame);
  tool->info_label = GTK_LABEL(gtk_label_new(""));
  gtk_label_set_xalign(tool->info_label, 0.0f);
  gtk_label_set_wrap(tool->info_label, TRUE);
  gtk_widget_set_margin_start(GTK_WIDGET(tool->info_label), 10);
  gtk_widget_set_margin_end(GTK_WIDGET(tool->info_label), 10);
  gtk_widget_set_margin_top(GTK_WIDGET(tool->info_label), 10);
  gtk_widget_set_margin_bottom(GTK_WIDGET(tool->info_label), 10);
  gtk_frame_set_child(GTK_FRAME(frame), GTK_WIDGET(tool->info_label));

  button = gtk_button_new_with_label("Copy to Editor");
  tool->copy_button = button;
  g_signal_connect(button, "clicked", G_CALLBACK(on_periodic_copy_clicked), tool);
  gtk_box_append(GTK_BOX(right), button);

  button = gtk_button_new_with_label("Apply to Selected Atom");
  tool->apply_button = button;
  g_signal_connect(button, "clicked", G_CALLBACK(on_periodic_apply_clicked), tool);
  gtk_box_append(GTK_BOX(right), button);

  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_destroy), window);
  gtk_box_append(GTK_BOX(right), button);

  g_signal_connect(window, "destroy", G_CALLBACK(on_periodic_table_tool_destroy), tool);

  self->periodic_table_tool = tool;
  gdis_gtk4_window_refresh_periodic_table_tool(self);
  gtk_window_present(GTK_WINDOW(window));
}

struct _GdisExecutableSpec
{
  const char *id;
  const char *label;
  const char *probe_a;
  const char *probe_b;
};

static const GdisExecutableSpec gdis_executable_specs[] = {
  {"gulp", "GULP", "gulp", NULL},
  {"gamess", "GAMESS", "rungms", "gamess"},
  {"monty", "Monty", "monty", NULL},
  {"qbox", "Qbox", "qbox", "qb"},
  {"siesta", "SIESTA", "siesta", NULL},
  {"vasp", "VASP", "vasp_std", "vasp"},
  {"uspex", "USPEX", "USPEX", "uspex"},
  {"povray", "POV-Ray", "povray", NULL},
  {"convert", "ImageMagick", "magick", "convert"}
};

static const GdisExecutableSpec *
gdis_executable_spec_for_id(const char *id)
{
  for (guint i = 0; i < G_N_ELEMENTS(gdis_executable_specs); i++)
    {
      if (g_strcmp0(gdis_executable_specs[i].id, id) == 0)
        return &gdis_executable_specs[i];
    }
  return NULL;
}

static gchar *
gdis_executable_detect_path(const GdisExecutableSpec *spec)
{
  gchar *path;
  const char *probes[3];

  g_return_val_if_fail(spec != NULL, NULL);

  path = g_find_program_in_path(spec->probe_a);
  if (!path && spec->probe_b)
    path = g_find_program_in_path(spec->probe_b);
  if (path)
    return path;

  probes[0] = spec->probe_a;
  probes[1] = spec->probe_b;
  probes[2] = NULL;
  for (guint i = 0; probes[i] != NULL && !path; i++)
    {
      const char *probe;

      probe = probes[i];
      if (!probe || probe[0] == '\0')
        continue;

      {
        g_autofree gchar *homebrew_candidate = NULL;
        g_autofree gchar *local_candidate = NULL;

        homebrew_candidate = g_build_filename("/opt/homebrew/bin", probe, NULL);
        if (g_file_test(homebrew_candidate, G_FILE_TEST_IS_EXECUTABLE))
          path = g_strdup(homebrew_candidate);

        if (!path)
          {
            local_candidate = g_build_filename("/usr/local/bin", probe, NULL);
            if (g_file_test(local_candidate, G_FILE_TEST_IS_EXECUTABLE))
              path = g_strdup(local_candidate);
          }
      }
    }
  return path;
}

static void
on_task_manager_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisTaskManagerTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->task_manager_tool = NULL;
  g_free(tool);
}

static void
gdis_gtk4_window_refresh_task_manager_tool(GdisGtk4Window *self)
{
  GdisTaskManagerTool *tool;
  GString *text;
  GdisIsoSurface *surface;
  GPtrArray *records;
  GPtrArray *undo_stack;

  g_return_if_fail(self != NULL);

  tool = self->task_manager_tool;
  if (!tool || !tool->buffer)
    return;

  text = g_string_new("");
  g_string_append_printf(text, "Loaded models: %u\n", self->models->len);
  g_string_append_printf(text, "Active model: %s\n\n",
                         self->active_model ? self->active_model->basename : "none");

  for (guint i = 0; i < self->models->len; i++)
    {
      GdisModel *model;

      model = g_ptr_array_index(self->models, i);
      g_string_append_printf(text,
                             "%c [%u] %s  |  %u atoms  |  %u bonds\n",
                             model == self->active_model ? '*' : '-',
                             i + 1,
                             model->basename,
                             model->atom_count,
                             model->bond_count);
    }

  g_string_append(text, "\nOpen GTK4 tools:\n");
  g_string_append_printf(text, "  Display properties: %s\n", self->display_tool ? "open" : "closed");
  g_string_append_printf(text, "  Editing: %s\n", self->edit_tool ? "open" : "closed");
  g_string_append_printf(text, "  Measurements: %s\n", self->measure_tool ? "open" : "closed");
  g_string_append_printf(text, "  Diffraction: %s\n", self->diffraction_tool ? "open" : "closed");
  g_string_append_printf(text, "  Surface builder: %s\n", self->surface_tool ? "open" : "closed");
  g_string_append_printf(text, "  Iso-surfaces: %s\n", self->isosurface_tool ? "open" : "closed");
  g_string_append_printf(text, "  Animation: %s\n", self->animation_tool ? "open" : "closed");
  g_string_append_printf(text, "  Qbox: %s\n", self->qbox_tool ? "open" : "closed");
  g_string_append_printf(text, "  Periodic table: %s\n", self->periodic_table_tool ? "open" : "closed");
  g_string_append_printf(text, "  Executable paths: %s\n", self->exec_paths_tool ? "open" : "closed");

  if (self->active_model)
    {
      records = gdis_gtk4_window_get_measurement_records(self, self->active_model, FALSE);
      undo_stack = gdis_gtk4_window_get_undo_stack(self, self->active_model, FALSE);
      surface = gdis_gtk4_window_get_isosurface(self, self->active_model, FALSE);
      g_string_append(text, "\nActive model state:\n");
      g_string_append_printf(text, "  Selection size: %u\n", self->selected_atoms ? self->selected_atoms->len : 0u);
      g_string_append_printf(text, "  Pick history size: %u\n", self->picked_atoms ? self->picked_atoms->len : 0u);
      g_string_append_printf(text, "  Undo snapshots: %u\n", undo_stack ? undo_stack->len : 0u);
      g_string_append_printf(text, "  Saved measurements: %u\n", records ? records->len : 0u);
      g_string_append_printf(text, "  Iso-surface triangles: %u\n",
                             surface && surface->triangles ? surface->triangles->len : 0u);
    }

  g_string_append_printf(text,
                         "\nConfigured executable paths: %u\n",
                         self->executable_paths ? g_hash_table_size(self->executable_paths) : 0u);
  if (self->qbox_last_summary && self->qbox_last_summary[0])
    g_string_append_printf(text, "Last Qbox activity: %s\n", self->qbox_last_summary);
  if (self->qbox_last_workdir)
    g_string_append_printf(text, "Last Qbox workdir: %s\n", self->qbox_last_workdir);
  if (self->qbox_last_save_path)
    g_string_append_printf(text, "Last Qbox XML: %s\n", self->qbox_last_save_path);
  gtk_text_buffer_set_text(tool->buffer, text->str, -1);
  g_string_free(text, TRUE);
}

static void
on_task_manager_action_clicked(GtkButton *button, gpointer user_data)
{
  GdisGtk4Window *self;
  const char *action;

  self = user_data;
  action = g_object_get_data(G_OBJECT(button), "task-action");
  if (!self || !action)
    return;

  if (g_strcmp0(action, "display") == 0)
    gdis_gtk4_window_present_display_tool(self);
  else if (g_strcmp0(action, "edit") == 0)
    gdis_gtk4_window_present_edit_tool(self);
  else if (g_strcmp0(action, "measure") == 0)
    gdis_gtk4_window_present_measure_tool(self);
  else if (g_strcmp0(action, "isosurface") == 0)
    gdis_gtk4_window_present_isosurface_tool(self);
  else if (g_strcmp0(action, "qbox") == 0)
    gdis_gtk4_window_present_qbox_tool(self);
  else if (g_strcmp0(action, "exec-paths") == 0)
    gdis_gtk4_window_present_executable_paths_tool(self, NULL);
  else if (g_strcmp0(action, "close-model") == 0)
    {
      if (self->active_model)
        gdis_gtk4_window_remove_model(self, self->active_model);
    }

  gdis_gtk4_window_refresh_task_manager_tool(self);
}

static void
gdis_gtk4_window_present_task_manager_tool(GdisGtk4Window *self)
{
  GdisTaskManagerTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *scroller;
  GtkWidget *text_view;
  GtkWidget *actions;
  GtkWidget *button;

  g_return_if_fail(self != NULL);

  if (self->task_manager_tool && GTK_IS_WINDOW(self->task_manager_tool->window))
    {
      gdis_gtk4_window_refresh_task_manager_tool(self);
      gtk_window_present(GTK_WINDOW(self->task_manager_tool->window));
      return;
    }

  tool = g_new0(GdisTaskManagerTool, 1);
  tool->owner = self;

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Task Manager");
  gtk_window_set_default_size(GTK_WINDOW(window), 720, 420);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(scroller, TRUE);
  gtk_widget_set_vexpand(scroller, TRUE);
  gtk_box_append(GTK_BOX(root), scroller);

  text_view = gtk_text_view_new();
  gtk_text_view_set_editable(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_monospace(GTK_TEXT_VIEW(text_view), TRUE);
  tool->buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  actions = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_box_append(GTK_BOX(root), actions);

  for (guint i = 0; i < 7; i++)
    {
      const char *label;
      const char *action;

      label = "";
      action = "";
      switch (i)
        {
        case 0:
          label = "Display";
          action = "display";
          break;
        case 1:
          label = "Edit";
          action = "edit";
          break;
        case 2:
          label = "Measure";
          action = "measure";
          break;
        case 3:
          label = "Iso-surfaces";
          action = "isosurface";
          break;
        case 4:
          label = "Exec Paths";
          action = "exec-paths";
          break;
        case 5:
          label = "Qbox";
          action = "qbox";
          break;
        case 6:
          label = "Close Model";
          action = "close-model";
          break;
        default:
          break;
        }

      button = gtk_button_new_with_label(label);
      g_object_set_data_full(G_OBJECT(button), "task-action", g_strdup(action), g_free);
      g_signal_connect(button, "clicked", G_CALLBACK(on_task_manager_action_clicked), self);
      gtk_box_append(GTK_BOX(actions), button);
    }

  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_destroy), window);
  gtk_box_append(GTK_BOX(actions), button);

  g_signal_connect(window, "destroy", G_CALLBACK(on_task_manager_tool_destroy), tool);

  self->task_manager_tool = tool;
  gdis_gtk4_window_refresh_task_manager_tool(self);
  gtk_window_present(GTK_WINDOW(window));
}

static void
on_exec_paths_tool_destroy(GtkWindow *window, gpointer user_data)
{
  GdisExecutablePathsTool *tool;

  (void) window;

  tool = user_data;
  if (!tool)
    return;

  if (tool->owner)
    tool->owner->exec_paths_tool = NULL;
  g_clear_pointer(&tool->entries, g_free);
  g_clear_pointer(&tool->focus_backend, g_free);
  g_free(tool);
}

static void
gdis_gtk4_window_refresh_executable_paths_tool(GdisGtk4Window *self)
{
  GdisExecutablePathsTool *tool;
  GString *preview;
  const GdisExecutableSpec *focus_spec;

  g_return_if_fail(self != NULL);

  tool = self->exec_paths_tool;
  if (!tool)
    return;

  for (guint i = 0; i < G_N_ELEMENTS(gdis_executable_specs); i++)
    {
      const char *stored_path;
      g_autofree gchar *detected_path = NULL;

      stored_path = self->executable_paths ?
        g_hash_table_lookup(self->executable_paths, gdis_executable_specs[i].id) :
        NULL;
      if (!stored_path || !stored_path[0])
        detected_path = gdis_executable_detect_path(&gdis_executable_specs[i]);
      gtk_editable_set_text(GTK_EDITABLE(tool->entries[i]),
                            stored_path && stored_path[0] ? stored_path :
                            detected_path ? detected_path : "");
    }

  preview = g_string_new("");
  focus_spec = gdis_executable_spec_for_id(tool->focus_backend);
  if (focus_spec)
    {
      const char *path_text;

      path_text = gtk_editable_get_text(GTK_EDITABLE(tool->entries[focus_spec - gdis_executable_specs]));
      if (g_strcmp0(focus_spec->id, "qbox") == 0)
        g_string_append_printf(preview,
                               "Focused backend: %s\n"
                               "Executable: %s\n"
                               "Active model: %s\n\n"
                               "Qbox now has a native GTK4 deck editor and run window. This executable-path page still matters because the Qbox tool reads the session path stored here when it launches.\n",
                               focus_spec->label,
                               (path_text && path_text[0]) ? path_text : "(not set)",
                               self->active_model ? self->active_model->basename : "none");
      else
        g_string_append_printf(preview,
                               "Focused backend: %s\n"
                               "Executable: %s\n"
                               "Active model: %s\n\n"
                               "The GTK4 rebuild now restores executable-path management here. Native diffraction and surface workflows already run directly in GTK4; backend-specific input-deck editors for %s are still a later port, but this window is now where the session executable path lives.\n",
                               focus_spec->label,
                               (path_text && path_text[0]) ? path_text : "(not set)",
                               self->active_model ? self->active_model->basename : "none",
                               focus_spec->label);
    }
  else
    {
      g_string_append(preview,
                      "Configure session executable paths for restored legacy backends here.\n"
                      "These paths are kept in the current GTK4 session and reused by the computation menus.");
    }
  gtk_text_buffer_set_text(tool->preview_buffer, preview->str, -1);
  g_string_free(preview, TRUE);
}

static void
on_exec_detect_clicked(GtkButton *button, gpointer user_data)
{
  GdisExecutablePathsTool *tool;
  guint index;
  g_autofree gchar *path = NULL;

  (void) button;

  tool = user_data;
  if (!tool)
    return;

  index = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(button), "spec-index"));
  if (index >= G_N_ELEMENTS(gdis_executable_specs))
    return;

  path = gdis_executable_detect_path(&gdis_executable_specs[index]);
  gtk_editable_set_text(GTK_EDITABLE(tool->entries[index]), path ? path : "");
  gdis_gtk4_window_refresh_executable_paths_tool(tool->owner);
}

static void
on_exec_apply_clicked(GtkButton *button, gpointer user_data)
{
  GdisExecutablePathsTool *tool;
  GdisGtk4Window *self;

  (void) button;

  tool = user_data;
  self = tool ? tool->owner : NULL;
  if (!self || !self->executable_paths)
    return;

  for (guint i = 0; i < G_N_ELEMENTS(gdis_executable_specs); i++)
    {
      const char *text;

      text = gtk_editable_get_text(GTK_EDITABLE(tool->entries[i]));
      if (text && text[0])
        g_hash_table_replace(self->executable_paths,
                             g_strdup(gdis_executable_specs[i].id),
                             g_strdup(text));
      else
        g_hash_table_remove(self->executable_paths, gdis_executable_specs[i].id);
    }

  gdis_gtk4_window_log(self, "Saved session executable-path settings.\n");
  gdis_gtk4_window_refresh_executable_paths_tool(self);
  gdis_gtk4_window_refresh_task_manager_tool(self);
}

static void
on_exec_detect_all_clicked(GtkButton *button, gpointer user_data)
{
  GdisExecutablePathsTool *tool;

  (void) button;

  tool = user_data;
  if (!tool)
    return;

  for (guint i = 0; i < G_N_ELEMENTS(gdis_executable_specs); i++)
    {
      g_autofree gchar *path = gdis_executable_detect_path(&gdis_executable_specs[i]);
      if (path)
        gtk_editable_set_text(GTK_EDITABLE(tool->entries[i]), path);
    }
  gdis_gtk4_window_refresh_executable_paths_tool(tool->owner);
}

static void
gdis_gtk4_window_present_executable_paths_tool(GdisGtk4Window *self,
                                               const char *focus_backend)
{
  GdisExecutablePathsTool *tool;
  GtkWidget *window;
  GtkWidget *root;
  GtkWidget *grid;
  GtkWidget *label;
  GtkWidget *button;
  GtkWidget *scroller;
  GtkWidget *text_view;
  GtkWidget *actions;

  g_return_if_fail(self != NULL);

  if (self->exec_paths_tool && GTK_IS_WINDOW(self->exec_paths_tool->window))
    {
      g_free(self->exec_paths_tool->focus_backend);
      self->exec_paths_tool->focus_backend = g_strdup(focus_backend);
      gdis_gtk4_window_refresh_executable_paths_tool(self);
      gtk_window_present(GTK_WINDOW(self->exec_paths_tool->window));
      return;
    }

  tool = g_new0(GdisExecutablePathsTool, 1);
  tool->owner = self;
  tool->entry_count = G_N_ELEMENTS(gdis_executable_specs);
  tool->entries = g_new0(GtkWidget *, tool->entry_count);
  tool->focus_backend = g_strdup(focus_backend);

  window = gtk_window_new();
  tool->window = window;
  gtk_window_set_application(GTK_WINDOW(window), self->app);
  gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(self->window));
  gtk_window_set_title(GTK_WINDOW(window), "Executable Paths");
  gtk_window_set_default_size(GTK_WINDOW(window), 780, 520);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
  gtk_widget_set_margin_start(root, 12);
  gtk_widget_set_margin_end(root, 12);
  gtk_widget_set_margin_top(root, 12);
  gtk_widget_set_margin_bottom(root, 12);
  gtk_window_set_child(GTK_WINDOW(window), root);

  grid = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
  gtk_box_append(GTK_BOX(root), grid);

  for (guint i = 0; i < G_N_ELEMENTS(gdis_executable_specs); i++)
    {
      label = gtk_label_new(gdis_executable_specs[i].label);
      gtk_widget_set_halign(label, GTK_ALIGN_START);
      gtk_grid_attach(GTK_GRID(grid), label, 0, (gint) i, 1, 1);

      tool->entries[i] = gtk_entry_new();
      gtk_widget_set_hexpand(tool->entries[i], TRUE);
      gtk_grid_attach(GTK_GRID(grid), tool->entries[i], 1, (gint) i, 1, 1);

      button = gtk_button_new_with_label("Detect");
      g_object_set_data(G_OBJECT(button), "spec-index", GUINT_TO_POINTER(i));
      g_signal_connect(button, "clicked", G_CALLBACK(on_exec_detect_clicked), tool);
      gtk_grid_attach(GTK_GRID(grid), button, 2, (gint) i, 1, 1);
    }

  scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(scroller, TRUE);
  gtk_widget_set_vexpand(scroller, TRUE);
  gtk_box_append(GTK_BOX(root), scroller);

  text_view = gtk_text_view_new();
  gtk_text_view_set_editable(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(text_view), FALSE);
  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_WORD_CHAR);
  tool->preview_buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), text_view);

  actions = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
  gtk_box_append(GTK_BOX(root), actions);

  button = gtk_button_new_with_label("Detect All");
  g_signal_connect(button, "clicked", G_CALLBACK(on_exec_detect_all_clicked), tool);
  gtk_box_append(GTK_BOX(actions), button);

  button = gtk_button_new_with_label("Apply");
  g_signal_connect(button, "clicked", G_CALLBACK(on_exec_apply_clicked), tool);
  gtk_box_append(GTK_BOX(actions), button);

  button = gtk_button_new_with_label("Close");
  g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_window_destroy), window);
  gtk_box_append(GTK_BOX(actions), button);

  g_signal_connect(window, "destroy", G_CALLBACK(on_exec_paths_tool_destroy), tool);

  self->exec_paths_tool = tool;
  gdis_gtk4_window_refresh_executable_paths_tool(self);
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
                             "No active model\nOpen a file to load a structure.");
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
                             "\nSelected atom:\n  Label: %s\n  Element: %s\n  FF Type: %s\n  Region: %d\n  Serial: %u\n  Position: %.4f %.4f %.4f\n",
                             selected_atom->label,
                             selected_atom->element,
                             selected_atom->ff_type ? selected_atom->ff_type : "",
                             selected_atom->region,
                             selected_atom->serial,
                             selected_atom->position[0],
                             selected_atom->position[1],
                             selected_atom->position[2]);
    }
  if (self->content_buffer)
    gdis_text_buffer_set(self->content_buffer, "%s", content->str);
  g_string_free(content, TRUE);

  editing = g_string_new(" Sel  Serial  Elem  FFType     Region  Label                X           Y           Z\n");
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
                             "  %s   %5u   %-3s   %-10s %6d  %-16s %10.4f %10.4f %10.4f\n",
                             marker,
                             atom->serial,
                             atom->element,
                             atom->ff_type ? atom->ff_type : "",
                             atom->region,
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

static gint
gdis_gtk4_window_compare_models(const GdisModel *left, const GdisModel *right)
{
  const char *left_name;
  const char *right_name;
  g_autofree gchar *left_key = NULL;
  g_autofree gchar *right_key = NULL;
  gint result;

  if (left == right)
    return 0;
  if (!left)
    return -1;
  if (!right)
    return 1;

  left_name = (left->basename && left->basename[0]) ? left->basename : left->path;
  right_name = (right->basename && right->basename[0]) ? right->basename : right->path;
  if (!left_name)
    left_name = "";
  if (!right_name)
    right_name = "";

  left_key = g_utf8_collate_key_for_filename(left_name, -1);
  right_key = g_utf8_collate_key_for_filename(right_name, -1);
  result = g_strcmp0(left_key, right_key);
  if (result != 0)
    return result;

  return g_strcmp0(left->path, right->path);
}

static gint
gdis_gtk4_window_get_model_index(GdisGtk4Window *self, GdisModel *model)
{
  g_return_val_if_fail(self != NULL, -1);

  if (!model)
    return -1;

  for (guint i = 0; i < self->models->len; i++)
    {
      if (g_ptr_array_index(self->models, i) == model)
        return (gint) i;
    }

  return -1;
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
  gdis_gtk4_window_sync_display_tool(self);

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
  gdis_gtk4_window_refresh_surface_tool(self);
  gdis_gtk4_window_refresh_isosurface_tool(self);
  gdis_gtk4_window_refresh_animation_tool(self);
  gdis_gtk4_window_refresh_qbox_tool(self);
  gdis_gtk4_window_refresh_periodic_table_tool(self);
  gdis_gtk4_window_refresh_task_manager_tool(self);
  gdis_gtk4_window_update_undo_action(self);
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
  GtkWidget *previous_button;
  gchar *label;
  guint insert_index;

  g_return_val_if_fail(self != NULL, NULL);
  g_return_val_if_fail(model != NULL, NULL);

  insert_index = self->models->len;
  for (guint i = 0; i < self->models->len; i++)
    {
      GdisModel *existing_model;

      existing_model = g_ptr_array_index(self->models, i);
      if (gdis_gtk4_window_compare_models(model, existing_model) < 0)
        {
          insert_index = i;
          break;
        }
    }
  g_ptr_array_insert(self->models, insert_index, model);

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

  previous_button = NULL;
  if (insert_index > 0)
    {
      GdisModel *previous_model;

      previous_model = g_ptr_array_index(self->models, insert_index - 1);
      previous_button = gdis_gtk4_window_find_model_button(self, previous_model->path);
    }
  gtk_box_insert_child_after(GTK_BOX(self->model_list), button, previous_button);
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

  if (self->measurement_records)
    g_hash_table_remove(self->measurement_records, model);
  if (self->undo_stacks)
    g_hash_table_remove(self->undo_stacks, model);
  if (self->iso_surfaces)
    g_hash_table_remove(self->iso_surfaces, model);

  g_ptr_array_remove_index(self->models, index);
  gdis_gtk4_window_set_active_model(self, replacement);
  gdis_gtk4_window_refresh_model_buttons(self);
  gdis_gtk4_window_refresh_edit_tool(self);
  gdis_gtk4_window_refresh_measure_tool(self);
  gdis_gtk4_window_refresh_animation_tool(self);
  gdis_gtk4_window_refresh_qbox_tool(self);
  gdis_gtk4_window_refresh_task_manager_tool(self);
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
  const GdisElementInfo *info;

  info = gdis_element_lookup(element);
  return info ? info->draw_radius : 0.82;
}

static void
gdis_set_element_color(cairo_t *cr, const char *element, gdouble alpha)
{
  const GdisElementInfo *info;

  info = gdis_element_lookup(element);
  cairo_set_source_rgba(cr,
                        info ? info->color_rgb[0] : 0.58,
                        info ? info->color_rgb[1] : 0.74,
                        info ? info->color_rgb[2] : 0.86,
                        alpha);
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

static gint
gdis_projected_triangle_compare(gconstpointer left, gconstpointer right)
{
  const GdisProjectedTriangle *a;
  const GdisProjectedTriangle *b;

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
  else if (self->active_model->periodic)
    {
      guint periodic_dims;
      gdouble average_shift[3];
      gdouble periodic_center[3];

      periodic_dims = self->active_model->periodicity;
      if (periodic_dims == 0u)
        periodic_dims = 3u;
      periodic_dims = MIN(periodic_dims, 3u);

      average_shift[0] = ((gdouble) (image_positive[0] - 1 - image_negative[0])) * 0.5;
      average_shift[1] = ((gdouble) (image_positive[1] - 1 - image_negative[1])) * 0.5;
      average_shift[2] = ((gdouble) (image_positive[2] - 1 - image_negative[2])) * 0.5;

      periodic_center[0] = average_shift[0] * a_vec[0] +
                           average_shift[1] * b_vec[0] +
                           average_shift[2] * c_vec[0];
      periodic_center[1] = average_shift[0] * a_vec[1] +
                           average_shift[1] * b_vec[1] +
                           average_shift[2] * c_vec[1];
      periodic_center[2] = average_shift[0] * a_vec[2] +
                           average_shift[1] * b_vec[2] +
                           average_shift[2] * c_vec[2];

      if (periodic_dims > 0u)
        {
          periodic_center[0] += 0.5 * a_vec[0];
          periodic_center[1] += 0.5 * a_vec[1];
          periodic_center[2] += 0.5 * a_vec[2];
        }
      if (periodic_dims > 1u)
        {
          periodic_center[0] += 0.5 * b_vec[0];
          periodic_center[1] += 0.5 * b_vec[1];
          periodic_center[2] += 0.5 * b_vec[2];
        }
      if (periodic_dims > 2u)
        {
          periodic_center[0] += 0.5 * c_vec[0];
          periodic_center[1] += 0.5 * c_vec[1];
          periodic_center[2] += 0.5 * c_vec[2];
        }

      center_out[0] = periodic_center[0];
      center_out[1] = periodic_center[1];
      if (periodic_dims >= 3u)
        center_out[2] = periodic_center[2];
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

  line_one = "Use File > Open to load one or more structure files.";
  line_two = "Left-drag atom or Alt-drag anywhere to rotate, left-drag empty space to box select, right-drag anywhere to rotate, and scroll to zoom.";
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
  bottom_left = g_strdup_printf("Drag atom or Alt/Right-drag: rotate    Drag empty space: box select    Shift+Click/Box: add-toggle    Scroll: zoom    Zoom %.2fx",
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

  if (self->drag_mode == GDIS_DRAG_MODE_BOX_SELECT)
    {
      gdouble rect_x;
      gdouble rect_y;
      gdouble rect_w;
      gdouble rect_h;

      rect_x = MIN(self->press_x, self->drag_current_x);
      rect_y = MIN(self->press_y, self->drag_current_y);
      rect_w = fabs(self->drag_current_x - self->press_x);
      rect_h = fabs(self->drag_current_y - self->press_y);

      cairo_set_source_rgba(cr, 0.42, 0.74, 0.95, 0.16);
      cairo_rectangle(cr, rect_x, rect_y, rect_w, rect_h);
      cairo_fill_preserve(cr);
      cairo_set_source_rgba(cr, 0.42, 0.74, 0.95, 0.95);
      cairo_set_line_width(cr, 1.4);
      cairo_stroke(cr);
    }
  cairo_restore(cr);

  g_free(top_left);
  g_free(bottom_left);
}

static void
gdis_draw_isosurface(cairo_t *cr,
                     GdisGtk4Window *self,
                     const GdisIsoSurface *surface,
                     const gdouble center[3],
                     gdouble scale,
                     int width,
                     int height)
{
  GdisProjectedTriangle *triangles;

  g_return_if_fail(cr != NULL);
  g_return_if_fail(self != NULL);

  if (!surface || !surface->triangles || surface->triangles->len == 0)
    return;

  triangles = g_new(GdisProjectedTriangle, surface->triangles->len);
  for (guint i = 0; i < surface->triangles->len; i++)
    {
      const GdisIsoTriangle *triangle;
      gdouble average_normal[3];
      gdouble rotated_normal[3];
      gdouble accumulated_depth;

      triangle = &g_array_index(surface->triangles, GdisIsoTriangle, i);
      accumulated_depth = 0.0;
      for (guint v = 0; v < 3; v++)
        {
          gdouble depth;

          gdis_project_point(triangle->vertices[v].position,
                             center,
                             self->rotation_x,
                             self->rotation_y,
                             scale,
                             width,
                             height,
                             &triangles[i].screen[v][0],
                             &triangles[i].screen[v][1],
                             &depth);
          accumulated_depth += depth;
        }

      average_normal[0] = triangle->vertices[0].normal[0] +
                          triangle->vertices[1].normal[0] +
                          triangle->vertices[2].normal[0];
      average_normal[1] = triangle->vertices[0].normal[1] +
                          triangle->vertices[1].normal[1] +
                          triangle->vertices[2].normal[1];
      average_normal[2] = triangle->vertices[0].normal[2] +
                          triangle->vertices[1].normal[2] +
                          triangle->vertices[2].normal[2];
      gdis_rotate_point(average_normal,
                        self->rotation_x,
                        self->rotation_y,
                        rotated_normal);
      triangles[i].depth = accumulated_depth / 3.0;
      triangles[i].shade = CLAMP(0.20 + 0.80 * fabs(rotated_normal[2]), 0.15, 1.0);
    }

  qsort(triangles,
        surface->triangles->len,
        sizeof(*triangles),
        gdis_projected_triangle_compare);

  cairo_save(cr);
  cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
  for (guint i = 0; i < surface->triangles->len; i++)
    {
      const gdouble shade = triangles[i].shade;

      cairo_move_to(cr, triangles[i].screen[0][0], triangles[i].screen[0][1]);
      cairo_line_to(cr, triangles[i].screen[1][0], triangles[i].screen[1][1]);
      cairo_line_to(cr, triangles[i].screen[2][0], triangles[i].screen[2][1]);
      cairo_close_path(cr);
      cairo_set_source_rgba(cr,
                            surface->color_rgba[0] * shade,
                            surface->color_rgba[1] * shade,
                            surface->color_rgba[2] * shade,
                            surface->color_rgba[3]);
      cairo_fill_preserve(cr);
      cairo_set_source_rgba(cr,
                            surface->color_rgba[0] * CLAMP(shade + 0.10, 0.0, 1.0),
                            surface->color_rgba[1] * CLAMP(shade + 0.10, 0.0, 1.0),
                            surface->color_rgba[2] * CLAMP(shade + 0.10, 0.0, 1.0),
                            MIN(surface->color_rgba[3] + 0.12, 0.75));
      cairo_set_line_width(cr, 0.6);
      cairo_stroke(cr);
    }
  cairo_restore(cr);

  g_free(triangles);
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

  {
    GdisIsoSurface *surface;

    surface = gdis_gtk4_window_get_isosurface(self, self->active_model, FALSE);
    if (surface)
      gdis_draw_isosurface(cr, self, surface, center, scale, width, height);
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
          gboolean picked;
          gdouble ring_radius;

          atom = g_ptr_array_index(self->active_model->atoms, draw_order[i].atom_index);
          selected = (draw_order[i].atom_index == self->selected_atom_index);
          group_selected = gdis_gtk4_window_atom_array_contains(self->selected_atoms,
                                                                draw_order[i].atom_index);
          picked = gdis_gtk4_window_atom_array_contains(self->picked_atoms,
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

          if (picked)
            {
              ring_radius = draw_order[i].radius + (selected ? 5.5 : 2.8);
              cairo_arc(cr,
                        draw_order[i].screen_x,
                        draw_order[i].screen_y,
                        ring_radius,
                        0.0,
                        2.0 * G_PI);
              cairo_set_source_rgba(cr, 0.58, 0.95, 0.74, 0.92);
              cairo_set_line_width(cr, 1.5);
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

static gboolean
gdis_gtk4_window_export_view_png(GdisGtk4Window *self,
                                 const char *path,
                                 gint width,
                                 gint height,
                                 GError **error)
{
  cairo_surface_t *surface;
  cairo_t *cr;
  cairo_status_t status;

  g_return_val_if_fail(self != NULL, FALSE);
  g_return_val_if_fail(path != NULL, FALSE);

  if (width < 32 || height < 32)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "PNG export size must be at least 32x32.");
      return FALSE;
    }

  surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
  status = cairo_surface_status(surface);
  if (status != CAIRO_STATUS_SUCCESS)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Could not create the export surface: %s",
                  cairo_status_to_string(status));
      cairo_surface_destroy(surface);
      return FALSE;
    }

  cr = cairo_create(surface);
  viewer_draw(GTK_DRAWING_AREA(self->viewer_area), cr, width, height, self);
  cairo_destroy(cr);

  status = cairo_surface_write_to_png(surface, path);
  cairo_surface_destroy(surface);
  if (status != CAIRO_STATUS_SUCCESS)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not write '%s': %s",
                  path,
                  cairo_status_to_string(status));
      return FALSE;
    }

  return TRUE;
}

static void
gdis_gtk4_window_remove_export_dir(const char *dir_path)
{
  GDir *dir;
  const gchar *name;

  if (!dir_path || dir_path[0] == '\0')
    return;

  dir = g_dir_open(dir_path, 0, NULL);
  if (!dir)
    return;

  while ((name = g_dir_read_name(dir)) != NULL)
    {
      g_autofree gchar *path = NULL;

      path = g_build_filename(dir_path, name, NULL);
      g_remove(path);
    }

  g_dir_close(dir);
  g_rmdir(dir_path);
}

static gboolean
gdis_gtk4_window_run_argv(const gchar * const *argv, GError **error)
{
  gint wait_status;
  gchar *stderr_text;

  stderr_text = NULL;
  if (!g_spawn_sync(NULL,
                    (gchar **) argv,
                    NULL,
                    G_SPAWN_SEARCH_PATH,
                    NULL,
                    NULL,
                    NULL,
                    &stderr_text,
                    &wait_status,
                    error))
    {
      g_free(stderr_text);
      return FALSE;
    }

  if (!g_spawn_check_wait_status(wait_status, error))
    {
      if (stderr_text && stderr_text[0] != '\0')
        g_prefix_error(error, "%s", stderr_text);
      g_free(stderr_text);
      return FALSE;
    }

  g_free(stderr_text);
  return TRUE;
}

static gboolean
gdis_gtk4_window_export_animation_sequence(GdisGtk4Window *self,
                                           const char *output_dir,
                                           const char *prefix,
                                           gint width,
                                           gint height,
                                           guint *frame_count_out,
                                           GError **error)
{
  GdisAnimationSourceType source;
  GdisModel *saved_model;
  guint source_count;
  gint active_index;
  guint saved_frame_index;
  gint saved_model_index;
  gboolean success;

  g_return_val_if_fail(self != NULL, FALSE);

  source = gdis_gtk4_window_get_animation_source(self, &source_count, &active_index);
  if (source == GDIS_ANIMATION_SOURCE_NONE || source_count == 0u)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "There is nothing to export yet.");
      return FALSE;
    }
  if (!output_dir || output_dir[0] == '\0')
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Choose an output directory first.");
      return FALSE;
    }
  if (g_mkdir_with_parents(output_dir, 0755) != 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not create %s.",
                  output_dir);
      return FALSE;
    }

  saved_model = self->active_model;
  saved_frame_index = self->active_model ? gdis_model_get_current_frame_index(self->active_model) : 0u;
  saved_model_index = gdis_gtk4_window_get_model_index(self, self->active_model);
  success = TRUE;

  for (guint i = 0; i < source_count; i++)
    {
      g_autofree gchar *basename = NULL;
      g_autofree gchar *path = NULL;

      gdis_gtk4_window_animation_step_to_index(self, (gint) i, FALSE);
      basename = g_strdup_printf("%s_%04u.png",
                                 (prefix && prefix[0] != '\0') ? prefix : "gdis_capture",
                                 i + 1u);
      path = g_build_filename(output_dir, basename, NULL);
      if (!gdis_gtk4_window_export_view_png(self, path, width, height, error))
        {
          success = FALSE;
          break;
        }
    }

  if (source == GDIS_ANIMATION_SOURCE_MODEL_FRAMES && saved_model)
    {
      GError *restore_error;

      restore_error = NULL;
      gdis_model_set_frame_index(saved_model, saved_frame_index, &restore_error);
      g_clear_error(&restore_error);
      gdis_gtk4_window_update_details(self);
      gdis_gtk4_window_refresh_viewer(self);
      gdis_gtk4_window_refresh_measure_tool(self);
      gdis_gtk4_window_refresh_edit_tool(self);
      gdis_gtk4_window_refresh_zmatrix_tool(self);
    }
  else if (source == GDIS_ANIMATION_SOURCE_SESSION_MODELS && saved_model_index >= 0)
    {
      gdis_gtk4_window_animation_step_to_index(self, saved_model_index, FALSE);
    }

  if (frame_count_out)
    *frame_count_out = source_count;
  return success;
}

static gboolean
gdis_gtk4_window_export_movie(GdisGtk4Window *self,
                              const char *output_dir,
                              const char *prefix,
                              gint width,
                              gint height,
                              guint fps,
                              guint format_index,
                              gboolean keep_frames,
                              gchar **movie_path_out,
                              guint *frame_count_out,
                              GError **error)
{
  g_autofree gchar *ffmpeg_path = NULL;
  g_autofree gchar *frames_dir = NULL;
  g_autofree gchar *movie_path = NULL;
  g_autofree gchar *sequence_pattern = NULL;
  g_autofree gchar *movie_name = NULL;
  g_autofree gchar *fps_text = NULL;
  g_autofree gchar *palette_dir = NULL;
  g_autofree gchar *palette_path = NULL;
  guint frame_count;
  gboolean success;
  const gchar *basename;

  g_return_val_if_fail(self != NULL, FALSE);

  ffmpeg_path = g_find_program_in_path("ffmpeg");
  if (!ffmpeg_path)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "ffmpeg is required for movie export.");
      return FALSE;
    }
  if (!output_dir || output_dir[0] == '\0')
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_FAILED,
                  "Choose an output directory first.");
      return FALSE;
    }
  if (g_mkdir_with_parents(output_dir, 0755) != 0)
    {
      g_set_error(error,
                  GDIS_MODEL_ERROR,
                  GDIS_MODEL_ERROR_IO,
                  "Could not create %s.",
                  output_dir);
      return FALSE;
    }

  if (keep_frames)
    frames_dir = g_canonicalize_filename(output_dir, NULL);
  else
    frames_dir = g_dir_make_tmp("gdis-movie-XXXXXX", error);
  if (!frames_dir)
    return FALSE;

  if (!gdis_gtk4_window_export_animation_sequence(self,
                                                  frames_dir,
                                                  prefix,
                                                  width,
                                                  height,
                                                  &frame_count,
                                                  error))
    {
      if (!keep_frames)
        gdis_gtk4_window_remove_export_dir(frames_dir);
      return FALSE;
    }

  basename = (prefix && prefix[0] != '\0') ? prefix : "gdis_capture";
  movie_name = g_strdup_printf("%s.%s", basename, format_index == 1u ? "gif" : "mp4");
  movie_path = g_build_filename(output_dir, movie_name, NULL);
  {
    g_autofree gchar *pattern_name = NULL;

    pattern_name = g_strdup_printf("%s_%%04d.png", basename);
    sequence_pattern = g_build_filename(frames_dir, pattern_name, NULL);
  }
  fps_text = g_strdup_printf("%u", MAX(fps, 1u));
  success = FALSE;

  if (format_index == 1u)
    {
      const gchar * const palette_argv[] = {
        "ffmpeg", "-y", "-framerate", fps_text, "-i", sequence_pattern,
        "-vf", "palettegen", NULL, NULL
      };
      const gchar * const gif_argv[] = {
        "ffmpeg", "-y", "-framerate", fps_text, "-i", sequence_pattern,
        "-i", NULL, "-lavfi", "paletteuse", NULL, NULL
      };
      gchar **palette_argv_mut;
      gchar **gif_argv_mut;

      palette_dir = g_dir_make_tmp("gdis-gif-XXXXXX", error);
      if (!palette_dir)
        goto movie_cleanup;
      palette_path = g_build_filename(palette_dir, "palette.png", NULL);

      palette_argv_mut = g_new0(gchar *, G_N_ELEMENTS(palette_argv));
      gif_argv_mut = g_new0(gchar *, G_N_ELEMENTS(gif_argv));
      memcpy(palette_argv_mut, palette_argv, sizeof(palette_argv));
      memcpy(gif_argv_mut, gif_argv, sizeof(gif_argv));
      palette_argv_mut[8] = palette_path;
      gif_argv_mut[7] = palette_path;
      gif_argv_mut[10] = movie_path;

      if (!gdis_gtk4_window_run_argv((const gchar * const *) palette_argv_mut, error))
        {
          g_free(palette_argv_mut);
          g_free(gif_argv_mut);
          goto movie_cleanup;
        }
      if (!gdis_gtk4_window_run_argv((const gchar * const *) gif_argv_mut, error))
        {
          g_free(palette_argv_mut);
          g_free(gif_argv_mut);
          goto movie_cleanup;
        }

      g_free(palette_argv_mut);
      g_free(gif_argv_mut);
    }
  else
    {
      const gchar * const movie_argv[] = {
        "ffmpeg", "-y", "-framerate", fps_text, "-i", sequence_pattern,
        "-movflags", "+faststart",
        "-pix_fmt", "yuv420p",
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        "-c:v", "libx264",
        NULL, NULL
      };
      gchar **movie_argv_mut;

      movie_argv_mut = g_new0(gchar *, G_N_ELEMENTS(movie_argv));
      memcpy(movie_argv_mut, movie_argv, sizeof(movie_argv));
      movie_argv_mut[14] = movie_path;
      if (!gdis_gtk4_window_run_argv((const gchar * const *) movie_argv_mut, error))
        {
          g_free(movie_argv_mut);
          goto movie_cleanup;
        }
      g_free(movie_argv_mut);
    }

  success = TRUE;

movie_cleanup:
  if (palette_dir)
    gdis_gtk4_window_remove_export_dir(palette_dir);
  if (!keep_frames)
    gdis_gtk4_window_remove_export_dir(frames_dir);

  if (!success)
    return FALSE;

  if (movie_path_out)
    *movie_path_out = g_steal_pointer(&movie_path);
  if (frame_count_out)
    *frame_count_out = frame_count;
  return TRUE;
}

static void
gdis_gtk4_window_handle_completed_pick_mode(GdisGtk4Window *self)
{
  GError *error;
  guint atom_index_a;
  guint atom_index_b;

  g_return_if_fail(self != NULL);

  if (self->click_mode == GDIS_CLICK_MODE_SELECT ||
      !self->picked_atoms ||
      self->picked_atoms->len < 2)
    return;

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
              if (!gdis_gtk4_window_parse_entry_uint(self->edit_tool->bond_order_entry,
                                                    "Bond order",
                                                    &order,
                                                    &parse_error))
                {
                  g_clear_error(&parse_error);
                  order = 1;
                }
            }

          gdis_gtk4_window_push_undo_snapshot(self, NULL);
          if (!gdis_model_add_explicit_bond(self->active_model,
                                            atom_index_a,
                                            atom_index_b,
                                            (guint8) order,
                                            &error))
            {
              gdis_gtk4_window_discard_undo_snapshot(self);
              gdis_gtk4_window_log(self, "Pick Add Bond failed: %s\n",
                                   error ? error->message : "unknown error");
            }
          else
            {
              gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
              gdis_gtk4_window_log(self, "Pick Add Bond completed.\n");
            }
        }
      else
        {
          gdis_gtk4_window_push_undo_snapshot(self, NULL);
          if (!gdis_model_remove_bond(self->active_model,
                                      atom_index_a,
                                      atom_index_b,
                                      &error))
            {
              gdis_gtk4_window_discard_undo_snapshot(self);
              gdis_gtk4_window_log(self, "Pick Remove Bond failed: %s\n",
                                   error ? error->message : "unknown error");
            }
          else
            {
              gdis_gtk4_window_refresh_after_model_edit(self, FALSE);
              gdis_gtk4_window_log(self, "Pick Remove Bond completed.\n");
            }
        }
    }

  g_clear_error(&error);
  gdis_gtk4_window_set_click_mode(self, GDIS_CLICK_MODE_SELECT);
  gdis_gtk4_window_log(self, "Click mode returned to Select.\n");
}

static guint
gdis_gtk4_window_hit_test_atom_at(GdisGtk4Window *self,
                                  gdouble x,
                                  gdouble y,
                                  GdisProjectedAtom *projected_out)
{
  GdisProjectedAtom *projected;
  gdouble center[3];
  gdouble scale;
  guint atom_count;
  guint cell_count;
  guint best_index;
  guint best_projected_index;
  gdouble best_distance2;

  g_return_val_if_fail(self != NULL, INVALID_ATOM_INDEX);

  projected = NULL;
  atom_count = 0;
  cell_count = 0;
  best_index = INVALID_ATOM_INDEX;
  best_projected_index = INVALID_ATOM_INDEX;
  best_distance2 = G_MAXDOUBLE;

  if (!self->active_model || !self->viewer_area)
    return INVALID_ATOM_INDEX;

  if (!gdis_prepare_projection(self,
                               gtk_widget_get_width(self->viewer_area),
                               gtk_widget_get_height(self->viewer_area),
                               &projected,
                               &atom_count,
                               &cell_count,
                               center,
                               &scale))
    return INVALID_ATOM_INDEX;
  (void) cell_count;

  for (guint i = 0; i < atom_count; i++)
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

  if (best_index != INVALID_ATOM_INDEX && projected_out)
    *projected_out = projected[best_projected_index];

  g_free(projected);
  return best_index;
}

static gboolean
gdis_gtk4_window_apply_box_selection(GdisGtk4Window *self,
                                     gdouble x0,
                                     gdouble y0,
                                     gdouble x1,
                                     gdouble y1,
                                     gboolean toggle_existing)
{
  GdisProjectedAtom *projected;
  GArray *boxed_atoms;
  GArray *base_selection;
  gdouble center[3];
  gdouble scale;
  gdouble min_x;
  gdouble min_y;
  gdouble max_x;
  gdouble max_y;
  guint atom_count;
  guint cell_count;

  g_return_val_if_fail(self != NULL, FALSE);

  projected = NULL;
  boxed_atoms = NULL;
  base_selection = NULL;
  atom_count = 0;
  cell_count = 0;

  if (!self->active_model || !self->viewer_area)
    return FALSE;

  if (!gdis_prepare_projection(self,
                               gtk_widget_get_width(self->viewer_area),
                               gtk_widget_get_height(self->viewer_area),
                               &projected,
                               &atom_count,
                               &cell_count,
                               center,
                               &scale))
    return FALSE;
  (void) center;
  (void) scale;
  (void) cell_count;

  min_x = MIN(x0, x1);
  min_y = MIN(y0, y1);
  max_x = MAX(x0, x1);
  max_y = MAX(y0, y1);
  boxed_atoms = g_array_new(FALSE, FALSE, sizeof(guint));

  for (guint i = 0; i < atom_count; i++)
    {
      gdouble radius;

      radius = projected[i].radius;
      if (projected[i].screen_x + radius < min_x ||
          projected[i].screen_x - radius > max_x ||
          projected[i].screen_y + radius < min_y ||
          projected[i].screen_y - radius > max_y)
        continue;

      if (!gdis_gtk4_window_atom_array_contains(boxed_atoms, projected[i].atom_index))
        g_array_append_val(boxed_atoms, projected[i].atom_index);
    }

  g_free(projected);

  if (boxed_atoms->len == 0)
    {
      g_array_free(boxed_atoms, TRUE);
      if (!toggle_existing)
        {
          self->selected_atom_index = INVALID_ATOM_INDEX;
          gdis_gtk4_window_clear_selected_atoms(self);
          gdis_gtk4_window_finish_selection_change(self, "Selection cleared.\n");
        }
      return FALSE;
    }

  base_selection = toggle_existing ? gdis_gtk4_window_copy_atom_array(self->selected_atoms)
                                   : g_array_new(FALSE, FALSE, sizeof(guint));

  for (guint i = 0; i < boxed_atoms->len; i++)
    {
      guint atom_index;

      atom_index = g_array_index(boxed_atoms, guint, i);
      if (toggle_existing)
        gdis_gtk4_window_toggle_atom_in_array(base_selection, atom_index);
      else if (!gdis_gtk4_window_atom_array_contains(base_selection, atom_index))
        g_array_append_val(base_selection, atom_index);
    }

  gdis_gtk4_window_clear_selected_atoms(self);
  if (base_selection->len > 0)
    g_array_append_vals(self->selected_atoms, base_selection->data, base_selection->len);
  self->selected_atom_index = (self->selected_atoms && self->selected_atoms->len > 0)
                              ? g_array_index(self->selected_atoms, guint, 0)
                              : INVALID_ATOM_INDEX;
  gdis_gtk4_window_finish_selection_change(self,
                                           toggle_existing ?
                                           "Shift box-selection toggled atoms in the current selection.\n" :
                                           "Box-selected atoms in the viewer.\n");

  g_array_free(base_selection, TRUE);
  g_array_free(boxed_atoms, TRUE);
  return TRUE;
}

static void
gdis_gtk4_window_select_atom_at(GdisGtk4Window *self,
                                gdouble x,
                                gdouble y,
                                GdkModifierType modifiers)
{
  GdisProjectedAtom projected_hit;
  guint best_index;

  g_return_if_fail(self != NULL);

  best_index = gdis_gtk4_window_hit_test_atom_at(self, x, y, &projected_hit);
  if (best_index != INVALID_ATOM_INDEX)
    {
      const GdisAtom *atom;

      atom = g_ptr_array_index(self->active_model->atoms, best_index);
      if (gdis_gtk4_window_should_record_viewer_picks(self))
        {
          gdis_gtk4_window_remember_atom_pick(self, best_index);
          gdis_gtk4_window_log(self,
                               "Picked atom in %s: %s [%s] #%u at %.4f %.4f %.4f (image offset: %d %d %d, pick set: %u)\n",
                               self->active_model->basename,
                               atom->label,
                               atom->element,
                               atom->serial,
                               atom->position[0],
                               atom->position[1],
                               atom->position[2],
                               projected_hit.image_offset[0],
                               projected_hit.image_offset[1],
                               projected_hit.image_offset[2],
                               self->picked_atoms ? self->picked_atoms->len : 0);
          gdis_gtk4_window_handle_completed_pick_mode(self);
          gdis_gtk4_window_update_details(self);
          gdis_gtk4_window_refresh_viewer(self);
          gdis_gtk4_window_refresh_measure_tool(self);
          gdis_gtk4_window_refresh_edit_tool(self);
          return;
        }

      self->selected_atom_index = best_index;
      if ((modifiers & GDK_SHIFT_MASK) != 0)
        {
          GArray *previous_selection;
          GArray *clicked_group;

          previous_selection = gdis_gtk4_window_copy_atom_array(self->selected_atoms);
          gdis_gtk4_window_apply_selection_mode(self, best_index);
          clicked_group = gdis_gtk4_window_copy_atom_array(self->selected_atoms);
          gdis_gtk4_window_clear_selected_atoms(self);
          if (previous_selection->len > 0)
            g_array_append_vals(self->selected_atoms,
                                previous_selection->data,
                                previous_selection->len);
          for (guint i = 0; i < clicked_group->len; i++)
            gdis_gtk4_window_toggle_atom_in_array(self->selected_atoms,
                                                  g_array_index(clicked_group, guint, i));
          g_array_free(previous_selection, TRUE);
          g_array_free(clicked_group, TRUE);
          gdis_gtk4_window_finish_selection_change(self,
                                                   "Shift-toggled the clicked atoms into the current selection.\n");
        }
      else
        {
          g_autofree gchar *message = NULL;

          gdis_gtk4_window_apply_selection_mode(self, best_index);
          message = g_strdup_printf("Selected atom in %s: %s [%s] #%u at %.4f %.4f %.4f (image offset: %d %d %d, mode: %s, selected set: %u)\n",
                                    self->active_model->basename,
                                    atom->label,
                                    atom->element,
                                    atom->serial,
                                    atom->position[0],
                                    atom->position[1],
                                    atom->position[2],
                                    projected_hit.image_offset[0],
                                    projected_hit.image_offset[1],
                                    projected_hit.image_offset[2],
                                    gdis_gtk4_window_selection_mode_label(self->selection_mode),
                                    self->selected_atoms ? self->selected_atoms->len : 0);
          gdis_gtk4_window_finish_selection_change(self, message);
        }
    }
  else if ((modifiers & GDK_SHIFT_MASK) == 0 &&
           !gdis_gtk4_window_should_record_viewer_picks(self) &&
           self->selected_atom_index != INVALID_ATOM_INDEX)
    {
      self->selected_atom_index = INVALID_ATOM_INDEX;
      gdis_gtk4_window_clear_selected_atoms(self);
      gdis_gtk4_window_finish_selection_change(self, "Selection cleared.\n");
    }
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
  gdis_gtk4_window_sync_display_tool(self);
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
  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_set_image_limits(self->active_model, limits, &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
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
  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_set_image_limits(self->active_model, limits, &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
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
  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_confine_atoms_to_cell(self->active_model, &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
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
  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_confine_molecules_to_cell(self->active_model, &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
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
  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_force_p1(self->active_model, &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
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
  gdis_gtk4_window_push_undo_snapshot(self, NULL);
  if (!gdis_model_make_supercell(self->active_model, repeat_a, repeat_b, repeat_c, &error))
    {
      gdis_gtk4_window_discard_undo_snapshot(self);
      gdis_gtk4_window_log(self, "Make supercell failed: %s\n",
                           error ? error->message : "unknown error");
      g_clear_error(&error);
      return;
    }

  self->selected_atom_index = INVALID_ATOM_INDEX;
  self->fragment_anchor_index = INVALID_ATOM_INDEX;
  gdis_gtk4_window_clear_selected_atoms(self);
  gdis_gtk4_window_clear_atom_picks(self);
  gdis_gtk4_window_clear_saved_measurements(self,
                                            self->active_model,
                                            "Saved measurements were cleared because the supercell rebuild changed atom indices.\n");
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
  GdkModifierType modifiers;

  self = user_data;
  modifiers = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(gesture));
  self->press_modifiers = modifiers;
  self->press_x = start_x;
  self->press_y = start_y;
  self->drag_current_x = start_x;
  self->drag_current_y = start_y;
  self->drag_origin_x = self->rotation_x;
  self->drag_origin_y = self->rotation_y;
  if (gdis_gtk4_window_should_record_viewer_picks(self))
    {
      self->drag_mode = GDIS_DRAG_MODE_NONE;
      return;
    }

  if ((modifiers & GDK_ALT_MASK) != 0 ||
      gdis_gtk4_window_hit_test_atom_at(self, start_x, start_y, NULL) != INVALID_ATOM_INDEX)
    self->drag_mode = GDIS_DRAG_MODE_ROTATE;
  else
    self->drag_mode = GDIS_DRAG_MODE_BOX_SELECT;
}

static void
on_viewer_secondary_drag_begin(GtkGestureDrag *gesture,
                               gdouble start_x,
                               gdouble start_y,
                               gpointer user_data)
{
  GdisGtk4Window *self;

  self = user_data;
  self->press_modifiers = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(gesture));
  self->press_x = start_x;
  self->press_y = start_y;
  self->drag_current_x = start_x;
  self->drag_current_y = start_y;
  self->drag_origin_x = self->rotation_x;
  self->drag_origin_y = self->rotation_y;
  self->drag_mode = gdis_gtk4_window_should_record_viewer_picks(self)
                      ? GDIS_DRAG_MODE_NONE
                      : GDIS_DRAG_MODE_ROTATE;
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
  if (self->drag_mode == GDIS_DRAG_MODE_ROTATE)
    {
      self->rotation_y = self->drag_origin_y + offset_x * 0.010;
      self->rotation_x = self->drag_origin_x + offset_y * 0.010;
      gdis_gtk4_window_refresh_viewer(self);
    }
  else if (self->drag_mode == GDIS_DRAG_MODE_BOX_SELECT)
    {
      self->drag_current_x = self->press_x + offset_x;
      self->drag_current_y = self->press_y + offset_y;
      gdis_gtk4_window_refresh_viewer(self);
    }
}

static void
on_viewer_secondary_drag_update(GtkGestureDrag *gesture,
                                gdouble offset_x,
                                gdouble offset_y,
                                gpointer user_data)
{
  on_viewer_drag_update(gesture, offset_x, offset_y, user_data);
}

static void
on_viewer_drag_end(GtkGestureDrag *gesture,
                   gdouble offset_x,
                   gdouble offset_y,
                   gpointer user_data)
{
  GdisGtk4Window *self;

  (void) gesture;

  self = user_data;
  if (self->drag_mode == GDIS_DRAG_MODE_BOX_SELECT &&
      (fabs(offset_x) > 5.0 || fabs(offset_y) > 5.0))
    {
      self->drag_current_x = self->press_x + offset_x;
      self->drag_current_y = self->press_y + offset_y;
      gdis_gtk4_window_apply_box_selection(self,
                                           self->press_x,
                                           self->press_y,
                                           self->drag_current_x,
                                           self->drag_current_y,
                                           (self->press_modifiers & GDK_SHIFT_MASK) != 0);
    }

  self->drag_mode = GDIS_DRAG_MODE_NONE;
  self->drag_current_x = self->press_x;
  self->drag_current_y = self->press_y;
  gdis_gtk4_window_refresh_viewer(self);
}

static void
on_viewer_secondary_drag_end(GtkGestureDrag *gesture,
                             gdouble offset_x,
                             gdouble offset_y,
                             gpointer user_data)
{
  on_viewer_drag_end(gesture, offset_x, offset_y, user_data);
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
  self->drag_current_x = x;
  self->drag_current_y = y;
  self->press_modifiers = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(gesture));
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
    gdis_gtk4_window_select_atom_at(self, x, y, self->press_modifiers);
}

static void
on_open_multiple_dialog_complete(GObject *source_object,
                                 GAsyncResult *result,
                                 gpointer user_data)
{
  GdisGtk4Window *self;
  GtkFileDialog *dialog;
  GListModel *files;
  GError *error;
  guint n_files;
  guint loaded_count;

  dialog = GTK_FILE_DIALOG(source_object);
  self = user_data;

  error = NULL;
  files = gtk_file_dialog_open_multiple_finish(dialog, result, &error);
  if (!files)
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

  n_files = g_list_model_get_n_items(files);
  loaded_count = 0;
  for (guint i = 0; i < n_files; i++)
    {
      GFile *file;
      gchar *path;
      gboolean loaded;

      file = g_list_model_get_item(files, i);
      if (!file)
        continue;

      path = g_file_get_path(file);
      if (path)
        {
          loaded = gdis_gtk4_window_add_model_from_path(self, path, TRUE);
          if (loaded)
            loaded_count++;
          g_free(path);
        }

      g_object_unref(file);
    }

  if (loaded_count > 0)
    gdis_gtk4_window_log(self, "Loaded %u model file%s from selection.\n",
                         loaded_count,
                         loaded_count == 1 ? "" : "s");

  g_object_unref(files);
  g_object_unref(dialog);
}

static void
present_open_dialog(GdisGtk4Window *self)
{
  GtkFileDialog *dialog;

  g_return_if_fail(self != NULL);

  dialog = gtk_file_dialog_new();
  gtk_file_dialog_set_title(dialog, "Open model(s)");
  gtk_file_dialog_open_multiple(dialog,
                                GTK_WINDOW(self->window),
                                NULL,
                                on_open_multiple_dialog_complete,
                                self);
}

static gboolean
gdis_gtk4_window_model_path_supported(const char *path)
{
  g_return_val_if_fail(path != NULL, FALSE);

  return gdis_model_format_from_path(path) != GDIS_MODEL_FORMAT_UNKNOWN;
}

static guint
gdis_gtk4_window_open_models_from_directory(GdisGtk4Window *self,
                                            const char *dir_path,
                                            gboolean *open_failed_out)
{
  GDir *dir;
  GError *error;
  const gchar *name;
  guint loaded_count;
  guint skipped_unsupported_count;
  guint failed_count;

  g_return_val_if_fail(self != NULL, 0);
  g_return_val_if_fail(dir_path != NULL, 0);

  if (open_failed_out)
    *open_failed_out = FALSE;

  error = NULL;
  dir = g_dir_open(dir_path, 0, &error);
  if (!dir)
    {
      if (open_failed_out)
        *open_failed_out = TRUE;

      gdis_gtk4_window_log(self,
                           "Could not open directory: %s%s%s%s\n",
                           dir_path,
                           error ? " (" : "",
                           error ? error->message : "",
                           error ? ")" : "");
      if (error)
        {
          gdis_gtk4_window_log(self,
                               "macOS may block direct folder reads for app bundles in protected locations.\n"
                               "Use File > Open Models Folder... once and pick this folder to grant access.\n");
          g_error_free(error);
        }
      return 0;
    }

  loaded_count = 0;
  skipped_unsupported_count = 0;
  failed_count = 0;
  while ((name = g_dir_read_name(dir)) != NULL)
    {
      g_autofree gchar *candidate = NULL;
      gboolean loaded;

      if (name[0] == '.')
        continue;

      candidate = g_build_filename(dir_path, name, NULL);
      if (!g_file_test(candidate, G_FILE_TEST_IS_REGULAR))
        continue;
      if (!gdis_gtk4_window_model_path_supported(candidate))
        {
          skipped_unsupported_count++;
          continue;
        }

      loaded = gdis_gtk4_window_add_model_from_path(self, candidate, TRUE);
      if (loaded)
        loaded_count++;
      else
        failed_count++;
    }

  g_dir_close(dir);
  gdis_gtk4_window_log(self,
                       "Loaded %u model file%s from folder: %s (skipped unsupported: %u, failed to load: %u)\n",
                       loaded_count,
                       loaded_count == 1 ? "" : "s",
                       dir_path,
                       skipped_unsupported_count,
                       failed_count);
  return loaded_count;
}

static void
on_open_models_folder_dialog_complete(GObject *source_object,
                                      GAsyncResult *result,
                                      gpointer user_data)
{
  GdisGtk4Window *self;
  GtkFileDialog *dialog;
  GFile *folder;
  GError *error;
  gchar *path;

  dialog = GTK_FILE_DIALOG(source_object);
  self = user_data;

  error = NULL;
  folder = gtk_file_dialog_select_folder_finish(dialog, result, &error);
  if (!folder)
    {
      if (error)
        {
          if (!g_error_matches(error, GTK_DIALOG_ERROR, GTK_DIALOG_ERROR_CANCELLED) &&
              !g_error_matches(error, G_IO_ERROR, G_IO_ERROR_CANCELLED))
            gdis_gtk4_window_log(self, "Open folder failed: %s\n", error->message);
          g_error_free(error);
        }
      g_object_unref(dialog);
      return;
    }

  path = g_file_get_path(folder);
  if (path)
    {
      gdis_gtk4_window_open_models_from_directory(self, path, NULL);
      g_free(path);
    }

  g_object_unref(folder);
  g_object_unref(dialog);
}

static void
present_open_models_folder_dialog_at(GdisGtk4Window *self,
                                     const char *initial_dir)
{
  GtkFileDialog *dialog;

  g_return_if_fail(self != NULL);

  dialog = gtk_file_dialog_new();
  gtk_file_dialog_set_title(dialog, "Open all models in folder");
  if (initial_dir && initial_dir[0] != '\0')
    {
      GFile *initial_folder;

      initial_folder = g_file_new_for_path(initial_dir);
      gtk_file_dialog_set_initial_folder(dialog, initial_folder);
      g_object_unref(initial_folder);
    }
  gtk_file_dialog_select_folder(dialog,
                                GTK_WINDOW(self->window),
                                NULL,
                                on_open_models_folder_dialog_complete,
                                self);
}

static void
present_open_models_folder_dialog(GdisGtk4Window *self)
{
  present_open_models_folder_dialog_at(self, NULL);
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
                       "No active model\nOpen a file to load a structure.");

  gtk_frame_set_child(GTK_FRAME(frame), box);

  return frame;
}

static GtkWidget *
build_model_list(GdisGtk4Window *self)
{
  GtkWidget *list;

  list = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
  gtk_widget_set_hexpand(list, TRUE);

  self->model_list = list;

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
  GtkGesture *secondary_drag;
  GtkGesture *click;
  GtkEventController *scroll;

  frame = gtk_frame_new(NULL);
  gtk_widget_add_css_class(frame, "gdis-viewer-frame");
  gtk_widget_set_hexpand(frame, TRUE);
  gtk_widget_set_vexpand(frame, TRUE);
  gtk_widget_set_size_request(frame, GDIS_VIEWER_MIN_WIDTH, GDIS_VIEWER_MIN_HEIGHT);

  drawing_area = gtk_drawing_area_new();
  gtk_widget_set_hexpand(drawing_area, TRUE);
  gtk_widget_set_vexpand(drawing_area, TRUE);
  gtk_drawing_area_set_content_width(GTK_DRAWING_AREA(drawing_area), GDIS_VIEWER_MIN_WIDTH);
  gtk_drawing_area_set_content_height(GTK_DRAWING_AREA(drawing_area), GDIS_VIEWER_MIN_HEIGHT);
  gtk_widget_set_cursor_from_name(drawing_area, "crosshair");
  self->viewer_area = drawing_area;
  gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawing_area), viewer_draw, self, NULL);

  drag = gtk_gesture_drag_new();
  gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag), GDK_BUTTON_PRIMARY);
  g_signal_connect(drag, "drag-begin", G_CALLBACK(on_viewer_drag_begin), self);
  g_signal_connect(drag, "drag-update", G_CALLBACK(on_viewer_drag_update), self);
  g_signal_connect(drag, "drag-end", G_CALLBACK(on_viewer_drag_end), self);
  gtk_widget_add_controller(drawing_area, GTK_EVENT_CONTROLLER(drag));

  secondary_drag = gtk_gesture_drag_new();
  gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(secondary_drag), GDK_BUTTON_SECONDARY);
  g_signal_connect(secondary_drag, "drag-begin", G_CALLBACK(on_viewer_secondary_drag_begin), self);
  g_signal_connect(secondary_drag, "drag-update", G_CALLBACK(on_viewer_secondary_drag_update), self);
  g_signal_connect(secondary_drag, "drag-end", G_CALLBACK(on_viewer_secondary_drag_end), self);
  gtk_widget_add_controller(drawing_area, GTK_EVENT_CONTROLLER(secondary_drag));

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
  gtk_widget_set_size_request(scroller, 760, GDIS_STATUS_MIN_HEIGHT);
  gtk_scrolled_window_set_min_content_height(GTK_SCROLLED_WINDOW(scroller), GDIS_STATUS_MIN_HEIGHT);

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
  GtkWidget *sidebar;
  GtkWidget *sidebar_scroller;
  GtkWidget *viewer;
  GtkWidget *status;

  main_paned = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
  gtk_widget_set_hexpand(main_paned, TRUE);
  gtk_widget_set_vexpand(main_paned, TRUE);
  gtk_widget_set_size_request(main_paned, GDIS_MAIN_ROOT_WIDTH, GDIS_MAIN_ROOT_HEIGHT);
  gtk_paned_set_wide_handle(GTK_PANED(main_paned), TRUE);
  self->main_paned = main_paned;

  right_paned = gtk_paned_new(GTK_ORIENTATION_VERTICAL);
  gtk_widget_set_hexpand(right_paned, TRUE);
  gtk_widget_set_vexpand(right_paned, TRUE);
  gtk_widget_set_size_request(right_paned, GDIS_VIEWER_MIN_WIDTH, 640);
  gtk_paned_set_wide_handle(GTK_PANED(right_paned), TRUE);
  self->right_paned = right_paned;

  viewer = build_viewer_area(self);
  status = build_status_view(self);

  sidebar = build_sidebar(self);
  sidebar_scroller = gtk_scrolled_window_new();
  gtk_widget_set_hexpand(sidebar_scroller, FALSE);
  gtk_widget_set_vexpand(sidebar_scroller, TRUE);
  gtk_widget_set_size_request(sidebar_scroller, GDIS_MAIN_SIDEBAR_WIDTH, -1);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sidebar_scroller),
                                 GTK_POLICY_NEVER,
                                 GTK_POLICY_AUTOMATIC);
  gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(sidebar_scroller), sidebar);

  gtk_paned_set_start_child(GTK_PANED(main_paned), sidebar_scroller);
  gtk_paned_set_end_child(GTK_PANED(main_paned), right_paned);
  gtk_paned_set_resize_start_child(GTK_PANED(main_paned), FALSE);
  gtk_paned_set_shrink_start_child(GTK_PANED(main_paned), FALSE);
  gtk_paned_set_resize_end_child(GTK_PANED(main_paned), TRUE);
  gtk_paned_set_shrink_end_child(GTK_PANED(main_paned), FALSE);
  gtk_paned_set_position(GTK_PANED(main_paned), GDIS_MAIN_SIDEBAR_WIDTH);

  gtk_paned_set_start_child(GTK_PANED(right_paned), viewer);
  gtk_paned_set_end_child(GTK_PANED(right_paned), status);
  gtk_paned_set_resize_start_child(GTK_PANED(right_paned), TRUE);
  gtk_paned_set_shrink_start_child(GTK_PANED(right_paned), FALSE);
  gtk_paned_set_resize_end_child(GTK_PANED(right_paned), FALSE);
  gtk_paned_set_shrink_end_child(GTK_PANED(right_paned), FALSE);
  gtk_paned_set_position(GTK_PANED(right_paned), GDIS_RIGHT_PANED_POSITION);

  return main_paned;
}

static void
gdis_gtk4_window_restore_layout(GdisGtk4Window *self)
{
  g_return_if_fail(self != NULL);

  if (self->main_paned)
    gtk_paned_set_position(GTK_PANED(self->main_paned), GDIS_MAIN_SIDEBAR_WIDTH);

  if (self->right_paned)
    gtk_paned_set_position(GTK_PANED(self->right_paned), GDIS_RIGHT_PANED_POSITION);
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
    gtk_window_set_default_size(GTK_WINDOW(self->window),
                                GDIS_MAIN_WINDOW_WIDTH,
                                GDIS_MAIN_WINDOW_HEIGHT);

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
  GMenu *visualization_menu;
  GMenu *building_menu;
  GMenu *computation_menu;
  GMenu *analysis_menu;
  GMenu *view_menu;
  GMenu *view_general_section;
  GMenu *view_mode_section;
  GMenu *view_tools_section;
  GMenu *help_menu;

  root = g_menu_new();

  file_menu = g_menu_new();
  g_menu_append(file_menu, "New", "app.new-model");
  g_menu_append(file_menu, "Open", "app.open");
  g_menu_append(file_menu, "Open ./models/*", "app.open-models-glob");
  g_menu_append(file_menu, "Open Models Folder...", "app.open-models-folder");
  g_menu_append(file_menu, "Save", "app.save");
  g_menu_append(file_menu, "Save As", "app.save-as");
  g_menu_append(file_menu, "Close", "app.close-model");
  g_menu_append(file_menu, "Quit", "app.quit");
  g_menu_append_submenu(root, "File", G_MENU_MODEL(file_menu));

  edit_menu = g_menu_new();
  g_menu_append(edit_menu, "Undo", "app.undo");
  g_menu_append(edit_menu, "Select All", "app.select-all");
  g_menu_append(edit_menu, "Invert Selection", "app.invert-selection");
  g_menu_append(edit_menu, "Delete Selected", "app.delete-selected");
  g_menu_append_submenu(root, "Edit", G_MENU_MODEL(edit_menu));

  tools_menu = g_menu_new();
  visualization_menu = g_menu_new();
  g_menu_append(visualization_menu, "Animation...", "app.animation");
  g_menu_append(visualization_menu, "Iso-surfaces...", "app.isosurface");
  g_menu_append(visualization_menu, "Periodic table...", "app.periodic-table");
  g_menu_append_submenu(tools_menu, "Visualization", G_MENU_MODEL(visualization_menu));

  building_menu = g_menu_new();
  g_menu_append(building_menu, "Editing...", "app.edit");
  g_menu_append(building_menu, "Dislocations...", "app.dislocations");
  g_menu_append(building_menu, "Docking...", "app.docking");
  g_menu_append(building_menu, "Dynamics...", "app.building-dynamics");
  g_menu_append(building_menu, "Surfaces...", "app.surface");
  g_menu_append(building_menu, "Zmatrix...", "app.zmatrix");
  g_menu_append_submenu(tools_menu, "Building", G_MENU_MODEL(building_menu));

  computation_menu = g_menu_new();
  g_menu_append(computation_menu, "Diffraction...", "app.diffraction");
  g_menu_append(computation_menu, "GULP...", "app.gulp");
  g_menu_append(computation_menu, "GAMESS...", "app.gamess");
  g_menu_append(computation_menu, "Monty...", "app.monty");
  g_menu_append(computation_menu, "Qbox...", "app.qbox");
  g_menu_append(computation_menu, "SIESTA...", "app.siesta");
  g_menu_append(computation_menu, "VASP...", "app.vasp");
  g_menu_append(computation_menu, "USPEX...", "app.uspex");
  g_menu_append_submenu(tools_menu, "Computation", G_MENU_MODEL(computation_menu));

  analysis_menu = g_menu_new();
  g_menu_append(analysis_menu, "Dynamics...", "app.analysis-dynamics");
  g_menu_append(analysis_menu, "Measurements...", "app.measure");
  g_menu_append(analysis_menu, "Plots...", "app.plots");
  g_menu_append_submenu(tools_menu, "Analysis", G_MENU_MODEL(analysis_menu));

  g_menu_append_submenu(root, "Tools", G_MENU_MODEL(tools_menu));

  view_menu = g_menu_new();
  view_general_section = g_menu_new();
  g_menu_append(view_general_section, "Display properties...", "app.render");
  g_menu_append(view_general_section, "Reset model images", "app.reset-images");
  g_menu_append_section(view_menu, NULL, G_MENU_MODEL(view_general_section));

  view_mode_section = g_menu_new();
  g_menu_append(view_mode_section, "Normal mode", "app.normal-mode");
  g_menu_append(view_mode_section, "Recording mode", "app.recording-mode");
  g_menu_append_section(view_menu, NULL, G_MENU_MODEL(view_mode_section));

  view_tools_section = g_menu_new();
  g_menu_append(view_tools_section, "Task manager...", "app.task-manager");
  g_menu_append(view_tools_section, "Executable paths...", "app.executable-paths");
  g_menu_append_section(view_menu, NULL, G_MENU_MODEL(view_tools_section));
  g_menu_append_submenu(root, "View", G_MENU_MODEL(view_menu));

  help_menu = g_menu_new();
  g_menu_append(help_menu, "About...", "app.about");
  g_menu_append(help_menu, "Manual...", "app.manual");
  g_menu_append_submenu(root, "Help", G_MENU_MODEL(help_menu));

  g_object_unref(file_menu);
  g_object_unref(edit_menu);
  g_object_unref(tools_menu);
  g_object_unref(visualization_menu);
  g_object_unref(building_menu);
  g_object_unref(computation_menu);
  g_object_unref(analysis_menu);
  g_object_unref(view_menu);
  g_object_unref(view_general_section);
  g_object_unref(view_mode_section);
  g_object_unref(view_tools_section);
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

  if (g_strcmp0(name, "open-models-folder") == 0)
    {
      present_open_models_folder_dialog(self);
      return;
    }

  if (g_strcmp0(name, "open-models-glob") == 0)
    {
      g_autofree gchar *models_dir = NULL;
      gboolean open_failed;

      models_dir = gdis_gtk4_window_resolve_path("./models");
      if (!models_dir || !g_file_test(models_dir, G_FILE_TEST_IS_DIR))
        {
          gdis_gtk4_window_log(self, "Could not resolve ./models from this launch location.\n");
          return;
        }

      open_failed = FALSE;
      gdis_gtk4_window_open_models_from_directory(self, models_dir, &open_failed);
      if (open_failed)
        present_open_models_folder_dialog_at(self, models_dir);
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

  if (g_strcmp0(name, "undo") == 0)
    {
      if (!gdis_gtk4_window_perform_undo(self))
        gdis_gtk4_window_log(self, "Undo requested, but no saved model edit is available.\n");
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

  if (g_strcmp0(name, "reset-images") == 0)
    {
      on_reset_image_limits_clicked(NULL, self);
      return;
    }

  if (g_strcmp0(name, "about") == 0)
    {
      g_autofree char *report = NULL;

      report = g_strdup(
        "GDIS GTK4 Rebuild\n\n"
        "Current restored core:\n"
        "  - native GTK4 viewer with rotate / zoom / atom picking\n"
        "  - legacy-derived loaders for XYZ, PDB, ARC/CAR, and CIF, with multi-frame XYZ/ARC trajectories restored\n"
        "  - structure editing, bond editing, measurements, and undo\n"
        "  - periodic image controls, slab surface builder, and powder diffraction tool\n"
        "  - selected-fragment iso-surfaces, periodic table, and analytic electron-density approximation\n"
        "  - per-model animation playback with Linux-style Control/Processing/Rendering panels, PNG/movie recording export, Z-matrix row editing with geometry rebuild, native dislocation transforms, docking project generation, and executable-path management\n\n"
        "Important docs in this repo:\n"
        "  - gtk4_rebuild/gtk4_app/README.md\n"
        "  - gtk4_rebuild/RESTORATION_AUDIT.md\n"
        "  - gtk4_rebuild/legacy_snapshot/src/gui_main.c");
      gdis_gtk4_window_present_report(self, "About GDIS GTK4", report);
      gdis_gtk4_window_log(self, "Opened GTK4 rebuild summary.\n");
      return;
    }

  if (g_strcmp0(name, "manual") == 0)
    {
      g_autofree char *report = NULL;

      report = g_strdup(
        "GTK4 Usage Guide\n\n"
        "Basic workflow:\n"
        "  1. File > Open to load a structure.\n"
        "  2. Drag an atom to rotate, or use Alt-drag / right-drag anywhere to rotate, and scroll to zoom.\n"
        "  3. Click atoms to build the current selection and pick history.\n"
        "  4. Use Tools > Building > Editing for atom edits, add atom, and bond editing.\n"
        "  5. Use Tools > Analysis > Measurements for distance / angle / torsion.\n"
        "  6. Use Tools > Visualization > Periodic table to copy or apply an element.\n"
        "  7. Use Tools > Visualization > Iso-surfaces for molecular, promolecule, Hirshfeld-style, or analytic electron-density surfaces.\n"
        "  8. Use Tools > Building > Zmatrix to rebuild an internal-coordinate editor from the whole model or current selection, then Recompute Geometry to apply edits.\n"
        "  9. Use Tools > Building > Dislocations or Docking for the restored native builders.\n"
        " 10. Use Tools > Computation > Qbox to generate an editable starter deck, write it into a working directory, and launch the configured qbox/qb executable.\n"
        " 11. Use Tools > Visualization > Animation for the tabbed playback panel, then open Record to export PNG snapshots, sequences, MP4 movies, or animated GIFs.\n"
        " 12. Use View > Reset Model Images for periodic image cleanup.\n\n"
        "Current GTK4-native tools:\n"
        "  - Editing\n"
        "  - Measurements\n"
        "  - Surface Builder\n"
        "  - Powder Diffraction\n"
        "  - Periodic Table\n"
        "  - Animation (active-model frames when available, otherwise loaded-model playback)\n"
        "  - Recording export (PNG sequence / MP4 / GIF)\n"
        "  - Z-matrix editor with explicit geometry rebuild\n"
        "  - Dislocation builder\n"
        "  - Docking project generator\n"
        "  - Qbox deck editor / runner plus session executable-path integration\n"
        "  - Executable Paths / Task Manager, including GULP, GAMESS, Monty, Qbox, SIESTA, VASP, and USPEX session paths\n\n"
        "Detailed parity audit:\n"
        "  gtk4_rebuild/RESTORATION_AUDIT.md\n\n"
        "Legacy reference:\n"
        "  gtk4_rebuild/legacy_snapshot/src/gui_main.c");
      gdis_gtk4_window_present_report(self, "GTK4 Manual", report);
      gdis_gtk4_window_log(self, "Opened GTK4 usage guide.\n");
      return;
    }

  if (g_strcmp0(name, "animation") == 0)
    {
      gdis_gtk4_window_present_animation_tool(self);
      gdis_gtk4_window_log(self, "Opened animation tool.\n");
      return;
    }

  if (g_strcmp0(name, "periodic-table") == 0)
    {
      gdis_gtk4_window_present_periodic_table_tool(self);
      gdis_gtk4_window_log(self, "Opened periodic table.\n");
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
      gdis_gtk4_window_present_surface_tool(self);
      gdis_gtk4_window_log(self, "Opened surface builder.\n");
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

  if (g_strcmp0(name, "dislocations") == 0)
    {
      gdis_gtk4_window_present_dislocation_tool(self);
      gdis_gtk4_window_log(self, "Opened dislocation tool.\n");
      return;
    }

  if (g_strcmp0(name, "docking") == 0)
    {
      gdis_gtk4_window_present_docking_tool(self);
      gdis_gtk4_window_log(self, "Opened docking tool.\n");
      return;
    }

  if (g_strcmp0(name, "qbox") == 0)
    {
      gdis_gtk4_window_present_qbox_tool(self);
      gdis_gtk4_window_log(self, "Opened Qbox tool.\n");
      return;
    }

  if (g_strcmp0(name, "qbox-use-last-xml") == 0)
    {
      gdis_gtk4_window_present_qbox_tool(self);
      on_qbox_use_last_xml_clicked(NULL, self);
      gdis_gtk4_window_log(self, "Triggered Qbox: Use Last XML.\n");
      return;
    }

  if (g_strcmp0(name, "qbox-continue-last") == 0)
    {
      gdis_gtk4_window_present_qbox_tool(self);
      on_qbox_continue_clicked(NULL, self);
      gdis_gtk4_window_log(self, "Triggered Qbox: Continue Last.\n");
      return;
    }

  if (g_strcmp0(name, "qbox-import-result") == 0)
    {
      gdis_gtk4_window_present_qbox_tool(self);
      on_qbox_import_result_clicked(NULL, self);
      gdis_gtk4_window_log(self, "Triggered Qbox: Import Result.\n");
      return;
    }

  if (g_strcmp0(name, "qbox-results") == 0)
    {
      gdis_gtk4_window_present_qbox_tool(self);
      on_qbox_results_clicked(NULL, self);
      gdis_gtk4_window_log(self, "Triggered Qbox: Results.\n");
      return;
    }

  if (g_strcmp0(name, "building-dynamics") == 0 ||
      g_strcmp0(name, "analysis-dynamics") == 0)
    {
      gdis_gtk4_window_present_placeholder_feature(
        self,
        "Dynamics...",
        "The original dynamics builders and analysis dialogs are not ported yet in this GTK4 rebuild.");
      gdis_gtk4_window_log(self, "Opened dynamics status.\n");
      return;
    }

  if (g_strcmp0(name, "zmatrix") == 0)
    {
      gdis_gtk4_window_present_zmatrix_tool(self);
      gdis_gtk4_window_log(self, "Opened Z-matrix tool.\n");
      return;
    }

  if (g_strcmp0(name, "gulp") == 0 ||
      g_strcmp0(name, "gamess") == 0 ||
      g_strcmp0(name, "monty") == 0 ||
      g_strcmp0(name, "siesta") == 0 ||
      g_strcmp0(name, "vasp") == 0 ||
      g_strcmp0(name, "uspex") == 0)
    {
      gdis_gtk4_window_present_executable_paths_tool(self, name);
      gdis_gtk4_window_log(self, "Opened external computation setup for %s.\n", name);
      return;
    }

  if (g_strcmp0(name, "plots") == 0)
    {
      gdis_gtk4_window_present_placeholder_feature(
        self,
        "Plots...",
        "The original plotting dialogs are not restored yet in the GTK4 rebuild.");
      gdis_gtk4_window_log(self, "Opened plots status.\n");
      return;
    }

  if (g_strcmp0(name, "isosurface") == 0)
    {
      gdis_gtk4_window_present_isosurface_tool(self);
      gdis_gtk4_window_log(self, "Opened iso-surface tool.\n");
      return;
    }

  if (g_strcmp0(name, "render") == 0)
    {
      gdis_gtk4_window_present_display_tool(self);
      gdis_gtk4_window_log(self, "Opened display properties.\n");
      return;
    }

  if (g_strcmp0(name, "normal-mode") == 0)
    {
      gdis_gtk4_window_log(self, "Normal mode selected. The GTK4 rebuild already uses the standard viewer interaction mode.\n");
      return;
    }

  if (g_strcmp0(name, "recording-mode") == 0)
    {
      gdis_gtk4_window_present_recording_tool(self);
      gdis_gtk4_window_log(self, "Opened recording tool.\n");
      return;
    }

  if (g_strcmp0(name, "task-manager") == 0)
    {
      gdis_gtk4_window_present_task_manager_tool(self);
      gdis_gtk4_window_log(self, "Opened task manager.\n");
      return;
    }

  if (g_strcmp0(name, "executable-paths") == 0)
    {
      gdis_gtk4_window_present_executable_paths_tool(self, NULL);
      gdis_gtk4_window_log(self, "Opened executable paths.\n");
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

  if (g_strcmp0(name, "select-all") == 0)
    {
      gdis_gtk4_window_select_all_atoms(self);
      return;
    }

  if (g_strcmp0(name, "invert-selection") == 0)
    {
      gdis_gtk4_window_invert_selected_atoms(self);
      return;
    }

  if (g_strcmp0(name, "delete-selected") == 0)
    {
      if (!self->active_model)
        {
          gdis_gtk4_window_log(self, "Delete selected requested, but no active model is loaded.\n");
          return;
        }

      if (self->selected_atoms && self->selected_atoms->len > 0)
        on_delete_selected_group_clicked(NULL, self);
      else
        on_delete_selected_atom_clicked(NULL, self);
      return;
    }

  gdis_gtk4_window_log(self, "Activated action: %s\n", name);
}

static void
install_actions(GtkApplication *app)
{
  const GActionEntry entries[] = {
    {.name = "open", .activate = action_dispatch},
    {.name = "open-models-glob", .activate = action_dispatch},
    {.name = "open-models-folder", .activate = action_dispatch},
    {.name = "save", .activate = action_dispatch},
    {.name = "save-as", .activate = action_dispatch},
    {.name = "undo", .activate = action_dispatch},
    {.name = "close-model", .activate = action_dispatch},
    {.name = "quit", .activate = action_dispatch},
    {.name = "new-model", .activate = action_dispatch},
    {.name = "edit", .activate = action_dispatch},
    {.name = "select-all", .activate = action_dispatch},
    {.name = "invert-selection", .activate = action_dispatch},
    {.name = "delete-selected", .activate = action_dispatch},
    {.name = "render", .activate = action_dispatch},
    {.name = "animation", .activate = action_dispatch},
    {.name = "periodic-table", .activate = action_dispatch},
    {.name = "measure", .activate = action_dispatch},
    {.name = "isosurface", .activate = action_dispatch},
    {.name = "dislocations", .activate = action_dispatch},
    {.name = "docking", .activate = action_dispatch},
    {.name = "building-dynamics", .activate = action_dispatch},
    {.name = "zmatrix", .activate = action_dispatch},
    {.name = "surface", .activate = action_dispatch},
    {.name = "diffraction", .activate = action_dispatch},
    {.name = "gulp", .activate = action_dispatch},
    {.name = "gamess", .activate = action_dispatch},
    {.name = "monty", .activate = action_dispatch},
    {.name = "qbox", .activate = action_dispatch},
    {.name = "qbox-use-last-xml", .activate = action_dispatch},
    {.name = "qbox-continue-last", .activate = action_dispatch},
    {.name = "qbox-import-result", .activate = action_dispatch},
    {.name = "qbox-results", .activate = action_dispatch},
    {.name = "siesta", .activate = action_dispatch},
    {.name = "vasp", .activate = action_dispatch},
    {.name = "uspex", .activate = action_dispatch},
    {.name = "analysis-dynamics", .activate = action_dispatch},
    {.name = "plots", .activate = action_dispatch},
    {.name = "reset-images", .activate = action_dispatch},
    {.name = "reset-view", .activate = action_dispatch},
    {.name = "normal-mode", .activate = action_dispatch},
    {.name = "recording-mode", .activate = action_dispatch},
    {.name = "task-manager", .activate = action_dispatch},
    {.name = "executable-paths", .activate = action_dispatch},
    {.name = "about", .activate = action_dispatch},
    {.name = "manual", .activate = action_dispatch}
  };

  g_action_map_add_action_entries(G_ACTION_MAP(app),
                                  entries,
                                  G_N_ELEMENTS(entries),
                                  app);

  gtk_application_set_accels_for_action(app, "app.new-model",
                                        (const char *[]) {"<Primary>n", NULL});
  gtk_application_set_accels_for_action(app, "app.open",
                                        (const char *[]) {"<Primary>o", NULL});
  gtk_application_set_accels_for_action(app, "app.open-models-glob",
                                        (const char *[]) {"<Primary><Shift>o", NULL});
  gtk_application_set_accels_for_action(app, "app.open-models-folder",
                                        (const char *[]) {"<Primary><Alt>o", NULL});
  gtk_application_set_accels_for_action(app, "app.save",
                                        (const char *[]) {"<Primary>s", NULL});
  gtk_application_set_accels_for_action(app, "app.undo",
                                        (const char *[]) {"<Primary>z", NULL});
  gtk_application_set_accels_for_action(app, "app.close-model",
                                        (const char *[]) {"<Primary>w", NULL});
  gtk_application_set_accels_for_action(app, "app.quit",
                                        (const char *[]) {"<Primary>q", NULL});
  gtk_application_set_accels_for_action(app, "app.edit",
                                        (const char *[]) {"<Primary>e", NULL});
  gtk_application_set_accels_for_action(app, "app.select-all",
                                        (const char *[]) {"<Primary>a", NULL});
  gtk_application_set_accels_for_action(app, "app.invert-selection",
                                        (const char *[]) {"<Primary>i", NULL});
  gtk_application_set_accels_for_action(app, "app.render",
                                        (const char *[]) {"<Primary>d", NULL});
  gtk_application_set_accels_for_action(app, "app.reset-images",
                                        (const char *[]) {"<Primary>r", NULL});
  gtk_application_set_accels_for_action(app, "app.qbox-use-last-xml",
                                        (const char *[]) {"<Primary><Shift>u", NULL});
  gtk_application_set_accels_for_action(app, "app.qbox-continue-last",
                                        (const char *[]) {"<Primary><Shift>c", NULL});
  gtk_application_set_accels_for_action(app, "app.qbox-import-result",
                                        (const char *[]) {"<Primary><Shift>i", NULL});
  gtk_application_set_accels_for_action(app, "app.qbox-results",
                                        (const char *[]) {"<Primary><Shift>t", NULL});

  {
    GAction *undo_action;

    undo_action = g_action_map_lookup_action(G_ACTION_MAP(app), "undo");
    if (undo_action && G_IS_SIMPLE_ACTION(undo_action))
      g_simple_action_set_enabled(G_SIMPLE_ACTION(undo_action), FALSE);
  }
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
  self->measurement_records = g_hash_table_new_full(g_direct_hash,
                                                    g_direct_equal,
                                                    NULL,
                                                    gdis_measurement_record_array_free);
  self->undo_stacks = g_hash_table_new_full(g_direct_hash,
                                            g_direct_equal,
                                            NULL,
                                            (GDestroyNotify) g_ptr_array_unref);
  self->iso_surfaces = g_hash_table_new_full(g_direct_hash,
                                             g_direct_equal,
                                             NULL,
                                             (GDestroyNotify) gdis_isosurface_free);
  self->executable_paths = g_hash_table_new_full(g_str_hash,
                                                 g_str_equal,
                                                 g_free,
                                                 g_free);
  self->selected_atom_index = INVALID_ATOM_INDEX;
  self->selection_mode = GDIS_SELECTION_MODE_ATOMS;
  self->click_mode = GDIS_CLICK_MODE_SELECT;
  self->fragment_anchor_index = INVALID_ATOM_INDEX;
  self->drag_mode = GDIS_DRAG_MODE_NONE;
  gdis_gtk4_window_reset_view(self);

  install_actions(app);

  self->window = gtk_application_window_new(app);
  gtk_window_set_title(GTK_WINDOW(self->window), "GDIS GTK4");
  gtk_window_set_default_size(GTK_WINDOW(self->window),
                              GDIS_MAIN_WINDOW_WIDTH,
                              GDIS_MAIN_WINDOW_HEIGHT);
  gdis_gtk4_window_install_css(self->window);

  root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_widget_add_css_class(root, "gdis-root");
  gtk_widget_set_size_request(root, GDIS_MAIN_ROOT_WIDTH, GDIS_MAIN_ROOT_HEIGHT);
  gtk_window_set_child(GTK_WINDOW(self->window), root);

  menu_model = build_menu_bar_model();
#ifdef __APPLE__
  gtk_application_set_menubar(app, menu_model);
  gdis_macos_menu_install(app);
#endif
  menu_bar = gtk_popover_menu_bar_new_from_model(menu_model);
  g_object_unref(menu_model);

  gtk_box_append(GTK_BOX(root), menu_bar);
  gtk_box_append(GTK_BOX(root), build_toolbar(self));
  gtk_box_append(GTK_BOX(root), gtk_separator_new(GTK_ORIENTATION_HORIZONTAL));
  gtk_box_append(GTK_BOX(root), build_main_content(self));
  gdis_gtk4_window_restore_layout(self);

  g_object_set_data_full(G_OBJECT(self->window), WINDOW_DATA_KEY, self, gdis_gtk4_window_free);
  gdis_gtk4_window_update_undo_action(self);

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

  gtk_window_set_default_size(GTK_WINDOW(self->window),
                              GDIS_MAIN_WINDOW_WIDTH,
                              GDIS_MAIN_WINDOW_HEIGHT);
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
gdis_gtk4_window_activate_action(GdisGtk4Window *self, const char *action_name)
{
  g_return_if_fail(self != NULL);
  g_return_if_fail(self->app != NULL);
  g_return_if_fail(action_name != NULL);

  if (!action_name[0])
    return;

  g_action_group_activate_action(G_ACTION_GROUP(self->app), action_name, NULL);
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
