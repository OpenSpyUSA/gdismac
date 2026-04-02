/*
Copyright (C) 2003 by Sean David Fleming

sean@ivec.org

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

The GNU GPL can also be found at http://www.gnu.org
*/

/* irix */
#define _BSD_SIGNALS 1

#include <signal.h>
#include <stdio.h>
#include <locale.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#ifdef __APPLE__
#include <limits.h>
#include <mach-o/dyld.h>
#include "mac_activation.h"
#endif
#ifndef __WIN32
#include <sys/wait.h>
#endif

#include "gdis.h"
#include "command.h"
#include "file.h"
#include "parse.h"
#include "task.h"
#include "host.h"
#include "render.h"
#include "matrix.h"
#include "model.h"
#include "numeric.h"
#include "module.h"
#include "library.h"
#include "interface.h"
#include "gui_image.h"

/* main data structures */
struct sysenv_pak sysenv;
struct elem_pak elements[MAX_ELEMENTS];

#ifdef __APPLE__
/***********************************/
/* canonicalize a CLI path at launch */
/***********************************/
static gchar *mac_cli_canonicalize(const gchar *path)
{
gchar *cwd, *full;

if (!path || !*path)
  return(NULL);

cwd = g_get_current_dir();
full = g_canonicalize_filename(path, cwd);
g_free(cwd);

return(full);
}

/******************************************/
/* retain terminal cwd when launching mac app */
/******************************************/
static void mac_capture_launch_directory(gint argc, char *argv[])
{
gint i;
gchar *cwd;

if (g_getenv("GDIS_START_DIR"))
  return;

for (i=1 ; i<argc ; i++)
  {
  if (!argv[i] || !*argv[i])
    continue;
  if (argv[i][0] == '-')
    continue;

  cwd = g_get_current_dir();
  if (cwd && *cwd)
    g_setenv("GDIS_START_DIR", cwd, TRUE);
  g_free(cwd);
  return;
  }
}

/*************************************/
/* normalize argv for bundle relaunch */
/*************************************/
static gchar **mac_cli_normalize_argv(gint argc,
                                      char *argv[],
                                      const gchar *program_path)
{
gint i;
gchar **args;

args = g_new0(gchar *, argc + 1);
args[0] = g_strdup(program_path);

for (i=1 ; i<argc ; i++)
  {
  if (!argv[i])
    continue;

  if (argv[i][0] == '-' || g_path_is_absolute(argv[i]))
    args[i] = g_strdup(argv[i]);
  else
    args[i] = mac_cli_canonicalize(argv[i]);
  }

return(args);
}

/*********************************************/
/* use the app bundle even when raw bin is run */
/*********************************************/
static void mac_relaunch_bundle_if_needed(gint argc, char *argv[])
{
uint32_t exec_path_size;
gchar *exec_path, *resolved_exec;
gchar *exec_dir, *repo_root, *bundle_exec;
gchar **bundle_argv;

if (g_getenv("GDIS_RAW_RELAUNCHED"))
  return;
if (argc > 1 && argv[1] && g_ascii_strncasecmp(argv[1], "-c", 3) == 0)
  return;

exec_path_size = PATH_MAX;
exec_path = g_malloc(exec_path_size);
if (_NSGetExecutablePath(exec_path, &exec_path_size) != 0)
  {
  g_free(exec_path);
  exec_path = g_malloc(exec_path_size);
  if (_NSGetExecutablePath(exec_path, &exec_path_size) != 0)
    {
    g_free(exec_path);
    return;
    }
  }

resolved_exec = g_canonicalize_filename(exec_path, NULL);
g_free(exec_path);
if (!resolved_exec)
  return;

if (g_strrstr(resolved_exec, ".app/Contents/MacOS/"))
  {
  g_free(resolved_exec);
  return;
  }

exec_dir = g_path_get_dirname(resolved_exec);
repo_root = g_path_get_dirname(exec_dir);
bundle_exec = g_build_filename(repo_root, "dist", "GDIS.app",
                               "Contents", "MacOS", "GDIS", NULL);

if (g_file_test(bundle_exec, G_FILE_TEST_IS_EXECUTABLE))
  {
  mac_capture_launch_directory(argc, argv);
  g_setenv("GDIS_RAW_RELAUNCHED", "1", TRUE);
  bundle_argv = mac_cli_normalize_argv(argc, argv, bundle_exec);
  execv(bundle_exec, bundle_argv);
  g_strfreev(bundle_argv);
  }

g_free(bundle_exec);
g_free(repo_root);
g_free(exec_dir);
g_free(resolved_exec);
}
#endif

#ifdef WITH_GUI
/***********************************/
/* deferred command-line startup state */
/***********************************/
struct gui_cli_pak
{
GSList *paths;
gchar *directory;
gint last_width;
gint last_height;
gint stable_polls;
gboolean windows_presented;
};

/***********************************/
/* canonicalize a CLI path once at launch */
/***********************************/
static gchar *gui_cli_canonicalize(const gchar *path)
{
gchar *cwd, *full;

if (!path || !*path)
  return(NULL);

cwd = g_get_current_dir();
full = g_canonicalize_filename(path, cwd);
g_free(cwd);

return(full);
}

/***********************************/
/* does a path look like a loadable file? */
/***********************************/
static gboolean gui_cli_path_supported(const gchar *path)
{
if (!path || !*path)
  return(FALSE);

return(get_file_info((gpointer) path, BY_EXTENSION) != NULL);
}

/***********************************************/
/* tell the user when they passed a directory */
/***********************************************/
static void gui_cli_append_directory(struct gui_cli_pak *startup,
                                     const gchar *dirname)
{
gchar *msg;

g_return_if_fail(startup != NULL);

if (!startup->directory)
  startup->directory = gui_cli_canonicalize(dirname);

msg = g_strdup_printf("%s is a directory.\nOpening the file browser there instead.\nUse a real file path such as ./models/deoxy.pdb or ./models/* to load files directly.\n",
                      dirname);
gui_text_show(ERROR, msg);
g_printerr("%s", msg);
g_free(msg);
}

/*********************************************/
/* queue a command-line file or expand folder */
/*********************************************/
static void gui_cli_append_argument(struct gui_cli_pak *startup,
                                    const gchar *arg)
{
gchar *full, *msg;

g_return_if_fail(startup != NULL);

if (!arg || !*arg)
  return;
if (g_ascii_strncasecmp(arg, "-c", 3) == 0)
  return;

full = gui_cli_canonicalize(arg);
if (!full)
  return;

if (g_file_test(full, G_FILE_TEST_IS_DIR))
  {
  gui_cli_append_directory(startup, full);
  g_free(full);
  return;
  }

if (gui_cli_path_supported(full))
  {
  startup->paths = g_slist_append(startup->paths, full);
  return;
  }

msg = g_strdup_printf("Skipping unsupported input: %s\n", arg);
gui_text_show(ERROR, msg);
g_printerr("%s", msg);
g_free(msg);
g_free(full);
}

/******************************************/
/* free deferred command-line startup data */
/******************************************/
static void gui_cli_free(struct gui_cli_pak *startup)
{
GSList *list;

if (!startup)
  return;

for (list=startup->paths ; list ; list=g_slist_next(list))
  g_free(list->data);
g_slist_free(startup->paths);
g_free(startup->directory);
g_free(startup);
}

/******************************************/
/* load queued command-line files in the UI */
/******************************************/
static gboolean gui_cli_ui_ready(struct gui_cli_pak *startup)
{
#ifdef __APPLE__
gint width, height;

if (!startup)
  return(FALSE);
if (!sysenv.glarea || !GTK_IS_WIDGET(sysenv.glarea))
  return(FALSE);
if (!sysenv.display_box || !GTK_IS_WIDGET(sysenv.display_box))
  return(FALSE);
if (!GTK_WIDGET_REALIZED(sysenv.glarea) || !GTK_WIDGET_MAPPED(sysenv.glarea))
  return(FALSE);
if (sysenv.tpane && GTK_WIDGET_VISIBLE(sysenv.tpane))
  return(FALSE);
width = sysenv.glarea->allocation.width;
height = sysenv.glarea->allocation.height;
if (width < 200 || height < 200)
  return(FALSE);
if (sysenv.display_box->allocation.width != width ||
    sysenv.display_box->allocation.height != height)
  return(FALSE);
if (sysenv.width != width || sysenv.height != height)
  return(FALSE);
if (startup->last_width != width || startup->last_height != height)
  {
  startup->last_width = width;
  startup->last_height = height;
  startup->stable_polls = 0;
  return(FALSE);
  }
if (startup->stable_polls < 2)
  {
  startup->stable_polls++;
  return(FALSE);
  }
#endif

return(TRUE);
}

static gboolean gui_cli_post_load_sync(gpointer data)
{
(void) data;

if (sysenv.active_model)
  gui_view_default();
if (sysenv.glarea)
  {
  gtk_widget_queue_resize(sysenv.glarea);
  gtk_widget_queue_draw(sysenv.glarea);
  }

return(FALSE);
}

static gboolean gui_cli_load_files(gpointer data)
{
GSList *list;
struct gui_cli_pak *startup = data;

if (!startup)
  return(FALSE);

#ifdef __APPLE__
if (!startup->windows_presented)
  {
  startup->windows_presented = TRUE;
  gui_macos_present_windows();
  return(TRUE);
  }
if (!gui_cli_ui_ready(startup))
  {
  return(TRUE);
  }
#endif

if (startup->directory && !startup->paths)
  {
  g_free(sysenv.cwd);
  sysenv.cwd = g_strdup(startup->directory);
  g_setenv("GDIS_START_DIR", startup->directory, TRUE);
  (void) chdir(startup->directory);
  file_load_dialog();
  gui_cli_free(startup);
  return(FALSE);
  }

for (list=startup->paths ; list ; list=g_slist_next(list))
  file_load(list->data, NULL);

#ifdef __APPLE__
gui_cli_post_load_sync(NULL);
g_timeout_add(150, gui_cli_post_load_sync, NULL);
#endif

gui_cli_free(startup);

return(FALSE);
}

/***********************************/
/* defer CLI file loads until ready */
/***********************************/
static void gui_cli_schedule_files(gint argc, char *argv[])
{
gint i;
struct gui_cli_pak *startup;

startup = g_malloc0(sizeof(struct gui_cli_pak));

for (i=1 ; i<argc ; i++)
  gui_cli_append_argument(startup, argv[i]);

if (!startup->paths && !startup->directory)
  {
  gui_cli_free(startup);
  return;
  }

#ifdef __APPLE__
g_timeout_add(150, gui_cli_load_files, startup);
#else
g_idle_add(gui_cli_load_files, startup);
#endif
}
#endif

#ifdef __APPLE__
/************************************/
/* runtime config token replacement */
/************************************/
static gchar *mac_replace_token(const gchar *input, const gchar *token, const gchar *value)
{
gchar **parts;
gchar *result;

g_assert(input != NULL);
g_assert(token != NULL);
g_assert(value != NULL);

parts = g_strsplit(input, token, -1);
result = g_strjoinv(value, parts);
g_strfreev(parts);

return(result);
}

/*************************************/
/* write expanded runtime config file */
/*************************************/
static void mac_write_runtime_config(const gchar *template_path,
                                     const gchar *output_path,
                                     const gchar *app_root)
{
gchar *input, *expanded;
gsize length;

if (!g_file_get_contents(template_path, &input, &length, NULL))
  return;

expanded = mac_replace_token(input, "@APP_ROOT@", app_root);
g_file_set_contents(output_path, expanded, -1, NULL);

g_free(expanded);
g_free(input);
}

/**********************************/
/* configure a bundled mac launch */
/**********************************/
static void mac_bundle_bootstrap(void)
{
uint32_t exec_path_size;
gchar *exec_path, *resolved_exec;
gchar *macos_dir, *contents_dir, *app_root;
gchar *resources_dir, *etc_dir, *share_dir;
gchar *gtk_template, *pixbuf_template;
gchar *runtime_template, *runtime_dir;
gchar *gtk_output, *pixbuf_output;
const gchar *home;

exec_path_size = PATH_MAX;
exec_path = g_malloc(exec_path_size);
if (_NSGetExecutablePath(exec_path, &exec_path_size) != 0)
  {
  g_free(exec_path);
  exec_path = g_malloc(exec_path_size);
  if (_NSGetExecutablePath(exec_path, &exec_path_size) != 0)
    {
    g_free(exec_path);
    return;
    }
  }

resolved_exec = g_canonicalize_filename(exec_path, NULL);
g_free(exec_path);
if (resolved_exec)
  g_setenv("GDIS_EXEC_PATH", resolved_exec, TRUE);

macos_dir = g_path_get_dirname(resolved_exec);
contents_dir = g_path_get_dirname(macos_dir);
app_root = g_path_get_dirname(contents_dir);

resources_dir = g_build_filename(contents_dir, "Resources", NULL);
etc_dir = g_build_filename(resources_dir, "etc", NULL);
share_dir = g_build_filename(resources_dir, "share", NULL);
gtk_template = g_build_filename(etc_dir, "gtk-2.0", "gtk.immodules.in", NULL);
pixbuf_template = g_build_filename(etc_dir, "gdk-pixbuf.loaders.in", NULL);

if (g_file_test(gtk_template, G_FILE_TEST_EXISTS)
    && g_file_test(pixbuf_template, G_FILE_TEST_EXISTS))
  {
  runtime_template = g_build_filename(g_get_tmp_dir(), "gdis-XXXXXX", NULL);
  runtime_dir = g_mkdtemp(runtime_template);
  if (runtime_dir)
    {
    gtk_output = g_build_filename(runtime_dir, "gtk.immodules", NULL);
    pixbuf_output = g_build_filename(runtime_dir, "gdk-pixbuf.loaders", NULL);

    mac_write_runtime_config(gtk_template, gtk_output, app_root);
    mac_write_runtime_config(pixbuf_template, pixbuf_output, app_root);

    home = g_getenv("GDIS_START_DIR");
    if (!home || !*home)
      home = g_get_home_dir();

    if (home && *home)
      {
      g_setenv("GDIS_START_DIR", home, TRUE);
      (void) chdir(home);
      }
    else
      {
      g_setenv("GDIS_START_DIR", app_root, TRUE);
      }

    g_setenv("GTK_DATA_PREFIX", resources_dir, TRUE);
    g_setenv("GTK_EXE_PREFIX", resources_dir, TRUE);
    g_setenv("GTK_IM_MODULE", "quartz", TRUE);
    g_setenv("GTK_IM_MODULE_FILE", gtk_output, TRUE);
    g_setenv("GDK_PIXBUF_MODULE_FILE", pixbuf_output, TRUE);
    g_setenv("XDG_DATA_DIRS", share_dir, TRUE);

    g_free(gtk_output);
    g_free(pixbuf_output);
    }
  g_free(runtime_template);
  }

g_set_prgname("GDIS");
g_set_application_name("GDIS");

g_free(pixbuf_template);
g_free(gtk_template);
g_free(share_dir);
g_free(etc_dir);
g_free(resources_dir);
g_free(app_root);
g_free(contents_dir);
g_free(macos_dir);
g_free(resolved_exec);
}

#endif

/********************************/
/* sfc hash table value cleanup */
/********************************/
void free_sfc(gpointer sfc_list)
{
free_slist((GSList *) sfc_list);
}

/******************/
/* system cleanup */
/******************/
void sys_free(void)
{

file_free();

free_slist(sysenv.render.light_list);

g_hash_table_destroy(sysenv.sfc_table);
g_hash_table_destroy(sysenv.library);
g_hash_table_destroy(sysenv.manual);
g_hash_table_destroy(sysenv.image_table);
module_free();

#ifndef _WIN32
host_free_all();
#endif
}

/**********************************/
/* read the init file if possible */
/**********************************/
gint read_gdisrc(void)
{
gint num_tokens;
#ifdef UNUSED_BUT_SET
gint i, element_flag;
#endif
gchar line[LINELEN], **buff;
gdouble version;
FILE *fp;

/* attempt to open */
fp = fopen(sysenv.init_file, "r");

/* check for an old/bad gdisrc */
if (fp)
  {
/* scan the 1st line */
  buff = get_tokenized_line(fp, &num_tokens);
  if (g_ascii_strncasecmp(*buff,"gdis",4) == 0)
    {
/* test for right version */
    if (num_tokens < 2)
      return(1);
    version = str_to_float(*(buff+1));
    sscanf(*(buff+1),"%lf",&version);
/* 0.75 was when the new element scheme was implemented */
    if (version < 0.75)
      return(1);
    }
  else
    return(1);
  g_strfreev(buff);
  }
else
  return(1);

/* read until EOF, or (failsafe) reached elements[] array allocation */
#ifdef UNUSED_BUT_SET
element_flag=0;
i=0;
#endif
for (;;)
  {
/* NB: we need line[] to get the OpenGL font */
  if (fgetline(fp, line))
    break;
  buff = tokenize(line, &num_tokens);
  if (!buff)
    break;

/* decide what to read */
  if (g_ascii_strncasecmp("size",*buff,4) == 0)
    {
    if (num_tokens > 2)
      {
      sysenv.width = str_to_float(*(buff+1));
      sysenv.height = str_to_float(*(buff+2));
      if (sysenv.width < MIN_WIDTH)
        sysenv.width = MIN_WIDTH;
      if (sysenv.height < MIN_HEIGHT)
        sysenv.height = MIN_HEIGHT;
#ifdef __APPLE__
      /*
       * Older broken macOS builds persisted a 200x200 display pane, which
       * leaves the GL canvas effectively unusable on restart. Treat that
       * legacy minimum-sized value as invalid and fall back to sane defaults.
       */
      if (sysenv.width == MIN_WIDTH && sysenv.height == MIN_HEIGHT)
        {
        sysenv.width = START_WIDTH;
        sysenv.height = START_HEIGHT;
        }
#endif
      }
    } 

/* model pane width */
  if (g_ascii_strncasecmp("pane",*buff,4) == 0)
    {
    if (num_tokens > 2)
      {
      sysenv.tree_width = str_to_float(*(buff+1));
      sysenv.tray_height = str_to_float(*(buff+2));
      }
    }

/* model tree/property divider */
  if (g_ascii_strncasecmp("divide",*buff,6) == 0)
    {
    if (num_tokens > 1)
      sysenv.tree_divider = (gint) str_to_float(*(buff+1));
    }

/* povray executable */
  if (g_ascii_strncasecmp("povray_p", *buff, 8) == 0)
    sysenv.povray_path = parse_strip_newline(&line[12]);

  if (g_ascii_strncasecmp("povray_e", *buff, 8) == 0)
    {
    if (num_tokens > 1)
      {
      g_free(sysenv.povray_exe);
      sysenv.povray_exe = g_strdup(*(buff+1));
      }
    }

/* animation creation tool */
  if (g_ascii_strncasecmp("convert_p", *buff, 9) == 0)
    sysenv.convert_path = parse_strip_newline(&line[13]);

  if (g_ascii_strncasecmp("convert_e", *buff, 9) == 0)
    {
    if (num_tokens > 1)
      {
      g_free(sysenv.convert_exe);
      sysenv.convert_exe = g_strdup(*(buff+1));
      }
    }

/* image viewing */
  if (g_ascii_strncasecmp("viewer_p", *buff, 8) == 0)
    sysenv.viewer_path = parse_strip_newline(&line[12]);

  if (g_ascii_strncasecmp("viewer_e", *buff, 8) == 0)
    {
    if (num_tokens > 1)
      {
      g_free(sysenv.viewer_exe);
      sysenv.viewer_exe = g_strdup(*(buff+1));
      }
    }

/* GULP */
  if (g_ascii_strncasecmp("gulp_p", *buff, 6) == 0)
    sysenv.gulp_path = parse_strip_newline(&line[10]);

  if (g_ascii_strncasecmp("gulp_e", *buff, 6) == 0)
    {
    if (num_tokens > 1)
      {
      g_free(sysenv.gulp_exe);
      sysenv.gulp_exe = g_strdup(*(buff+1));
      }
    }
/* OpenGL drawing font */
  if (g_ascii_strncasecmp("gl_font",*buff,7) == 0)
    if (num_tokens > 1)
      strcpy(sysenv.gl_fontname, g_strstrip(&line[8]));

/* model tree box */
  if (g_ascii_strncasecmp("mtb",*buff,3) == 0)
    if (num_tokens > 1)
      sysenv.mtb_on = (gint) str_to_float(*(buff+1));

/* model properties box */
  if (g_ascii_strncasecmp("mpb",*buff,3) == 0)
    if (num_tokens > 1)
      sysenv.mpb_on = (gint) str_to_float(*(buff+1));

/* model symmetry box */
  if (g_ascii_strncasecmp("msb",*buff,3) == 0)
    if (num_tokens > 1)
      sysenv.msb_on = (gint) str_to_float(*(buff+1));

/* atom properties box */
  if (g_ascii_strncasecmp("apb",*buff,3) == 0)
    if (num_tokens > 1)
      sysenv.apb_on = (gint) str_to_float(*(buff+1));

/* halo type */
  if (g_ascii_strncasecmp("halo",*buff,4) == 0)
    if (num_tokens > 1)
      sysenv.render.halos = (gint) str_to_float(*(buff+1));

/* low quality rotation */
  if (g_ascii_strncasecmp("fast",*buff,4) == 0)
    if (num_tokens > 1)
      sysenv.render.fast_rotation = (gint) str_to_float(*(buff+1));

  if (g_ascii_strncasecmp(*buff, "colour_bg", 9) == 0)
    {
    if (num_tokens > 3)
      {
      sysenv.render.bg_colour[0] = str_to_float(*(buff+1));
      sysenv.render.bg_colour[1] = str_to_float(*(buff+2));
      sysenv.render.bg_colour[2] = str_to_float(*(buff+3));
      }
    }
  if (g_ascii_strncasecmp(*buff, "colour_morph", 11) == 0)
    {
    if (num_tokens > 3)
      {
      sysenv.render.morph_colour[0] = str_to_float(*(buff+1));
      sysenv.render.morph_colour[1] = str_to_float(*(buff+2));
      sysenv.render.morph_colour[2] = str_to_float(*(buff+3));
      }
    }
  if (g_ascii_strncasecmp("colour_rsurf", *buff, 12) == 0)
    {
    if (num_tokens > 3)
      {
      sysenv.render.rsurf_colour[0] = str_to_float(*(buff+1));
      sysenv.render.rsurf_colour[1] = str_to_float(*(buff+2));
      sysenv.render.rsurf_colour[2] = str_to_float(*(buff+3));
      }
    }

/* cleanup */
  g_strfreev(buff);
  }

/* parse for elements */
rewind(fp);
read_elem_data(fp, MODIFIED);

return(0);
}

/*********************************************/
/* write setup & elements to the gdisrc file */
/*********************************************/
gint write_gdisrc(void)
{
FILE *fp;

fp = fopen(sysenv.init_file,"w");
if (!fp)
  {
  printf("Error: unable to create %s\n", sysenv.init_file);
  return(1);
  }

/* save final dimensions */
if (sysenv.mpane)
  if (GTK_IS_WIDGET(sysenv.mpane))
    sysenv.tree_width = GTK_WIDGET(sysenv.mpane)->allocation.width;
if (sysenv.tpane)
  if (GTK_IS_WIDGET(sysenv.tpane))
    sysenv.tray_height = GTK_WIDGET(sysenv.tpane)->allocation.height;

fprintf(fp,"gdis %f\n", VERSION);
fprintf(fp,"canvas %d\n", sysenv.canvas);
fprintf(fp,"size %d %d\n", sysenv.width,sysenv.height);
fprintf(fp,"pane %d %d\n", sysenv.tree_width, sysenv.tray_height);
fprintf(fp,"divider %d\n", sysenv.tree_divider);
fprintf(fp,"gl_font %s\n", sysenv.gl_fontname);
fprintf(fp,"mtb %d\n", sysenv.mtb_on);
fprintf(fp,"mpb %d\n", sysenv.mpb_on);
fprintf(fp,"msb %d\n", sysenv.msb_on);
fprintf(fp,"apb %d\n", sysenv.apb_on);
fprintf(fp,"halos %d\n", sysenv.render.halos);
fprintf(fp,"fast %d\n", sysenv.render.fast_rotation);

fprintf(fp,"colour_bg %f %f %f\n", sysenv.render.bg_colour[0],
                                   sysenv.render.bg_colour[1],
                                   sysenv.render.bg_colour[2]);
fprintf(fp,"colour_morph %f %f %f\n", sysenv.render.morph_colour[0],
                                      sysenv.render.morph_colour[1],
                                      sysenv.render.morph_colour[2]);
fprintf(fp,"colour_rsurf %f %f %f\n", sysenv.render.rsurf_colour[0],
                                      sysenv.render.rsurf_colour[1],
                                      sysenv.render.rsurf_colour[2]);

if (sysenv.babel_path)
  fprintf(fp,"babel_path %s\n", sysenv.babel_path);
if (sysenv.convert_path)
  fprintf(fp,"convert_path %s\n", sysenv.convert_path);
if (sysenv.gulp_path)
  fprintf(fp,"gulp_path %s\n", sysenv.gulp_path);
if (sysenv.gamess_path)
  fprintf(fp,"gamess_path %s\n", sysenv.gamess_path);
if (sysenv.vasp_path)
  fprintf(fp,"vasp_path %s\n", sysenv.vasp_path);
if (sysenv.uspex_path)
  fprintf(fp,"uspex_path %s\n", sysenv.uspex_path);
if (sysenv.mpirun_path)
  fprintf(fp,"mpirun_path %s\n", sysenv.mpirun_path);
if (sysenv.povray_path)
  fprintf(fp,"povray_path %s\n", sysenv.povray_path);
if (sysenv.viewer_path)
  fprintf(fp,"viewer_path %s\n", sysenv.viewer_path);


fprintf(fp,"babel_exe %s\n", sysenv.babel_exe);
fprintf(fp,"convert_exe %s\n", sysenv.convert_exe);
fprintf(fp,"gulp_exe %s\n", sysenv.gulp_exe);
fprintf(fp,"gamess_exe %s\n", sysenv.gamess_exe);
fprintf(fp,"vasp_exe %s\n", sysenv.vasp_exe);
fprintf(fp,"uspex_exe %s\n", sysenv.uspex_exe);
fprintf(fp,"mpirun_exe %s\n", sysenv.mpirun_exe);
fprintf(fp,"povray_exe %s\n", sysenv.povray_exe);
fprintf(fp,"viewer_exe %s\n", sysenv.viewer_exe);



/* write the non-default element data */
write_elem_data(fp);

fclose(fp);
return(0);
}

/**************/
/* main setup */
/**************/
#define DEBUG_SYS_INIT 0
void sys_init(gint argc, gchar *argv[])
{
gchar *temp;
const gchar *ctemp;
struct light_pak *light;
FILE *fp;

/* top level structure initialization */
sysenv.max_threads = -1;
sysenv.task_list = NULL;
sysenv.host_list = NULL;
sysenv.dialog_list = NULL;
sysenv.glconfig = NULL;
sysenv.stereo = FALSE;
sysenv.stereo_windowed = FALSE;
sysenv.stereo_fullscreen = FALSE;
sysenv.canvas_list = NULL;
sysenv.canvas_rows = 1;
sysenv.canvas_cols = 1;
sysenv.width = START_WIDTH;
sysenv.height = START_HEIGHT;
sysenv.snapshot = FALSE;
sysenv.snapshot_filename = NULL;
sysenv.tree_width = START_WIDTH/4;
sysenv.tree_divider = -1;
sysenv.tray_height = 60;
sysenv.mpane = NULL;
sysenv.tpane = NULL;
sysenv.fps = 40;
sysenv.moving = FALSE;
sysenv.roving = FALSE;
sysenv.refresh_canvas = FALSE;
sysenv.refresh_dialog = FALSE;
sysenv.refresh_tree = FALSE;
sysenv.refresh_properties = FALSE;
sysenv.refresh_text = FALSE;
sysenv.mtb_on = TRUE;
sysenv.mpb_on = TRUE;
sysenv.msb_on = TRUE;
sysenv.apb_on = TRUE;
sysenv.lmb_on = FALSE;
sysenv.pib_on = FALSE;

/*edit_dialog*/
sysenv.cedit.edit_anim_n=0;
sysenv.cedit.edit_chirality[0]=6.;
sysenv.cedit.edit_chirality[1]=6.;
sysenv.cedit.edit_length=1.44;
sysenv.cedit.edit_basis[0]=NULL;
sysenv.cedit.edit_basis[1]=NULL;
sysenv.cedit.edit_spatial_colour[0]=1.0;
sysenv.cedit.edit_spatial_colour[1]=0.0;
sysenv.cedit.edit_spatial_colour[2]=0.0;
sysenv.cedit.edit_label_colour[0]=1.0;
sysenv.cedit.edit_label_colour[1]=1.0;
sysenv.cedit.edit_label_colour[2]=1.0;
sysenv.cedit.spatial_list=NULL;
sysenv.cedit.spatial_tree=NULL;
sysenv.cedit.spatial_selected=NULL;
sysenv.cedit.diffract_model=NULL;
/*edit_widget*/
sysenv.cedit.apd_data=NULL;
sysenv.cedit.apd_core=NULL;

/* default to single model display */
sysenv.mal = NULL;
sysenv.displayed[0] = -1;
sysenv.active_model = NULL;
sysenv.canvas = TRUE;
sysenv.select_mode = CORE;
sysenv.calc_pbonds = TRUE;

/* cairo --OVHPA*/
sysenv.cairo_surface=NULL;

/* default masks */
sysenv.file_type = DATA;
sysenv.babel_type = AUTO;
sysenv.num_elements = sizeof(elements) / sizeof(struct elem_pak);
sysenv.elements = NULL;

/* rendering setup */
sysenv.render.width = 600;
sysenv.render.height = 600;
sysenv.render.vp_dist = 5.0;
sysenv.render.zone_size = 10.0;
sysenv.render.show_energy = FALSE;

/* stereo defaults */
sysenv.render.stereo_quadbuffer = FALSE;
sysenv.render.stereo_use_frustum = TRUE;
sysenv.render.stereo_eye_offset = 1.0;
sysenv.render.stereo_parallax = 1.0;
sysenv.render.stereo_left = TRUE;
sysenv.render.stereo_right = TRUE;

/* default light */
light = g_malloc(sizeof(struct light_pak));
light->type = DIRECTIONAL;
VEC3SET(light->x, -1.0, 1.0, -1.0);
VEC3SET(light->colour, 1.0, 1.0, 1.0);
light->ambient = 0.2;
light->diffuse = 0.8;
light->specular = 0.7;
/*** JJM new values for better (IMHO) lighting ***/
light->diffuse = 0.5;
//light->specular = 0.3;
sysenv.render.light_list = g_slist_append(sysenv.render.light_list, light);
sysenv.render.num_lights = 1;

sysenv.render.perspective = FALSE;
sysenv.render.antialias = FALSE;
sysenv.render.wire_show_hidden = FALSE;
sysenv.render.fog = FALSE;
sysenv.render.wire_model = FALSE;
sysenv.render.wire_surface = FALSE;
sysenv.render.shadowless = FALSE;
sysenv.render.animate = FALSE;
sysenv.render.no_povray_exec = FALSE;
sysenv.render.no_keep_tempfiles = TRUE;
sysenv.render.animate_type = ANIM_GIF;
sysenv.render.animate_file = g_strdup("movie");
sysenv.render.atype = FALSE;
sysenv.render.axes = FALSE;
sysenv.render.morph_finish = g_strdup("F_Glass4");
sysenv.render.ref_index = 2.5;
sysenv.render.delay = 20.0;
sysenv.render.mpeg_quality = 95.0;
sysenv.render.sphere_quality = 3.0;
sysenv.render.cylinder_quality = 9.0;
sysenv.render.ribbon_quality = 10.0;
sysenv.render.ms_grid_size = 0.5;
sysenv.render.auto_quality = FALSE;
sysenv.render.fast_rotation = FALSE;
sysenv.render.halos = FALSE;
sysenv.render.ambience = 0.2;
sysenv.render.diffuse = 0.9;
sysenv.render.specular = 0.9;
sysenv.render.transmit = 1.0;
sysenv.render.ghost_opacity = 0.5;
sysenv.render.ball_radius = 0.3;
sysenv.render.stick_radius = 0.1;
sysenv.render.stick_thickness = GTKGL_LINE_WIDTH;
sysenv.render.line_thickness = GTKGL_LINE_WIDTH;
sysenv.render.frame_thickness = GTKGL_LINE_WIDTH;
sysenv.render.geom_line_width = 3.0;
sysenv.render.cpk_scale = 1.0;
sysenv.render.fog_density = 0.35;
sysenv.render.fog_start = 0.5;
sysenv.render.ribbon_curvature = 0.5;
sysenv.render.ribbon_thickness = 1.0;
sysenv.render.phonon_scaling = 4.0;
sysenv.render.ahl_strength = 0.7;
sysenv.render.ahl_size = 20;
sysenv.render.shl_strength = 0.8;
sysenv.render.shl_size = 20;
VEC3SET(sysenv.render.fg_colour, 1.0, 1.0, 1.0);
VEC3SET(sysenv.render.bg_colour, 0.0, 0.0, 0.0);
VEC3SET(sysenv.render.morph_colour, 0.1, 0.1, 0.8);
VEC3SET(sysenv.render.rsurf_colour, 0.0, 0.3, 0.6);
VEC3SET(sysenv.render.label_colour, 1.0, 1.0, 0.0);
VEC3SET(sysenv.render.title_colour, 0.0, 1.0, 1.0);
VEC3SET(sysenv.render.ribbon_colour, 0.0, 0.0, 1.0);

/* setup directory and file pointers */
sysenv.cwd = g_get_current_dir();
const gchar *envdir = g_getenv("GDIS_START_DIR");
if (envdir)
  {
  g_free(sysenv.cwd);
  sysenv.cwd = g_strdup(envdir);
  }

#if DEBUG_SYS_INIT
printf("cwd: [%s]\n", sysenv.cwd);
#endif

/* generate element file full pathname */
/* sometimes this returns the program name, and sometimes it doesn't */
#ifdef __APPLE__
{
const gchar *resolved_exec = g_getenv("GDIS_EXEC_PATH");
uint32_t exec_path_size = PATH_MAX;
gchar *exec_path;

if (resolved_exec && *resolved_exec)
  {
  temp = g_strdup(resolved_exec);
  }
else
  {
  exec_path = g_malloc(exec_path_size);
  if (_NSGetExecutablePath(exec_path, &exec_path_size) != 0)
    {
    g_free(exec_path);
    exec_path = g_malloc(exec_path_size);
    if (_NSGetExecutablePath(exec_path, &exec_path_size) != 0)
      {
      g_free(exec_path);
      exec_path = NULL;
      }
    }

  if (exec_path)
    {
    temp = g_canonicalize_filename(exec_path, NULL);
    g_free(exec_path);
    }
  else
    temp = NULL;
  }
}
#else
temp = g_find_program_in_path(argv[0]);
if (!temp && argv[0] && *argv[0])
  {
  if (g_path_is_absolute(argv[0]))
    temp = g_strdup(argv[0]);
  else
    temp = g_canonicalize_filename(argv[0], NULL);
  }
#endif
/* remove program name (if attached) */
if (temp && g_file_test(temp, G_FILE_TEST_IS_DIR))
  sysenv.gdis_path = temp;
else if (temp)
  {
  sysenv.gdis_path = g_path_get_dirname(temp);
  g_free(temp);
  }

#if DEBUG_SYS_INIT
printf("%s path: [%s]\n", argv[0], sysenv.gdis_path);
#endif

if (sysenv.gdis_path)
  sysenv.elem_file = g_build_filename(sysenv.gdis_path, ELEM_FILE, NULL);
else
  {
  printf("WARNING: gdis directory not found.\n");
  sysenv.elem_file = g_build_filename(sysenv.cwd, ELEM_FILE, NULL);
  }

/* generate gdisrc full pathname */
ctemp = g_get_home_dir();
if (ctemp)
  sysenv.init_file = g_build_filename(ctemp, INIT_FILE, NULL);
else
  {
  printf("WARNING: home directory not found.\n");
  if (sysenv.gdis_path)
    sysenv.init_file = g_build_filename(sysenv.gdis_path, INIT_FILE, NULL);
  else
    {
    if (sysenv.cwd)
      sysenv.init_file = g_build_filename(sysenv.cwd, INIT_FILE, NULL);
    else
      {
      printf("FATAL: current directory not found.\n");
      exit(-1);
      }
    }
  }

/* defaults */
#if _WIN32
sysenv.babel_exe = g_strdup("babel.exe");
sysenv.povray_exe = g_strdup("povray.exe");
sysenv.convert_exe = g_strdup("convert.exe");
sysenv.viewer_exe = g_strdup("display.exe");
sysenv.gulp_exe = g_strdup("gulp.exe");
sysenv.gamess_exe = g_strdup("wingamess");
#else
sysenv.babel_exe = g_strdup("babel");
sysenv.povray_exe = g_strdup("povray");
sysenv.convert_exe = g_strdup("convert");
sysenv.viewer_exe = g_strdup("display");
sysenv.gulp_exe = g_strdup("gulp");
sysenv.gamess_exe = g_strdup("run_gms_for_gdis");
sysenv.vasp_exe = g_strdup("vasp");
sysenv.uspex_exe = g_strdup("USPEX");
sysenv.mpirun_exe = g_strdup("mpirun");
#endif


strcpy(sysenv.gl_fontname, GL_FONT);

/* atomic scattering factor coefficients */
sysenv.sfc_table = g_hash_table_new_full(&g_str_hash,
                                         &hash_strcmp,
                                         &g_free,
                                         &free_sfc);

/* IMPORTANT this must be done before gdisrc is parsed as */
/* setup the element data before any possible modifications */
/* can be read in from the gdisrc file */
printf("scanning: %-50s ", sysenv.elem_file);
fp = fopen(sysenv.elem_file, "rt");
if (fp)
  {
  read_elem_data(fp, DEFAULT);
  fclose(fp);
  printf("[ok]\n");
  }
else
  {
/* missing default element info is fatal */
  printf("[not found]\n");
  exit(1);
  }

/* if failed to read gdisrc - write a new one */
/* TODO - if overwrite old version, save/rename first? */
sysenv.babel_path = NULL;
sysenv.convert_path = NULL;
sysenv.gulp_path = NULL;
sysenv.gamess_path = NULL;
sysenv.vasp_path = NULL;
sysenv.uspex_path = NULL;
sysenv.mpirun_path = NULL;
sysenv.povray_path = NULL;
sysenv.viewer_path = NULL;
if (read_gdisrc())
  {
  printf("creating: %-50s ", sysenv.init_file);
  if (write_gdisrc())
    printf("[error]\n");
  else
    printf("[ok]\n");
  }
else
  printf("scanning: %-50s [ok]\n", sysenv.init_file);

/* get executable paths */
if (!sysenv.babel_path)
  sysenv.babel_path = g_find_program_in_path(sysenv.babel_exe);
if (!sysenv.convert_path)
  sysenv.convert_path = g_find_program_in_path(sysenv.convert_exe);
if (!sysenv.gulp_path)
  sysenv.gulp_path = g_find_program_in_path(sysenv.gulp_exe);
if (!sysenv.vasp_path)
  sysenv.vasp_path = g_find_program_in_path(sysenv.vasp_exe);
if (!sysenv.uspex_path)
  sysenv.uspex_path = g_find_program_in_path(sysenv.uspex_exe);
if (!sysenv.mpirun_path)
  sysenv.mpirun_path = g_find_program_in_path(sysenv.mpirun_exe);
if (!sysenv.povray_path)
  sysenv.povray_path = g_find_program_in_path(sysenv.povray_exe);

/* display program */
if (!sysenv.viewer_path)
  sysenv.viewer_path = g_find_program_in_path(sysenv.viewer_exe);

/* alternative (FIXME - starting to get ugly - redo) */
if (!sysenv.viewer_path)
  sysenv.viewer_path = g_find_program_in_path("eog");

/* HACK FIX - GAMESS path contains the path, not the path + file */
temp = g_find_program_in_path(sysenv.gamess_exe);

if (temp)
  sysenv.gamess_path = g_path_get_dirname(temp);
else
#if __WIN32
  sysenv.gamess_path = g_strdup("c:\\wingamess");
#else
  sysenv.gamess_path = NULL;
#endif

g_free(temp);

/* tree init */
sysenv.num_trees = 0;
/* copied selection */
sysenv.select_source = NULL;
sysenv.module_list = NULL;
sysenv.projects = NULL;
sysenv.manual = NULL;
sysenv.image_table = NULL;
sysenv.surfaces = NULL;

/* check for eps format support by pixbuf */
#ifdef CAIRO_HAS_PS_SURFACE
sysenv.have_eps=TRUE;
#else
/* check for eps format support by pixbuf <- very unlikely: remove?*/
sysenv.have_eps=FALSE;
GSList *formats = gdk_pixbuf_get_formats ();
GSList *format=formats;
for(;format;format=g_slist_next(format)){
	GdkPixbufFormat *this_format=format->data;
	if (strcmp("eps", gdk_pixbuf_format_get_name(this_format)) == 0){
		if(gdk_pixbuf_format_is_writable(this_format)) sysenv.have_eps=TRUE;
	}
}
g_slist_free(formats);
#endif

/* module */
file_init();
#ifdef WITH_GUI
help_init();
#endif
init_math();
library_init();
}

/********/
/* MAIN */
/********/
gint main(int argc, char *argv[])
{
/* NBL read this 1st as it affects canvas type, but command */
/* line arguments should be able to overide */
#ifdef __APPLE__
mac_capture_launch_directory(argc, argv);
mac_relaunch_bundle_if_needed(argc, argv);
mac_bundle_bootstrap();
#endif
sys_init(argc, argv);
module_setup();
sysenv.write_gdisrc = FALSE;

/* command line? */
if (argc > 1)
  if (g_ascii_strncasecmp(argv[1],"-c",3) == 0)
    sysenv.canvas = FALSE;


#ifdef WITH_GUI
/* set up main window and event handlers (or start non-gui mode) */
if (sysenv.canvas)
  {
  gint i;

/* initialize threads and set up the queue */
  task_queue_init();

/* CURRENT - absorb these into gui_init() */
/*
  gtk_init(&argc, &argv);
  gtk_gl_init(&argc, &argv);
*/

/* main interface */
  gui_init(argc, argv);
#ifdef __APPLE__
  mac_activate_foreground_app_now();
  g_idle_add(mac_activate_foreground_app, NULL);
#endif
/* internationalization: set locale to C to avoid _BUG_*/
setlocale(LC_ALL, "C");

/* task update timer */
/* TODO - put this in gui_widget_handler? */
  g_timeout_add(1000, (GSourceFunc) &update_task_info, NULL);

/* screen redraw timer */
  g_timeout_add(25, (GSourceFunc) &gui_canvas_handler, NULL);

/* CURRENT - gui updates done through a timer */
/* this was done so that threads can safely schedule a gui update */
  g_timeout_add(200, (GSourceFunc) &gui_widget_handler, NULL);

/* process arguments as files to load */
#ifdef WITH_GUI
  gui_cli_schedule_files(argc, argv);
#else
  for (i=1 ; i<argc ; i++)
    file_load(argv[i], NULL);
#endif

/* thread-safe handling */
  gdk_threads_enter();
  gtk_main();
  gdk_threads_leave();
  }
else
#else
sysenv.canvas = FALSE;
#endif
  command_main_loop(argc, argv);

return(0);
}

/* CURRENT */
/* routines that are not cleanly separable when we build with no GUI */
#ifndef WITH_GUI
void gui_text_show(gint type, gchar *msg)
{
printf("%s", msg);
}

void gui_refresh(gint dummy)
{
}

void tree_select_delete(void)
{
}
void tree_select_active(void)
{
}
void tree_select_model(struct model_pak *m)
{
}
void tree_model_add(struct model_pak *m)
{
}
void tree_model_refresh(struct model_pak *m)
{
}

void dialog_destroy_type(gint dummy)
{
}

gpointer graph_new(const gchar *dummy, struct model_pak *m)
{
return(NULL);
}

void graph_add_data(gint a, gdouble *b, gdouble c, gdouble d, gpointer e)
{
}

void graph_set_yticks(gint a, gint b, gpointer c)
{
}

void graph_free_list(struct model_pak *m)
{
}

void meas_graft_model(struct model_pak *m)
{
}

void meas_prune_model(struct model_pak *m)
{
}

/* FIXME - this should actually be replaced by gui_refresh() */
void redraw_canvas(gint dummy)
{
}
#endif
