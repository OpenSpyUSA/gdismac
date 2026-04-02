#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import plistlib
import re
import shlex
import shutil
import stat
import subprocess
import sys
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = ROOT / "src"
BIN_DIR = ROOT / "bin"
DIST_DIR = ROOT / "dist"
BUILD_DIR = DIST_DIR / "build"
APP_NAME = "GDIS"
APP_BUNDLE = DIST_DIR / f"{APP_NAME}.app"
CONTENTS_DIR = APP_BUNDLE / "Contents"
MACOS_DIR = CONTENTS_DIR / "MacOS"
FRAMEWORKS_DIR = CONTENTS_DIR / "Frameworks"
RESOURCES_DIR = CONTENTS_DIR / "Resources"
ETC_DIR = RESOURCES_DIR / "etc"
SHARE_DIR = RESOURCES_DIR / "share"


def run(command: list[str], cwd: Path | None = None, capture: bool = False) -> str:
    print("+", shlex.join(command))
    if capture:
        return subprocess.check_output(command, cwd=cwd, text=True)
    subprocess.run(command, cwd=cwd, check=True)
    return ""


def pkg_config(variable: str, package: str) -> str:
    return run(["pkg-config", f"--variable={variable}", package], capture=True).strip()


def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise SystemExit(f"Missing required file: {path}")


def make_writable(path: Path) -> None:
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IWUSR)


def copy_file(source: Path, destination: Path) -> None:
    ensure_exists(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        make_writable(destination)
    shutil.copy2(source, destination)


def otool_dependencies(path: Path) -> list[str]:
    output = run(["otool", "-L", str(path)], capture=True)
    dependencies = []
    for line in output.splitlines()[1:]:
        line = line.strip()
        if not line:
            continue
        dependencies.append(line.split(" (", 1)[0])
    return dependencies


def should_bundle_dependency(path: str) -> bool:
    if path.startswith("@"):
        return False
    if path.startswith("/System/Library/"):
        return False
    if path.startswith("/usr/lib/"):
        return False
    return path.startswith("/")


def install_name_for(path: Path) -> str:
    return f"@executable_path/../Frameworks/{path.name}"


def bundle_dependency(source_name: str) -> Path:
    source = Path(source_name)
    ensure_exists(source)
    source = source.resolve()
    destination = FRAMEWORKS_DIR / source.name
    if not destination.exists():
        copy_file(source, destination)
        make_writable(destination)
    return destination


def rewrite_binary_links(path: Path) -> None:
    queue = [path]
    seen: set[Path] = set()

    while queue:
        current = queue.pop()
        current = current.resolve()
        if current in seen:
            continue
        seen.add(current)

        if FRAMEWORKS_DIR in current.parents:
            make_writable(current)
            run(["install_name_tool", "-id", install_name_for(current), str(current)])

        changes: list[str] = []
        for dependency in otool_dependencies(current):
            if not should_bundle_dependency(dependency):
                continue
            bundled = bundle_dependency(dependency)
            changes.extend(["-change", dependency, install_name_for(bundled)])
            queue.append(bundled)

        if changes:
            make_writable(current)
            run(["install_name_tool", *changes, str(current)])


def resource_relative_path(source: Path) -> Path:
    parts = source.parts
    if "lib" in parts:
        index = parts.index("lib")
        return Path("lib").joinpath(*parts[index + 1 :])
    if "share" in parts:
        index = parts.index("share")
        return Path("share").joinpath(*parts[index + 1 :])
    return Path(source.name)


def copy_resource_module(source_name: str) -> Path:
    source = Path(source_name)
    ensure_exists(source)
    source = source.resolve()
    destination = RESOURCES_DIR / resource_relative_path(source)
    copy_file(source, destination)
    return destination


def rewrite_runtime_config(text: str) -> str:
    def replacer(match: re.Match[str]) -> str:
        original = match.group(0)
        if not original.startswith("/"):
            return original
        original_path = Path(original)
        if "lib" in original_path.parts or "share" in original_path.parts:
            relative_path = resource_relative_path(original_path)
            return f"@APP_ROOT@/Contents/Resources/{relative_path.as_posix()}"
        return original

    return re.sub(r"/[^\" \n:]+", replacer, text)


def query_module_paths(text: str) -> list[str]:
    return re.findall(r'^"([^"]+)"\s*$', text, flags=re.MULTILINE)


def write_text(path: Path, content: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    if executable:
        mode = path.stat().st_mode
        path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def create_info_plist() -> None:
    info = {
        "CFBundleDevelopmentRegion": "English",
        "CFBundleDisplayName": APP_NAME,
        "CFBundleExecutable": APP_NAME,
        "CFBundleIconFile": "GDIS.icns",
        "CFBundleIdentifier": "org.sean.gdis",
        "CFBundleInfoDictionaryVersion": "6.0",
        "CFBundleName": APP_NAME,
        "CFBundlePackageType": "APPL",
        "CFBundleShortVersionString": "1.0",
        "CFBundleVersion": "1.0",
        "LSMinimumSystemVersion": "12.0",
        "NSPrincipalClass": "NSApplication",
        # GTK2 + gtkglext renders incorrectly on modern Retina displays unless
        # the app opts out of HiDPI. Without this, the GL scene ends up scaled
        # into the lower-left portion of the macOS window.
        "NSHighResolutionCapable": False,
    }
    with (CONTENTS_DIR / "Info.plist").open("wb") as handle:
        plistlib.dump(info, handle)


def create_launcher() -> None:
    source = BUILD_DIR / "gdis_launcher.c"
    launcher = textwrap.dedent(
        r"""
        #include <errno.h>
        #include <libgen.h>
        #include <limits.h>
        #include <mach-o/dyld.h>
        #include <stdbool.h>
        #include <stdio.h>
        #include <stdlib.h>
        #include <string.h>
        #include <unistd.h>

        static void die_message(const char *message)
        {
          fprintf(stderr, "%s\n", message);
          exit(EXIT_FAILURE);
        }

        static void die_errno(const char *message)
        {
          perror(message);
          exit(EXIT_FAILURE);
        }

        static char *xstrdup(const char *value)
        {
          char *copy = strdup(value);
          if (!copy)
            die_errno("strdup");
          return copy;
        }

        static char *dirname_copy(const char *path)
        {
          char *scratch, *dirpath, *copy;

          scratch = xstrdup(path);
          dirpath = dirname(scratch);
          copy = xstrdup(dirpath);
          free(scratch);

          return copy;
        }

        static char *join_path(const char *lhs, const char *rhs)
        {
          size_t lhs_len, rhs_len, extra;
          char *path;

          lhs_len = strlen(lhs);
          rhs_len = strlen(rhs);
          extra = (lhs_len > 0 && lhs[lhs_len - 1] != '/') ? 1 : 0;

          path = malloc(lhs_len + rhs_len + extra + 1);
          if (!path)
            die_errno("malloc");

          snprintf(path, lhs_len + rhs_len + extra + 1, "%s%s%s", lhs, extra ? "/" : "", rhs);

          return path;
        }

        static char *read_file(const char *path, size_t *length)
        {
          FILE *stream;
          long size;
          char *buffer;

          stream = fopen(path, "rb");
          if (!stream)
            die_errno(path);

          if (fseek(stream, 0, SEEK_END) != 0)
            die_errno(path);
          size = ftell(stream);
          if (size < 0)
            die_errno(path);
          if (fseek(stream, 0, SEEK_SET) != 0)
            die_errno(path);

          buffer = malloc((size_t) size + 1);
          if (!buffer)
            die_errno("malloc");

          if (fread(buffer, 1, (size_t) size, stream) != (size_t) size)
            die_errno(path);
          buffer[size] = '\0';

          if (fclose(stream) != 0)
            die_errno(path);

          *length = (size_t) size;
          return buffer;
        }

        static void write_file(const char *path, const char *buffer, size_t length)
        {
          FILE *stream;

          stream = fopen(path, "wb");
          if (!stream)
            die_errno(path);

          if (fwrite(buffer, 1, length, stream) != length)
            die_errno(path);

          if (fclose(stream) != 0)
            die_errno(path);
        }

        static char *replace_all(const char *input, const char *needle, const char *replacement)
        {
          const char *cursor, *match;
          size_t input_len, needle_len, replacement_len, matches, output_len;
          char *output, *writer;

          input_len = strlen(input);
          needle_len = strlen(needle);
          replacement_len = strlen(replacement);
          matches = 0;

          if (needle_len == 0)
            return xstrdup(input);

          cursor = input;
          while ((match = strstr(cursor, needle)))
            {
              matches++;
              cursor = match + needle_len;
            }

          output_len = input_len + matches * (replacement_len - needle_len);
          output = malloc(output_len + 1);
          if (!output)
            die_errno("malloc");

          cursor = input;
          writer = output;
          while ((match = strstr(cursor, needle)))
            {
              size_t segment_len = (size_t) (match - cursor);

              memcpy(writer, cursor, segment_len);
              writer += segment_len;
              memcpy(writer, replacement, replacement_len);
              writer += replacement_len;
              cursor = match + needle_len;
            }

          strcpy(writer, cursor);
          return output;
        }

        static void write_runtime_config(const char *template_path, const char *output_path, const char *app_root)
        {
          size_t length;
          char *input, *expanded;

          input = read_file(template_path, &length);
          (void) length;
          expanded = replace_all(input, "@APP_ROOT@", app_root);
          write_file(output_path, expanded, strlen(expanded));

          free(expanded);
          free(input);
        }

        static char *create_temp_dir(void)
        {
          const char *tmp_root;
          char *template;

          tmp_root = getenv("TMPDIR");
          if (!tmp_root || !tmp_root[0])
            tmp_root = "/tmp";

          template = join_path(tmp_root, "gdis-XXXXXX");
          if (!mkdtemp(template))
            die_errno("mkdtemp");

          return template;
        }

        int main(int argc, char *argv[])
        {
          uint32_t exec_path_size;
          char *exec_path;
          char resolved_exec_path[PATH_MAX];
          const char *resolved_exec;
          char *macos_dir, *contents_dir, *app_root;
          char *resources_dir, *etc_dir, *share_dir;
          char *gtk_template_dir, *gtk_template;
          char *pixbuf_template, *gtk_output, *pixbuf_output;
          char *runtime_dir, *gdis_bin;
          char **child_argv;
          const char *home;
          int i;

          exec_path_size = PATH_MAX;
          exec_path = malloc(exec_path_size);
          if (!exec_path)
            die_errno("malloc");

          if (_NSGetExecutablePath(exec_path, &exec_path_size) != 0)
            {
              free(exec_path);
              exec_path = malloc(exec_path_size);
              if (!exec_path)
                die_errno("malloc");
              if (_NSGetExecutablePath(exec_path, &exec_path_size) != 0)
                die_message("Unable to determine executable path.");
            }

          resolved_exec = realpath(exec_path, resolved_exec_path);
          if (!resolved_exec)
            resolved_exec = exec_path;

          macos_dir = dirname_copy(resolved_exec);
          contents_dir = dirname_copy(macos_dir);
          app_root = dirname_copy(contents_dir);

          resources_dir = join_path(contents_dir, "Resources");
          etc_dir = join_path(resources_dir, "etc");
          share_dir = join_path(resources_dir, "share");
          gtk_template_dir = join_path(etc_dir, "gtk-2.0");
          gtk_template = join_path(gtk_template_dir, "gtk.immodules.in");
          pixbuf_template = join_path(etc_dir, "gdk-pixbuf.loaders.in");
          runtime_dir = create_temp_dir();
          gtk_output = join_path(runtime_dir, "gtk.immodules");
          pixbuf_output = join_path(runtime_dir, "gdk-pixbuf.loaders");
          gdis_bin = join_path(macos_dir, "gdis-bin");

          write_runtime_config(gtk_template, gtk_output, app_root);
          write_runtime_config(pixbuf_template, pixbuf_output, app_root);

          home = getenv("HOME");
          if (home && home[0])
            {
              setenv("GDIS_START_DIR", home, 1);
              (void) chdir(home);
            }
          else
            {
              setenv("GDIS_START_DIR", app_root, 1);
            }

          setenv("GTK_DATA_PREFIX", resources_dir, 1);
          setenv("GTK_EXE_PREFIX", resources_dir, 1);
          setenv("GTK_IM_MODULE", "quartz", 1);
          setenv("GTK_IM_MODULE_FILE", gtk_output, 1);
          setenv("GDK_PIXBUF_MODULE_FILE", pixbuf_output, 1);
          setenv("XDG_DATA_DIRS", share_dir, 1);

          child_argv = calloc((size_t) argc + 1, sizeof(char *));
          if (!child_argv)
            die_errno("calloc");

          child_argv[0] = gdis_bin;
          for (i = 1; i < argc; i++)
            child_argv[i] = argv[i];

          execv(gdis_bin, child_argv);
          die_errno(gdis_bin);
          return EXIT_FAILURE;
        }
        """
    ).strip() + "\n"
    write_text(source, launcher)
    run(
        [
            "clang",
            "-O2",
            "-Wall",
            "-Wextra",
            "-mmacosx-version-min=12.0",
            str(source),
            "-o",
            str(MACOS_DIR / APP_NAME),
        ]
    )


def create_bundle_skeleton() -> None:
    binary_path = BIN_DIR / "gdis-bin"
    if not binary_path.exists():
        binary_path = BIN_DIR / "gdis"

    if APP_BUNDLE.exists():
        shutil.rmtree(APP_BUNDLE)
    MACOS_DIR.mkdir(parents=True, exist_ok=True)
    FRAMEWORKS_DIR.mkdir(parents=True, exist_ok=True)
    ETC_DIR.mkdir(parents=True, exist_ok=True)
    SHARE_DIR.mkdir(parents=True, exist_ok=True)

    create_info_plist()
    copy_file(SRC_DIR / "GDIS.icns", RESOURCES_DIR / "GDIS.icns")
    copy_file(binary_path, MACOS_DIR / APP_NAME)
    copy_file(BIN_DIR / "gdis.elements", MACOS_DIR / "gdis.elements")
    copy_file(BIN_DIR / "gdis.library", MACOS_DIR / "gdis.library")
    copy_file(BIN_DIR / "gdis.manual", MACOS_DIR / "gdis.manual")


def copy_optional_tree(source: Path, destination: Path) -> None:
    if not source.exists():
        return
    shutil.copytree(source, destination, dirs_exist_ok=True)


def bundle_runtime_support() -> None:
    gtk_prefix = Path(pkg_config("prefix", "gtk+-2.0"))

    copy_optional_tree(gtk_prefix / "share" / "locale", SHARE_DIR / "locale")
    copy_optional_tree(gtk_prefix / "share" / "themes", SHARE_DIR / "themes")

    gtk_query = run(["gtk-query-immodules-2.0"], capture=True)
    for module_path in query_module_paths(gtk_query):
        bundled_module = copy_resource_module(module_path)
        rewrite_binary_links(bundled_module)
    write_text(ETC_DIR / "gtk-2.0" / "gtk.immodules.in", rewrite_runtime_config(gtk_query))

    pixbuf_query = run(["gdk-pixbuf-query-loaders"], capture=True)
    for loader_path in query_module_paths(pixbuf_query):
        bundled_loader = copy_resource_module(loader_path)
        rewrite_binary_links(bundled_loader)
    write_text(ETC_DIR / "gdk-pixbuf.loaders.in", rewrite_runtime_config(pixbuf_query))


def iter_nested_code() -> list[Path]:
    nested = sorted(FRAMEWORKS_DIR.glob("*.dylib"))
    nested.extend(sorted(RESOURCES_DIR.rglob("*.so")))
    macos_files = sorted(path for path in MACOS_DIR.iterdir() if path.is_file())
    nested.extend(path for path in macos_files if path.name != APP_NAME)
    nested.extend(path for path in macos_files if path.name == APP_NAME)
    return nested


def ad_hoc_sign() -> None:
    for path in iter_nested_code():
        run(["codesign", "--force", "--sign", "-", str(path)])
        run(["codesign", "--verify", "--strict", str(path)])

    run(["codesign", "--force", "--sign", "-", str(APP_BUNDLE)])
    run(["codesign", "--verify", "--deep", "--strict", str(APP_BUNDLE)])


def create_zip_archive() -> None:
    archive = DIST_DIR / "GDIS-macos-arm64.zip"
    if archive.exists():
        archive.unlink()
    run(
        [
            "ditto",
            "-c",
            "-k",
            "--keepParent",
            "--sequesterRsrc",
            str(APP_BUNDLE),
            str(archive),
        ]
    )


def create_dmg_archive() -> None:
    staging_dir = DIST_DIR / "dmg-root"
    archive = DIST_DIR / "GDIS-macos-arm64.dmg"
    if archive.exists():
        archive.unlink()
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(APP_BUNDLE, staging_dir / APP_BUNDLE.name)
    os.symlink("/Applications", staging_dir / "Applications")
    try:
        run(
            [
                "hdiutil",
                "create",
                "-volname",
                APP_NAME,
                "-srcfolder",
                str(staging_dir),
                "-ov",
                "-format",
                "UDZO",
                str(archive),
            ]
        )
    finally:
        shutil.rmtree(staging_dir, ignore_errors=True)


def build_artifacts(skip_build: bool, create_dmg: bool, create_zip: bool) -> None:
    DIST_DIR.mkdir(parents=True, exist_ok=True)

    if not skip_build:
        run(["./install", "default"], cwd=ROOT)

    binary_path = BIN_DIR / "gdis-bin"
    if not binary_path.exists():
        binary_path = BIN_DIR / "gdis"

    ensure_exists(binary_path)
    ensure_exists(BIN_DIR / "gdis.elements")
    ensure_exists(BIN_DIR / "gdis.library")
    ensure_exists(BIN_DIR / "gdis.manual")

    create_bundle_skeleton()
    rewrite_binary_links(MACOS_DIR / APP_NAME)
    bundle_runtime_support()
    ad_hoc_sign()

    if create_zip:
        create_zip_archive()
    if create_dmg:
        create_dmg_archive()


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a portable macOS bundle for GDIS.")
    parser.add_argument("--skip-build", action="store_true", help="Reuse the existing binary in bin/.")
    parser.add_argument("--no-dmg", action="store_true", help="Skip DMG creation.")
    parser.add_argument("--no-zip", action="store_true", help="Skip ZIP creation.")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    build_artifacts(
        skip_build=args.skip_build,
        create_dmg=not args.no_dmg,
        create_zip=not args.no_zip,
    )
    print(f"Portable bundle created at {APP_BUNDLE}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
