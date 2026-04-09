#!/bin/sh
set -eu

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
CONTENTS_DIR="$(dirname "$SCRIPT_DIR")"
APP_ROOT="$(dirname "$CONTENTS_DIR")"
RESOURCES_DIR="${CONTENTS_DIR}/Resources"
FRAMEWORKS_DIR="${CONTENTS_DIR}/Frameworks"
BIN_PATH="${SCRIPT_DIR}/gdis-gtk4-bin"
RUNTIME_DIR=""

if [ ! -x "$BIN_PATH" ]; then
  /usr/bin/osascript -e 'display alert "gdis-gtk4 binary not found in app bundle" message "Rebuild the app with: make -C gtk4_rebuild/gtk4_app dist-app"' >/dev/null 2>&1 || true
  exit 1
fi

cleanup()
{
  if [ -n "$RUNTIME_DIR" ] && [ -d "$RUNTIME_DIR" ]; then
    rm -rf "$RUNTIME_DIR"
  fi
}

trap cleanup EXIT HUP INT TERM

ensure_runtime_dir()
{
  if [ -z "$RUNTIME_DIR" ]; then
    RUNTIME_DIR="$(mktemp -d "${TMPDIR:-/tmp}/gdis-gtk4.XXXXXX")"
  fi
}

rewrite_template()
{
  src_path="$1"
  dest_path="$2"

  if [ ! -f "$src_path" ]; then
    return 0
  fi

  mkdir -p "$(dirname "$dest_path")"
  sed "s#@APP_ROOT@#${APP_ROOT}#g" "$src_path" > "$dest_path"
}

if [ -d "$FRAMEWORKS_DIR" ]; then
  export DYLD_FALLBACK_LIBRARY_PATH="${FRAMEWORKS_DIR}${DYLD_FALLBACK_LIBRARY_PATH:+:${DYLD_FALLBACK_LIBRARY_PATH}}"
fi

export PATH="${SCRIPT_DIR}${PATH:+:${PATH}}"
export XDG_DATA_DIRS="${RESOURCES_DIR}/share${XDG_DATA_DIRS:+:${XDG_DATA_DIRS}}"
export GSETTINGS_SCHEMA_DIR="${RESOURCES_DIR}/share/glib-2.0/schemas"

if [ -f "${RESOURCES_DIR}/etc/fonts/fonts.conf" ]; then
  export FONTCONFIG_PATH="${RESOURCES_DIR}/etc/fonts"
  export FONTCONFIG_FILE="${RESOURCES_DIR}/etc/fonts/fonts.conf"
  ensure_runtime_dir
  mkdir -p "${RUNTIME_DIR}/cache"
  export XDG_CACHE_HOME="${RUNTIME_DIR}/cache"
fi

if [ -d "${RESOURCES_DIR}/openmpi" ]; then
  export OMPI_MCA_component_path="${RESOURCES_DIR}/openmpi"
fi

if [ -d "${RESOURCES_DIR}/qbox-pseudos" ]; then
  export GDIS_GTK4_QBOX_BUNDLED_PSEUDO_DIR="${RESOURCES_DIR}/qbox-pseudos"
fi

if [ -x "${SCRIPT_DIR}/qbox" ]; then
  export GDIS_GTK4_QBOX_BUNDLED_EXEC="${SCRIPT_DIR}/qbox"
fi

if [ -f "${RESOURCES_DIR}/etc/gdk-pixbuf.loaders.in" ]; then
  ensure_runtime_dir
  rewrite_template "${RESOURCES_DIR}/etc/gdk-pixbuf.loaders.in" "${RUNTIME_DIR}/gdk-pixbuf.loaders"
  export GDK_PIXBUF_MODULE_FILE="${RUNTIME_DIR}/gdk-pixbuf.loaders"
fi

CURRENT_DIR="$(pwd -P 2>/dev/null || printf '%s' '/')"
if [ "$CURRENT_DIR" = "/" ] && [ -n "${HOME:-}" ] && [ -d "${HOME}" ]; then
  # Finder launches commonly start in "/", which makes file dialogs awkward.
  # Use the user's home directory as a safe default, but keep terminal cwd intact.
  cd "${HOME}"
fi

exec "$BIN_PATH" "$@"
