#!/bin/sh
set -eu

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
BIN_PATH="${SCRIPT_DIR}/gdis-gtk4-bin"

if [ ! -x "$BIN_PATH" ]; then
  /usr/bin/osascript -e 'display alert "gdis-gtk4 binary not found in app bundle" message "Rebuild the app with: make -C gtk4_rebuild/gtk4_app dist-app"' >/dev/null 2>&1 || true
  exit 1
fi

# If this app sits inside the repository tree, switch to the repo root so
# legacy relative sample paths (examples/, models/) resolve consistently.
SEARCH_DIR="$SCRIPT_DIR"
while [ "$SEARCH_DIR" != "/" ]; do
  if [ -d "$SEARCH_DIR/examples" ] && [ -d "$SEARCH_DIR/models" ]; then
    cd "$SEARCH_DIR"
    break
  fi
  SEARCH_DIR="$(dirname "$SEARCH_DIR")"
done

exec "$BIN_PATH" "$@"
