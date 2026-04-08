#!/bin/sh
set -eu

if [ "$#" -lt 1 ]; then
  echo "Usage: scripts/git-stage-selected.sh <path> [<path> ...]"
  echo "Example:"
  echo "  scripts/git-stage-selected.sh .gitignore gtk4_rebuild/gtk4_app"
  exit 1
fi

# Keep working tree changes, but clear the staging area first.
git restore --staged :/
git add -- "$@"

echo
echo "Staged paths:"
git diff --cached --name-only
