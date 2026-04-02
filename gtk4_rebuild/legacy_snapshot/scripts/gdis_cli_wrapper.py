#!/usr/bin/env python3

from __future__ import annotations

import os
import sys


def script_dir() -> str:
    return os.path.dirname(os.path.abspath(__file__))


def repo_root() -> str:
    return os.path.dirname(script_dir())


def bundle_executable(app_bundle: str) -> str:
    return os.path.join(app_bundle, "Contents", "MacOS", "GDIS")


def normalize_args(argv: list[str]) -> list[str]:
    normalized: list[str] = []
    for arg in argv:
        if arg.startswith("-") or os.path.isabs(arg):
            normalized.append(arg)
        else:
            normalized.append(os.path.abspath(arg))
    return normalized


def main() -> int:
    root = repo_root()
    directory = script_dir()
    raw_binary = os.path.join(directory, "gdis-bin")
    if not os.path.exists(raw_binary):
        raw_binary = os.path.join(root, "bin", "gdis-bin")
    app_bundle = os.path.join(root, "dist", "GDIS.app")
    bundle_binary = bundle_executable(app_bundle)

    if len(sys.argv) > 1 and sys.argv[1] == "-c":
        os.execv(raw_binary, [raw_binary, *sys.argv[1:]])

    if os.path.exists(bundle_binary):
        args = normalize_args(sys.argv[1:])
        env = os.environ.copy()
        env.setdefault("GDIS_START_DIR", os.getcwd())
        os.execve(bundle_binary, [bundle_binary, *args], env)

    os.execv(raw_binary, [raw_binary, *sys.argv[1:]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
