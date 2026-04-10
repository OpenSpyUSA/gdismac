#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import stat
import subprocess
import sys
from pathlib import Path


APP_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = APP_DIR.parent.parent
BUILD_DIR = APP_DIR / "build"
MACOS_SRC_DIR = APP_DIR / "macos"
DIST_DIR = REPO_ROOT / "dist"

APP_NAME = "GDIS"
APP_BUNDLE = DIST_DIR / f"{APP_NAME}.app"
CONTENTS_DIR = APP_BUNDLE / "Contents"
MACOS_DIR = CONTENTS_DIR / "MacOS"
FRAMEWORKS_DIR = CONTENTS_DIR / "Frameworks"
RESOURCES_DIR = CONTENTS_DIR / "Resources"
ETC_DIR = RESOURCES_DIR / "etc"
SHARE_DIR = RESOURCES_DIR / "share"
LIB_DIR = RESOURCES_DIR / "lib"
OPENMPI_DIR = RESOURCES_DIR / "openmpi"
QBOX_PSEUDOS_DIR = RESOURCES_DIR / "qbox-pseudos"

APP_BINARY_NAME = "gdis-gtk4-bin"
APP_LAUNCHER_NAME = "gdis-gtk4"
APP_ZIP = DIST_DIR / "GDIS-macos-arm64.zip"
APP_DMG = DIST_DIR / "GDIS-macos-arm64.dmg"

RELEASE_DIR = DIST_DIR / "GDIS Portable"
RELEASE_ZIP = DIST_DIR / "GDIS-Portable-macos-arm64.zip"
RELEASE_DMG = DIST_DIR / "GDIS-Portable-macos-arm64.dmg"

STARTER_PSEUDOS = {
    "H": "H_ONCV_PBE-1.2.xml",
    "C": "C_ONCV_PBE-1.2.xml",
    "O": "O_ONCV_PBE-1.2.xml",
}
STARTER_PSEUDO_BASE_URL = "http://quantum-simulation.org/potentials/sg15_oncv/xml"

SAMPLE_FILES = [
    REPO_ROOT / "examples" / "water.xyz",
    REPO_ROOT / "examples" / "benzene.xyz",
    REPO_ROOT / "examples" / "rocksalt_demo.cif",
    REPO_ROOT / "examples" / "water_motion.ani",
    REPO_ROOT / "models" / "deoxy.pdb",
]


def run(command: list[str], cwd: Path | None = None, capture: bool = False) -> str:
    print("+", shlex.join(command))
    if capture:
        return subprocess.check_output(command, cwd=cwd, text=True)
    subprocess.run(command, cwd=cwd, check=True)
    return ""


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


def copy_tree(source: Path, destination: Path) -> None:
    if not source.exists():
        return
    shutil.copytree(source, destination, dirs_exist_ok=True)


def write_text(path: Path, content: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    if executable:
        mode = path.stat().st_mode
        path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def pkg_config(variable: str, package: str) -> str:
    return run(["pkg-config", f"--variable={variable}", package], capture=True).strip()


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
    source = Path(source_name).resolve()
    ensure_exists(source)
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
    if "etc" in parts:
        index = parts.index("etc")
        return Path("etc").joinpath(*parts[index + 1 :])
    return Path(source.name)


def copy_resource_module(source_name: str) -> Path:
    source = Path(source_name).resolve()
    ensure_exists(source)
    destination = RESOURCES_DIR / resource_relative_path(source)
    copy_file(source, destination)
    return destination


def rewrite_runtime_config(text: str) -> str:
    def replacer(match: re.Match[str]) -> str:
        original = match.group(0)
        if not original.startswith("/"):
            return original
        original_path = Path(original)
        parts = original_path.parts
        if "lib" in parts or "share" in parts or "etc" in parts:
            relative = resource_relative_path(original_path)
            return f"@APP_ROOT@/Contents/Resources/{relative.as_posix()}"
        return original

    return re.sub(r"/[^\" \n:]+", replacer, text)


def query_module_paths(text: str) -> list[str]:
    return re.findall(r'^"([^"]+)"\s*$', text, flags=re.MULTILINE)


def macho_files() -> list[Path]:
    files: list[Path] = []
    files.extend(sorted(FRAMEWORKS_DIR.glob("*.dylib")))
    files.extend(sorted(path for path in MACOS_DIR.iterdir() if path.is_file() and path.name != APP_LAUNCHER_NAME))
    files.extend(sorted(RESOURCES_DIR.rglob("*.so")))
    files.extend(sorted(RESOURCES_DIR.rglob("*.dylib")))
    return [path for path in files if not path.is_symlink()]


def ad_hoc_sign() -> None:
    for path in macho_files():
        run(["codesign", "--force", "--sign", "-", str(path)])
        run(["codesign", "--verify", "--strict", str(path)])

    run(["codesign", "--force", "--deep", "--sign", "-", str(APP_BUNDLE)])
    run(["codesign", "--verify", "--deep", "--strict", str(APP_BUNDLE)])


def find_qbox_binary() -> Path | None:
    candidates = [
        os.environ.get("GDIS_QBOX_BINARY"),
        str(Path.home() / "qbox-public" / "src" / "qb"),
        shutil.which("qbox"),
        shutil.which("qb"),
    ]
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser()
        if path.exists() and os.access(path, os.X_OK):
            return path.resolve()
    return None


def find_openmpi_component_dir() -> Path | None:
    try:
        output = run(["ompi_info", "--path", "pkglibdir"], capture=True).strip()
    except Exception:
        output = ""

    if ":" in output:
        candidate = output.split(":", 1)[1].strip()
        path = Path(candidate)
        if path.exists():
            return path.resolve()

    fallback = Path("/opt/homebrew/opt/open-mpi/lib/openmpi")
    if fallback.exists():
        return fallback.resolve()
    return None


def download_file(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    run(["curl", "-fL", url, "-o", str(destination)])


def find_local_starter_pseudo(filename: str) -> Path | None:
    candidates = [
        APP_DIR / "qbox_jobs" / "water" / filename,
        RELEASE_DIR / f"{APP_NAME}.app" / "Contents" / "Resources" / "qbox-pseudos" / filename,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def prepare_starter_pseudos() -> None:
    if QBOX_PSEUDOS_DIR.exists():
        shutil.rmtree(QBOX_PSEUDOS_DIR)
    QBOX_PSEUDOS_DIR.mkdir(parents=True, exist_ok=True)

    for symbol, filename in STARTER_PSEUDOS.items():
        destination = QBOX_PSEUDOS_DIR / filename
        local_source = find_local_starter_pseudo(filename)
        if local_source:
            copy_file(local_source, destination)
            continue
        download_file(f"{STARTER_PSEUDO_BASE_URL}/{filename}", destination)


def sanitize_fontconfig() -> None:
    fonts_conf = ETC_DIR / "fonts" / "fonts.conf"
    conf_readme = ETC_DIR / "fonts" / "conf.d" / "README"

    if fonts_conf.exists():
        text = fonts_conf.read_text(encoding="utf-8")
        text = re.sub(
            r"^\s*<cachedir>/opt/homebrew/var/cache/fontconfig</cachedir>\n?",
            "",
            text,
            flags=re.MULTILINE,
        )
        fonts_conf.write_text(text, encoding="utf-8")

    if conf_readme.exists():
        conf_readme.unlink()


def bundle_runtime_support(include_qbox: bool) -> None:
    homebrew_prefix = Path("/opt/homebrew")

    copy_tree(homebrew_prefix / "share" / "gtk-4.0", SHARE_DIR / "gtk-4.0")
    copy_tree(homebrew_prefix / "share" / "icons", SHARE_DIR / "icons")
    copy_tree(homebrew_prefix / "share" / "themes", SHARE_DIR / "themes")
    copy_tree(homebrew_prefix / "share" / "glib-2.0", SHARE_DIR / "glib-2.0")
    copy_tree(homebrew_prefix / "share" / "fontconfig", SHARE_DIR / "fontconfig")
    copy_tree(homebrew_prefix / "etc" / "fonts", ETC_DIR / "fonts")
    sanitize_fontconfig()

    schemas_dir = SHARE_DIR / "glib-2.0" / "schemas"
    if schemas_dir.exists():
        run(["glib-compile-schemas", str(schemas_dir)])

    pixbuf_query = run(["gdk-pixbuf-query-loaders"], capture=True)
    for loader_path in query_module_paths(pixbuf_query):
        bundled_loader = copy_resource_module(loader_path)
        rewrite_binary_links(bundled_loader)
    write_text(ETC_DIR / "gdk-pixbuf.loaders.in", rewrite_runtime_config(pixbuf_query))

    if include_qbox:
        component_dir = find_openmpi_component_dir()
        if component_dir:
            copy_tree(component_dir, OPENMPI_DIR)
            for module in sorted(OPENMPI_DIR.rglob("*.so")):
                rewrite_binary_links(module)


def create_bundle_skeleton(skip_qbox: bool) -> Path | None:
    app_binary = BUILD_DIR / "gdis-gtk4"
    ensure_exists(app_binary)

    qbox_binary = None if skip_qbox else find_qbox_binary()

    if APP_BUNDLE.exists():
        shutil.rmtree(APP_BUNDLE)

    MACOS_DIR.mkdir(parents=True, exist_ok=True)
    FRAMEWORKS_DIR.mkdir(parents=True, exist_ok=True)
    ETC_DIR.mkdir(parents=True, exist_ok=True)
    SHARE_DIR.mkdir(parents=True, exist_ok=True)
    LIB_DIR.mkdir(parents=True, exist_ok=True)

    copy_file(MACOS_SRC_DIR / "Info.plist", CONTENTS_DIR / "Info.plist")
    copy_file(REPO_ROOT / "src" / "GDIS.icns", RESOURCES_DIR / "GDIS.icns")
    copy_file(MACOS_SRC_DIR / "gdis-gtk4-launcher.sh", MACOS_DIR / APP_LAUNCHER_NAME)
    launcher_path = MACOS_DIR / APP_LAUNCHER_NAME
    launcher_path.chmod(launcher_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    copy_file(app_binary, MACOS_DIR / APP_BINARY_NAME)
    make_writable(MACOS_DIR / APP_BINARY_NAME)

    if qbox_binary:
        copy_file(qbox_binary, MACOS_DIR / "qbox")
        make_writable(MACOS_DIR / "qbox")

    return qbox_binary


def build_bundle(skip_build: bool, skip_qbox: bool) -> None:
    if not skip_build:
        run(["make", "-C", str(APP_DIR), "all"])

    qbox_binary = create_bundle_skeleton(skip_qbox)
    rewrite_binary_links(MACOS_DIR / APP_BINARY_NAME)
    if qbox_binary:
        rewrite_binary_links(MACOS_DIR / "qbox")
        prepare_starter_pseudos()
    bundle_runtime_support(include_qbox=qbox_binary is not None)
    ad_hoc_sign()


def build_sample_models_folder(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for sample in SAMPLE_FILES:
        copy_file(sample, destination / sample.name)

    examples_readme = REPO_ROOT / "examples" / "README.md"
    if examples_readme.exists():
        copy_file(examples_readme, destination / "README.md")


def write_release_docs(destination: Path, include_sample_models: bool = True) -> None:
    contains_lines = [
        "- GDIS.app",
        "- models",
    ]
    if include_sample_models:
        contains_lines.insert(1, "- Sample Models")
    contains_lines.extend(
        [
            "- README.txt",
            "- Qbox Setup Instructions.txt",
        ]
    )

    if include_sample_models:
        quick_start_line = '2. Open files from the "Sample Models" folder or the full "models" folder.'
        extra_readme = """
Suggested first files:
- water.xyz
- benzene.xyz
- rocksalt_demo.cif
- deoxy.pdb
- water_motion.ani

About the models folders:
- "Sample Models" contains a small curated set of good starter files.
- "models" contains the full repo models directory for broader testing.
"""
        qbox_model_line = "2. Load a model such as water.xyz or deoxy.pdb."
        qbox_intro = "- That means the included water.xyz, benzene.xyz, and deoxy.pdb samples are the safest first Qbox tests."
    else:
        quick_start_line = '2. Open files from the "models" folder.'
        extra_readme = """
Note:
- This ZIP intentionally omits the separate "Sample Models" folder.
- Use the full "models" folder included here for testing and sharing.
"""
        qbox_model_line = "2. Load a model from the bundled models folder, such as deoxy.pdb."
        qbox_intro = "- That means models using only H, C, and O are the safest first Qbox tests."

    readme = f"""GDIS Portable for macOS

Thank you for trying the modern GTK4 macOS build of GDIS.

This package contains:
{chr(10).join(contains_lines)}

Recommended first launch:
1. Open GDIS.app directly from this folder, or move it to Applications first.
{quick_start_line}
3. If macOS blocks the app on first launch, right-click GDIS.app and choose Open.

Package notes:
- This build is intended for Apple Silicon (arm64) Macs.
- The app is ad-hoc signed for local distribution and is not notarized.
- The included app bundle is self-contained and intended to be easier to run than the thin development build.
- Qbox is bundled for local single-process use inside the app.
{extra_readme}
"""

    qbox_notes = f"""Qbox Setup Instructions

What is already included:
- GDIS.app contains a bundled Qbox executable.
- Starter pseudopotentials for H, C, and O are included.
{qbox_intro}

Quick Qbox workflow:
1. Open GDIS.app.
{qbox_model_line}
3. Open Tools > Computation > Qbox...
4. Click Detect if you want to confirm the bundled Qbox path.
5. Click Regenerate to prepare the inputs, then click Run.

If your model uses other elements:
- Use "Setup Local Pseudos" in the Qbox window.
- GDIS may need internet access the first time so it can download missing pseudo XML files.

Current scope:
- This package is aimed at local single-process Qbox use through the bundled executable.
- More advanced MPI launch setups may still need machine-specific tuning.
"""

    write_text(destination / "README.txt", readme)
    write_text(destination / "Qbox Setup Instructions.txt", qbox_notes)


def create_release_folder() -> None:
    if RELEASE_DIR.exists():
        shutil.rmtree(RELEASE_DIR)
    RELEASE_DIR.mkdir(parents=True, exist_ok=True)

    shutil.copytree(APP_BUNDLE, RELEASE_DIR / APP_BUNDLE.name)
    build_sample_models_folder(RELEASE_DIR / "Sample Models")
    copy_tree(REPO_ROOT / "models", RELEASE_DIR / "models")
    write_release_docs(RELEASE_DIR)


def create_release_zip_staging() -> Path:
    staging_root = DIST_DIR / "zip-root"
    staged_release = staging_root / RELEASE_DIR.name

    if staging_root.exists():
        shutil.rmtree(staging_root)

    shutil.copytree(RELEASE_DIR, staged_release)
    shutil.rmtree(staged_release / "Sample Models", ignore_errors=True)
    write_release_docs(staged_release, include_sample_models=False)
    return staged_release


def create_zip_archive(source: Path, archive: Path) -> None:
    if archive.exists():
        archive.unlink()
    run(["ditto", "-c", "-k", "--keepParent", "--sequesterRsrc", str(source), str(archive)])


def create_dmg_archive(source: Path, archive: Path, volume_name: str) -> None:
    staging_dir = DIST_DIR / "dmg-root"
    if archive.exists():
        archive.unlink()
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, staging_dir / source.name)
    os.symlink("/Applications", staging_dir / "Applications")
    try:
        run(
            [
                "hdiutil",
                "create",
                "-volname",
                volume_name,
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


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a portable macOS GTK4 GDIS release.")
    parser.add_argument("--skip-build", action="store_true", help="Reuse the existing build/gdis-gtk4 binary.")
    parser.add_argument("--skip-qbox", action="store_true", help="Do not bundle Qbox into the app.")
    parser.add_argument("--no-zip", action="store_true", help="Skip ZIP archive creation.")
    parser.add_argument("--no-dmg", action="store_true", help="Skip DMG creation.")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    DIST_DIR.mkdir(parents=True, exist_ok=True)
    build_bundle(skip_build=args.skip_build, skip_qbox=args.skip_qbox)
    create_release_folder()

    create_zip_archive(APP_BUNDLE, APP_ZIP)
    create_dmg_archive(APP_BUNDLE, APP_DMG, APP_NAME)

    if not args.no_zip:
        staged_release = create_release_zip_staging()
        try:
            create_zip_archive(staged_release, RELEASE_ZIP)
        finally:
            shutil.rmtree(staged_release.parent, ignore_errors=True)
    if not args.no_dmg:
        create_dmg_archive(RELEASE_DIR, RELEASE_DMG, "GDIS Portable")

    print(f"Portable app bundle created at {APP_BUNDLE}")
    print(f"Portable release folder created at {RELEASE_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
