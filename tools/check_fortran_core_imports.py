#!/usr/bin/env python3
"""Check runtime Fortran files for forbidden CATChem core ownership imports."""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path


LOGGER = logging.getLogger("check_fortran_core_imports")

FORBIDDEN_MODULES = {
    "chemspeciesutils_mod",
    "chemstate_mod",
    "configmanager_mod",
    "constants",
    "diagnosticinterface_mod",
    "diagnosticmanager_mod",
    "catchem_emis_data_mod",
    "gridgeometry_mod",
    "metstate_mod",
    "precision_mod",
    "species_mod",
    "statemanager_mod",
    "timestate_mod",
    "unitconversion_mod",
}

CURRENTLY_ALLOWED_IMPORTS = {
    "src/api/CATChem_API.F90": {
        "precision_mod",
    },
    "drivers/nuopc/catchem_nuopc_interface.F90": {
        "constants",
        "catchem_emis_data_mod",
        "precision_mod",
    },
    "drivers/nuopc/catchem_emis_mod.F90": {
        "constants",
        "catchem_emis_data_mod",
        "precision_mod",
    },
}

USE_RE = re.compile(r"^\s*use\s*,?\s*(?:non_intrinsic\s*::\s*)?([a-z0-9_]+)", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed checker options.
    """
    parser = argparse.ArgumentParser(
        description="Fail when runtime Fortran files import forbidden core ownership modules."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Repository root to scan.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat the current migration allowlist as forbidden too.",
    )
    parser.add_argument(
        "paths",
        nargs="*",
        default=("src/api", "drivers/nuopc"),
        help="Runtime files or directories to scan, relative to --root.",
    )
    return parser.parse_args()


def iter_fortran_files(root: Path, paths: list[str]) -> list[Path]:
    """Collect Fortran files under the requested runtime paths.

    Parameters
    ----------
    root
        Repository root directory.
    paths
        Files or directories to scan, relative to ``root`` unless absolute.

    Returns
    -------
    list[Path]
        Sorted Fortran source files.
    """
    files: set[Path] = set()
    for raw_path in paths:
        path = Path(raw_path)
        target = path if path.is_absolute() else root / path
        if target.is_file() and target.suffix.lower() in {".f90", ".f", ".for", ".ftn"}:
            files.add(target)
        elif target.is_dir():
            for suffix in ("*.F90", "*.f90", "*.F", "*.f", "*.for", "*.ftn"):
                files.update(target.rglob(suffix))
    return sorted(files)


def find_forbidden_imports(root: Path, files: list[Path], strict: bool) -> list[str]:
    """Find forbidden imports in runtime Fortran files.

    Parameters
    ----------
    root
        Repository root directory.
    files
        Fortran files to scan.
    strict
        If true, ignore the temporary migration allowlist.

    Returns
    -------
    list[str]
        Human-readable violation messages.
    """
    violations: list[str] = []
    for file_path in files:
        relative = file_path.relative_to(root).as_posix()
        allowed = set() if strict else CURRENTLY_ALLOWED_IMPORTS.get(relative, set())
        for line_number, line in enumerate(file_path.read_text(encoding="utf-8").splitlines(), start=1):
            match = USE_RE.match(line)
            if match is None:
                continue
            module_name = match.group(1).lower()
            if module_name in FORBIDDEN_MODULES and module_name not in allowed:
                violations.append(f"{relative}:{line_number}: forbidden import of {module_name}")
    return violations


def main() -> int:
    """Run the static import check.

    Returns
    -------
    int
        Zero when no violations are found, non-zero otherwise.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    root = args.root.resolve()
    files = iter_fortran_files(root, list(args.paths))
    violations = find_forbidden_imports(root, files, args.strict)
    if violations:
        LOGGER.error("FATAL ERROR: forbidden CATChem Fortran core imports detected")
        for violation in violations:
            LOGGER.error(violation)
        return 1
    LOGGER.info("checked %d Fortran runtime files for forbidden core imports", len(files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
