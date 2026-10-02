#!/usr/bin/env python3
"""Verify Dockerfile base image argument compliance."""

import argparse
import logging
import os
import re
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

BUILD_ARG = "BASE_IMAGE"


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Verify that Dockerfile requires the BASE_IMAGE build argument without default values."
    )
    parser.add_argument(
        "--dockerfile",
        default="docker/Dockerfile",
        help="Path to Dockerfile (default: docker/Dockerfile).",
    )
    parser.add_argument(
        "--target",
        default=os.environ.get("TARGET") or os.environ.get("GITHUB_BASE_REF") or os.environ.get("GITHUB_REF_NAME", ""),
        help="Target branch being evaluated (e.g. develop, main).",
    )
    return parser.parse_args()


def verify_dockerfile(dockerfile_path: Path, build_arg: str = BUILD_ARG) -> None:
    """Verify that Dockerfile requires the build argument without default and consumes it in FROM."""
    if not dockerfile_path.is_file():
        raise FileNotFoundError(f"Dockerfile not found at {dockerfile_path}")

    content = dockerfile_path.read_text(encoding="utf-8")
    lines = content.splitlines()

    arg_declared = False
    arg_has_default = False
    default_val = ""
    from_uses_arg = False

    for line in lines:
        stripped = re.sub(r"#.*$", "", line).strip()
        if not stripped:
            continue

        arg_match = re.match(r"^ARG\s+([A-Za-z_][A-Za-z0-9_]*)(?:\s*=\s*(.*))?$", stripped)
        if arg_match:
            name, default = arg_match.group(1), arg_match.group(2)
            if name == build_arg:
                arg_declared = True
                if default is not None:
                    arg_has_default = True
                    default_val = default.strip()

        from_match = re.match(r"^FROM\s+([^\s]+)", stripped)
        if from_match:
            from_image = from_match.group(1)
            expected_patterns = [
                f"${{{build_arg}}}",
                f"${build_arg}",
            ]
            if any(from_image == p or from_image.startswith(p) for p in expected_patterns):
                from_uses_arg = True
            break

    if not arg_declared:
        raise ValueError(
            f"Build argument '{build_arg}' is not declared in {dockerfile_path}. "
            f"Expected 'ARG {build_arg}' before FROM instruction."
        )

    if arg_has_default:
        raise ValueError(
            f"Build argument '{build_arg}' has default value '{default_val}' in {dockerfile_path}. "
            f"The Dockerfile should require a build arg with the image (use 'ARG {build_arg}' without a default)."
        )

    if not from_uses_arg:
        raise ValueError(
            f"FROM instruction in {dockerfile_path} does not consume '${{{build_arg}}}'. "
            f"Expected 'FROM ${{{build_arg}}}'."
        )


def main() -> int:
    """Run verification checks."""
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    dockerfile_path = Path(args.dockerfile)
    target = args.target.strip()

    logger.info("Evaluating Dockerfile: %s", dockerfile_path)
    if target:
        logger.info("Target branch context: %s", target)
    logger.info("Required build argument: %s", BUILD_ARG)

    try:
        verify_dockerfile(dockerfile_path, BUILD_ARG)
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    logger.info("Success: Dockerfile requires '%s' build argument without default.", BUILD_ARG)
    logger.info("Success: FROM instruction consumes '${%s}'.", BUILD_ARG)
    logger.info("Success: Dockerfile complies with external base image build argument requirements.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
