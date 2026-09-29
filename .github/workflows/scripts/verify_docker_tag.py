#!/usr/bin/env python3
"""Verify Dockerfile base image argument and target branch tag compliance."""

import argparse
import logging
import os
import re
import sys
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Verify that Dockerfile requires and appropriately accepts the base image "
            "build argument prescribed by CI."
        )
    )
    parser.add_argument(
        "--target",
        default=os.environ.get("TARGET") or os.environ.get("GITHUB_BASE_REF") or os.environ.get("GITHUB_REF_NAME", ""),
        help="Target branch to evaluate (e.g. develop, main). Defaults to TARGET or GitHub Actions ref.",
    )
    parser.add_argument(
        "--image",
        default=os.environ.get("BASE_IMAGE") or os.environ.get("IMAGE", ""),
        help="Prescribed base image (e.g. noaaepic/ufschem-spack-base-ubuntu-gcc-13-dev:latest).",
    )
    parser.add_argument(
        "--build-arg",
        default=os.environ.get("BUILD_ARG", "BASE_IMAGE"),
        help="Name of the required build argument in Dockerfile (default: BASE_IMAGE).",
    )
    parser.add_argument(
        "--docker-org",
        default=os.environ.get("DOCKER_ORG", "noaaepic"),
        help="Docker organization namespace (default: noaaepic or DOCKER_ORG env var).",
    )
    parser.add_argument(
        "--dockerfile",
        default="docker/Dockerfile",
        help="Path to Dockerfile (default: docker/Dockerfile).",
    )
    parser.add_argument(
        "--log-level",
        default=os.environ.get("LOG_LEVEL", "INFO").upper(),
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Set the logging level (default: INFO).",
    )
    return parser.parse_args()


def verify_dockerfile(dockerfile_path: Path, build_arg: str) -> tuple[bool, str]:
    """Verify that Dockerfile requires the build argument without default and consumes it in FROM."""
    if not dockerfile_path.is_file():
        return False, f"Dockerfile not found at {dockerfile_path}"

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
        return (
            False,
            f"Build argument '{build_arg}' is not declared in {dockerfile_path}. "
            f"Expected 'ARG {build_arg}' before FROM instruction.",
        )

    if arg_has_default:
        return (
            False,
            f"Build argument '{build_arg}' has default value '{default_val}' in {dockerfile_path}. "
            f"The Dockerfile should require a build arg with the image (use 'ARG {build_arg}' without a default).",
        )

    if not from_uses_arg:
        return (
            False,
            f"FROM instruction in {dockerfile_path} does not consume '${{{build_arg}}}'. "
            f"Expected 'FROM ${{{build_arg}}}'.",
        )

    return True, ""


def resolve_prescribed_image(target: str, image_arg: str, docker_org: str) -> str:
    """Resolve the prescribed base image for the given target branch."""
    if image_arg:
        return image_arg.strip()
    org = docker_org.strip() or "noaaepic"
    if target == "main":
        return f"{org}/ufschem-spack-base-ubuntu-gcc-13:latest"
    return f"{org}/ufschem-spack-base-ubuntu-gcc-13-dev:latest"


def verify_image_tag(target: str, image: str) -> tuple[bool, str]:
    """Verify that the prescribed base image tag complies with target branch requirements."""
    if target == "develop":
        if not (image.endswith("-dev:latest") or ":-dev:latest" in image):
            return (
                False,
                f'Invalid prescribed base image "{image}" for branch "develop". '
                f'Merges to develop must prescribe a base image ending with "-dev:latest".',
            )
        return True, f'Prescribed base image "{image}" correctly uses "-dev:latest" for develop branch.'
    if target == "main":
        if not image.endswith(":latest") or "-dev" in image:
            return (
                False,
                f'Invalid prescribed base image "{image}" for branch "main". '
                f'Merges to main must prescribe a base image ending with ":latest" without "-dev".',
            )
        return (
            True,
            f'Prescribed base image "{image}" correctly targets production ":latest" without "-dev" for main branch.',
        )
    return True, f'Target branch "{target}" is neither develop nor main; skipping branch tag enforcement.'


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level, logging.INFO),
        format="%(levelname)s: %(message)s",
    )

    target = args.target.strip()
    dockerfile_path = Path(args.dockerfile)
    build_arg = args.build_arg.strip()
    docker_org = args.docker_org.strip()

    prescribed_image = resolve_prescribed_image(target, args.image, docker_org)

    logger.info("Evaluating Dockerfile: %s", dockerfile_path)
    logger.info("Target branch: %s", target or "(none)")
    logger.info("Required build argument: %s", build_arg)
    logger.info("Prescribed base image from CI: %s", prescribed_image)

    # 1. Verify Dockerfile requires the build arg and uses it in FROM
    df_ok, df_err = verify_dockerfile(dockerfile_path, build_arg)
    if not df_ok:
        logger.error("%s", df_err)
        return 1
    logger.info("Success: Dockerfile requires '%s' build argument without default.", build_arg)
    logger.info("Success: FROM instruction consumes '${%s}'.", build_arg)

    # 2. Verify prescribed base image complies with target branch requirements
    if target:
        tag_ok, tag_msg = verify_image_tag(target, prescribed_image)
        if not tag_ok:
            logger.error("%s", tag_msg)
            return 1
        logger.info("Success: %s", tag_msg)
    else:
        logger.warning("No target branch specified; skipping branch-specific tag enforcement.")

    logger.info(
        "Success: Dockerfile will appropriately accept the prescribed base image argument (%s=%s).",
        build_arg,
        prescribed_image,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
