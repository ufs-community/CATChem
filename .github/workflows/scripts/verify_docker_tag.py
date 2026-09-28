#!/usr/bin/env python3
"""Verify Dockerfile base image tag compliance for target branches."""

import argparse
import os
import re
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify Dockerfile base image tag complies with target branch requirements."
    )
    parser.add_argument(
        "--target",
        default=os.environ.get("TARGET") or os.environ.get("GITHUB_BASE_REF") or os.environ.get("GITHUB_REF_NAME", ""),
        help="Target branch to evaluate (e.g. develop, main). Defaults to TARGET or GitHub Actions ref.",
    )
    parser.add_argument(
        "--dockerfile",
        default="docker/Dockerfile",
        help="Path to Dockerfile (default: docker/Dockerfile).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    target = args.target.strip()
    dockerfile_path = Path(args.dockerfile)

    if not dockerfile_path.is_file():
        print(f"::error file={dockerfile_path}::Dockerfile not found at {dockerfile_path}")
        return 1

    content = dockerfile_path.read_text(encoding="utf-8")
    match = re.search(r"^\s*(?:ARG\s+BASE_IMAGE\s*=\s*|FROM\s+)([^\s]+)", content, re.MULTILINE)
    if not match:
        print(f"::error file={dockerfile_path}::Could not find BASE_IMAGE or FROM instruction in {dockerfile_path}")
        return 1

    base_image = match.group(1).strip()
    print(f"Detected base image in {dockerfile_path}: {base_image}")
    print(f"Evaluating for target branch: {target}")

    if target == "develop":
        if not (base_image.endswith("-dev:latest") or ":-dev:latest" in base_image):
            print(
                f'::error file={dockerfile_path}::Invalid base image "{base_image}" for branch "develop". '
                f'Merges to develop must use "*-dev:latest" in {dockerfile_path}.'
            )
            return 1
        print(f'Success: Base image "{base_image}" correctly uses "-dev:latest" for develop branch.')
    elif target == "main":
        if not base_image.endswith(":latest") or "-dev" in base_image:
            print(
                f'::error file={dockerfile_path}::Invalid base image "{base_image}" for branch "main". '
                f'Merges to main must use ":latest" without "-dev" in {dockerfile_path}.'
            )
            return 1
        print(f'Success: Base image "{base_image}" correctly targets production ":latest" without "-dev" for main branch.')
    else:
        print(f'Notice: Target branch "{target}" is neither develop nor main; skipping tag enforcement.')

    return 0


if __name__ == "__main__":
    sys.exit(main())
