#!/usr/bin/env python3
"""Prepare Git repository branch state for prospective semantic release evaluation."""

import argparse
import logging
import os
import subprocess
import sys

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Prepare branch state with a synthetic squash commit for preview release evaluation."
    )
    parser.add_argument(
        "--target",
        default=os.environ.get("TARGET") or os.environ.get("GITHUB_BASE_REF") or os.environ.get("GITHUB_REF_NAME", ""),
        help="Target branch for release evaluation (e.g. develop, main).",
    )
    parser.add_argument(
        "--pr-title",
        default=os.environ.get("PR_TITLE", ""),
        help="Pull request title used for synthetic squash commit message.",
    )
    return parser.parse_args()


def run_git(args: list[str], check: bool = True) -> subprocess.CompletedProcess[str]:
    """Execute a git command with error handling."""
    cmd = ["git"] + args
    return subprocess.run(cmd, capture_output=True, text=True, check=check)


def prepare_preview(target: str, pr_title: str) -> None:
    """Prepare repository branch and synthetic squash merge commit if PR title is present."""
    if not target:
        raise ValueError("Target branch must be specified for release preview.")

    logger.info("Evaluating prospective release against target branch: %s", target)
    run_git(["config", "user.name", "github-actions[bot]"])
    run_git(["config", "user.email", "github-actions[bot]@users.noreply.github.com"])

    if pr_title.strip():
        pr_head = run_git(["rev-parse", "HEAD"]).stdout.strip()
        logger.info("PR head SHA: %s", pr_head)
        run_git(["checkout", "-B", target, f"origin/{target}"])

        merge_res = run_git(["merge", "--squash", pr_head], check=False)
        if merge_res.returncode != 0:
            logger.warning("Squash merge encountered conflicts or errors; aborting merge.")
            run_git(["merge", "--abort"], check=False)

        logger.info("Applying synthetic squash-merge commit from PR title: %s", pr_title)
        run_git(["commit", "--allow-empty", "-m", pr_title])
    else:
        logger.info("No PR title provided; checking out target branch: %s", target)
        run_git(["checkout", "-B", target])


def main() -> int:
    """Main execution entry point."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    try:
        prepare_preview(args.target.strip(), args.pr_title.strip())
    except subprocess.CalledProcessError as exc:
        logger.error("Git command failed (exit %d): %s\nStderr: %s", exc.returncode, exc.cmd, exc.stderr)
        return exc.returncode
    except ValueError as exc:
        logger.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
