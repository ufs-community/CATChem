#!/usr/bin/env python3
"""Evaluate release outputs and report semantic release actions to GitHub Actions Step Summary."""

import argparse
import logging
import os
import subprocess
import sys

logger = logging.getLogger(__name__)


def str_to_bool(val: str | bool) -> bool:
    """Convert string or boolean value to boolean."""
    if isinstance(val, bool):
        return val
    return str(val).strip().lower() in ("true", "1", "yes")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Evaluate release outputs and generate summary report.")
    parser.add_argument(
        "--target",
        default=os.environ.get("TARGET") or os.environ.get("GITHUB_BASE_REF") or os.environ.get("GITHUB_REF_NAME", ""),
        help="Target branch evaluated.",
    )
    parser.add_argument(
        "--version",
        default=os.environ.get("VERSION", ""),
        help="Release version produced by PSR.",
    )
    parser.add_argument(
        "--tag",
        default=os.environ.get("TAG", ""),
        help="Git tag produced by PSR.",
    )
    parser.add_argument(
        "--released",
        nargs="?",
        const="true",
        default=os.environ.get("RELEASED", "false"),
        help="Whether a release was triggered (true/false).",
    )
    parser.add_argument(
        "--dry-run",
        nargs="?",
        const="true",
        default=os.environ.get("DRY_RUN", "false"),
        help="Whether the run was in dry-run / preview mode (true/false).",
    )
    return parser.parse_args()


def is_prerelease(version: str, target: str) -> bool:
    """Determine whether the release is a pre-release."""
    return "-" in version or target == "develop"


def write_github_output(outputs: dict[str, str]) -> None:
    """Write key-value pairs to GITHUB_OUTPUT environment file."""
    output_path = os.environ.get("GITHUB_OUTPUT")
    if not output_path:
        return
    with open(output_path, "a", encoding="utf-8") as f:
        for k, v in outputs.items():
            f.write(f"{k}={v}\n")


def get_git_diff() -> str:
    """Capture prospective git diff from dry run file modifications."""
    try:
        res = subprocess.run(["git", "diff", "HEAD"], capture_output=True, text=True, check=True)
        return res.stdout.strip()
    except Exception as exc:
        logger.warning("Failed to capture git diff: %s", exc)
        return ""


def generate_report(target: str, version: str, tag: str, released: bool, dry_run: bool, prerelease: bool) -> str:
    """Build Markdown report content."""
    lines = ["### 🚀 Semantic Release Plan"]
    if dry_run:
        lines.append("**Mode**: Preview (Dry Run / No-op) — No tags or releases created in this run.\n")
    else:
        lines.append("**Mode**: Production Execution\n")

    lines.append("| Parameter | Value |")
    lines.append("|---|---|")
    lines.append(f"| Target Branch | `{target}` |")
    lines.append(f"| Next Version | `{version or 'None'}` |")
    lines.append(f"| Git Tag | `{tag or 'None'}` |")
    lines.append(f"| Will Release? | `{str(released).lower()}` |")
    lines.append(f"| Is Prerelease? | `{str(prerelease).lower()}` |\n")

    if released:
        lines.append(f"#### 📦 Actions on Merge to `{target}`:")
        lines.append(f"- Git tag `{tag}` will be created.")
        lines.append("- Release notes and `CHANGELOG.md` will be published.")
        if prerelease:
            lines.append(f"- Container `catchem-ubuntu-gcc-13-dev:{version}` will be built and pushed.")
            lines.append("- Docker tag `:latest` will be updated on `catchem-ubuntu-gcc-13-dev`.")
        else:
            lines.append(f"- Container `catchem-ubuntu-gcc-13:{version}` will be built and pushed.")
            lines.append("- Docker tag `:latest` will be updated on `catchem-ubuntu-gcc-13`.")

        if dry_run:
            diff_text = get_git_diff()
            if diff_text:
                lines.append("\n#### 📝 Projected Repository Diff")
                lines.append("<details open>")
                lines.append("<summary>Click to collapse projected file changes</summary>\n")
                lines.append("```diff")
                lines.append(diff_text)
                lines.append("```")
                lines.append("</details>")
    else:
        lines.append("#### ℹ️ No release will be triggered on merge based on current commit history.")

    return "\n".join(lines) + "\n"


def main() -> int:
    """Main execution entry point."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()

    target = args.target.strip()
    version = args.version.strip()
    tag = args.tag.strip()
    released = str_to_bool(args.released)
    dry_run = str_to_bool(args.dry_run)

    prerelease = is_prerelease(version, target)

    # 1. Export outputs for downstream jobs
    write_github_output(
        {
            "is_prerelease": "true" if prerelease else "false",
            "target_branch": target,
        }
    )

    # 2. Build and publish report
    report_content = generate_report(target, version, tag, released, dry_run, prerelease)
    print(report_content)

    step_summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if step_summary_path:
        with open(step_summary_path, "a", encoding="utf-8") as f:
            f.write(report_content)

    return 0


if __name__ == "__main__":
    sys.exit(main())
