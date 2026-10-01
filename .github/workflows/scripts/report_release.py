#!/usr/bin/env python3
"""Evaluate release outputs and report semantic release actions to GitHub Actions Step Summary."""

import argparse
import logging
import os
import re
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
    parser.add_argument(
        "--current-version",
        default=os.environ.get("CURRENT_VERSION", ""),
        help="Current release version being evaluated against.",
    )
    return parser.parse_args()


def is_prerelease(version: str, target: str) -> bool:
    """Determine whether the release is a pre-release."""
    return "-" in version or target == "develop"


def write_github_output(outputs: dict[str, str]) -> None:
    """Write key-value pairs to GITHUB_OUTPUT environment file."""
    if out := os.environ.get("GITHUB_OUTPUT"):
        with open(out, "a", encoding="utf-8") as f:
            f.writelines(f"{k}={v}\n" for k, v in outputs.items())


def get_git_diff() -> str:
    """Capture prospective git diff from dry run file modifications."""
    try:
        return subprocess.check_output(["git", "diff", "HEAD"], text=True, stderr=subprocess.PIPE).strip()
    except subprocess.CalledProcessError:
        return ""


def extract_version_from_toml(text: str) -> str:
    """Extract project version from TOML content using regex."""
    match = re.search(r'^\s*version\s*=\s*["\']([^"\']+)["\']', text, re.MULTILINE)
    return match.group(1).strip().lstrip("v") if match else ""


def determine_current_version(target: str, released: bool, dry_run: bool) -> str:
    """Determine the current version before this release evaluation."""
    refs = (
        ["HEAD~1", "HEAD~2"]
        if (not dry_run and released)
        else ([f"origin/{target}", target] if target else []) + ["HEAD~1", "HEAD"]
    )
    for ref in refs:
        try:
            return (
                subprocess.check_output(
                    ["git", "describe", "--tags", "--abbrev=0", ref],
                    text=True,
                    stderr=subprocess.DEVNULL,
                )
                .strip()
                .lstrip("v")
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass
        try:
            content = subprocess.check_output(
                ["git", "show", f"{ref}:pyproject.toml"],
                text=True,
                stderr=subprocess.DEVNULL,
            )
            if ver := extract_version_from_toml(content):
                return ver
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass

    try:
        with open("pyproject.toml", encoding="utf-8") as f:
            if ver := extract_version_from_toml(f.read()):
                return ver
    except OSError:
        pass

    return "None"


def generate_report(
    target: str,
    current_version: str,
    version: str,
    tag: str,
    released: bool,
    dry_run: bool,
    prerelease: bool,
) -> str:
    """Build Markdown report content."""
    lines = ["### 🚀 Semantic Release Plan"]
    if dry_run:
        lines.append("**Mode**: Preview (Dry Run / No-op) — No tags or releases created in this run.\n")
    else:
        lines.append("**Mode**: Production Execution\n")

    table_rows = [
        ("Target Branch", target),
        ("Current Version", current_version or "None"),
        ("Next Version", version or "None"),
        ("Git Tag", tag or "None"),
        ("Will Release?", str(released).lower()),
        ("Is Prerelease?", str(prerelease).lower()),
    ]
    lines.extend(["| Parameter | Value |", "|---|---|"] + [f"| {p} | `{v}` |" for p, v in table_rows] + [""])

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
    current_version = args.current_version.strip() or determine_current_version(target, released, dry_run)

    prerelease = is_prerelease(version, target)

    # 1. Export outputs for downstream jobs
    write_github_output(
        {
            "is_prerelease": "true" if prerelease else "false",
            "target_branch": target,
            "current_version": current_version,
        }
    )

    # 2. Build and publish report
    try:
        report_content = generate_report(target, current_version, version, tag, released, dry_run, prerelease)
    except subprocess.CalledProcessError as exc:
        logger.error("Git command failed (exit %d): %s\nStderr: %s", exc.returncode, exc.cmd, exc.stderr)
        return exc.returncode

    print(report_content)

    if step_summary_path := os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(step_summary_path, "a", encoding="utf-8") as f:
            f.write(report_content)

    return 0


if __name__ == "__main__":
    sys.exit(main())
