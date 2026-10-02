#!/usr/bin/env python3
"""Evaluate release outputs and report semantic release actions to GitHub Actions Step Summary."""

import argparse
import logging
import os
import re
import subprocess
import sys

logger = logging.getLogger(__name__)


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
        help="Release version produced by PSR in production.",
    )
    parser.add_argument(
        "--tag",
        default=os.environ.get("TAG", ""),
        help="Git tag produced by PSR in production.",
    )
    parser.add_argument(
        "--released",
        nargs="?",
        const="true",
        default=os.environ.get("RELEASED", "false"),
        help="Whether a release was triggered in production (true/false).",
    )
    parser.add_argument(
        "--version-merge",
        default=os.environ.get("VERSION_MERGE", ""),
        help="Release version produced by PSR for merge strategy preview.",
    )
    parser.add_argument(
        "--tag-merge",
        default=os.environ.get("TAG_MERGE", ""),
        help="Git tag produced by PSR for merge strategy preview.",
    )
    parser.add_argument(
        "--released-merge",
        nargs="?",
        const="true",
        default=os.environ.get("RELEASED_MERGE", ""),
        help="Whether a release was triggered for merge strategy preview (true/false).",
    )
    parser.add_argument(
        "--version-squash",
        default=os.environ.get("VERSION_SQUASH", ""),
        help="Release version produced by PSR for squash strategy preview.",
    )
    parser.add_argument(
        "--tag-squash",
        default=os.environ.get("TAG_SQUASH", ""),
        help="Git tag produced by PSR for squash strategy preview.",
    )
    parser.add_argument(
        "--released-squash",
        nargs="?",
        const="true",
        default=os.environ.get("RELEASED_SQUASH", ""),
        help="Whether a release was triggered for squash strategy preview (true/false).",
    )
    parser.add_argument(
        "--diff-merge-file",
        default=os.environ.get("DIFF_MERGE_FILE", ""),
        help="Path to file containing projected git diff for merge strategy.",
    )
    parser.add_argument(
        "--diff-squash-file",
        default=os.environ.get("DIFF_SQUASH_FILE", ""),
        help="Path to file containing projected git diff for squash strategy.",
    )
    parser.add_argument(
        "--dry-run",
        nargs="?",
        const="true",
        default=os.environ.get("DRY_RUN", "false"),
        help="Whether the run was in dry-run / preview mode (true/false).",
    )
    parser.add_argument(
        "--org",
        default=os.environ.get("DOCKER_ORG", "noaaepic"),
        help="Target Docker Hub organization namespace (default: noaaepic).",
    )
    parser.add_argument(
        "--image-name",
        default=os.environ.get("IMAGE_NAME", "catchem-ubuntu-gcc-13"),
        help="Base container image name (default: catchem-ubuntu-gcc-13).",
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

    raise RuntimeError(
        "Unable to determine current repository version from git tags or pyproject.toml."
    )


def read_diff_file(path: str) -> str:
    """Read diff content from file path if it exists."""
    if not path:
        return ""

    if not os.path.isfile(path):
        logger.info("Diff file '%s' does not exist; no difference to report.", path)
        return ""

    try:
        with open(path, encoding="utf-8") as f:
            content = f.read().strip()
            if not content:
                logger.info("Diff file '%s' contains no changes (zero difference).", path)
            return content
    except OSError as exc:
        logger.warning("Failed to read diff file '%s': %s", path, exc)
        return ""


def describe_actions(
    strategy_label: str,
    released: bool,
    version: str,
    tag: str,
    prerelease: bool,
    org: str,
    image_name: str,
) -> list[str]:
    """Generate action items for a given release strategy."""
    prefix = f"{org}/" if org else ""
    header = f"#### 📦 Actions if merged via {strategy_label}:" if strategy_label else "#### 📦 Actions on Merge:"
    lines = [header]
    if released:
        lines.append(f"- Git tag `{tag}` will be created.")
        lines.append("- Release notes and `CHANGELOG.md` will be published.")
        if prerelease:
            dev_image = f"{prefix}{image_name}-dev"
            lines.append(f"- Container `{dev_image}:{version}` will be built and pushed.")
            lines.append(f"- Docker tag `:latest` will be updated on `{dev_image}`.")
        else:
            prod_image = f"{prefix}{image_name}"
            lines.append(f"- Container `{prod_image}:{version}` will be built and pushed.")
            lines.append(
                f"- Docker tag `:latest` will be updated to `{version}` on `{prod_image}`."
            )
    else:
        lines.append("- ℹ️ No release will be triggered based on conventional commit rules.")
    return lines


def generate_report(
    target: str,
    current_version: str,
    dry_run: bool,
    org: str,
    image_name: str,
    version: str = "",
    tag: str = "",
    released: bool = False,
    version_merge: str = "",
    tag_merge: str = "",
    released_merge: bool = False,
    version_squash: str = "",
    tag_squash: str = "",
    released_squash: bool = False,
    diff_merge: str = "",
    diff_squash: str = "",
) -> str:
    """Build Markdown report content."""
    is_dual = dry_run
    lines = ["### 🚀 Semantic Release Plan"]
    if dry_run:
        lines.append("**Mode**: Preview (Dry Run / No-op) — No tags or releases created in this run.\n")
    else:
        lines.append("**Mode**: Production Execution\n")

    if is_dual:
        prerelease_merge = is_prerelease(version_merge, target)
        prerelease_squash = is_prerelease(version_squash, target)
        table_rows = [
            ("Target Branch", target, target),
            ("Current Version", current_version or "None", current_version or "None"),
            ("Next Version", version_merge or "None", version_squash or "None"),
            ("Git Tag", tag_merge or "None", tag_squash or "None"),
            ("Will Release?", str(released_merge).lower(), str(released_squash).lower()),
            ("Is Prerelease?", str(prerelease_merge).lower(), str(prerelease_squash).lower()),
        ]
        lines.extend(
            ["| Parameter | Merge Commit | Squash Merge |", "|---|---|---|"]
            + [f"| {p} | `{vm}` | `{vs}` |" for p, vm, vs in table_rows]
            + [""]
        )
        lines.extend(
            describe_actions(
                "Merge Commit", released_merge, version_merge, tag_merge, prerelease_merge, org, image_name
            )
        )
        lines.append("")
        lines.extend(
            describe_actions(
                "Squash Merge", released_squash, version_squash, tag_squash, prerelease_squash, org, image_name
            )
        )

        if diff_merge:
            lines.append("\n#### 📝 Projected Repository Diff (Merge Commit)")
            lines.append("<details open>")
            lines.append("<summary>Click to collapse projected file changes</summary>\n")
            lines.append("```diff")
            lines.append(diff_merge)
            lines.append("```")
            lines.append("</details>")
        if diff_squash:
            lines.append("\n#### 📝 Projected Repository Diff (Squash Merge)")
            lines.append("<details open>")
            lines.append("<summary>Click to collapse projected file changes</summary>\n")
            lines.append("```diff")
            lines.append(diff_squash)
            lines.append("```")
            lines.append("</details>")
    else:
        eff_version = version or version_merge
        eff_tag = tag or tag_merge
        eff_released = released or released_merge
        prerelease = is_prerelease(eff_version, target)
        table_rows_single = [
            ("Target Branch", target),
            ("Current Version", current_version or "None"),
            ("Next Version", eff_version or "None"),
            ("Git Tag", eff_tag or "None"),
            ("Will Release?", str(eff_released).lower()),
            ("Is Prerelease?", str(prerelease).lower()),
        ]
        lines.extend(["| Parameter | Value |", "|---|---|"] + [f"| {p} | `{v}` |" for p, v in table_rows_single] + [""])
        lines.extend(describe_actions("", eff_released, eff_version, eff_tag, prerelease, org, image_name))
        diff_text = diff_merge or get_git_diff()
        if dry_run and diff_text:
            lines.append("\n#### 📝 Projected Repository Diff")
            lines.append("<details open>")
            lines.append("<summary>Click to collapse projected file changes</summary>\n")
            lines.append("```diff")
            lines.append(diff_text)
            lines.append("```")
            lines.append("</details>")

    return "\n".join(lines) + "\n"


def main() -> int:
    """Main execution entry point."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()

    target = args.target.strip()
    version = args.version.strip()
    tag = args.tag.strip()
    released = str(args.released).lower() in ("true", "1", "yes")
    dry_run = str(args.dry_run).lower() in ("true", "1", "yes")
    org = args.org.strip()
    image_name = args.image_name.strip()

    version_merge = args.version_merge.strip()
    tag_merge = args.tag_merge.strip()
    released_merge = str(args.released_merge).lower() in ("true", "1", "yes")

    version_squash = args.version_squash.strip()
    tag_squash = args.tag_squash.strip()
    released_squash = str(args.released_squash).lower() in ("true", "1", "yes")

    diff_merge = read_diff_file(args.diff_merge_file)
    diff_squash = read_diff_file(args.diff_squash_file)

    try:
        current_version = args.current_version.strip() or determine_current_version(
            target, released or released_merge, dry_run
        )
    except RuntimeError as exc:
        logger.error("%s", exc)
        return 1

    primary_version = version if not dry_run else (version_merge or version_squash or version)
    primary_released = released if not dry_run else (released_merge or released_squash or released)
    primary_tag = tag if not dry_run else (tag_merge or tag_squash or tag)
    prerelease = is_prerelease(primary_version, target)

    # 1. Export outputs for downstream jobs
    write_github_output(
        {
            "is_prerelease": "true" if prerelease else "false",
            "target_branch": target,
            "current_version": current_version,
            "new_release_version": primary_version,
            "new_release_published": "true" if primary_released else "false",
            "tag": primary_tag,
        }
    )

    # 2. Build and publish report
    report_content = generate_report(
        target=target,
        current_version=current_version,
        dry_run=dry_run,
        org=org,
        image_name=image_name,
        version=version,
        tag=tag,
        released=released,
        version_merge=version_merge,
        tag_merge=tag_merge,
        released_merge=released_merge,
        version_squash=version_squash,
        tag_squash=tag_squash,
        released_squash=released_squash,
        diff_merge=diff_merge,
        diff_squash=diff_squash,
    )

    print(report_content)

    if step_summary_path := os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(step_summary_path, "a", encoding="utf-8") as f:
            f.write(report_content)

    return 0


if __name__ == "__main__":
    sys.exit(main())
