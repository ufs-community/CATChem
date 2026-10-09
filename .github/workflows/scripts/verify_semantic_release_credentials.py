#!/usr/bin/env python3
"""Validate GitHub App credentials, repository permissions, and ruleset bypass for Semantic Release."""

import argparse
import fnmatch
import json
import logging
import os
import subprocess
import sys

logger = logging.getLogger(__name__)


def check_secrets(app_id: str, private_key: str) -> int:
    """Validate presence and structure of Semantic Release GitHub App secrets."""
    missing = [
        name for name, val in [("SEMVER_APP_ID", app_id), ("SEMVER_APP_PRIVATE_KEY", private_key)] if not val.strip()
    ]
    if missing:
        logger.error("Missing required Semantic Release repository secrets: %s", " ".join(missing))
        return 1

    key = private_key.strip()
    checks = [
        (
            key.endswith(".pem") or key.startswith(("/", "~")),
            "SEMVER_APP_PRIVATE_KEY is a file path. Paste file contents instead.",
        ),
        (
            "-----BEGIN" not in key or "-----END" not in key,
            "SEMVER_APP_PRIVATE_KEY missing BEGIN/END RSA PRIVATE KEY markers.",
        ),
        (
            (key.startswith('"') and key.endswith('"')) or (key.startswith("'") and key.endswith("'")),
            "SEMVER_APP_PRIVATE_KEY wrapped in quotes. Remove quotes in Secrets.",
        ),
    ]
    for failed, err in checks:
        if failed:
            logger.error("%s", err)
            return 1

    logger.info("All required Semantic Release secrets are present and formatted properly.")
    return 0


def run_gh_api(endpoint: str, method: str = "GET", fields: dict[str, str] | None = None) -> tuple[int, str]:
    """Execute a GitHub API command using gh CLI."""
    cmd = ["gh", "api", endpoint, "-X", method]
    if fields:
        for k, v in fields.items():
            cmd.extend(["-f", f"{k}={v}"])
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    return res.returncode, res.stdout or res.stderr


def get_default_branch(repo: str) -> str:
    """Return the repository's default branch name."""
    code, out = run_gh_api(f"repos/{repo}")
    if code != 0:
        raise RuntimeError(f"Failed to query repository {repo}: {out.strip()}")
    return str(json.loads(out).get("default_branch", ""))


def ruleset_targets_branch(detail: dict, branch: str, default_branch: str) -> bool:
    """Return True when a branch ruleset's ref_name conditions apply to ``branch``."""
    ref = f"refs/heads/{branch}"
    cond = detail.get("conditions", {}).get("ref_name", {})

    def matches(pattern: str) -> bool:
        if pattern == "~ALL":
            return True
        if pattern == "~DEFAULT_BRANCH":
            return branch == default_branch
        return fnmatch.fnmatchcase(ref, pattern)

    included = any(matches(p) for p in cond.get("include", []))
    excluded = any(matches(p) for p in cond.get("exclude", []))
    return included and not excluded


def check_permissions(repo: str, sha: str, run_id: str, branches: list[str], app_id: str, app_slug: str) -> int:
    """Verify repository access, ref mutation permissions, and branch ruleset bypass."""
    logger.info("Verifying GitHub App permissions for %s...", repo)

    code, out = run_gh_api("/installation/repositories")
    if code != 0:
        logger.error("Failed to query installation repositories: %s", out.strip())
        return 1

    try:
        data = json.loads(out)
        repo_names = [r.get("full_name") for r in data.get("repositories", [])]
        if repo not in repo_names:
            logger.error("GitHub App is not installed on repository %s.", repo)
            return 1
    except json.JSONDecodeError as exc:
        logger.error("Invalid JSON returned from installation repositories API: %s\nOutput: %s", exc, out.strip())
        return 1

    logger.info("GitHub App repository access confirmed for %s.", repo)

    test_ref = f"tags/ci-perm-check-{run_id}"
    logger.info("Testing Git ref creation (%s)...", test_ref)
    code, out = run_gh_api(f"repos/{repo}/git/refs", method="POST", fields={"ref": f"refs/{test_ref}", "sha": sha})
    if code != 0:
        logger.error("Failed to create test Git ref on %s. Verify App has 'Contents: Read and write'.\n%s", repo, out)
        return 1

    logger.info("Git ref creation succeeded. Cleaning up test ref...")
    code, out = run_gh_api(f"repos/{repo}/git/refs/{test_ref}", method="DELETE")
    if code != 0:
        logger.error(
            "Failed to delete test Git ref refs/%s on %s; delete it manually. "
            "Verify App has 'Contents: Read and write'.\n%s",
            test_ref,
            repo,
            out,
        )
        return 1
    logger.info("Test ref refs/%s deleted.", test_ref)

    # Check ruleset bypass configuration
    code, out = run_gh_api(f"repos/{repo}/rulesets")
    if code != 0:
        logger.warning("Failed to query repository rulesets for %s (exit %d): %s", repo, code, out.strip())
    else:
        try:
            rulesets_data = json.loads(out)
        except json.JSONDecodeError as exc:
            logger.error("Failed to parse JSON response from rulesets API: %s\nOutput: %s", exc, out.strip())
            return 1

        try:
            default_branch = get_default_branch(repo)
        except RuntimeError as exc:
            logger.error("%s", exc)
            return 1

        missing_bypass: list[str] = []
        for rs in rulesets_data:
            if not (rs_id := rs.get("id")):
                continue
            d_code, d_out = run_gh_api(f"repos/{repo}/rulesets/{rs_id}")
            if d_code != 0:
                logger.warning("Failed to query detail for ruleset %s (exit %d): %s", rs_id, d_code, d_out.strip())
                continue
            try:
                detail = json.loads(d_out)
            except json.JSONDecodeError as exc:
                logger.error("Failed to parse JSON for ruleset %s: %s\nOutput: %s", rs_id, exc, d_out.strip())
                continue

            if not any(r.get("type") == "pull_request" for r in detail.get("rules", [])):
                continue

            rs_name = detail.get("name", str(rs_id))
            governed = [b for b in branches if ruleset_targets_branch(detail, b, default_branch)]
            if not governed:
                logger.info("Ruleset '%s' does not govern any of %s; skipping.", rs_name, branches)
                continue

            can_bypass = detail.get("current_user_can_bypass") == "always"
            has_bypass = can_bypass or any(
                a.get("actor_type") == "Integration" and a.get("bypass_mode") == "always"
                for a in detail.get("bypass_actors", [])
            )

            actor_desc = app_slug or f"ID {app_id}"
            for branch in governed:
                if has_bypass:
                    logger.info(
                        "Ruleset '%s' includes App (%s) in bypass list for %s.", rs_name, actor_desc, branch
                    )
                else:
                    logger.error(
                        "Ruleset '%s' requires PRs on %s; App (%s) not in bypass list.",
                        rs_name,
                        branch,
                        actor_desc,
                    )
                    missing_bypass.append(f"{rs_name}:{branch}")

        if missing_bypass:
            logger.error(
                "Semantic release cannot push to protected branches without bypass: %s", ", ".join(missing_bypass)
            )
            return 1

    logger.info("Semantic release credential verification completed successfully.")
    return 0


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Verify Semantic Release secrets and repository permissions.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: check-secrets
    secrets_parser = subparsers.add_parser("check-secrets", help="Check presence and formatting of secrets.")
    secrets_parser.add_argument(
        "--app-id",
        default=os.environ.get("SEMVER_APP_ID", ""),
        help="Semantic release GitHub App ID.",
    )
    secrets_parser.add_argument(
        "--private-key",
        default=os.environ.get("SEMVER_APP_PRIVATE_KEY", ""),
        help="Semantic release GitHub App RSA private key.",
    )

    # Subcommand: check-permissions
    perm_parser = subparsers.add_parser("check-permissions", help="Verify repository permissions and ruleset bypass.")
    perm_parser.add_argument(
        "--repo",
        default=os.environ.get("REPO") or os.environ.get("GITHUB_REPOSITORY", ""),
        help="GitHub repository (owner/repo).",
    )
    perm_parser.add_argument(
        "--sha",
        default=os.environ.get("SHA") or os.environ.get("GITHUB_SHA", ""),
        help="Git commit SHA for test ref creation.",
    )
    perm_parser.add_argument(
        "--run-id",
        default=os.environ.get("RUN_ID") or os.environ.get("GITHUB_RUN_ID", ""),
        help="GitHub Actions run ID for unique ref names.",
    )
    perm_parser.add_argument(
        "--branches",
        nargs="+",
        default=["develop", "main"],
        help="Branches to evaluate for ruleset bypass.",
    )
    perm_parser.add_argument(
        "--app-id",
        default=os.environ.get("APP_ID") or os.environ.get("SEMVER_APP_ID", ""),
        help="GitHub App ID.",
    )
    perm_parser.add_argument(
        "--app-slug",
        default=os.environ.get("APP_SLUG", ""),
        help="GitHub App slug.",
    )

    return parser.parse_args()


def main() -> int:
    """Main execution entry point."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()

    if args.command == "check-secrets":
        return check_secrets(args.app_id, args.private_key)

    if args.command == "check-permissions":
        branches = [item for b in args.branches for item in b.replace(",", " ").split() if item]
        return check_permissions(args.repo, args.sha, args.run_id, branches, args.app_id, args.app_slug)

    return 1


if __name__ == "__main__":
    sys.exit(main())
