#!/usr/bin/env python3
"""Validate GitHub App credentials, repository permissions, and ruleset bypass for Semantic Release."""

import argparse
import json
import logging
import os
import subprocess
import sys

logger = logging.getLogger(__name__)


def write_github_output(outputs: dict[str, str]) -> None:
    """Write key-value pairs to GITHUB_OUTPUT environment file."""
    if out := os.environ.get("GITHUB_OUTPUT"):
        with open(out, "a", encoding="utf-8") as f:
            f.writelines(f"{k}={v}\n" for k, v in outputs.items())


def check_secrets(app_id: str, private_key: str, allow_missing: bool) -> int:
    """Validate presence and structure of Semantic Release GitHub App secrets."""
    missing = [
        name for name, val in [("SEMVER_APP_ID", app_id), ("SEMVER_APP_PRIVATE_KEY", private_key)] if not val.strip()
    ]
    if missing:
        missing_str = " ".join(missing)
        if allow_missing:
            logger.info("Semantic Release secrets not provided (%s).", missing_str)
            logger.info("Skipping credential verification for preview run.")
            write_github_output({"skip_verification": "true"})
            return 0
        print(f"::error::Missing required Semantic Release repository secrets: {missing_str}")
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
            print(f"::error::{err}")
            return 1

    write_github_output({"skip_verification": "false"})
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


def check_permissions(repo: str, sha: str, run_id: str, branches: list[str], app_id: str, app_slug: str) -> int:
    """Verify repository access, ref mutation permissions, and branch ruleset bypass."""
    logger.info("Verifying GitHub App permissions for %s...", repo)

    code, out = run_gh_api("/installation/repositories")
    if code != 0:
        print(f"::error::Failed to query installation repositories: {out.strip()}")
        return 1

    try:
        data = json.loads(out)
        repo_names = [r.get("full_name") for r in data.get("repositories", [])]
        if repo not in repo_names:
            print(f"::error::GitHub App is not installed on repository {repo}.")
            return 1
    except json.JSONDecodeError:
        print(f"::error::Invalid JSON returned from installation repositories API: {out.strip()}")
        return 1

    logger.info("GitHub App repository access confirmed for %s.", repo)

    test_ref = f"tags/ci-perm-check-{run_id}"
    logger.info("Testing Git ref creation (%s)...", test_ref)
    code, out = run_gh_api(f"repos/{repo}/git/refs", method="POST", fields={"ref": f"refs/{test_ref}", "sha": sha})
    if code != 0:
        print(f"::error::Failed to create test Git ref on {repo}. Verify App has 'Contents: Read and write'.\n{out}")
        return 1

    logger.info("Git ref creation succeeded. Cleaning up test ref...")
    run_gh_api(f"repos/{repo}/git/refs/{test_ref}", method="DELETE")

    # Check ruleset bypass configuration
    code, out = run_gh_api(f"repos/{repo}/rulesets")
    if code == 0:
        try:
            for rs in json.loads(out):
                if not (rs_id := rs.get("id")):
                    continue
                d_code, d_out = run_gh_api(f"repos/{repo}/rulesets/{rs_id}")
                if d_code != 0:
                    continue
                detail = json.loads(d_out)
                if not any(r.get("type") == "pull_request" for r in detail.get("rules", [])):
                    continue

                rs_name = detail.get("name", str(rs_id))
                can_bypass = detail.get("current_user_can_bypass") == "always"
                has_bypass = can_bypass or any(
                    a.get("actor_type") == "Integration" and a.get("bypass_mode") == "always"
                    for a in detail.get("bypass_actors", [])
                )

                actor_desc = app_slug or f"ID {app_id}"
                for branch in branches:
                    if has_bypass:
                        logger.info(
                            "Ruleset '%s' includes App (%s) in bypass list for %s.", rs_name, actor_desc, branch
                        )
                    else:
                        print(
                            f"::warning::Ruleset '{rs_name}' requires PRs on {branch}; "
                            f"App ({actor_desc}) not in bypass list."
                        )
        except json.JSONDecodeError:
            pass

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
    secrets_parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Gracefully skip if secrets are missing (e.g. preview dry run).",
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
        return check_secrets(args.app_id, args.private_key, args.allow_missing)

    if args.command == "check-permissions":
        branches = [item for b in args.branches for item in b.replace(",", " ").split() if item]
        return check_permissions(args.repo, args.sha, args.run_id, branches, args.app_id, args.app_slug)

    return 1


if __name__ == "__main__":
    sys.exit(main())
