#!/usr/bin/env python3
"""Validate Docker Hub credentials and registry push authorization for CI workflows."""

import argparse
import base64
import json
import logging
import os
import sys
import urllib.error
import urllib.request

logger = logging.getLogger(__name__)


def check_secrets(username: str, token: str) -> int:
    """Validate presence of Docker repository secrets.

    Returns:
        int: 0 if valid, 1 if missing required secrets.
    """
    missing = [
        name for name, val in [("DOCKER_USERNAME", username), ("DOCKERHUB_TOKEN", token)] if not val.strip()
    ]
    if missing:
        logger.error("Missing required Docker repository secrets: %s", " ".join(missing))
        return 1

    logger.info("All required Docker secrets are present: DOCKER_USERNAME=%s", username)
    return 0


def get_docker_auth_token(repo: str, username: str, token: str) -> str:
    """Obtain a scoped bearer token from Docker authentication service."""
    url = f"https://auth.docker.io/token?service=registry.docker.io&scope=repository:{repo}:push,pull"
    basic_auth = base64.b64encode(f"{username}:{token}".encode()).decode("utf-8")
    req = urllib.request.Request(
        url,
        headers={"Authorization": f"Basic {basic_auth}"},
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            bearer_token = str(data.get("token") or data.get("access_token") or "")
            if not bearer_token:
                raise ValueError(f"No token returned in auth response for {repo}: {data}")
            return bearer_token
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise ValueError(f"HTTP {exc.code} fetching Docker auth token for {repo}: {body}") from exc
    except urllib.error.URLError as exc:
        raise ValueError(f"Network error fetching Docker auth token for {repo}: {exc.reason}") from exc


def check_single_repository_push(repo: str, username: str, token: str) -> bool:
    """Verify push permissions for a single Docker Hub repository.

    Returns:
        bool: True if authorized, False otherwise.
    """
    logger.info("Checking push permissions for repository: %s...", repo)
    try:
        bearer_token = get_docker_auth_token(repo, username, token)
    except ValueError as exc:
        logger.error("%s", exc)
        return False

    upload_url = f"https://registry-1.docker.io/v2/{repo}/blobs/uploads/"
    req = urllib.request.Request(
        upload_url,
        headers={"Authorization": f"Bearer {bearer_token}"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            if resp.status == 202:
                logger.info("Push access successfully verified for %s (HTTP 202 Accepted).", repo)
                location: str | None = resp.headers.get("Location")
                if location:
                    del_url = location if location.startswith("http") else f"https://registry-1.docker.io{location}"
                    del_req = urllib.request.Request(
                        del_url,
                        headers={"Authorization": f"Bearer {bearer_token}"},
                        method="DELETE",
                    )
                    try:
                        with urllib.request.urlopen(del_req, timeout=10):
                            pass
                    except (urllib.error.URLError, TimeoutError, OSError) as exc:
                        logger.warning(
                            "Could not cancel upload session for %s (%s); Docker Hub expires it automatically.",
                            repo,
                            exc,
                        )
                return True
            logger.error("Unexpected status %s checking push for %s.", resp.status, repo)
            return False
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        logger.error("Registry push check failed for %s (HTTP %s): %s", repo, exc.code, body)
        return False
    except urllib.error.URLError as exc:
        logger.error("Network error checking registry push for %s: %s", repo, exc.reason)
        return False


def check_push(repositories: list[str], username: str, token: str) -> int:
    """Verify push permissions for a list of repositories."""
    valid_repos = [r.strip() for r in repositories if r.strip()]
    if not valid_repos:
        logger.error("No repositories specified for push access verification.")
        return 1

    results = [check_single_repository_push(r, username, token) for r in valid_repos]
    if all(results):
        logger.info("All specified Docker Hub repositories verified successfully.")
        return 0
    return 1


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Verify Docker Hub secrets and push authorization.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: check-secrets
    secrets_parser = subparsers.add_parser("check-secrets", help="Check presence of required secrets.")
    secrets_parser.add_argument(
        "--username",
        default=os.environ.get("DOCKER_USERNAME", ""),
        help="Docker Hub username.",
    )
    secrets_parser.add_argument(
        "--token",
        default=os.environ.get("DOCKERHUB_TOKEN", ""),
        help="Docker Hub token.",
    )

    # Subcommand: check-push
    push_parser = subparsers.add_parser("check-push", help="Verify push permissions for repositories.")
    push_parser.add_argument(
        "--repositories",
        nargs="+",
        default=[],
        help="List of repository names to test (e.g. noaaepic/catchem-ubuntu-gcc-13).",
    )
    push_parser.add_argument(
        "--username",
        default=os.environ.get("DOCKER_USERNAME", ""),
        help="Docker Hub username.",
    )
    push_parser.add_argument(
        "--token",
        default=os.environ.get("DOCKERHUB_TOKEN", ""),
        help="Docker Hub token.",
    )

    return parser.parse_args()


def main() -> int:
    """Main execution entry point."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()

    if args.command == "check-secrets":
        return check_secrets(args.username, args.token)

    if args.command == "check-push":
        repos = [item for r in args.repositories for item in r.replace(",", " ").split() if item]
        return check_push(repos, args.username, args.token)

    return 1


if __name__ == "__main__":
    sys.exit(main())
