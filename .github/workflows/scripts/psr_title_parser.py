#!/usr/bin/env python3
"""Python Semantic Release commit parser that evaluates only the commit title.

The squash-merge body on develop/main carries every PR commit message verbatim, so
release decisions must come from the first line alone: type and scope from
``type(scope):`` and a breaking change only from ``!``. Body footers such as
``BREAKING CHANGE:`` are ignored.
"""

from semantic_release.commit_parser.conventional import ConventionalCommitParser

# Referenced by name from pyproject.toml ([tool.semantic_release] commit_parser).
__all__ = ["TitleOnlyCommitParser"]


class TitleOnlyCommitParser(ConventionalCommitParser):
    """Conventional commit parser restricted to the first line of the message."""

    def parse_message(self, message: str):
        return super().parse_message(message.split("\n", 1)[0].strip())
