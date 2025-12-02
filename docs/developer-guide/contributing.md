# Contributing to CATChem

Thank you for your interest in contributing to CATChem! We welcome contributions from everyone, and we are excited to have you as part of our community. This guide will walk you through the process of contributing to the project.

## Code of Conduct

All members of the CATChem community are expected to abide by our [Code of Conduct](https://www.contributor-covenant.org/). Please make sure you have read and understood it.

## Ways to Contribute

There are many ways to contribute to CATChem:

- **Code Contributions**: Implement new features, fix bugs, or improve performance.
- **Documentation**: Write tutorials, improve the developer guides, or clarify the API documentation.
- **Issue Reporting**: Report bugs or request new features by opening an issue.
- **Code Reviews**: Help us maintain the quality of the codebase by reviewing pull requests.
- **Community Support**: Help other users and developers on our GitHub Discussions forum.

## Your First Contribution

If you are looking for a good place to start, check out the issues labeled ["good first issue"](https://github.com/UFS-Community/CATChem/labels/good%20first%20issue) on our GitHub issue tracker. These are issues that are well-suited for new contributors.

## Development Workflow

We follow an issue-driven development workflow. Here is the process for making a contribution:

### 1. Set Up Your Development Environment

- **Fork and Clone**: Fork the [CATChem repository](https://github.com/UFS-Community/CATChem) on GitHub, and then clone your fork to your local machine.
- **Create a Branch**: Create a new branch for your feature or bugfix. The branch name should be descriptive, e.g., `feature/new-settling-scheme` or `bugfix/fix-memory-leak`.
- **Build the Code**: Follow the instructions in the [User Guide](../user-guide/index.md) to build CATChem. For development, it is recommended to build in debug mode.

### 2. Make Your Changes

- **Write Your Code**: Make your changes, following the [Coding Standards](coding-standards.md).
- **Add Tests**: All new features and bug fixes must be accompanied by tests. See the [Testing Guide](testing.md) for more information.
- **Update Documentation**: If your changes affect the user-facing API or behavior of the model, please update the documentation accordingly.

### 3. Commit Your Changes

We use the [Conventional Commits](https://www.conventionalcommits.org/) specification for our commit messages. This helps us automatically generate changelogs and understand the history of the project.

Each commit message should have the following format:

```
<type>(<scope>): <subject>

<body>

<footer>
```

- **type**: `feat` (new feature), `fix` (bug fix), `docs` (documentation), `style`, `refactor`, `test`, `chore`.
- **scope** (optional): The part of the codebase that is affected (e.g., `chemistry`, `emissions`, `docs`).
- **subject**: A short, imperative-tense description of the change.

**Example:**

```
feat(chemistry): add support for aqueous-phase reactions

- Implement a new solver for aqueous-phase chemistry.
- Add a database of Henry's Law constants.
- Include unit tests for the new solver.

Closes #123
```

### 4. Create a Pull Request

When you are ready to submit your contribution, push your branch to your fork and open a pull request to the `develop` branch of the main CATChem repository.

- **Fill out the template**: Your pull request should include a clear description of the changes you have made and why. Please fill out the pull request template completely.
- **Pass the checks**: Your pull request will be automatically tested by our continuous integration (CI) system. All checks must pass before your pull request can be merged.
- **Respond to feedback**: The maintainers will review your pull request and may request changes. Please be responsive to their feedback.

## Getting Help

If you have any questions or need help with your contribution, please don't hesitate to reach out to us on our [GitHub Discussions](https://github.com/UFS-Community/CATChem/discussions) forum.

Thank you for contributing to CATChem!
