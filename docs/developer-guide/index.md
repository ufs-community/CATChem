# CATChem Developer Guide

Welcome to the CATChem Developer Guide! Whether you are a new contributor or an experienced developer, this guide provides the information you need to understand, develop, and extend the CATChem framework.

## Getting Started

If you are new to CATChem, here are the first steps to get you started:

1.  **[Contributor Guide](contributing.md)**: This is the essential starting point for all contributors. It covers the basics of how to contribute to the project, including our development workflow and pull request process.
2.  **[Coding Standards](coding-standards.md)**: Before writing any code, please review our coding standards. This will help you write code that is consistent, maintainable, and performant.
3.  **[Architecture Overview](architecture.md)**: For a high-level overview of the software architecture and design principles of CATChem, please refer to this guide.

## Core Concepts

CATChem is built on a set of modern software engineering principles and patterns. Understanding these core concepts is key to developing effectively within the framework. For a deep dive into the foundational infrastructure of CATChem, including state management, configuration, diagnostics, and error handling, please see the **[Core Concepts](../core-concepts/index.md)** section.

## Development Guides

Ready to dive in? These guides provide detailed information on specific areas of development:

<div class="grid cards" markdown>

-   [:material-cog: **Core Concepts**](../core-concepts/index.md)

    ---

    A deep dive into the foundational infrastructure of CATChem, including state management, configuration, diagnostics, and error handling.

-   [:material-puzzle: **Process Development**](processes/index.md)

    ---

    Learn how to create new atmospheric processes and schemes, and how to integrate them into the CATChem framework.

-   [:material-merge: **Integration**](integration/index.md)

    ---

    Guides for integrating CATChem with various host models, such as FV3, and other frameworks like CCPP and NUOPC.

-   [:material-test-tube: **Testing**](testing.md)

    ---

    An overview of the testing framework, including how to write unit tests, integration tests, and validation benchmarks.

-   [:material-book-open-variant-outline: **Documentation**](documentation.md)

    ---

    Our standards and best practices for writing clear, comprehensive documentation for your code.

</div>

## Development Workflow

Our development workflow is based on Git and GitHub. Here is a typical cycle:

1.  **Set up your environment**: Fork the repository, clone it, and create a development build using CMake.
2.  **Create a feature branch**: All new work should be done in a feature branch.
3.  **Develop, test, and document**: Write your code, add corresponding tests, and update the documentation.
4.  **Commit and push**: Use conventional commit messages and push your branch to your fork.
5.  **Create a pull request**: Open a pull request to the main CATChem repository for review.

For detailed instructions, please see the [Contributor Guide](contributing.md).

## Getting Help

We are here to help! If you have questions or need support, please use our community resources:

-   **GitHub Discussions**: For technical discussions, questions, and answers.
-   **Issue Tracker**: For bug reports and feature requests.
-   **Developer Meetings**: We hold regular virtual meetings for developers. Please contact the development team for details.
