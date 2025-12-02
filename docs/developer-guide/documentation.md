# Documentation Guide

Good documentation is essential for a successful software project. This guide explains how to write, build, and contribute to the documentation for CATChem.

## Documentation Philosophy

Our goal is to have documentation that is:

- **Clear and Concise**: Easy to understand for both users and developers.
- **Comprehensive**: Covers all aspects of the model, from installation to advanced development.
- **Up-to-date**: Kept in sync with the latest changes to the codebase.
- **Accessible**: Easy to navigate and search.

## The Documentation Stack

CATChem's documentation is built using a modern, open-source toolchain:

- **MkDocs**: A fast and simple static site generator for building project documentation.
- **Material for MkDocs**: A beautiful and highly configurable theme for MkDocs.
- **Doxygen**: A tool for generating API documentation directly from comments in the source code.
- **MkDoxy**: A plugin that integrates Doxygen-generated API documentation into the MkDocs site.

The documentation is written in a combination of Markdown (for narrative guides) and Doxygen-style comments in the Fortran source code (for the API reference).

## Documentation Structure

The documentation source files are located in the `docs/` directory.

```
docs/
├── index.md                  # Main landing page
├── user-guide/               # End-user documentation
├── processes/                # Process-specific documentation used in the user-guide
├── developer-guide/          # Technical documentation for developers
├── api/                      # API reference
├── ufschem/                  # How CATChem is integrated into the UFS to create UFS-Chem
├── community/                # How to get involved?
├── evaluation/               # How CATChem and UFS-Chem are evaluated
├── assets/                   # Images, stylesheets, scripts
└── ../mkdocs.yml             # Documentation configuration
```

## Writing Documentation

### Markdown Guides

The narrative documentation, such as the user and developer guides, is written in Markdown. We follow the [CommonMark](https://commonmark.org/) specification, with some extensions provided by the Material for MkDocs theme.

**Example Markdown:**

```markdown
# My Awesome Feature

This is a guide to my awesome new feature.

!!! note
    This is a note. Use admonitions to highlight important information.

## How it Works

Here is a code block with syntax highlighting:

```fortran
subroutine my_awesome_subroutine()
  print *, "Hello, world!"
end subroutine
```
```

For more information on the available Markdown extensions, please see the [Material for MkDocs documentation](https://squidfunk.github.io/mkdocs-material/reference/).

### API Documentation with Doxygen

The API documentation is generated automatically from Doxygen-style comments in the Fortran source code. All public procedures, types, and modules must be documented.

**Example Doxygen Comments:**

```fortran
!> @brief A brief, one-line description of the subroutine.
!>
!> A more detailed description of what the subroutine does,
!> its parameters, and any important notes for the user.
!>
!> @param[in]  input_arg   Description of the input argument.
!> @param[out] output_arg  Description of the output argument.
!> @param[out] rc          The return code (CC_SUCCESS on success).
subroutine my_subroutine(input_arg, output_arg, rc)
  implicit none
  real, intent(in) :: input_arg
  real, intent(out) :: output_arg
  integer, intent(out) :: rc
  ! ...
end subroutine
```

Key Doxygen tags to use:

- `@brief`: A short, one-line description.
- `@param[in|out|inout]`: Describes a parameter.
- `@return`: Describes the return value of a function.
- `@note`: An optional note with additional information.
- `@warning`: An important warning.

## Building the Documentation

To build and preview the documentation locally, you will need to set up a Python environment.

### Setup

1.  **Create a Conda Environment**: We recommend using Conda to manage your documentation dependencies.
    ```bash
    conda env create -f docs/environment-docs.yml
    conda activate catchem-docs
    ```
2.  **Install Dependencies with Pip (Alternative)**: If you are not using Conda, you can install the dependencies with pip.
    ```bash
    pip install -r docs/requirements.txt
    ```

### Local Development Server

The easiest way to work on the documentation is to use the MkDocs built-in development server. It will automatically rebuild the documentation whenever you save a file and reload your browser.

```bash
mkdocs serve
```

You can then view the documentation at `http://127.0.0.1:8000`.

### Building the Site

To build the static HTML site, run the following command:

```bash
mkdocs build
```

The generated site will be in the `site/` directory.

## Best Practices for Documentation

- **Write for Your Audience**: Keep your target audience in mind. Is this a guide for users, or for developers?
- **Be Clear and Direct**: Use simple language and avoid jargon where possible.
- **Provide Examples**: Code examples are often the best way to explain how something works.
- **Keep it Up-to-date**: When you make a change to the code, make sure you also update the corresponding documentation.
- **Cross-Reference**: Link to other relevant parts of the documentation.
- **Review Your Writing**: Read your documentation out loud to catch awkward phrasing and typos.
