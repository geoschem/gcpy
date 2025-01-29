# Contributing Guidelines

Thank you for looking into contributing to GCPy! GCPy is a grass-roots package that relies on contributions from community members like you. Whether you're new to GCPy or a longtime user, you're a valued member of the community, and we want you to feel empowered to contribute.

## We use GitHub and ReadTheDocs

We use GitHub to host the GCPy source code, to track issues, user questions, and feature requests, and to accept pull requests: [https://github.com/geoschem/gcpy](https://github.com/geoschem/gcpy). Please help out as you can in response to issues and user questions.

GCPy documentation can be found at [gcpy.readthedocs.io](https://gcpy.readthedocs.io).

## When should I submit updates?

Submit bug fixes right away, as these will be given the highest priority.  Please see **[Support Guidelines](https://gcpy.readthedocs.io/en/stable/reference/SUPPORT.html)** for more information.

The practical aspects of submitting code updates are listed below.

## How can I submit updates?

We use **GitHub Flow**, so all changes happen through [pull requests](https://help.github.com/articles/creating-a-pull-request/). This workflow is [described here](https://docs.github.com/en/get-started/using-github/github-flow).

As the author you are responsible for:
- Testing your changes
- Updating the user documentation (if applicable)
- Supporting issues and questions related to your changes

### Process for submitting code updates

  1. Create or log into your [GitHub](https://github.com/) account.
  2. [Fork the GCPy repository](https://help.github.com/articles/fork-a-repo/) into your Github account.
  3. Clone your fork of the GCPy repositories to your computer system.
  4. Add your modifications into a [new branch](https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell) off the **main** branch.
  5. Add a sentence to the `CHANGELOG.md` file describing your update.
  6. Test your update thoroughly and make sure that it works.
  7. Review the coding conventions and checklists for code and data updates listed below.
  8. Create a [pull request in GitHub](https://help.github.com/articles/creating-a-pull-request/).
  9. The [GEOS-Chem Support Team](https://wiki.geos-chem.org/GEOS-Chem_Support_Team) will add your updates into the development branch for an upcoming GCPy version.
  10. If the benchmark simulations reveal a problem with your update, the GCST will request that you take further corrective action.

### Coding conventions

We recommend that GCPy developers adhere to the [PEP-8 Python style guide](https://peps.python.org/pep-0008/).  You can run `pylint` on all source code files that you modify to ensure adherence to PEP-8 style conventions.

### Checklist for submitting code updates

  1. Include thorough comments in all submitted code.
  2. Include full citations for references at the top of relevant source code modules.
  3. Remove extraneous code updates (e.g. testing options, other science).
  4. Check that you have updated the `CHANGELOG.md` file.
  5. Run `pylint` on each source code file that you have modified to check for bugs and conformance to the PEP-8 style conventions.

## How can I request a new feature?

We accept feature requests through issues on GitHub. To request a new feature, **[open a new issue](https://github.com/geoschem/gcpy/issues/new/choose)** and select the feature request template. Please include all the information that migth be relevant, including the motivation for the feature.

## How can I report a bug?

Please see **[Support Guidelines](https://gcpy.readthedocs.io/en/stable/reference/SUPPORT.html)**.

## Where can I ask for help?

Please see **[Support Guidelines](https://gcpy.readthedocs.io/en/stable/reference/SUPPORT.html)**
