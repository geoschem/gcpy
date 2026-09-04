# Security Policy

## Supported Versions

GCPy does not maintain long-term support branches. Security fixes are only made against the latest release and the `main`/`dev` development branches.

| Version         | Supported          |
| ---------------- | ------------------ |
| Latest release    | :white_check_mark: |
| Older releases     | :x:                 |

## Reporting a Vulnerability

GCPy is an analysis toolkit for GEOS-Chem model output; it does not run as a network service and does not handle authentication, secrets, or untrusted network input directly. Even so, if you discover a security vulnerability (for example, unsafe deserialization, arbitrary code execution when reading a data/config file, or a vulnerable dependency), please report it privately rather than opening a public GitHub issue by using [GitHub's private vulnerability reporting](https://github.com/geoschem/gcpy/security/advisories/new) for this repository.

Please do not disclose the issue publicly until the GCST has had a chance to investigate and, if needed, issue a fix.

## What to Expect

The GEOS-Chem Support Team will acknowledge your report, investigate, and work with you on a fix and disclosure timeline. Since GCPy is maintained by a small team, response times may vary, but security reports are treated as a high priority ahead of feature work. 

For non-security bugs and user questions, please see [SUPPORT.md](SUPPORT.md) and use [GitHub issues](https://github.com/geoschem/gcpy/issues/new/choose) instead.
