# Security Policy

## Supported Versions

KPP is a source-code generator run locally as a development-time tool;
there is no hosted service to secure. Security fixes are applied to
the latest release on the `main` branch. We do not backport fixes to
older released versions.

## Reporting a Vulnerability

If you discover a security vulnerability in KPP (for example, in the
scanner/parser in `src/`, which processes user-supplied `.kpp`/`.def`/
`.eqn`/`.spc` input files), please report it privately rather than
opening a public GitHub issue.

**Preferred method:** Use GitHub's private vulnerability reporting
feature for this repository:
[https://github.com/KineticPreProcessor/KPP/security/advisories/new](https://github.com/KineticPreProcessor/KPP/security/advisories/new)

**Alternative:** Contact a maintainer listed in
[`AUTHORS.txt`](AUTHORS.txt) directly.

Please include:
- A description of the vulnerability and its potential impact
- Steps to reproduce it (a minimal `.kpp`/mechanism file that triggers
  the issue is ideal)
- The KPP version or commit affected

## What to Expect

We will acknowledge receipt of your report, investigate, and work with
you to understand and address the issue. We will credit reporters
(unless they prefer to remain anonymous) when a fix is released.
