# Pinned LNCS class

`llncs.cls` is the unmodified class from Springer's current
[LaTeX2e proceedings template](https://link.springer.com/series/558/information-for-authors-and-editors),
retrieved on 2026-09-17. Only line endings were normalized for the source tree.

- Publisher archive: https://cms-resources.apps.public.k8s.springernature.io/springer-cms/rest/v1/content/27851904/data/LaTeX2e%20Proceedings%20Template%20ZIP
- Class identification: `2026/09/03 v2.25`.
- SHA-256 of the LF-normalized local file:
  `04ec99046ad9867f6cd54ffd1ef8c47c767b327d1baec77536db0528cd492fa2`.

The previous local class identified itself as `2025/02/25 v2.26`.
The apparently lower version number above is the publisher's current file,
not a locally renamed or modified class. Build from `paper/` so LaTeX uses
this pinned copy; the log must identify `2026/09/03 v2.25`.
The publisher class also writes `DescriptionTexts.txt` to the build-output
directory. This generated file is not a research result or source artifact.

Both builds retain `\documentclass{llncs}` with default typography and
`\pagestyle{plain}`. Do not patch the class to save pages.
