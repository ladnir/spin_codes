# Artifact validation: 2026-09-07

This records local validation, not a published release or a fresh complete
numerical proof replay.

## Checks performed

- `reproduce.py quick`: passed. It checks plot data, five finite ledgers,
  selected-map algebra and spectra, manuscript timing values, and all 31
  asymptotic manifest entries.
- Eight artifact regression tests: passed. One constructs an isolated compact
  source tree in a temporary directory, with no Git metadata, bulk receipts,
  or references back to the development checkout, and runs the quick checks.
  Other tests cover path containment, modified manifests, documentation links,
  missing evidence, and refusal to overwrite an existing archive.
- `reproduce.py evidence`: 788 authenticated BCH/performance pins, zero
  missing dependencies and zero mismatches in the development checkout.
  Of these, 515 are outside Git and belong in the evidence archive.
- Repository hygiene check: passed for tracked files. Newly restored compact
  inputs are under 5 MiB; the evidence ZIP stays under ignored `output/`.
- A Windows/Linux Python 3.11 CI workflow has been added for quick checks and
  regression tests. It has not yet run on GitHub; local validation used Windows
  and Python 3.14.

## Restored compact evidence

The asymptotic manifest and `golay_ba3_rm2sub_joint_interval_d11.json` were
restored from the earlier paper worktree, with their recorded hashes checked.
The historical selected-map receipt is copied under `artifact/data/asymptotic/`
with the original hash `e2658fefdb787d5946c4c0a1dc4c8c530414fae1f6ef8bcf9c40298e6e13f117`.
Twenty-four existing files needed CRLF-to-LF correction to match the manifest.
Every corrected byte sequence matches its original hash; no mathematical or
algorithmic source content was changed. Git attributes now preserve these
conventions across checkouts. Frozen BCH pins still pass unchanged.

## Local evidence archive

The paper was rebuilt with clickable repository references and visually
checked at the edited passages and parameter/theorem pages. The final TeX
log has no warnings or unresolved references. PDF: 48 pages, SHA-256
`6b648fd31e9f68964f49f60f4b70da06246861125f8588220d33d96618f06167`.

File: `output/artifact/bch-evidence-2026-09-07.zip`.

- 788 pinned files plus `evidence-manifest.json`.
- Uncompressed pinned inputs: 242,557,416 bytes.
- Archive size: 66,545,176 bytes (about 63.5 MiB).
- SHA-256: `c1f51575155893203df70c33b01d1004a3b95c6e35069c5855956c9c6e7f1eb6`.
- Every archived file was reopened and checked against its pinned SHA-256.

This archive is a local release candidate, not a public download. It contains
retained BCH evidence, not the small-BCH raw grid or a self-contained numerical
runtime. Source code in the companion Git revision remains necessary.

## Next release step

Review and commit the compact artifact and manuscript, then test the numerical
replay paths in a fresh environment. Publish the evidence archive separately
and pin the final Git revision and archive checksum in the paper. The current
paper hyperlinks use `main` and will become useful after integration; they
are not yet immutable release references.
