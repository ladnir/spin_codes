# Expander-code enumerator project

This directory contains the new, unified project on enumerator bounds for
Expand--Accumulate (EA) and Expand--Convolute (EC) codes.

The manuscript separates two cases.  The binary case targets bit-valued
applications such as Silent OT.  The finite-field case studies genuine field
mixing and emphasizes field orders near `2^128`, as needed for Silent VOLE.

## Layout

- `paper/` contains the unified submission manuscript.
- `notes/previous_draft/` preserves the manuscript sources and research notes
  formerly stored in `enumerator_paper/`.  It is source material, not the
  active manuscript.
- `scripts/` contains the diagnostic and outward-rounded verification code.
- `results/` contains certificates and recorded experiment outputs.
- `certificate_manifest.json` maps each finite claim to its verifier, inputs,
  hashes, precision, and expected result.
- `SUBMISSION.md` records the venue-independent submission and artifact
  checklist.
- `HANDOFF_ASTRA.md` records the mathematical, artifact, implementation, and
  release context for the independent completion review.

Install the tested Python packages and check the artifact package with:

```powershell
python -m pip install -r requirements.txt
python scripts/verify_manifest.py --strict-versions --list
```

Run `python scripts/verify_manifest.py --run headline` for the principal
certificates or use `--run all` for every rigorous numerical claim.  The
driver executes commands serially.

The large-field certificates include two rate-one-half EC results over
`F_(2^128)` at 98.415% of the field's GV distance.  Left/right degrees 24/12
with convolution memory 5 give failure probability below `2^-45.03`.
Degrees 26/13 with memory 4 give failure probability below `2^-135.20`.
These are labeled finite-field ensembles, not binary matrices reinterpreted
over `F_(2^128)`.

At a fixed output-weight cutoff, the degree-24/12 certificate retains its
45.03-bit diagnostic margin as the field order changes near `2^128`.  If the
cutoff instead tracks a fixed fraction of each field's GV distance, the
requested distance increases with the field order and the margin declines.
Outward-rounded certificates give 45.84, 45.03, 44.24, and 38.78 failure bits
at field orders `2^127`, `2^128`, `2^129`, and `2^136`, respectively.  The
last row uses locally retuned markers; retuning gains only 0.093 bits.  The
full frozen-marker curve through `2^144` is recorded in
`results/finite_field_biregular_ec_gf128_rate_half_d24_m5_field_size_diagnostic.json`.
