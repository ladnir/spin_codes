# Reproducing the paper

Historical author-side instructions, not artifact requirements. The artifact
now covers only the [core implementation](README.md). Commands below may
require unpublished local research inputs; no complete replay bundle is promised.

Run commands from the repository root unless stated otherwise. Use `python3`
in place of `python` if needed. The quick path uses Python 3.11+ and no
third-party Python packages. It is intended for both Windows and Linux;
Windows is the locally tested platform for this revision.

On Windows, enable long paths for the clone because some frozen source paths
are long: `git clone -c core.longpaths=true https://github.com/ladnir/permute_conv.git`.
The Windows CI job enables this setting before checkout as well.

Automatic **Paper artifact source checks** run packaging, missing-evidence
rejection, exact-union adapter, and terminology tests on Windows and Linux.
They do not validate the paper's numerical certificates: those checks require
the separately supplied pinned evidence and use the commands below. A green
source-check job is not a certificate-replay result.

## 1. Check retained results

The selected finite results now use IMT. Supply their pinned evidence before
running these checks; a source-only checkout does not contain the bulk receipts.

```sh
python -B artifact/imt_reproduce.py check
python -B artifact/reproduce.py quick
python -B -m unittest discover -s artifact -p 'test_*.py'
```

The first command reports `SELECTED_FINITE_IMT_INTEGRATION_PASSED`: seven
certificate targets, four timing cells, 57 map words, and exact union checks.
The historical `quick` wrapper also checks the 130 retained RM2Sub geometries
and 31 imported asymptotic manifest entries. Its old progress labels are not
an IMT inventory. Neither command performs parameter search, interval replay,
or a timing experiment. The current finite check also authenticates the 130-cell
IMT Q1 grid and checks its three figures plus five matched Q1/full anchors.

Check the current asymptotic IMT theorem separately with
`python -B paper/check_imt_integration.py`; this requires the numerical
dependencies below and the separate local asymptotic IMT receipts.

The asymptotic checker uses the original manifest and one explicit relocation:
the historical `receipts/min_state/s19_rm2sub_selection.json` input is retained
as [data/asymptotic/s19_rm2sub_selection.json](data/asymptotic/s19_rm2sub_selection.json).
Its original hash is unchanged. Git attributes retain the manifest's LF bytes
on Windows; other frozen BCH sources keep their existing byte conventions.

## 2. Rebuild plots and paper

Install TeX Live with latexmk, BibTeX, PGFPlots, and placeins. Then run:

```sh
python -B artifact/imt_reproduce.py figures
python -B artifact/reproduce.py paper
```

The first command rewrites only the four generated IMT TeX figure inputs.
The second runs quick checks and builds `output/pdf/spin_codes_draft.pdf`.
The PDF uses the build date, so its bytes need not match an earlier build.
No plotting library is needed: the figures are native vector PGFPlots.

The three parameter slices come from `PARAMETER_NO_CONSTANT_Q1_v1.json` and
its authenticated replay receipt. The replay checked all 130 geometries and
39 exact map pairs, including seven separate log-domain checks. Regeneration
authenticates that evidence; it does not perform a fresh numerical evaluation.
The five BCH-256 Q1/full pairs come from the selected full-certificate ledgers
and use identical maps within each pair. The old `build_parameter_figures.py`
and its rounded RM2Sub tables remain available only for historical reproduction.

## 3. Inventory and authenticate evidence

For the selected finite IMT results:

```sh
python -B artifact/imt_reproduce.py inventory
python -B artifact/imt_reproduce.py inventory --output output/artifact/imt-inventory.json
python -B artifact/imt_reproduce.py pack --output output/artifact/imt-evidence.zip
```

Use fresh output names. Inventory fails on missing or mismatched files;
an unreadable root receipt also marks discovery incomplete. The accepted
local set has 753 files and 640,779,033 uncompressed bytes. The ZIP contains
only the eight accepted root receipts, their declared source pins, and an
inventory manifest. Packing streams the files, reopens the archive, and
verifies its member list and every SHA-256. It never overwrites an archive.
No production IMT archive has been created or published by this step.

This scope covers the selected finite certificates and timing bindings.
It does not include a runtime, installed dependencies, the separate asymptotic
IMT evidence, or the diagnostic Q1 grid by default. To include the current
parameter-figure inputs, add `--include-q1` to `inventory` or `pack`. The
expanded local inventory has ten root receipts, 766 files, and 641,495,560
uncompressed bytes, with no missing or mismatched files. It still does not
certify every grid cell. File authentication is not an interval replay or an
independent review of the proof.

The commands below retain their **historical BCH/RM2Sub** scope. They are
not substitutes for the current IMT inventory or package:

```sh
python -B artifact/reproduce.py inventory
python -B artifact/reproduce.py evidence
```

`inventory` emits JSON with every dependency's expected SHA-256, local byte
size, Git status, and missing/mismatch status. Missing files are reported
without making inventory fail; hash mismatches fail. `evidence` requires the
complete local set and fails if any dependency is absent or altered.

This inventory covers the retained BCH ladder audit and performance pins,
not every experiment in the repository. The current audit has 788 pins;
515 were outside Git before artifact packaging. A clean source checkout is
expected to lack those bulk dependencies. Do not treat a source-only check
as a full numerical replay.

When the complete evidence is available, package it without sweeping up
unrelated experiments:

```sh
python -B artifact/reproduce.py pack-evidence --output output/artifact/bch-evidence.zip
```

The ZIP contains only the 788 pinned files and an inventory manifest, using
repository-relative paths. The command verifies every archived SHA-256 and
prints the archive checksum. It does not overwrite an existing archive.
If interrupted, discard the incomplete output only after inspecting it and
retry with a fresh name. Archive timestamps can change the ZIP checksum.
Before restoring into a checkout, compare existing files with the inventory;
do not blindly overwrite a modified working tree.

## 4. Replay numerical proofs

Install numerical dependencies in a separate virtual environment. The local
finite-proof environment uses Python 3.14, python-flint 0.9.0, NumPy 2.4.4,
SciPy 1.18.0, and mpmath 1.3.0. These are recorded working-environment versions,
not a claim that every replay command has been retested in a fresh environment.

For current IMT replay, start with the
[finite results ledger](../workstreams/inner_design/finite_migration/PAPER_RESULTS.md)
and [asymptotic IMT guide](../workstreams/inner_design/imt_asymptotic/README.md).
They identify the accepted inputs and separate authentication from arithmetic
replay. The procedures in the next paragraphs describe historical RM2Sub proofs.

For historical asymptotic Structured SPIN, start with the [certificate README](../workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/README.md).
It lists the outward-arithmetic producers and their order. Use new output
paths; preserve the accepted receipts and their hashes. Some historical
producers still resolve inputs through old repository locations. The quick
manifest check is portable; a fully relocated end-to-end producer workflow
is a remaining release task.

For BCH-256, first supply all files required by `evidence`. Replay instructions
are size-specific:

- [K=2^16](../workstreams/bch_rm2sub_bridge/T128_S19_M16_CLOSURE.md).
- [K=2^18](../workstreams/bch_rm2sub_bridge/T128_S19_M18_CLOSURE.md).
- [K=2^20](../workstreams/bch_rm2sub_bridge/T128_S19_M20_CLOSURE.md).
- [K=2^22 and 2^24](../workstreams/bch_rm2sub_bridge/T128_S19_M22_M24_CLOSURE.md).

Replace a historical absolute Python path with your environment's interpreter.
Use fresh output directories. Search budgets may yield partial coverage;
only a complete accepted ledger establishes the stated finite result.
The [completion auditor](../workstreams/bch_rm2sub_bridge/audit_ladder_completion.py)
freshly reevaluates dense boxes at 768 bits, but authenticates existing sparse
512-bit replays rather than rerunning them. There is no one-command full
proof regeneration claim in this artifact revision.

## 5. Encoder correctness and performance

The current IMT measurements and exact map bindings are listed in the
[finite results ledger](../workstreams/inner_design/finite_migration/PAPER_RESULTS.md).
Use `python -B paper/build_imt_comparison.py --check` to validate the current
comparison table. The build below and the original campaign generator are
retained historical RM2Sub procedures, not commands to reproduce IMT timings.

See the [encoder README](../workstreams/bare_bch_rm2sub/README.md) for the API,
build commands, ISA requirements, and measurement procedure. A representative
GCC/Linux build is:

```sh
cmake -S workstreams/bare_bch_rm2sub -B out/bare-spin -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPIN_ARCH=native
cmake --build out/bare-spin -j 3
ctest --test-dir out/bare-spin --output-on-failure
```

The x86 implementation needs the SIMD/carryless-multiply instructions enabled
by its CMake flags; it has no runtime ISA dispatch. `native` alone does not
make an older CPU support explicitly enabled instructions. Published timings
use a Ryzen 9 7950X and GCC 15.2.0, not an arbitrary CI runner.

Benchmarks are opt-in and must run **one at a time**, with no competing
benchmark process. Use Linux for the published affinity-controlled procedure.
Do not interpret a timing difference on another machine as a failed proof.
The artifact quick command only checks transcription and retained records.

For the new same-host comparison against chosen BAA, Expand--Convolute,
and binary RAA, follow the
[comparison guide](../workstreams/transposed_comparison/README.md).
Its pinned dependency build and serial runner are separate from the frozen
SPIN benchmark. Validate its retained samples and generated paper table with
`python -B paper/build_imt_comparison.py --check`.
This check performs no measurements and is not yet included in `quick`.

## Release preparation

Before an external artifact submission:

1. Commit and review the compact artifact, then record its immutable revision
   in the paper. The current GitHub link is a repository entry point, not a
   versioned release claim.
2. Publish the selected IMT numerical evidence separately, with a file manifest,
   sizes, checksums, and an archive checksum. Do not commit raw worker trees.
3. Include the current Q1 inputs with `--include-q1` and test their producer
   and replay commands, not just figure regeneration. Restore the older small-BCH
   raw exports only if also releasing the historical RM2Sub sweep.
4. Test numerical replay in a fresh environment, resolving historical paths
   through explicit adapters instead of modifying frozen producers.
5. Record measured runtime, peak memory, machine, dependencies, and expected
   final results for each replay level. Re-run benchmarks separately.

No archive upload, Git push, or fresh benchmark was performed by this cleanup.
