# Reproducing the paper

Run commands from the repository root unless stated otherwise. Use `python3`
in place of `python` if needed. The quick path uses Python 3.11+ and no
third-party Python packages. It is intended for both Windows and Linux;
Windows is the locally tested platform for this revision.

On Windows, enable long paths for the clone because some frozen source paths
are long: `git clone -c core.longpaths=true https://github.com/ladnir/permute_conv.git`.
The Windows CI job enables this setting before checkout as well.

## 1. Check retained results

```sh
python -B artifact/reproduce.py quick
python -B -m unittest discover -s artifact -p 'test_*.py'
```

Expected: `QUICK CHECK PASSED`, including all 130 small-BCH geometries,
five finite certificate margins, 19 selected generators, all 524,288 inner
states, the exact kernel spectrum, six timing entries, and 31 asymptotic
manifest entries. This should be a short check on an ordinary desktop;
it performs no parameter search, interval replay, or timing experiment.

The asymptotic checker uses the original manifest and one explicit relocation:
the historical `receipts/min_state/s19_rm2sub_selection.json` input is retained
as [data/asymptotic/s19_rm2sub_selection.json](data/asymptotic/s19_rm2sub_selection.json).
Its original hash is unchanged. Git attributes retain the manifest's LF bytes
on Windows; other frozen BCH sources keep their existing byte conventions.

## 2. Rebuild plots and paper

Install TeX Live with latexmk, BibTeX, PGFPlots, and placeins. Then run:

```sh
python -B artifact/reproduce.py figures
python -B artifact/reproduce.py paper
```

The first command rewrites only the three generated TeX figure inputs.
The second runs quick checks and builds `output/pdf/spin_codes_draft.pdf`.
The PDF uses the build date, so its bytes need not match an earlier build.
No plotting library is needed: the figures are native vector PGFPlots.

Figures 1--3 are regenerated from the tracked rounded results table. Full
margins retain six decimal places; Q1 overlays add the reported loss only
where both values are present. Rebuilding them is not a fresh evaluation of
the original grid. The table records fingerprints of the raw exports, which
are not included in this checkout. The five BCH-256 curve points instead
come from exact rational ledgers checked against the manuscript.

## 3. Inventory and authenticate BCH evidence

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

For asymptotic Structured SPIN, start with the [certificate README](../workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/README.md).
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
`python -B workstreams/transposed_comparison/report.py --check`.
This check performs no measurements and is not yet included in `quick`.

## Release preparation

Before an external artifact submission:

1. Commit and review the compact artifact, then record its immutable revision
   in the paper. The current GitHub link is a repository entry point, not a
   versioned release claim.
2. Publish the missing BCH numerical evidence separately, with a file manifest,
   sizes, checksums, and an archive checksum. Do not commit raw worker trees.
3. Restore the small-BCH raw grid inputs and exports if claiming reproduction
   of the underlying sweep, rather than reproduction of its plotted results.
4. Test numerical replay in a fresh environment, resolving historical paths
   through explicit adapters instead of modifying frozen producers.
5. Record measured runtime, peak memory, machine, dependencies, and expected
   final results for each replay level. Re-run benchmarks separately.

No archive upload, Git push, or fresh benchmark was performed by this cleanup.
