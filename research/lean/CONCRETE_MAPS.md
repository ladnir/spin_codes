# Concrete IMT maps and spectrum replay

> **Historical snapshot — archived 2026-09-28.** Status statements and task recommendations below describe an earlier stage.
> See [STATUS.md](STATUS.md) and [NATIVE_RESULT.md](NATIVE_RESULT.md) for the closed native distance/rate result and verification scope.

The target remains the full distance theorem. T3a proves a contraction for a
numerical matrix; this work connects that matrix to the maps in
`research/paper/structured_imt_appendix.tex`, Table `tab:imt-maps`.

## Map definitions and proved structure

`ConcreteMapData.lean` contains the table's exact 19 basis images of `A` and
`C` transpose. Bit zero is the least significant bit. `cRows` transposes the
listed feedback rows into the 128 images of input basis vectors.

`PackedMap.lean` proves that the evaluator preserves XOR. It proves the
composition and inverse rules, a bijection between packed words and finite
supports, and equality between packed Hamming weight and support cardinality.
All arguments apply to every input; enumeration is not used to assume
linearity.

`ConcreteMaps.lean` defines the actual maps on packed words and on the finite
supports used by the transfer development. It proves:

- injectivity of `A` and `Ctranspose`, using checked left inverses;
- surjectivity of `C`, using a checked right inverse;
- preservation of XOR and of the empty support;
- weight five, nonzero value, and distinctness of the 128 actual images
  `C (inputBasis i)`.

`ConcreteMapPin.lean` supplements the original statement pins. The 13 finite
checks include row bounds, inverse basis identities, the nonzero-column test
for the expansion map, and the feedback column checks.

## Concrete Fourier and fiber bounds

`ConcretePairing.lean` proves that the actual `Cset` and `CtransposeSet` are
adjoint for the binary character pairing. Their basis incidences are the
transpose of the same checked table. It identifies the orthogonal code with
the actual kernel of `Cset` and identifies the abstract fibers with actual
syndrome fibers.

`ConcreteFourier.lean` instantiates Fourier inversion and Parseval for this
map. `ConcreteWeightSums.lean` groups the resulting sums by output weight
and equates the support counts with the packed counts used by the spectrum
replay. `ConcreteVariance.lean` proves that the squared fiber sizes sum to
the equal-syndrome pair count, then applies Cauchy–Schwarz over all `2^19-1`
nonzero syndromes, including empty fibers.

`ConcreteKernel.lean` proves that every nonzero kernel word has weight at
least four. Each feedback column has odd weight five, so kernel words have
even weight; distinct columns exclude weight two. The resulting packing
bounds apply both to a fiber and to its complements.

`ConcreteFiberBounds.lean` provides sufficient integer checks for the four
bounds in the paper: the non-kernel layer size, absolute Fourier sum,
constant-weight packing, and variance. `FiberCertificate.lean` justifies
checking the variance quadratic at `cap + 1`; this avoids trusting a rounded
square-root computation. `ConcreteCounts.lean` now connects the frozen
numerical caps to these theorems about the actual map.

Run `python -X utf8 -B scripts/check-map-fibers.py` for the sequential module,
statement-pin, and axiom audit. Its report is `scripts/map_data/fiber_verification.json`.
Only a `PASS` report denotes a completed audit.

`scripts/map_fiber_preview.py` prepares the next finite certificates from the
proposed spectrum totals. Its `UNTRUSTED_PREVIEW` report reproduces all 129
frozen kernel counts and finds a sufficient integer check for every frozen
cap: 125 Fourier checks, one packing check on each side, and two complement
checks at the endpoints. `fiber_candidates.json` contains the proposed pair
counts and absolute Fourier sums. The preview remains an untrusted artifact;
the separate Lean replay of all 129 rows has now passed.

`FiberNumericsDefs.lean` implements the exact Krawtchouk sum with the already
proved executable binomial coefficients. `FiberNumerics.lean` proves equality
with the mathematical sums, and `FiberNumericsBridge.lean` applies checked
rows to the actual kernel and fibers under an explicit spectrum identity.
`FiberFrozenData.lean` and `FiberFrozen.lean` identify the proposed arrays with
the arrays used by the numerical occupation matrix.

The separate `check-map-fiber-data.py` replay checks all 129 weights: each
module checks its Krawtchouk row, signed sum, squared sum, absolute sum, and
selected integer cap inequality. Its `kernel_fiber_data.json` report is `PASS`.
`FiberNumericsAll.lean` proves the universal count and cap theorems with an
explicit spectrum premise. `ConcreteCounts.lean` discharges that premise
using the actual transpose spectrum. Both modules now compile; the eleven
axiom audits report only `propext`, `Classical.choice`, and `Quot.sound`.
Run `python scripts/check-concrete-counts.py` to repeat the final assembly
and audit using the previously checked dependencies. Its report is
`scripts/map_data/concrete_counts_verification.json`.

The assembly uses finite elimination on an abstract predicate and separate
lemmas for the 129 caps. The endpoint arithmetic is rewritten in a general
lemma before the concrete indices are substituted. This avoids the kernel
memory failure of the earlier assembly. The local final check passed with
one worker and a 10,000 MB Lean memory limit.

`ConcreteShells.lean` defines the five actual weight shells of `Aset`. Their
positive weights exclude zero, their distinct weights make them disjoint,
and their five cardinalities imply coverage of all nonzero states. The
cardinalities are supplied by the complete `A` spectrum in `ConcreteCounts`.
Its `actualShellSystem` covers all nonzero states and has exactly the shell
sizes used by the numerical matrix.

The transfer induction must carry nonnegative budgets and keep diffuse mass
off zero. `TRANSFER_INVARIANT.md` explains the two gaps in the older interface
and the proved `LiveDominates` replacement used by the new sparse-moment
bridge. This bridge still requires the actual one-step law.

`check-map-connection.py` compiles the supplemental connection pins and
audits 16 proof closures. It records the current semantic-module hashes in
`connection_audit.json` without rebuilding numerical worker dependencies.

## Exact spectra

The generated spectrum has 512 blocks, each covering 1,024 consecutive packed
states. Both maps are evaluated on every state, including zero. The low ten
input bits use a shared kernel-checked table of images; each block supplies a
checked high image. `PackedMap.eval_split` proves this decomposition for the
original evaluator. A faster weight recurrence is proved equal to the
original bit-based weight.

`MapSpectrum.lean` proves that each checked block histogram counts its actual
inputs, that the blocks cover the entire input range once, and that the final
histogram counts the corresponding `Fin (2^19)` vectors. `MapSpectrumSum.lean`
justifies combining the block histograms.

The full replay of all 512 blocks has passed; see
`scripts/map_data/kernel_blocks.json`. The generator's Python totals alone
are not Lean evidence.

After all blocks pass, `map_spectrum_assemble.py` checks source hashes and
emits `MapSpectrumTotals.lean` and `MapSpectrumBridge.lean`. Both have passed
Lean checking and combine the numerical totals with the actual-map counting proof.

## Replay

Finish dependency builds before starting numerical workers. Do not run Lake
or rebuild an imported numerical module while the workers are active.

The Fourier-row workers additionally depend on `PolyIdentityDefs`,
`FiberNumericsDefs`, and `FiberNumericsData/Tables`; keep their object files
unchanged until that replay finishes. Its process/session locator is
`scripts/map_data/fiber_replay_session.json`.

```text
python -X utf8 -B scripts/concrete_maps.py
python -X utf8 -B scripts/check-concrete-maps.py
python -X utf8 -B scripts/map_spectrum_data.py
python -X utf8 -B scripts/check-map-spectrum.py --basis --jobs 1
python -X utf8 -B scripts/check-map-spectrum.py --jobs 3
python -X utf8 -B scripts/map_spectrum_assemble.py
```

The structural report is `scripts/map_data/structural_verification.json`.
Numerical reports record source hashes, dependency-object hashes, exit codes,
timings, and individual logs. Preserve a live replay across goal turns; poll
its existing process/session instead of starting a duplicate.

The expensive block checks use `decide +kernel`: Lean checks their proof in
the kernel without first repeating the calculation in the elaborator. The
sample axiom audit reports only `propext`. The earlier plain-`decide` run was
deliberately stopped after profiling this change; its partial report is
retained as `kernel_blocks_plain_decide.json`. Those earlier source hashes
are not reused as passes for the regenerated blocks. `--resume` reuses only
successful checks with matching source, output-object, and dependency hashes.

## Peach verification workspace

The user-authorized remote workspace is `/tmp/spin-lean-peach/project`.
The toolchain is `/tmp/spin-lean-peach/toolchain/lean-4.34.0-linux/bin`.
`prepare-peach-counts.py` packages the local checked dependency closure and
records source/object hashes in `peach_dependencies.json`. Project proof
objects are reused; mathlib comes from the pinned Linux cache/build.
`check-peach-counts.py` verifies those hashes, recompiles the final assembly,
and audits the eleven axiom closures. It does not rerun every certificate.

The first completed remote check passed: 6.00 seconds for `FiberNumericsAll`,
4.00 seconds for `ConcreteCounts`, and 3.75 seconds for the axiom audit. Peak
resident memory was 6.98 GiB. The source hashes match the local checked
sources. The downloaded report is `scripts/map_data/peach_counts_verification.json`.

Connect with the configured `connect-peach.ps1` wrapper. In the remote
project directory, put the toolchain above on `PATH` and run
`python3 scripts/check-peach-counts.py`. The installed WinSCP rejected the
configured SHA256 host-key syntax, so transfers used binary streams over
Plink with the exact pinned host-key fingerprint.

## Connections remaining at this checkpoint

The actual one-step domination law now covers every input weight. Both
special cancellation tables are identified with the maps; see
`CONCRETE_STEP.md` and `LOW_CANCELLATION.md`. The Bernoulli mixture is
identified with the direct input kernel and connected to the sparse iterate
bound. At this checkpoint, the concrete family, route linkage, selection
estimates, and occupation asymptotics were still required for the full
distance theorem.
