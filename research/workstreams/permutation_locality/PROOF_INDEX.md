# Packet-Route Proof Index

The current packed GL32 **R4** construction has a **complete >10% distance /
>52.05-bit whole-code certificate** at K=2^20, rate one-half. Its optimized
precomputed transposed encoder now takes **5.334629 ms** on 128-bit XOR
elements (about 5.32--5.35 ms in the latest matched confirmation).
Fresh 384-bit replay covers all q=1..32 sparse occupancies and an exhaustive
337-interval dense cover for q=33..2048. Their exact sum gives
52.052884355478 bits of setup-failure margin. See the
[R4 closure record](packed_mixing/R4_CLOSURE.md) and
[matched GFNI inner performance report](gfni_inner_REPORT.md).

The [hill-climb ledger](packed_mixing/HILL_CLIMB.md) separates counting
improvements from the change to four actual inner updates. Exact-map fusion
then removes 8.01% of sequential R4's measured time without changing its
distribution. Subsequent wide stores and compact GL32 coefficients remove
4.91% in a separate matched comparison and halve coefficient storage.
Keeping the same 19-state update packed for GFNI removes another 9.41%
against its matched scalar-fused control, without changing the encoder.
The separate [16-bit-state construction](packed_mixing/s16_closure/CLOSURE.md)
is also **closed at >10% distance / >55.08 bits of whole-code margin**.
It uses 64-bit physical steps with independent uniform GL16 updates.
Fresh 384-bit replay covers all q=1..32 sparse occupancies and 106 dense
intervals for q=33..2048; the aggregate margin is 55.089909760705 bits.
Matched confirmation measures 5.344505 ms versus 5.324117 ms for the
retained R4 code: essentially equal performance, not a speedup claim.
The earlier 128-bit-step screens remain in its
[progress ledger](packed_mixing/s16_closure/PROGRESS.md). The new proof
does not replace the retained 19-bit-state certificate or implementation.
The [encoder cost profile](encoder_cost_REPORT.md) attributes about 59%
of t64 time to packed GL32/BCH work. A matched footprint diagnostic finds
both substantial compute/layout work and a cache-related penalty; it
changes neither construction nor certificate.
The old **R2** construction remains certified at **>9.5% /
>36.68 bits**, with its retained **6.203224 ms** measurement; see its
[closure record](packed_mixing/FIRST_CLOSURE.md). R2's 10% sparse bound is
still partial. No production default or paper claim changes here.

The completed lower-distance investigation uses `shared4-gf16-r2`, whose measured
precomputed transpose time is 5.863 ms. The user accepts a 20-bit whole-code
setup-failure margin. See the [shared-route follow-on](gf16_packets/SHARED_SHUFFLE.md#lower-distance-search-2026-09-29).
Fresh 384-bit whole-code replays now certify **>7% relative distance** with
**>44.61 bits of margin**, and **>5% with >73.38 bits** for the same encoder.
Both cover q=1..32 directly and q=33..2048 through exhaustive dense intervals
(303 and 68, respectively). Separate scope/hash/aggregation audits pass.
The earlier >6% / >40.68-bit receipt is retained, but is dominated by the
new 7% certificate. Higher-distance 7.5%/8% screens remain partial.
The earlier higher-distance and four-update records remain separate below.

A [two-route exploration](TWO_ROUTE_EXPLORATION.md) on 2026-09-30 tested
fused BCH loads and a randomized local pair mixer. Fused loads did not
improve the certified independent-row implementation. The new mixer takes
7.315 ms with two updates; two dense diagnostics pass at 9.5% but fail at
10%. It has no complete certificate and does not replace any result below.

The [BCH implementation comparison](bch_compare_REPORT.md) measures the same
outer map with XOR circuits, GFNI tiles, and a new polynomial prototype.
GFNI remains fastest on the tested 128-bit-element workload. These are
implementation comparisons, not new distance certificates; the report
separately labels the alternate-message-basis polynomial experiment.

The [packed-mixing exploration](packed_mixing/README.md) keeps data in the
GFNI layout for new local random maps. Full GL32 block mixing replaces the
old packet label randomizers. Its retained R2 [minimal optimized driver](packed_driver_REPORT.md)
measures **6.203 ms**, versus **5.805 ms** for the same-build legacy
encoder. Routing-callback outlining explained a substantial part of the
earlier build sensitivity. This is a 6.86% cost for a stronger construction,
not a speedup over optimized legacy; the earlier timings remain documented.
Exact local distributions, canonical-block shortening counts, expected-CDF
transport, and tighter pointwise expected-shell bounds are implemented and
tested. The new sparse checker includes the full support range down to five,
without changing the old checker's domain. Fresh 384-bit replay closes
q=1..32 with 37.95969 bits of aggregate sparse margin and q=33..2048 with
37.45010 bits across **240 intervals**. Their exact sum gives
**36.68251 bits for the whole code**. Six narrow intervals needed one further
subdivision; all their children passed without changing the encoder. The
[packed-mixing progress ledger](packed_mixing/PROGRESS.md) records the exact
scope and historical receipts. This completed the earlier 9.5% / 20-bit
target. R2's previous point bound at mean 0.032 still fails at 10%; the
new whole-code 10% result above belongs to R4.

This index separates the retained constructions from the current experiment.
No old proof is superseded by changing the coordinate-shuffle distribution.
The sources, historical notes, and selected numerical witnesses are preserved;
the local snapshot below also captures the uncommitted sources.

The shared-GF16 closure investigation has a [progress ledger](gf16_packets/PROGRESS.md)
with quantitative checkpoints and an explicit plateau criterion. That earlier
four-update 10% investigation stopped at its plateau on 2026-09-29; its target
remains open. It is distinct from the completed lower-distance result above.

The follow-on is [two independently shuffled row pairs](gf16_packets/PAIRWISE_SHUFFLE.md).
It changes the coordinate-shuffle distribution; it does not inherit an old
distance certificate. Its separate ledger records exact counting results,
coverage, and matched implementation timings. It stopped at a documented
plateau on 2026-09-29; the full >10% / 40-bit target remains open.

The pairwise track's union-based comparison now verifies an interval
around a previously failing dense diagnostic at 384-bit precision. Its
full dense cover is still incomplete; the restricted sparse margin below
must not be read as a whole-code guarantee.
The latest cover retains 95 verified cells. Wider probes expose failed
singleton bounds at comparison means .128 and .192. Wider regional
tuning, pooled comparison selection, tighter nearby shell counts, and
finer variance partitions did not resolve them. A fresh 384-bit replay
confirms the failed mean-.128 bound numerically.

## Scope and Status

All packet rows below use BCH[256,128], K=2^20, N=2^21, IMT(128,19),
zero initial state, and no terminal flush. R is the number of IMT updates
per step. A group contains either two or four BCH rows. Occupancy q counts
nonzero groups, not nonzero individual rows. Every stated margin is a
bound on bad-setup probability, over the specified ideal setup distribution
sampled once and used for all messages. It is not an estimate of actual
minimum distance.

**Complete** means all nonzero messages are covered. **Partial** means only
the stated message classes are covered; its margin is not a whole-code
guarantee. A relative-distance claim is strict where marked `>`.

| Stable ID | Coordinate shuffle and packet mixing | R | Retained result | Status / entry point |
|---|---|---:|---|---|
| `canonical-gl32-width8-shared4-r4` | Uniform GL32 per canonical four-row/eight-column block; shared group shuffle; no GF16 stage | 4 | >10% distance / >52.05-bit margin; 5.854159 ms | **Complete, fresh 384-bit whole replay:** [R4 closure](packed_mixing/R4_CLOSURE.md) |
| `canonical-gl32-width8-shared4-r2` | Uniform GL32 per canonical four-row/eight-column block; shared group shuffle; no GF16 stage | 2 | >9.5% distance / >36.68-bit margin; 6.203224 ms | **Complete, fresh 384-bit whole replay:** [closure record](packed_mixing/FIRST_CLOSURE.md) |
| `shared4-lanes-r1` | Shared within four-row group; lane permutation | 1 | 10% target, q=1–13; >48.15-bit restricted margin | Partial: [ONE_COLUMN](ONE_COLUMN.md) |
| `shared4-lanes-r2` | Shared within four-row group; lane permutation | 2 | 10% target, q=1–14; >41.83-bit restricted margin | Partial: [TWO_UPDATES](TWO_UPDATES.md), [joint cancellation](JOINT_CANCELLATION.md) |
| `independent4-lanes-r2` | Independent per row; lane permutation | 2 | >0.5% distance; >43.74-bit margin | Complete: [first closure](independent_rows/dense_closure/FIRST_CLOSURE.md) |
| `independent2-lanes-r2` | Independent per row; two-bit packets | 2 | >9.25% distance; >49.11-bit margin | Complete: [first closure](two_bit/FIRST_CLOSURE.md), [sparse coverage](two_bit/COVER_STATUS.md) |
| `independent4-gf16-r2` | Independent per row; independent nonzero GF16 multipliers | 2 | >9% distance; >46.85-bit margin | Complete: [GF16 certificates](gf16_packets/README.md#complete-certificates-2026-09-29) |
| `independent4-gf16-r3` | Same independent-row GF16 ensemble | 3 | 9.9% target, q=1–48; selected further occupancies | Partial: [GF16 history](gf16_packets/README.md) |
| `independent4-gf16-r4` | Same independent-row GF16 ensemble | 4 | >10% distance; >48.65-bit margin | Complete: [GF16 record](gf16_packets/README.md) |
| `shared4-gf16-r2` | Shared within four-row group; independent nonzero GF16 multipliers | 2 | >7% distance / >44.61-bit margin; also >5% / >73.38 bits; 5.863 ms retained timing | **Complete, fresh 384-bit replays:** [lower-distance closure](gf16_packets/SHARED_SHUFFLE.md#lower-distance-search-2026-09-29) |
| `shared4-gf16-r4` | Shared within four-row group; independent nonzero GF16 multipliers | 4 | 10% target, q=1–22; >44.20-bit best restricted margin, >44.12 with independently regenerated q=22 | **Partial; search plateau:** [progress ledger](gf16_packets/PROGRESS.md) |
| `pairwise4-gf16-r4` | Independent shuffle per row pair; independent nonzero GF16 multipliers | 4 | >10% target: q=1–48 at 384-bit precision, q=49 at 256 bits; combined restricted margin >48.401 bits | **Partial; search plateau:** [construction and ledger](gf16_packets/PAIRWISE_SHUFFLE.md). Occupancies 50–2048 are not fully covered. |

The production one-bit route and its finite-size certificates are separate
from these packet experiments; see the [paper's finite certificates](../../paper/finite_certificates.tex)
and [finite appendix](../../paper/finite_appendix.tex). Neither production
defaults nor paper claims change with this exploration.

The two-bit 9.5% search and independent four-bit lane-only attempts at
higher distances remain partial. Their failed upper bounds do not prove
bad codewords exist. The earlier [two-track status](PROOF_TRACKS.md) is
retained as a historical checkpoint, not the current index.

## Why the New Proof Is Separate

For `shared4-gf16-r4`, a fixed four-row message tuple has a union support
of size u. The shared uniform coordinate shuffle places that support
uniformly among the u-subsets of 256 regions. Independent nonzero GF16
multipliers make its nonzero packet values independent and uniform.
Thus the conditional inner analysis transfers, but the number of message
tuples at each union support must be bounded afresh. The independent-row
outer mixture and its completed dense certificate do **not** transfer.

The shared-shuffle implementation measures 6.532 ms versus 7.506 ms for the certified
independent-row four-update control in the same precomputed transpose run.
The approximately 5.44 ms shared lane-only control is a different, still
partially proved construction. See the shared-shuffle record for timing scope.

## Retained Numerical Witnesses

Paths in this section are repository-relative. Numerical witnesses remain
local and outside version control. These are selected full witnesses and
their assembly receipts, not all exploratory output.

| Construction | Principal witness or receipt under `tmp/` | Verification entry point |
|---|---|---|
| Canonical GL32 R4, >10% / >52.05 bits | `packed-hill/r4-d10-whole-parallel-p384/whole.json`; search `packed-hill/dense-r4-ekr-d10-complete-p256.json` and `packed-hill/sparse-r4-ekr-d10-p256.json` | `packed_mixing/whole_hill_replay.py`; [complete replay command and hashes](packed_mixing/R4_CLOSURE.md) |
| Canonical GL32 R2, >9.5% / >36.68 bits | `packed-closure/gl32-d095-whole-p384/whole.json` | `packed_mixing/whole_replay.py`; [retained replay command and hashes](packed_mixing/FIRST_CLOSURE.md) |
| Independent four-bit lanes, >0.5% | `four-bit-dense-adaptive-half-final.json` | `independent_rows/dense_closure/assemble.py` |
| Independent two-bit, >9.25% | `two-bit-affine-dense-401-d0925.json` | `two_bit/closure.py` |
| Independent GF16 R2, >9% | `gf16-r2-d09-independent-p384-complete.json`; dense `gf16-birth-classes-d09-tilt3over16-q97.json` | `gf16_packets/assemble.py` |
| Independent GF16 R4, >10% | `gf16-r4-d10-assembly-p384-complete.json`; dense `gf16-lowthreads-r4-d10-q49.json` | `gf16_packets/assemble.py` |
| Independent GF16 R4, separate sparse replay | `gf16-r4-d10-q1-48-prefix-p384.json` | `gf16_packets/assemble.py:regenerate_sparse` |
| Shared GF16 R2, >5% / >73.38 bits | `shared-relaxed/shared-r2-d05-complete-p384.json`; dense `shared-relaxed/dense-d05-parallel.json` | `gf16_packets/shared_relaxed_strategy.py`; scope/sum check `gf16_packets/shared_relaxed_audit.py` |
| Shared GF16 R2, >6% / >40.68 bits | `shared-relaxed/shared-r2-d06-complete-p384.json`; dense `shared-relaxed/dense-pruned-d06-parallel.json` | Same fresh assembler and independent scope/sum check; [reproduction](gf16_packets/SHARED_SHUFFLE.md#lower-distance-search-2026-09-29) |
| Shared GF16 R2, >7% / >44.61 bits | `shared-relaxed/shared-r2-d07-complete-p384.json`; dense `shared-relaxed/dense-monotone-d07-tight.json`; fresh sparse `shared-relaxed/shared-r2-d07-complete-p384.sparse.json` | Same fresh assembler and independent scope/sum check; all 303 dense intervals and q=1..32 freshly replayed at 384 bits |
| Shared GF16 R4, partial q=1–22 | `shared-gf16-q1-p384.json`, `shared-gf16-r4-d10-q2-16-fill.json`, `shared-gf16-r4-d10-return-q16-32.json`, `shared-goal/q19-p384.json`, `shared-goal/joint-counts-q20-22-p384.json`, `shared-goal/wide-counts-q22-32.json`, `shared-goal/wide-counts-q22-23-p384.json` | `gf16_packets/shared_support.py`, `gf16_packets/shared_sparse.py`; [progress and reproduction](gf16_packets/PROGRESS.md) |
| Pairwise GF16 R4, partial q=1–49 and dense cells | `pairwise-goal/prefix1-48-p384.json`, `pairwise-goal/shell-cdf-sparse-probes.json`, `pairwise-goal/dense-cover-unioncost-parallel.json`; obstruction replay `pairwise-goal/regional-retune-0128-replay-p384.json` | `gf16_packets/pairwise_sparse.py`, `gf16_packets/pairwise_cover.py`, `gf16_packets/pairwise_retune.py`; [plateau ledger](gf16_packets/PAIRWISE_SHUFFLE.md) |

The lane-only assembly drivers freshly replay their dense witnesses and
combine them with the documented earlier sparse lemmas. They do not rerun
every sparse cover. Their linked records give the separate sparse commands.
The GF16 full assembly regenerates its sparse prefix as well as replaying
the dense partition; use the complete commands in its README.

Quick integrity/scope check of the retained 10% GF16 certificate:

```text
python -B research/workstreams/permutation_locality/gf16_packets/audit_complete.py tmp/gf16-r4-d10-assembly-p384-complete.json tmp/gf16-lowthreads-r4-d10-q49.json --distance .1 --bits 48 --prefix tmp/gf16-r4-d10-q1-48-prefix-p384.json
```

This check passed on 2026-09-29: all 427 dense intervals form a disjoint
exhaustive partition, the independently regenerated sparse endpoint matches,
and the exact dyadic aggregate is below 2^-48. This is a receipt audit,
not a fresh numerical replay of the inequalities.

## Preserved Snapshots

The immutable local ZIP is:

`research/.proof-archives/packet-proofs-before-shared-gf16-20260929.zip`

SHA256:
`6a40992c3de7176bb7f00a241aa89e5f41d69d56944c9c9a2366bf09e624d0dd`

It contains 4,113 files (194,455,975 compressed bytes): research/library
sources and notes, authenticated BCH premises, the principal witnesses
for the completed tracks above, shared-route receipts through the first
q=1–16 check, saved two-bit sparse partitions, and the matched timing logs.
The later conditional-return experiment is not in this earlier snapshot.
`MANIFEST.json` records each file's size and SHA256, the source Git commit,
the presence of uncommitted sources, and the Python/package versions.
Every member was read back and hash-checked after creation. This is a
local preservation copy, not a substitute for an off-machine backup.

```text
python -B research/workstreams/permutation_locality/proof_archive.py --verify research/.proof-archives/packet-proofs-before-shared-gf16-20260929.zip
```

The archive helper refuses to overwrite an existing snapshot. Its `--verify`
mode does not extract or overwrite workspace files. The archive directory
is ignored by Git, so the numerical data is not added to the repository.
Future snapshots should use new names. Keep old sources and receipts;
record new claims under their own ensemble ID and output paths.

The later shared-shuffle plateau checkpoint is preserved separately in:

`research/.proof-archives/packet-proofs-shared-gf16-plateau-20260929.zip`

SHA256:
`f5da6b5909d78dcb4ff374732d850621c76a6a676449a6b1f6e143f9560977f0`

It contains 4,196 files (195,295,905 compressed bytes). In addition to the
old source and witness scope, it preserves the conditional-return prefix,
all JSON receipts in `tmp/shared-goal/`, the goal's logs, and the updated
proof sources and plateau ledger. Its manifest and every archived member
were hash-checked after creation. The original archive is unchanged.
Both archives remain local, ignored by Git, and separate from a remote backup.

```text
python -B research/workstreams/permutation_locality/proof_archive.py --verify research/.proof-archives/packet-proofs-shared-gf16-plateau-20260929.zip
```

The pairwise plateau checkpoint is preserved in:

`research/.proof-archives/packet-proofs-pairwise-gf16-plateau-20260929.zip`

SHA256:
`41d3aaab47c0fe13c54c0df04420f3e309ff00e549debc4f8eeb40c50b6c6f12`

It contains 4,391 files (199,650,650 compressed bytes), including the
pairwise sources, partial witnesses, failed diagnostics, fresh obstruction
replay, test logs, matched performance logs, and plateau ledger. Its
manifest and every member were read back and hash-checked. Both earlier
archives were reverified and remain unchanged. This archive is also local
and ignored by Git; no raw experiment data was committed or pushed.
This index entry was added after snapshot creation.

```text
python -B research/workstreams/permutation_locality/proof_archive.py --verify research/.proof-archives/packet-proofs-pairwise-gf16-plateau-20260929.zip
```

The completed lower-distance shared-R2 investigation is preserved in:

`research/.proof-archives/packet-proofs-shared-r2-d07-complete-20260930.zip`

SHA256:
`3d4295e4e2c49e134443bd0ac8ce9fa1f7d74867e75284bd9f5ca94e47793bd3`

It contains 4,571 files (207,340,796 compressed bytes). It preserves the
completed 5%, 6%, and 7% shared-R2 receipts, all 303 final dense witnesses,
the fresh sparse replay, count witnesses, failed higher-distance probes,
proof sources, test logs, and matched timing records. The manifest and
every member were read back and hash-checked. All three earlier archives
were reverified and remain unchanged. This snapshot is local and ignored
by Git; no numerical data was committed or pushed. This index entry was
added after snapshot creation.

```text
python -B research/workstreams/permutation_locality/proof_archive.py --verify research/.proof-archives/packet-proofs-shared-r2-d07-complete-20260930.zip
```
