# Finite RM2Sub grid findings

Start with [Current engineering results](CURRENT_ENGINEERING_RESULTS.md) for
the consolidated numbers, trends, evidence classes, and 40-bit candidates.
[The numerical appendix](CURRENT_RESULTS_TABLES.md) includes all 130 BCH
geometries and the matched family references.

The grid evaluates 3,108 native parameter tuples. Write Q for the number
of nonzero outer rows and L for the total number of outer rows. Every
tuple has a Q1 calculation, Q2..64 bounds, composition-preserving
Q2/Q3/Q4 bounds, and a typed cover of Q65 through L. Thus the occupation
coverage is complete. Most full bounds remain too loose to select a code.
All numerical results use nearest binary64 arithmetic; the study produces
no new outward certificate.

The final snapshot contains 263,262 observations. All 72 tests pass.
The strict completion audit verifies 13,352 registered dependencies and
203,739 selected union components, with no unregistered completed batch.
The source-content and checkout-line-ending checks also pass.

The primary constituents are the four complete exact BCH spectra at
lengths 8, 32, 64, 128 and the four complete exact RM spectra at lengths
8, 32, 128, 512. The random reference uses a uniform full-rank rate-half
subspace at lengths 8 through 1024, in powers of two. One constituent is
reused across rows. For random references, conditional shell bounds use
one authenticated simultaneous spectrum event and pay its failure budget
once in a full union.

Each family is evaluated at k=2^16, 2^18, 2^20, 2^22, 2^24 and
t=64,128,256, with log2(t)+1 <= s <= 20. Maps are nested prefixes of
the frozen generator chains. Twelve length-1024 random tuples have fewer
rows than one epoch and are explicitly inapplicable. BCH-256 and the
partial-spectrum RM(5,11) do not enter this primary grid.

## What the complete dense sweep establishes

The deterministic typed cover gives 3,108 dense-range observations. Of
these, 1,574 select the typed coefficient bound and 1,534 select the
explicit all-message counting fallback. None of these coarse dense-range
observations has a positive margin. This records the limitations of this
finite witness search; it does not establish that the codes fail the
distance target. Selected exact-region composition refinements and the
earlier balanced occupation sweep give stronger full unions.

`complete_grid_unions.csv` records the selected full and partial margins.
Its JSON companion identifies every selected receipt and any setup event.
`complete_landscape_audit.json` checks the complete cover, source hashes,
conditional-event containment, union sums, database, and compressed
snapshots. `certified_results` remains empty because binary64 diagnostics
are not outward certificates.

The full typed cover is reproducible with
`python run_typed_range_grid.py --output-dir activation_occupation_grid_typed_all_v1`.
Completed batches are verified and skipped. The selected small-state RM
comparison uses
`python run_threeband_grid.py --output-dir activation_occupation_grid_rm49_boxes_tradeoff_v1 --steps 64 128 256 --states 11 --message-exponents 16 --blocks 512 --exact-only`.
Never run these numerical producers concurrently. Their receipts bind the
implementation and input files by SHA-256; use a new version and output
directory for changes to a bound or its witness search.

## Choosing t and s

The four full positive diagnostics are
RM(4,9) at k=2^16, t64/s18, with margin 42.9694638927 bits, and random
[512,256] at k=2^16, t64/s12, s18, and s20, with margins near 60 bits.
The random margins include the shared setup-event failure term; a larger
conditional margin cannot remove that term. These are candidates for
outward replay.

The targeted s=11 composition refinement does not close the full bound
for any of t=64,128,256. At t64, its unresolved tail starts at Q248;
larger t also retains middle-occupation slack. The full union selects the
best compatible old or new receipt for each range, so a weaker new bound
does not replace a stronger existing screen. Within this broad ledger, the smallest
observed full 40-bit RM choice is s=18 at t64. The separate outward
RM(4,9) certificate at k=2^16, t64/s14 supplies a smaller-state reference;
it is not imported into this ledger. See the current overview for both scopes.

Use `tested_parameter_frontiers.csv` to compare thresholds of 0, 20, 40,
and 60 bits at each message length. Its coverage column is part of the
claim: a Q1 or Q1..64 passing state is only a partial screen. The smallest
observed passing state is relative to this fixed family of maps and
witnesses; it does not prove global optimality.

For RM(4,9) at k=2^16, the Q1 screen first reaches 40 bits at s=11 for
each of the three epoch sizes. The corresponding implementation proxies
are:

| t | Q1 margin (bits) | Full diagnostic margin (bits) | Independent XORs per output bit | State updates per output bit |
|---|---:|---:|---:|---:|
| 64 | 40.2391 | -4,539.38 | 11.328125 | 1/64 |
| 128 | 40.3287 | -15,754.21 | 11.664063 | 1/128 |
| 256 | 40.2283 | -21,303.45 | 11.769531 | 1/256 |

Larger t reduces the frequency of state updates while slightly increasing
this independent-XOR proxy. All three choices remain on the proxy Pareto
frontier. The counts include the literal A and B linear maps and the
output XOR. They exclude field arithmetic, outer encoding, permutations,
SIMD, and shared XOR circuits. They support a benchmark shortlist, not a
runtime ranking.

The Q1..4 screen for RM(4,9) requires s=14 at k=2^18 for all three t
values to reach 40 bits. No tested s through 20 reaches that threshold
in this screen at k >= 2^20. The exact BCH [128,64,22] Q1..4 screen does
not reach 40 bits at any of the five message lengths. The lower thresholds
remain in the tables; a missed 40-bit threshold does not discard a tuple.

## Scaling beyond complete spectra

The 473 projections in `exact_family_q1_projections.csv` are engineering
estimates of Q1, fitted only to complete exact spectra. They do not project
the full occupation union. BCH lengths 256 and 512 use fitted distance
and minimum-shell-count proxies; these are not asserted properties of
future BCH codes. RM length 2048 uses exact minimum-distance and
minimum-shell-count formulas, with a fitted transfer cost. Random lengths
through 1024 are computed directly and are not extrapolated.

Holding out the largest known constituent gives median absolute Q1 errors
of about 13.66 bits for BCH and 7.71 bits for RM; maximum errors are about
19.02 and 14.09 bits. These errors are large relative to a narrow 40-bit
target margin. The estimates are useful for choosing the next constituent
size and witness search, but should not justify a near-threshold choice.
The fitting-window ranges are sensitivity checks, not confidence bounds.
See `EXACT_FAMILY_PROJECTIONS.md` for the model and authenticated training
rows.

The separate BCH-256 result is recorded in
`secondary_bch256_comparison.json`. Its reported full margin is about
50.44 bits at k=2^20, t64/s20, using a different selected inner map and
a one-bit larger cutoff. It is a qualitative comparison, excluded from
the fits and the primary grid. Its outward calculation was not replayed
in this worktree.

## Next implementation decision

Replay the chosen positive full RM diagnostic with outward arithmetic,
then benchmark its t/s alternatives sequentially using the actual field
and permutation implementation. For larger messages, prioritize tighter
dense witnesses on promising exact-spectrum or directly computed random
tuples. The existing grid supplies reproducible candidates and exposes
where the current proof relaxation, rather than a measured runtime,
limits selection.
