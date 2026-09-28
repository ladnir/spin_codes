# Rank-four diagnostics and the move to two columns

2026-09-24. This investigation led to the two-column candidate in TWO_COLUMNS.md.
The four-column refinements below remain useful diagnostics, not certificates.
They did not close rank four and do not refute the four-column construction.

## Count subspaces and their extensions together

For a rank-four row tuple, let U be its four-dimensional span. The even
linear combinations of the four rows form a three-dimensional subcode H.
If their support sizes are u and v, respectively, then exactly u-v columns
are all ones. For each fixed pair (U,H), 1344 ordered row bases induce this H.
The eight words outside H have total weight 8u-4v, so one has weight at most
floor(u-v/2).

Let H_v upper-bound the number of three-dimensional BCH subcodes with support
exactly v, and let A_<=w upper-bound the number of nonzero words up to weight w.
A valid cumulative bound on the corresponding tuples is

    1344 H_v A_<=floor(u-v/2).

The initial screen bounds H_v by the existing rank-three tuple CDF at v,
divided by 2520. This is an upper shell count, not a difference of CDF caps.

There is also an exact-exterior-weight bound. For a fixed H with support S
of size v, restriction outside S has fibers of size 2^dim(C_S). Write D(v)
for a verified upper bound on this shortened dimension. An extension with
union support exactly u has exterior weight u-v. Its number of candidate
cosets modulo H is at most

    (binom(256-v,u-v) 2^D(v) - 8 [u=v]) / 8.

The subtraction removes H itself when the exterior is zero. Each other coset
contains eight words with the same exterior. Multiply this bound by 1344 H_v
for a tuple-shell bound. Exact small-code tests check the coset identity,
the 1344 factor, and the restriction-fiber count.

These counts alone did not solve the inner-bound problem. The coarse screen's
untruncated log2 upper bound was +51.85; exact exterior shells reduced it
to +45.66. Both exceed one and are therefore vacuous probability bounds.

## Retain total weight and the first feedback distribution

Let W be the sum of the four row weights. The ordinary BCH spectrum bounds
the coefficient of x^W in the fourth power of its nonzero enumerator.
For union support u, first discard individual word weights above u.
This bounds joint support/total-weight classes, even after discarding
the independence condition on the four words.

Also, each of the v non-all-one occupied columns has weight between one
and three. Hence

    4u-3v <= W <= 4u-v.

This restricts which H-support counts contribute to each (u,W) class.
The screen uses width-16 weight intervals and pays for their upper endpoint
when removing a factor rho^W, for 0<rho<=1. Exact small-code enumeration
checks the resulting combined count caps.

The earlier inner envelope discarded the distribution of Bx after activation.
It retained only that the state was nonzero. Exact pair counts show substantial
loss from this relaxation. For independent local draws of shapes
(1,1,1,1) and (4,4,4,4), the feedback vectors never coincide. The old
arbitrary-state cancellation bound allows probability 1/8 for the second draw.
Among all 4761 ordered shape pairs, 2654 have zero collision probability;
286 of those have the same feedback parity, so parity alone cannot explain them.

The new diagnostic adds one coordinate for each normalized, nonzero Bx
distribution. It retains that distribution through empty lazy steps.
A uniform refresh returns to the old density representation; another active
step returns to the coarse representation, but uses the exact pair-collision
count for its lazy zero-state transition. This is memory in the analysis,
not a change to the encoder.

The implementation uses sparse block structure instead of dense 76-by-76
matrices. Independent dense calculations check its region recurrence and
backward actions. It maximizes complete shape actions on a continuation
vector. Taking entrywise maxima of the expanded matrices would sum bounds for
mutually exclusive outgoing memory classes, greatly loosening the result.

## Joint-support moments from the existing dual-distance proof

The existing BCH certificate gives dual distance at least 30. Its split-proof
replay passed again during this investigation. Thus projections onto at most
29 coordinates are uniform. For h independent uniform BCH words, the union
support has the first 29 moments of Binomial(256,1-2^-h).

Write Q=2^h. The Q-ary Krawtchouk polynomials K_j are orthogonal for this
binomial measure, with squared norms binom(256,j)(Q-1)^j. The degree-14
reproducing-kernel argument gives the point-mass bound

    Pr[union support = u]
      <= 1 / [sum_{j=0}^{14} K_j(u)^2 / (binom(256,j)(Q-1)^j)].

Multiply by 2^(128h) to upper-bound ordered tuples. For h=3, divide the
independent-tuple count upper bound by 168 to upper-bound three-spaces.
Dependent tuples only enlarge the sampled population, so they do not
invalidate these upper bounds. The script uses exact rational arithmetic
and tests its recurrence, orthogonality, and small-code inequalities.

The OA bound requires the dual-distance premise. `joint_oa_caps.py` states
that premise explicitly; it does not silently treat its algebra tests as
a proof of BCH dual distance. The successful two-column verifier does not
use this additional OA component.

## What the screens established

All entries below are binary64 diagnostics. Positive log2 values are vacuous,
not positive security margins. Counts use exact integers, but numerical moment
calculations need outward replay before they can certify anything.

| Four-column screen | Untruncated log2 union upper |
|---|---:|
| Total weight, old inner envelope | +87.65 |
| Total weight, first-orbit memory, flag counts, joint OA caps | +34.30 |
| Same, uniformly shuffle all sixteen positions in each bundle | +31.55 |

The full-window shuffle is a different distribution. It offers only a modest
improvement in these screens and has not been benchmarked. Do not transfer
the four-column implementation timing or partial certificate to it.

Reducing the bundle to two columns instead gave a -61.30 diagnostic for
rank four with the coarse inner envelope and refined counts. More importantly,
the simpler outward verifier then closed every single-group rank without these
refinements. Its rank-four margin is 51.18 bits, and its combined margin is
47.38 bits. TWO_COLUMNS.md records that verified result and the matched timings.

## Reproduction

The main diagnostic entry points are:

    python -B research/workstreams/permutation_locality/rank_four_flags.py --coarse
    python -B research/workstreams/permutation_locality/pair_shape_collision.py
    python -B research/workstreams/permutation_locality/rank_four_weight.py --memory --flags --oa
    python -B research/workstreams/permutation_locality/rank_four_weight.py --memory --flags --oa --full-window
    python -B research/workstreams/permutation_locality/rank_four_weight.py --columns 2 --flags --oa
    python -B research/bch_spectrum_work/bch_spectrum_codex_bundle/code/certify_bch_shift_rank_split.py --verify

For the diagnostic matrix operations, use one BLAS thread. No encoder benchmarks
run concurrently. These screens write no raw data files and change no production
defaults. The useful next proof target is multiple active groups for the
two-column distribution, rather than more tuning of the four-column relaxation.
