# Keeping the exact empty quiet-state entry

This control changes the proof calculation, not the G16 encoder. It uses
the saved paired nonzero-amplitude regional matrices at r=18,24,32 and
alpha=1/5. All numerical slopes remain fixed. The result is a FLOAT
screen, not a directed certificate or evidence of actual low-weight words.

## Positive empty/nonempty decomposition

Let R00(m) be the regional weight-generating cap for entering and leaving
the inner in state zero, when m outer groups occupy that region. Its
empty value is exactly R00(0)=1. For the r=18 witness, its log2 value is
approximately -28.89 at m=2 and -42.05 at m=18. The previous single
exponential envelope, forced through m=0, loses about 28 bits on the
nonempty part.

Fix B>0 and choose C<=1 so that R00(m)<=C B^m for every feasible m>0.
Then

    R00(m) <= C B^m + (1-C) 1[m=0].

Both terms are nonnegative. The second term requires every group mask
to avoid that region. The other three regional state transitions retain
their existing exponential envelopes.

Expand this identity only at the 0-to-0 positions of each chronological
sixteen-region state path. If that path has n00 such positions, retaining
k of them as allowed coordinates has coefficient

    C(n00,k) C^k (1-C)^(n00-k).

The other n00-k coordinates must be absent from every group mask. The
fixed-type mask laws depend only on cardinality, so the particular choice
of allowed 0-to-0 positions need not be enumerated. The elementary
symmetric polynomial instead uses k copies of the 0-to-0 coordinate
marker, plus all markers for the other state transitions.

The expansion preserves ordered state endpoints and path multiplicities.
Its cut coefficients remain outside the sixteen-column Holder norm and
the subsequent Minkowski bound. Type weights, constituent-row sets, and
support markers remain fixed across all paths and columns. Their original
sixteen-coordinate moments are not replaced by moments on a smaller
ambient space. The helper accepts both the older nonempty-column floors
and the new joint pre-H weight laws.

## Bounded control and remaining obstruction

The driver tests three existing floor templates and five joint pre-H
templates at each anchor: twenty-four controls, with no parameter fitting
and no new regional banks. Unlike the high/low split, this control bounds
all arities at each of the three selected r values. It does not cover a
continuous range of r.

| Active groups | Previous selected log2 bound | Empty-separated best |
|---:|---:|---:|
| 18 | +6209.338583 | +2700.251233 |
| 24 | +7998.822805 | +3789.765729 |
| 32 | +10360.850231 | +6039.627030 |

Every bound still fails. The calculation took 1.094 seconds. At r=18,
the selected bound is still dominated by admitted all-empty comparison
mass. At r=24 and r=32, the dominant column term instead alternates eight
0-to-1 transitions with eight 1-to-0 transitions and selects arity eight.
The empty/nonempty fix therefore removes a large identified loss but
exposes a different one. It does not justify another search over the
same slopes.

The eight helper tests compare against explicit labelled masks and
ordered matrix products, check the C=0 and C=1 endpoints, verify the
forbidden coefficients, and reproduce the earlier type laws when no
coordinate is forbidden. A separate reviewer checked the cut weights,
counting factors, and fixed-type chronology. The complete lightweight
regression passes 220 tests.

Artifact: `theory/hd_two_sum_gf16_emptycut_anchors_float.json`, SHA-256
`acee5349066ac7ad76528dd29bea1eab2c9dec1367425981bc6c5d5683fee3a6`.
Helper: `theory/hd_two_sum_gf16_state_emptycut.py`, SHA-256
`5cca0a50e5f0517a0a3066e09e50a98833c8deb6bb8e31925c4c966051f9c84f`.
Driver: `theory/hd_two_sum_gf16_emptycut_control.py`, SHA-256
`cd1990bc0c5ae23b76de4077fd6d6a46c42fda915c40c0a8cd0234d1f6d1b8b7`.

Next, inspect the dependence discarded in repeated activation/return
paths. Two exact constraints are available for the next bounded control:

- Both off-diagonal regional entries vanish at empty occupancy, whereas
  their current positive exponential envelopes do not. The alternating
  comparison path can therefore admit impossible empty transitions.
- A fixed source support has fixed adjacent-pair counts under the two
  matchings. Splitting by their maximum F restricts each regional full-pair
  count to f<=F. The current regional maximum discards that restriction.
  Any such source stratification must pay its source count only once,
  outside alpha, and retain the same source class across all regions.

Neither improvement has a measured gain yet. The complete certified
implementation and its certificate remain the fallback; none of these
partial results is promoted.
