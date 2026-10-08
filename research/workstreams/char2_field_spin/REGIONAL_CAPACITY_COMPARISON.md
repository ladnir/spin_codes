# More outer parity repairs the tested capacity deficit

The next candidate should spend work on outer redundancy rather than only
increase the inner state. GRS[32,22] with four state symbols currently passes
the strongest tested window-event upper bound and the one-row bound. Its
rate is 11/16, rather than the binary baseline's 1/2. GRS[34,24], at rate
12/17, is a second candidate with less proof headroom.

Neither has a complete distance certificate or a new encoding time. These
results select the next proof problem; they do not approve implementation.
The [cost comparison](PARITY_STATE_COST.md) separates arithmetic estimates
from runtime. The independent [dense-class theorem](REGIONAL_GLOBAL_ROUTE.md)
covers another part of the proof without using the inner's state analysis.

## Construction held fixed

Use a fixed systematic GRS[b,k] outer over GF256 on L source rows. Each
outer coordinate has its own region, containing one coordinate from every
row. Independently permute the L row labels in each region. These regional
permutations form one global routing stage. Sample independent nonzero
GF256 multipliers for all coordinates.

The inner has width sixteen and four state symbols. It uses the period-64
translation maps from [the scalar reference](REGION64_CANDIDATE.md), with
independent uniform GL(4,256) orientations at region boundaries. There is
no per-packet random matrix. Initial state is zero; final state is discarded.
The outer, diagonals, route, and orientations have the independence specified
in the earlier construction. The previously rejected [32,24] profile is not
being reused as a certificate for these new outers.

For requested K=2^ell, round L up to a multiple of 64 so that kL>=K. Put
N=bL and e=ceil(N/10)-1. The target excludes output weights through e.
All reported bounds concern ideal setup randomness over GF256. Scalar
extension preserves the resulting code distance over the compatible fields.

## One-row comparison

The exact return-count argument in `rate_design_region_refresh.py` gives:

| Outer | Rate | State | K20 one-row margin |
|---|---:|---:|---:|
| [32,24] | 3/4 | 5 | 35.45 bits |
| [32,24] | 3/4 | 6 | 51.45 bits |
| [33,24] | 8/11 | 4 | 39.30 bits |
| [32,23] | 23/32 | 4 | 49.45 bits |
| [34,24] | 12/17 | 4 | 46.68 bits |
| [32,22] | 11/16 | 4 | 70.64 bits |

Increasing the state does not monotonically improve this bound: it reduces
the minimum weight emitted by a nonzero state. Adding a region can also
allow one additional intermediate region inside the output-weight budget.
These discrete changes explain the weak [33,24] and state-five entries.
A failed upper bound does not prove that either construction has bad distance.

## The window event and why equal lengths were insufficient

Fix h outer coordinates, shared by a selected set of a source rows. Write
p=b-k. Shortening each row to these coordinates leaves h-p dimensions.
If all selected row positions fit inside one window in each selected region,
the source space has dimension a(h-p). Requiring zero state after each
window imposes at most sh constraints for state size s. Thus

    a(h-p)>sh

forces a nonzero word confined to the selected windows. If the last window
is a suffix in the final region, its terminal state need not vanish; replace
sh by s(h-1). These are sufficient conditions for a bad code, not necessary
conditions for every bad word.

An equal-window test sums over every h and every possible starting position,
but fixes all h lengths to floor(e/h), rounded down to whole packets. That
test gives about 65.78 bits for [32,23] and 58.49 bits for [33,24] at K20.
Allowing unequal lengths makes the corresponding bound substantially larger.

For a fixed a-subset of row labels, the number of routed position subsets
with span r is

    (L-r+1) binom(r-2,a-2),   a<=r<=L.

Here span includes both endpoints. Dividing by binom(L,a) gives its exact
probability. Replacing L-r+1 by L and extending the range gives a positive
generating-function upper bound

    L z^a / ((1-z)^(a-1) binom(L,a)).

For a suffix, use its extent from the first selected position to the final
position; the corresponding upper series is

    z^a / ((1-z)^a binom(L,a)).

Multiply h such series, sum coefficients through total length e, and union
over coordinate sets and row subsets. For each h, the all-return family has
probability at most

    binom(b,h) L^h binom(e-h,h(a-1)) / binom(L,a)^(h-1),
    a = floor(sh/(h-p))+1.

The final-suffix family has probability at most

    binom(b-1,h-1) L^(h-1) binom(e-h+1,h(a-1)+1) / binom(L,a)^(h-1),
    a = floor(s(h-1)/(h-p))+1.

Union both bounds over p<h<=b. This covers all starts and unequal lengths
with total budget at most e, and relaxes packet alignment. It still covers
only the stated kernel events: the rows share an outer support, with one
window per selected region. Different row supports, multiple windows in a
region, and low-weight words without a dimension-forced kernel remain open.

## Exact comparison across sizes

The following numbers are minus the binary logarithm of the rational
event-family upper bound. They are not complete code margins.

| Outer, state | K16 | K18 | K20 |
|---|---:|---:|---:|
| [32,23], 4 | 88.52 | 57.51 | 27.17 |
| [33,24], 4 | 81.53 | 50.39 | 19.84 |
| [34,24], 4 | 127.74 | 102.37 | 74.44 |
| [32,22], 4 | 142.25 | 117.14 | 90.85 |

The [32,22] one-row margins are 74.63/72.63/70.64 bits. Its rounded
(K,N) pairs are (66176,96256), (263296,382976), and (1048960,1525760).
The [34,24] one-row margins are 50.67/48.68/46.68 bits.

The dense-class theorem covers inverse inputs of weight at least .24N for
[32,22] and .25N for [34,24]. Its exact K16/K18/K20 margin lower bounds are
578/3381/14745 and 527/3244/14318 bits, respectively. These large numbers
do not describe the sparse class or a complete setup-failure guarantee.

The verified [32,22] outer needs 41 fixed GF256 products per four source
rows. Including the regional inner, its leading processed-byte arithmetic
volume is 177/88 per source symbol, compared with 385/192 for [34,24].
Those estimates differ by less than 0.31%; neither includes measured
stride or cache effects. Proof headroom, rather than this tiny arithmetic
difference, currently favors [32,22].

This comparison favors [32,22] as a proof candidate. It does not establish
that [32,23] fails forty bits: 27.17 is the margin of an upper bound, unlike
the previous actual failure lower bound for [32,24], state four.

## A separate dense packet-map experiment

A fixed Cauchy matrix on the sixteen inputs and state symbols gives an
exact local enumerator from submatrix ranks. It deterministically excludes
more low-input, low-output branches than identity feedthrough. The new
`cauchy_inner_transfer.py` evaluates that law with positive integer counts;
small-field exhaustive checks and mass identities pass.

Its law requires independent GL state refresh before each packet, unlike
the regional candidate above. It is a different construction. The tested
ordinary source-word upper bounds still fail for sparse and intermediate
row counts. The elementary implementation budget is also much larger.
Preserve this exact enumerator as a reference, but do not start a native
Cauchy kernel on the strength of the local proof alone.

## Next proof obligation and replay

The dense-class theorem handles all inverse-input weights above its stated
cutoff. The one-row theorem handles only one active source row. The window
bound handles only the specified dimension-forced event family. Their union
does not cover the complete code. The next task is a covering argument for
the remaining sparse multirow words, with conditional counts clipped before
averaging rare routes.

```sh
python -B research/workstreams/char2_field_spin/theory/regional_capacity_screen.py --outer 32 22 --state 4 --unequal-upper
python -B research/workstreams/char2_field_spin/theory/regional_global_route.py --self-test
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p regional_capacity_tests.py -v
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_cauchy_inner_transfer.py -v
```

No encoder benchmarks, production changes, commits, or pushes accompany
this checkpoint.
