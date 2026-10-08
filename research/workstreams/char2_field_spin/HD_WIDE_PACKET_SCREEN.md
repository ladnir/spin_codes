# Bounded screen of wide, unpunctured Cauchy packets

Larger state removes the previously observed fixed-window obstruction for
the geometries below. The established two-state rank cover still gives
negative sampled middle-range margins. No new positive proof candidate
or additional distance counterexample emerged from this bounded screen.

## Construction and scope

Apply a GRS [b,b/2] outer to independent source rows, with K=2^20 and N=2K.
Sample an independent uniform permutation within each outer-coordinate
region. After this sole route, a causal full-Cauchy packet consumes w inputs,
retains s state symbols, and emits all w outputs. There is no puncture.

The fixed packet maps are over GF256 embedded in GF65536. At each packet,
the input and state diagonal entries are independent uniform nonzero
GF65536 elements. The state block of the square Cauchy matrix is normalized
to the identity, as in `high_distance_fields.py`. The selected-output and
all-state rank identities are unchanged by increasing w and s. Disjoint
GF256 Cauchy node sets exist when w+s<=128, including w=s=64.

The new helper extends only the validation range of the existing local
rank formula. It does not propose a fast implementation of these dense maps.
No structured Cauchy multiplication cost was measured or established.

## The fixed-window obstruction no longer rejects these parameters

For a complete-packet prefix of W inputs per region, let T count the
coordinates of one outer row landing in those windows. Then
T has distribution Binomial(b,W/R), where R is the number of outer rows.
The expected shortened source dimension is

    R * E[(T-b/2)_+].

Zeroing the inner state after every regional window costs at most b*s
equations. The exact window test from `high_distance_window_obstruction.py`
therefore proves failure only when its expectation exceeds b*s. Here W is
the largest multiple of w strictly below 0.4R.

| Outer | w | s | R | W | Expected shortened dimension | Exit equations |
|---|---:|---:|---:|---:|---:|---:|
| [64,32] | 32 | 64 | 32768 | 13088 | 2780.911 | 4096 |
| [64,32] | 64 | 64 | 32768 | 13056 | 2674.841 | 4096 |
| [128,64] | 32 or 64 | 16 | 16384 | 6528 | 309.639 | 2048 |
| [128,64] | 32 or 64 | 32 | 16384 | 6528 | 309.639 | 4096 |

Every listed geometry survives this particular obstruction. A zero failure
lower bound is inconclusive; it does not prove a low probability of failure.
The expectations above are displayed as decimals but computed from exact
rational binomial sums.

## Sparse and middle-range screens

The local operator records whether the retained state subspace is zero or
nonzero. It uses the same selected-zero expansion, rank potential, and flat
cover as the existing helper. Computing the operator is polynomial in w
and s; it does not enumerate all 2^w selected-output sets individually.
The exact rank roots are the expensive step. The new helper caches them
at a small fixed grid of rank marks tau.

The initial grid was tau in {4,16,64,256}. For each tau, only the remaining
Chernoff marks were optimized. Searches used 20 or 25 iterations. A final
check used {4096,65535} and 35 iterations on [128,64], w=64, s=32.
Most middle-range optimizations reached their iteration limit, so the
table gives valid floating witness evaluations, not certified optima.

The sparse class fixes 16 active source rows. The middle classes fix total
outer-output weights 62915 and 838861, approximately 3% and 40% of N.
Margins are the negative base-two logarithm of the sampled upper bound:

| Outer | w | s | Active rows 16 | Input density .03 | Input density .4 |
|---|---:|---:|---:|---:|---:|
| [128,64] | 32 | 16 | +894.15 | -73514.17 | -625141.36 |
| [128,64] | 32 | 32 | +1360.33 | -73536.67 | -625140.46 |
| [64,32] | 32 | 64 | +541.25 | -76941.98 | -642707.24 |
| [128,64] | 64 | 32 | +1360.45 | -78207.55 | -625139.52 |

For the last geometry, the larger tau grid improves the sparse witness to
+2925.80 bits. The middle witnesses worsen to -297822.42 and -3070769.31 bits.

The positive sparse margins cover only one sampled active-row class.
They are not a 40-bit whole-code certificate. Negative middle margins mean
that the tested upper bound exceeds one; they do not exhibit bad codewords.
All evaluations use ordinary floating arithmetic after the exact rank-root
calculation. No interval enclosure or complete class cover was performed.

At these witnesses, increasing state from 16 to 32 barely changes the
middle result. The live-state transition maximizes over possible retained
ranks, which can discard the benefit of additional state. This observation
motivates a rank-resolved analysis if the construction remains interesting;
it does not prove that such an analysis will succeed.

## Reproduction and next decision

`theory/hd_wide_packet_screens.py` contains all 15 fixed witnesses. Its
`--check` mode compares 45 local operators with the established formula.
Those comparisons pass for widths through 32 and states through 32.
Its `--replay` mode evaluates the frozen witnesses without optimization.

```text
python -B research/workstreams/char2_field_spin/theory/hd_wide_packet_screens.py --check --replay
python -B research/workstreams/char2_field_spin/theory/hd_wide_packet_screens.py --window --active 16 --density .03 .4 --outer 128 64 --width 32 --state 32 --iterations 25
```

The next turn should not implement a dense wide packet on these results
alone. A useful next step is either a rank-resolved proof of the middle
range, or a different low-cost inner with an independently justified local
rank law. The proposed two additive states with a varying observation row
is such a different model; the present diagonal-refresh bounds cannot simply
be assigned to it without checking the joint conditional law.

Only this new note and the new proof helper were written for this screen.
No native code, build, performance benchmark, or remote command was run.
