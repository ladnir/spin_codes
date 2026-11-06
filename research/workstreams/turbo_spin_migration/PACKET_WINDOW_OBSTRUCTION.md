# A packet-window obstruction at the paper's growth constant

## Statement and scope

Consider the growing, rate-half Turbo family in
[TURBO_FAMILY_BRIDGE.md](TURBO_FAMILY_BRIDGE.md). Each local outer map is
linear from `2R` bits to `4R` bits, its output is partitioned into `R`
four-bit packets, and every outer position sends one packet to each
region. Each region independently applies a uniform permutation to its
`L` packet slots. The inner is linear, has a 22-bit state, starts in
state zero, and outputs zero when both its current input and state are
zero. Put `b=4R`, `N=bL`, and choose admissible lengths with

```
b = (39/4) log_2 N + O(1).
```

These properties alone imply that, with probability tending to one
over the regional permutations, there is a nonzero message whose
encoded output has relative weight at most `0.10+o(1)`. This holds for
*every* rate-half linear outer realization, including the selected
Golay--BAA outer, and for every realization of the inner maps. Thus the
four-bit packet family at the paper's `39/4` growth constant **cannot**
have a `>11%` asymptotic minimum-distance theorem. The obstruction is
about this growing family; it says nothing adverse about its finite
RS8 operating points.

## A concrete 20% window

Take `L` divisible by 80, so the final fifth of every region consists
of whole 64-bit inner steps. This includes the native paper sequence
`L=128m` along `m` divisible by five. Call an outer position *good* if its
packet slot is in that final fifth in at least `ceil(13R/25)` regions.
For each fixed position, the independent regional permutations give

```
A_i ~ Binomial(R,1/5).
```

Let `p_R=Pr[A_i>=ceil(13R/25)]`. Stirling's bound at the single count
`ceil(13R/25)` gives

```
p_R >= exp(-R D(13/25 || 1/5)-O(log R)).
```

Since `b=4R` and `log L/R -> (16/39) ln 2`, the strict inequality

```
D(13/25 || 1/5) = 0.2516696520... < (16/39) ln 2 = 0.2843680740...
```

implies `L p_R -> infinity`. The outward-rounded check is in
`verify_packet_window_obstruction.py`. Good positions are not
independent within a region, but this causes no loss: for any two
different positions and any pair of their length-`R` indicator strings,
the joint probability is at most `(L/(L-1))^R` times the independent
Bernoulli probability. Therefore the number of good positions has
variance at most `Lp_R+O(R/L)(Lp_R)^2`, and it is at least 300 with
probability `1-o(1)`.

Fix 300 good positions after sampling the route. For each one, require
that the outer word's packets be zero in all regions where that
position's slot falls outside the final fifth. The packet-to-region
map is a permutation, so these are at most
`4(R-ceil(13R/25))` linear constraints on `2R` message bits. Each
position retains a message subspace of dimension at least

```
2R - 4(R-ceil(13R/25)) >= (2/25)R.
```

Across 300 positions, the direct-sum dimension is at least `24R`.
Require the 22-bit inner state to be zero at the end of every one of
the `R` regions. These are at most `22R` additional homogeneous linear
constraints, leaving a subspace of dimension at least `2R`. Its state
enters every region at zero. No input packet appears before the final
fifth, so all output bits before that window are zero.

The output of this subspace is supported on at most `N/5` bit
coordinates. If the encoder is noninjective, a nonzero message already
has zero output. Otherwise the output subspace has dimension at least
`2R`. Averaging the weight of its nonzero codewords gives a word of
weight at most

```
(N/5) / (2(1-2^(-2R))) = (0.10+o(1)) N.
```

The argument uses linearity only. It does not assume the outer is
random, that packet labels are uniform, or that a feedback cancellation
has any particular probability.

## Stronger limiting upper bound

For a window fraction `f<1/2` and any `a>1/2`, a position with at least
`aR` slots in the window retains at least `4(a-1/2)R` message
dimensions. Choose any fixed number `Q>22/[4(a-1/2)]` of such
positions. The same boundary-state argument gives distance at most
`f/2+o(1)` whenever

```
D(a || f) < 4 ln 2 / (39/4).
```

Let `a` decrease to `1/2`. The critical window fraction solves
`D(1/2 || f_*)=4 ln 2/(39/4)`, hence

```
f_* = (1-sqrt(1-2^(-32/39)))/2 = 0.1706979467... .
```

For every `epsilon>0`, a fixed `Q` and an admissible length sequence
give distance at most `(f_*/2+epsilon)N` with probability `1-o(1)`.
Thus the asymptotic relative-distance upper bound along these lengths
is `f_*/2 = 0.0853489733...`. The explicit 20% window above is a
shorter standalone counterexample to the 11% target.
The limiting bound does not depend on the constant state width: for
any fixed width `s`, choose the fixed number of active positions
`Q>s/[4(a-1/2)]`. Increasing 22 to another fixed state width alone
cannot repair the `39/4` four-bit packet family.

## Design consequence

The one-, two-, and three-position first-moment bounds close, but
their success cannot extend to a full 11% proof for this exact growing
family. The missing occupation classes contain real low-weight words,
not merely a loose transfer inequality.

The dimension argument extends to packet width `p` when `b=pR` and
the local outer still has input dimension `b/2`. A position with at
least `aR` packet slots in a window retains at least
`p(a-1/2)R` message dimensions. A fixed number
`Q>s/[p(a-1/2)]` of such positions can force all `sR` boundary-state
coordinates to zero. At outer-growth constant `c`, the critical
window and distance upper bound are

```
f_*(p,c) = (1-sqrt(1-2^(-2p/c)))/2,
delta_upper(p,c) = f_*(p,c)/2.
```

At `c=39/4`, this gives distance upper bounds of approximately
`15.90%`, `12.56%`, `10.27%`, and `8.535%` for packet widths
`p=1,2,3,4`, respectively. Width two is the largest integer packet
width for which this argument leaves room for an 11% theorem at that
growth constant. This is a necessary screen, not a distance proof.

Two changes remove *this particular* obstruction without by
themselves proving 11% distance:

* Use two-bit packets in the growing member. At the same growth
  constant, the analogous window argument gives only a 12.56%
  distance upper bound, above the target.
* Keep four-bit packets but increase the outer growth constant to at
  least `4 ln 2 / D(1/2 || 0.22) = 14.73626...`. This is necessary to
  evade the window argument at 11%, not sufficient for a proof.

The finite four-bit Turbo encoder and its certificates remain intact.
Before changing the paper, decide whether the growing theorem should
use a tunable packet width or a larger growth constant, then redo the
full occupation analysis for that explicit construction.
