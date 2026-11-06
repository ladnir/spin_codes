# A packet-window obstruction at the paper's growth constant

## Statement and scope

Consider the growing, rate-half Turbo family in
[TURBO_FAMILY_BRIDGE.md](TURBO_FAMILY_BRIDGE.md). Each local outer map is
linear from `2R` bits to `4R` bits, its output is partitioned into `R`
four-bit packets, and every outer position sends one packet to each
region. The route may place those packets in the region's `L` slots by
*any* permutation; independence and uniformity are unnecessary. The
inner is linear, has a 22-bit state, starts in
state zero, and outputs zero when both its current input and state are
zero. Put `b=4R`, `N=bL`, and choose admissible lengths with

```
b = (39/4) log_2 N + O(1).
```

These properties alone imply that, for **every** such route, there is
a nonzero message whose
encoded output has relative weight at most `0.10+o(1)`. This holds for
*every* rate-half linear outer realization, including the selected
Golay--BAA outer, and for every realization of the inner maps. Thus the
four-bit packet family at the paper's `39/4` growth constant **cannot**
have a `>11%` asymptotic minimum-distance theorem. The obstruction is
about this growing family; it says nothing adverse about its finite
RS8 operating points.

## A concrete 20% window

Take `L` divisible by 80, so each region partitions into five
consecutive windows of whole 64-bit inner steps. This includes the
native paper sequence `L=128m` along `m` divisible by five. Fix an
arbitrary route. Choose one of the five windows independently and
uniformly in each region. Call an outer position *good* if its packet
slot lies in the chosen window in at least `ceil(13R/25)` regions.
For each fixed position, randomness of the **window choices**, with
the route held fixed, gives

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
`verify_packet_window_obstruction.py`. By linearity of expectation over
window choices, the expected number of good positions is exactly
`L p_R`, regardless of dependence between positions. Hence **some**
choice of one window per region has at least `L p_R` good positions,
and therefore at least 300 for all sufficiently large native lengths.
This is a deterministic existence argument for every route.

Fix such windows and 300 good positions. For each one, require
that the outer word's packets be zero in all regions where that
position's slot falls outside the final fifth. The packet-to-region
map is a permutation, so these are at most
`4(R-ceil(13R/25))` linear constraints on `2R` message bits. Each
position retains a message subspace of dimension at least

```
2R - 4(R-ceil(13R/25)) >= (2/25)R.
```

Across 300 positions, the direct-sum dimension is at least `24R`.
Require the 22-bit inner state to be zero at the **end of the chosen
window** in each of the `R` regions. These are at most `22R` additional
homogeneous linear constraints, leaving a subspace of dimension at
least `2R`. The state starts at zero. After each chosen window it is
zero and remains zero until the next window, since the input there is
zero. Thus input and output vanish outside the chosen windows.

The output of this subspace is supported on at most `N/5` bit
coordinates. If the encoder is noninjective, a nonzero message already
has zero output. Otherwise the output subspace has dimension at least
`2R`. Averaging the weight of its nonzero codewords gives a word of
weight at most

```
(N/5) / (2(1-2^(-2R))) = (0.10+o(1)) N.
```

The argument uses linearity only. It does not assume a random route or
outer, uniform packet labels, or any particular probability of
feedback cancellation.

## Stronger limiting upper bound

For a window fraction `f<1/2`, choose in each region a uniform cyclic
window of `floor(f L/16)` consecutive 64-bit steps. Each fixed packet
slot is included with probability `f+o(1)`, independently across
regions, for **any** fixed route. A cyclic window has at most two
ordinary interval components, so forcing the state to zero at the end
of each component costs at most `44R` linear constraints. For any
`a>1/2`, a position included in at least `aR` windows retains at
least `4(a-1/2)R` message dimensions. Choose a fixed
`Q>44/[4(a-1/2)]` of these positions. The same expectation and
dimension argument gives distance at most `f/2+o(1)` whenever

```
D(a || f) < 4 ln 2 / (39/4).
```

Let `a` decrease to `1/2`. The critical window fraction solves
`D(1/2 || f_*)=4 ln 2/(39/4)`, hence

```
f_* = (1-sqrt(1-2^(-32/39)))/2 = 0.1706979467... .
```

For every `epsilon>0`, a fixed `Q` gives distance at most
`(f_*/2+epsilon)N` for **every** route at all sufficiently large
admissible lengths. Thus the asymptotic relative-distance upper bound
is `f_*/2 = 0.0853489733...`. The explicit 20% window above is a
shorter standalone counterexample to the 11% target.
The limiting bound does not depend on the constant state width: for
any fixed width `s`, choose the fixed number of active positions
`Q>2s/[4(a-1/2)]`. Increasing 22 to another fixed state width alone
cannot repair the `39/4` four-bit packet family. The explicit 20%
window proof also survives any state width `s(R)=exp(o(R))`: its
good-position count grows as `exp(0.0326... R-o(R))`, whereas the
number of positions needed to overcome the state constraints is only
`O(s(R))`. In particular, a logarithmically growing state alone does
not remove the obstruction.

## Design consequence

The one-, two-, and three-position first-moment bounds close, but
their success cannot extend to a full 11% proof for this exact growing
family. The missing occupation classes contain real low-weight words,
not merely a loose transfer inequality.

The dimension argument extends to packet width `p` when `b=pR` and
the local outer still has input dimension `b/2`. A position with at
least `aR` packet slots in a window retains at least
`p(a-1/2)R` message dimensions. A fixed number
`Q>2s/[p(a-1/2)]` of such positions can force the at most `2sR`
window-exit state coordinates to zero. At outer-growth constant `c`, the critical
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
Changing only the distribution of the region permutations cannot fix
this asymptotic family: the windows are chosen *after* any route is
fixed. Even if each outer position visits the five bins in a perfectly
balanced way, a different bin can be chosen in every region. Averaging
over those `5^R` bin patterns still yields a pattern with many good
positions. A repair must change a property used by the argument, such as
packet width, the number of regions, outer-growth constant, or the
amount of state information crossing each window. Then the full
occupation analysis must be redone for that explicit construction.
