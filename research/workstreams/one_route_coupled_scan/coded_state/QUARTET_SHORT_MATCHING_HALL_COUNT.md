# Exact dual Hall counting and its negative-witness entropy

2026-11-01. The half-row flow equals a matching with nine candidate inputs
per interior source. This note gives an exact count of compressed dual
Hall witnesses. It also proves that the unfiltered count is too large
for the proposed 64-window moment at marker `2^(1/8)`. No numerical jobs
or probability certificate are asserted here.

The flow identity and the proof that one temporal edge suffices are in
[QUARTET_FIXED_ROW_FLOW.md](QUARTET_FIXED_ROW_FLOW.md). The concrete
remaining moment target is in
[QUARTET_FIXED_ROW_FLOW_TAIL.md](QUARTET_FIXED_ROW_FLOW_TAIL.md).

## Dual witnesses and quartet saturation

Fix the quartet partitions and all routed row labels in a window of
`T` updates and width `n`. Let `J_t` be the non-`H` inputs at update `t`.
Assume exact total half-row balance, so `2*sum_t |J_t|=n*T`. Sources have
capacity one and inputs in `J` have capacity two. A source can reach its
own input and the quartet inputs one or two updates later.

For an input set `C subset J`, let `N^-(C)` be its source predecessors.
The ordinary capacitated Hall identity, applied on the input side, gives
the exact integer deficit

```
delta = max_(C subset J) [2|C|-|N^-(C)|].                    (1)
```

The balance assumption makes this the same deficit as unmatched source
units. Without balance, the two deficits differ by the difference in
total capacities.

For `t>=2`, put `U_t=Gamma_t(C_t)`, the union of fresh quartets meeting
`C_t`. Enlarge `C_t` to `U_t intersect J_t`. Each added input contributes
two capacity units and adds at most its own source to `N^-(C)`. Its
predecessors at earlier updates were already included by `U_t`.
Consequently the score in (1) does not decrease. The sets `U_t` remain
unchanged. Every selected quartet must meet `J_t`; quartets containing
only `H` inputs cannot be selected.

At update one there are no preceding sources. Adding any missing input
of `J_1` increases the score by at least one. We may therefore require
`C_1=J_1`, without defining or counting a set `U_1`.

Set `U_(T+1)=U_(T+2)=empty` and define

```
B_t=U_(t+1) union U_(t+2).
```

For these saturated witnesses,

```
N^-_t(C)=C_t union B_t,
score(C)=sum_t [2|C_t|-|C_t union B_t|].              (2)
```

Thus (1) is unchanged when restricted to the saturated witnesses above.
The argument selects a maximizing witness; it does not claim uniqueness.

## A positive reverse-time kernel

Fix a marker `y>1` for the score. A variable `x` marks the number of
non-`H` input labels. Immediately before processing update `t>=2`, retain

```
A=U_(t+1),              B=U_(t+1) union U_(t+2),
a=|A|,                  b=|B|.
```

The three old coordinate labels are `A`, `E=B\A`, and `Z=B^c`.
The new labels are `A'=U_t`, `E'=A\U_t`, and `Z'=(U_t union A)^c`.
Use the same letters as formal variables in the following homogeneous
polynomial, with primes denoting new-label variables:

```
P_x = (1+x)^4 [Z Z' + y^-1 A E' + y^-1 E Z']^4
    + (A')^4 {
        [(1+xy)Z + (y^-1+xy)(A+E)]^4
        - [Z+y^-1(A+E)]^4
      }.                                               (3)
```

The first term leaves the quartet unselected. Its score contribution
is minus the number of old `B` positions. The second selects it. A
selected non-`H` coordinate contributes one to the score, while a
selected `H` coordinate contributes minus one exactly when it lies in
old `B`. The subtraction removes the all-`H` selected case. Its removal
leaves nonnegative coefficients in `x` and the old/new label variables;
coefficients are positive Laurent expressions in `y`.

Let `m=n/4`. For old state `(a,b)` and new state `(a',b')`, take the
coefficient of

```
A^a E^(b-a) Z^(n-b) (A')^a' (E')^(b'-a') (Z')^(n-b')
```

in `P_x^m`, and divide by `multinomial(n;a,b-a,n-b)`. Call the resulting
polynomial entry `K_x[(a,b),(a',b')]`. This normalization averages the
fresh uniform partition exactly: the old labels are uniformly assigned
to the ordered quartet positions. New selected quartets are counted,
and the non-`H` label assignment is enumerated by `x`.

A safe state set is

```
0<=a<=b<=n,             a divisible by four.
```

At `n=512` it has 33,153 states. Some states are unreachable; keeping
them does not alter an exact computation from the specified start.
No pointwise maximum over quartet types has been taken.

At update one, the forced choice `C_1=J_1` gives the terminal vector

```
q_x(a,b)=(y^-1+xy)^b (1+xy)^(n-b).                  (4)
```

The reverse-time initial state is `(a,b)=(0,0)`. Formula (4) retains the
free time boundary. It must not be replaced by a stationary bulk row.

## Exact regional quotas and the row-set union

Return to `n=512`, 16 updates per region, and `h=4096`. Write
`C0=binom(8192,4096)`. Define the coefficient operators

```
B=[x^4096] K_x^16,
b_end=[x^4096] K_x^15 q_x.                           (5)
```

For `r` complete regions, let `Z_y` be the sum of `y^score(C)` over all
saturated witnesses. Independent exact quotas and fresh partitions give

```
E Z_y = C0^(-r) e_0 B^(r-1) b_end.                  (6)
```

The same formal variable is used throughout each region before taking
its coefficient. A new coefficient is taken in the next region.
For the full encoder, `r=256`. Since `y^delta<=Z_y`, a single union over
the eligible row set gives the valid finite bound

```
Pr[some obstruction with h=4096]
 <= y^(-104860) C0^(-255) e_0 B^255 b_end.            (7)
```

This contains no union over state supports or rank-profile parameters.
It is nevertheless a witness union, rather than the moment of the
maximum in (1). That distinction is quantitatively substantial below.

A dense matrix on 33,153 states has over one billion entries. Equations
(3)--(7) specify a finite exact operator, not a practical dense-matrix
algorithm. A useful implementation would need coefficient transforms,
matrix-free actions, or an independently proved positive comparison.

## What exact quotas cancel

For a bulk region staying in state `(0,0)`, the path weight is
`(1+x)^8192`. Its quota coefficient divided by `C0` is exactly one.

For a bulk region staying in state `(n,n)` and selecting every quartet,
the path weight is

```
y^-8192 * [((1+xy^2)^4-1)]^2048.
```

At coefficient `x^4096`, every power of `y` cancels. The normalized
coefficient is at most one; it is the quota probability that every
quartet meets `J`. These are bulk statements. The initial and terminal
factors still matter.

A common Bernoulli tilt before regional composition can destroy this
cancellation: the empty phase prefers `x=1`, whereas the full phase
prefers `x=y^-2`. Exact quota extraction is safer than selecting one
tilt and then taking a maximum over both phases.

One positive decomposition can preserve both pure paths. Let `D_m` be
diagonal, retaining only their `m`-step stay weights, and put
`R_m=K_x^m-D_m`. The remainder counts all other paths and is nonnegative.
For doubling,

```
R_(2m)=R_m^2+R_m D_m+D_m R_m.
```

Thus their exact quota diagonals can be restored after separately
bounding the positive remainder. This removes a false phase penalty;
it does not remove the witness entropy proved next.

## A deterministic obstruction to the proposed raw moment count

Consider a 64-update window with four exact half-row quotas, so
`sum_t |J_t|=16384`. Updates three through 64 contain at least 15,360
non-`H` slots. One residue class modulo three therefore contains at
least 5,120 such slots, spread over at least 1,280 nonempty quartets.

Allow `U_t` to contain any subset of those quartets and set every other
`U_t` empty. Keep `C_1=J_1`. A selected quartet with `j>=1` non-`H`
positions contributes `2j` capacity units, at most `j` own predecessors,
and four predecessors at each of two earlier updates. Its incremental
score is at least `j-8>=-7`. Different candidate quartets have disjoint
predecessor triples in time, or disjoint quartets at a common time.
Possible overlap with `C_1` only increases the score.

Consequently, for every realization of the labels and partitions,

```
Z_y >= product_over_candidates (1+y^(j-8))
    >= (1+y^-7)^1280.                               (8)
```

At the proposed marker `y=2^(1/8)`, this is larger than
`(3/2)^1280>2^640`. In particular, the exact raw count cannot establish
the target `E[y^delta]<=2^75` by bounding `E Z_y`. This conclusion needs
no averaging, numeric optimization, or asymptotic estimate.

The factor in (8) comes from optional negative-score witnesses. It does
not lower-bound the maximum deficit and does not disprove the desired
moment target. It rules out this unfiltered witness sum at that marker.
It also explains why exact quotas and pure-phase restoration alone are
insufficient. No impossibility claim is made for every marker `y`.

Larger markers remain possible. For example, at `y=2`, a sufficient
64-window moment exponent is 1,509.8125 bits, using the upper bound
`binom(8192,4096)<2^8192` and allocating forty failure bits. At `y=4`,
the corresponding exponent is 3,148.25 bits.
The sparse-family floor (8) is small at those markers. Thus (8) does not
justify abandoning the exact count without evaluating the relevant
larger-marker bounds. The following normalization also removes a
provable amount of harmless witness multiplicity.

## Positive normalization around a maximizing witness

Fix the labels and partitions. There is one binary selection variable
for each quartet at updates `t>=2` that contains `j>=1` non-`H` inputs.
Its addition contributes `2j` capacity and introduces between zero and
`j+8` source predecessors. Consequently, uniformly over all other bits,

```
score(after addition)-score(before addition) in [j-8,2j].   (9)
```

The bound remains valid near the time boundary, where fewer predecessors
are possible. Put `M_j=max(8-j,2j)`. Starting from a maximizing witness
and flipping any chosen subset of bits loses at most the sum of their
`M_j`. Summing this inequality over all flip subsets proves

```
y^delta <= Z_y / product_over_nonempty_quartets(1+y^-M_j).  (10)
```

This quotient is pointwise, before averaging any randomness. It neither
changes the code nor normalizes field amplitudes. It is a consequence
of bounded changes in the Hall score on the witness cube.

An asymmetric version retains both bounds in (9). Define

```
a_j=8-j,             b_j=2j,
A_j=y^-a_j,          B_j=y^-b_j,
w_(0,j)=(1-A_j)/(1-A_j*B_j),
w_(1,j)=(1-B_j)/(1-A_j*B_j).                         (11)
```

Both weights are positive for `y>1`, and

```
w_(0,j)+A_j*w_(1,j)=1,
B_j*w_(0,j)+w_(1,j)=1.                              (12)
```

For a reference witness `u`, a flip from absent to present loses at most
`a_j`; a flip from present to absent loses at most `b_j`. Telescoping
these bounds along any flip subset and summing with product weights
gives

```
y^score(u) <= sum_v y^score(v) product_i w_(v_i,j_i).
```

Taking `u` to maximize the score yields the positive majorant

```
y^delta <= Z_y^weighted
 := sum_v y^score(v) product_i w_(v_i,j_i).           (13)
```

This needs no multiplicity-one rule and no separate union over actual
state supports. It handles arbitrary interactions between quartet
choices, because (9) holds after every intermediate flip. The symmetric
weights from (10) and the asymmetric weights in (11) give different
valid comparisons; neither is asserted to dominate the other for every
realization.

The weighted kernel is explicit. Split (3) into its unselected and
selected polynomials `P0_x` and `P1_x`, respectively. Replace it by

```
P_x^weighted = sum_(j=0)^4 x^j {
   w_(0,j) [x^j]P0_x + w_(1,j) [x^j]P1_x
}.
```

For `j=0`, use unselected weight one; the selected coefficient is zero.
Use either both weights `1/(1+y^-M_j)` for (10), or (11) for (13).
All coefficients remain nonnegative. The state space, terminal vector,
regional coefficient extraction, and row-set prefactor in (4)--(7)
remain unchanged. The resulting expression bounds `E[y^delta]` rather
than merely discarding the normalization after averaging.

Exact quota cancellation, a tractable action of this operator, and a
certified moment remain open. The next finite control should preserve
the regional coefficients and compare relevant markers, including
`y=2` or `4`. An alternative is a direct count of maximizing or positive
connected witnesses. No scalar local maximum or new construction is
required by these interfaces.

## The asymmetric count is a random-cover moment

For a quartet with `j>=1` non-`H` inputs, define

```
g_j=y^(-2j)*w_(0,j)=1-w_(1,j),
g_0=1.
```

Conditional on the labels and partitions, independently select each
nonempty quartet at updates `t>=2` with probability `w_(1,j)`. Keep
`C_1=J_1`. The second equality in (12), together with exact total balance,
gives an alternative expression for (13):

```
Z_y^weighted = E_selection[y^(nT-|N^-(C)|)].          (16)
```

Indeed, moving the factor `y^(-2j)` onto each unselected quartet changes
`2|C|` into `2|J|=nT`. The remaining selection weights sum to one. The
exponent in (16) counts sources uncovered by the sampled input set.
An individual uncovered count need not bound the matching deficit.
The positive moment after averaging the selections is the majorant.

This representation removes witness entropy algebraically. It does not
permit the regional quotas to be discarded without a large loss. To
see this, extend the same random cover to independent half-density
input labels and put

```
g_bar=2^-4 sum_(j=0)^4 binom(4,j) g_j.
```

The event that every quartet at updates two through `T` is unselected
has probability `g_bar^(128*(T-1))`. On that event only `C_1=J_1` can
cover sources, so the uncovered count is `n(T-1)+|H_1|`. Therefore the
iid-label cover moment is at least

```
(g_bar*y^4)^(128*(T-1)) * ((1+y)/2)^n.               (17)
```

At `y=2`, direct exact values are

```
g_1=127/511,    g_2=21/341,    g_3=31/2047,    g_4=1/273.
```

They give `16*g_bar>2`. Thus (17) already exceeds `2^8064` for a
64-update window, before its terminal factor. The needed exponent is
only 1,509.8125. Conditioning this iid moment on regional quotas is a
valid inequality, but cannot certify the target. The same return from
the all-unselected state prevents a pointwise local comparison from
repairing that iid relaxation.

Under exact quotas, the factors `g_j=y^(-2j)w_(0,j)` instead cancel the
total capacity contribution. Keeping that cancellation remains essential
even after the random-cover reformulation.

A separate exact-quota audit excludes `y=sqrt(2)` for this asymmetric
majorant. Its mean uncovered count is at least `38144/25` in a 64-update
window. Jensen's inequality gives a moment exponent at least `762.88`,
above the sufficient target `690.59375`. The proof uses only own-`H`
sources whose two future quartets are unselected, plus terminal `H`
sources; it preserves the exact quota law through a nine-label sampling
comparison. This limitation concerns the asymmetric cover majorant at
that marker. It does not exclude `y=2`, `y=4`, or the true deficit
moment. The detailed moment targets are tracked in
[the flow-tail note](QUARTET_FIXED_ROW_FLOW_TAIL.md).

## A factorization for a possible matrix-free action

The local polynomial need not be raised as a six-variable polynomial.
For each `r=0,...,4`, let `S_r(x)` be the weighted selected action for a
quartet containing `r` old `B` coordinates. Its old labels are fixed;
the action sums its current non-`H` labels. Define

```
A(x)=sum_(j=0)^4 binom(4,j) x^j w_(0,j),
S_r(x)=y^-r sum_(j=1)^4 w_(1,j) x^j y^j
                 * sum_l binom(r,l) binom(4-r,j-l) y^l,
L_(k,l)=[z^l](sum_(r=0)^4 binom(4,r) S_r(x) z^r)^k
                    / binom(4k,l).                         (14)
```

Use zero for invalid binomial coefficients and `L_(0,0)=1`. The
unselected action is `y^-r A(x)`. Choose `k` selected quartets. Their
union is a uniform set of `4k` coordinates, and, conditional on that
union, its quartet partition is uniform. If the union contains `i`
old `A` positions and `l` old `B` positions, the new state is
`(4k,4k+a-i)`. Therefore

```
K_x[(a,b),(4k,4k+a-i)]
 = binom(m,k)/binom(n,4k) * binom(a,i) * A(x)^(m-k) * y^-b
   * sum_l binom(b-a,l-i) binom(n-b,4k-l) y^l L_(k,l).       (15)
```

Equation (15) is the same exact kernel as (3), including normalized
weights. It uses at most 129 powers of one univariate quartic, followed
by hypergeometric transforms. It does not by itself make the sixteen
compositions in (5) cheap. It supplies a structured alternative to
materializing a matrix with over one billion entries.

Independent mathematical review passed the unweighted kernel,
regional and terminal normalization, the sparse-witness floor, both
positive normalizations, and the factorization (14)--(15). An exact
toy checker is being prepared separately. No numerical evaluation of
this operator has been run here.

## A block comparison without a dense matrix

For a fixed positive vector `phi(a,b)`, define the normalized regional
operator and final regional vector by

```
T_h=[x^h]K_x^16/C0,
v_end=[x^h]K_x^15 q_x/C0.
```

All coefficients are nonnegative. Therefore, pointwise in the physical
state `(a,b)`,

```
(T_h phi)(a,b)
 <= inf_(x>0) x^-h (K_x^16 phi)(a,b)/C0.             (18)
```

The value of `x` may differ between physical rows after all sixteen
actions have been composed. Taking the minimum between transformed
basis coefficients, or before each separate step, is not justified by
(18). The same comparison applies to `v_end` with fifteen actions and
the terminal vector retained.

If these pointwise upper bounds prove `T_h phi<=lambda phi` and
`v_end<=c phi`, then positivity gives

```
E[Z_y^weighted] <= c*lambda^255*phi(0,0).
```

The row-set union consequently gives the failure bound
`C0*y^(-104860)*c*lambda^255*phi(0,0)`. This is a block comparison;
it retains the sixteen-update quota before taking a row maximum.

The matrix-free action has a second explicit form. For an arbitrary
vector `f`, put `p=4k`, `q=n-4k`, and define

```
P_k(V,Z)=sum_l binom(4k,l)L_(k,l)y^l V^l Z^(4k-l),
Q_k(X,W)=sum_d binom(n-4k,d)f(4k,4k+d)X^d W^(n-4k-d),
G_f(X,Y,Z)=sum_k binom(m,k)A(x)^(m-k)
                      *P_k(X+Y,Z)*Q_k(X,Y+Z).
```

Then the exact action is

```
(K_x f)(a,b)
 = y^-b * [X^a Y^(b-a) Z^(n-b)]G_f
                    /multinomial(n;a,b-a,n-b).       (19)
```

In particular, the factor `y^-b` is essential. Expanding the two
polynomials recovers the hypergeometric coefficient in (15), with
`d=a-i`. Equivalently it can be incorporated by substituting `X/y,Y/y`
before extracting coefficients.

Dehomogenizing `Z=1` leaves 129 two-dimensional polynomial products of
total degree at most 512. Their Fourier products can be summed before
one inverse transform. At fixed `x,y`, the transforms of `P_k` can be
cached; a side length near 540 would require roughly 300 MB for real
transforms or 600 MB for full complex transforms. Each later action
would transform 129 new `Q_k` polynomials and perform one inverse.
These are storage and operation counts, not measured runtimes.

An ordinary double-precision FFT is not a certified implementation of
(19). The old multinomial coefficients alone span roughly 800 binary
orders of magnitude. An absolute transform error that is harmless at
central coefficients can dominate a tail coefficient, then become a
large row error after division. Clipping negative outputs does not
produce an upper bound.

Possible certified implementations include multiple exponential gauges
with explicit transform error bounds and positive treatment of tail
rows, or upward fixed-point rounding with exact NTT/CRT convolutions.
Either must account for the further range introduced by markers and
the potential. The latter would need many machine-word moduli even
before those factors. No such implementation is specified or claimed
practical here. A bounded exact check of (19) at small width is the next
appropriate validation step; the dense production matrix is unnecessary.

Read-only source audit of
[quartet_dual_hall_checks.py](quartet_dual_hall_checks.py) passed the
physical predecessor indices, saturation, balanced and unbalanced
deficit identities, both normalizations, the random-cover powers,
positive local coefficients, and the forced first-update endpoint.
That checker does not yet implement the large action (19).

The parent ran its exact checks serially in 3.58 seconds. They passed:
208,099 input subsets; 1,004 dual/matching maxima; 3,873 canonical
witnesses; 3,959 single-flip and 20,159 multiple-flip bounds; 4,016
normalized majorants; 2,000 balanced cover identities and eight
unbalanced controls; 810 local coefficients; 152 endpoint coefficients;
16 weight identities; and 12 quota coefficients. Markers were `y=2,4`.
An earlier launcher used the wrong main-function key and stopped before
the checks; the corrected invocation passed without a source change.

## A 129-state bound using consecutive hole intersections

The random-cover moment (16) has a smaller positive upper bound. Let
`V_t=U_t^c` be the union of unselected quartets at updates `t>=2`.
Set both formal sets `V_(T+1),V_(T+2)` to the full coordinate set.
The uncovered sources at time `t` are exactly

```
C_t^c intersect V_(t+1) intersect V_(t+2).
```

Dropping the first condition gives the deterministic bound

```
number uncovered
 <= n + |V_T| + sum_(u=2)^(T-1) |V_u intersect V_(u+1)|.    (20)
```

Here `n` bounds the last layer, `|V_T|` bounds the preceding layer,
and the sum handles all earlier layers. The endpoint `n` is charged
once for the full 4096-update sequence, not once per region.

Use the asymmetric probabilities `g_j` from (16), and define

```
G(x)=sum_(j=0)^4 binom(4,j)g_j x^j,
P(x)=(1+x)^4-G(x),
w_b(x)=binom(128,b)G(x)^b P(x)^(128-b),
H_(a,b)(y)=E[y^|A intersect V|],
K_x^(pair)[a,b]=H_(a,b)(y)w_b(x),                    (21)
```

where `|A|=4a` is fixed and `V` is a uniform subset of size `4b`.
Thus `H_(a,b)` is an exact hypergeometric moment. Both indices range
from zero through 128. The coefficients of `G,P` are nonnegative;
`g_0=1` makes the constant coefficient of `P` zero.

Conditional on the number of unselected quartets and the current
non-`H` count, their union is a uniform set of its size. This follows
from the uniform partition and invariance under coordinate permutations.
Consequently (21) averages each consecutive intersection exactly for
the majorant (20), while `x` retains the regional label quota.

Start at state `a=128`, representing `V_(T+1)` full. The first step
charges `|V_T|`; the following steps charge the consecutive overlaps
down through `V_2`. Define

```
B_pair=[x^4096](K_x^(pair))^16,
q_pair=[x^4096](1+x)^512 (K_x^(pair))^15 1.           (22)
```

The last factor leaves input-one labels unconstrained before quota
extraction. There is no selection coin or partition choice for `U_1`.
Equations (16) and (20)--(22) give

```
E[y^delta] <= y^512 C0^-256 e_128 B_pair^255 q_pair,
Pr[some obstruction with h=4096]
 <= y^(512-104860) C0^-255 e_128 B_pair^255 q_pair.   (23)
```

Each entry of `B_pair` may be upper-bounded by
`inf_(x>0) x^-4096 (K_x^(pair))^16`, and the endpoint has the analogous
fifteen-step bound. Taking minima after composition is valid because
all values refer to the same physical population states. A finite grid
of positive `x` values gives a weaker but still valid upper bound.
No dense 33,153-state matrix or Fourier extraction is needed.

Raw hypergeometric moments can reach `y^512`, and the label-generating
weights can be larger. Numerical evaluation therefore needs explicit
scales or logarithmic construction. Scales must be restored to the
common physical basis before entrywise minima across different `x`.
The positive algebra avoids cancellation, but a proof certificate still
requires justified rounding bounds. No production value of (23) is
asserted in this note.

Read-only source audit of
[quartet_cover_holes.py](quartet_cover_holes.py) passed its overlap Gram
formula, positive quartet probabilities, column weights, regional and
terminal compositions, restored-scale minima, and global normalization.
It reverses both matrix axes and the endpoint vector once, so the
reused power routine's initial index zero represents full holes.
The final expression charges `512*log2(y)` once, subtracts 256 quota
normalizers, and adds exactly one row-set union. Its numerical output
is explicitly a floating diagnostic; parent-owned checks and controls
are separate from this static audit.

The parent subsequently ran the exact source checks successfully and
two production controls, serially. At `y=2`, the conditional moment
upper exponent was `233328.245660`, giving tail exponent
`+136653.419868`. At `y=4`, those exponents were `681103.083957` and
`+479568.258165`. Neither control closes the primary row class, and
neither is an optimization claim.

These values bound the pair-cover comparison after regional Cauchy
minima. They do not evaluate the exact maximum-flow moment, the exact
normalized Hall moment, or even the exact pair-cover quota coefficient.
The losses from independent flip normalization, dropping own-input
absence in (20), and coefficient comparison remain distinct. A further
argument should identify or avoid one of these losses, rather than
interpret these upper bounds as lower bounds on the desired moment.

The source
[quartet_cover_holes_quota.py](quartet_cover_holes_quota.py) subsequently
added one positive polynomial action for the exact regional quota.
Read-only audit passed its physical degree offsets, structural zeros,
positive convolution followed by overlap mixing, common coefficient
scale, terminal-label step, and Collatz composition. The parent control
returned tail exponent `+128423.718` in 5.324 seconds. Its trial
potential had exact-quota row-ratio exponents from approximately
`8964.642` through `9061.593`, compared with the Cauchy comparison's
`9093.937`. This remains a positive floating control, with no optimality
claim; it reduces the coefficient-comparison loss without closing (23).

The pair comparison at `y=2` is now excluded by an exact rational
certificate, not only by positive floating controls. Fix 96 hole
quartets at every real layer. A single overlap term and explicit
regional label histograms give integer lower exponents `291` for each
of the 4,094 diagonal overlaps, `4070` for each of 255 bulk quota
coefficients, and `4322` for the endpoint coefficient. With quota upper
exponent `8186`, the comparison-tail expression is at least

```
2^[896-104860+4094*291+255*4070+4322-255*8186]
 = 2^42132.
```

The first-overlap and external endpoint exponents sum to `896`.
Independent audit passed the formulas and scope. The parent verified
the rational witness in
[quartet_cover_holes_floor.py](quartet_cover_holes_floor.py), together
with 120 scalar coefficient, 60 regional entry, 60 full-path, and
60 tilt checks. This excludes the exact pair-cover comparison at that
marker even with optimal potentials and perfect quota extraction. It
does not lower-bound the actual failure probability or refute a sharper
cover moment.

The next count retains population transitions and actual own-input
coverage through
[QUARTET_ADAPTIVE_HALL_NORMALIZATION.md](QUARTET_ADAPTIVE_HALL_NORMALIZATION.md).
That note proves the adaptive normalization, records its passing exact
checks, and separates the unimplemented population action from another
scalar family that has also been excluded exactly. The quantitative
low-death route to a full amplitude bound is now stated in
[QUARTET_FLOW_FIELD_BRIDGE.md](QUARTET_FLOW_FIELD_BRIDGE.md); its required
flow tails remain unproved.
