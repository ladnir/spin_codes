# A composable band of quartet histories at 45%

2026-10-26. Follow-up to [the 224/32 certificate](QUARTET_224_TWO_UPDATE.md).
The outer, route, quartet distribution, and acceptance rule are unchanged.

## Result

The proof now sums time-varying weights and turnover, rather than certifying
another prescribed schedule. It also covers every possible first activation
time. This is a restricted minimum-distance statement, not a whole-code proof.

Use the retained construction over J=GF(2^32), with K=2^20, N=2^21,
8192 RS[256,128] outer rows, 256 regional permutations, and 512 state
coordinates over 4096 updates. The recurrence, independently sampled
field-quartet maps, and two-matrix superregularity test are those of the
preceding certificate. All setup randomness is fixed for every message.

Let E be the event that some nonzero J-message has output weight at most
W=floor(0.45N)=943718 and satisfies the following conditions:

- Its exact number r of inactive outer rows is in [1850,2100].
- There is a first nonzero state Q_j. Every state Q_j,...,Q_4096 has weight
  in [216,232]. Earlier states are zero.
- At each transition after activation, the number of coordinates leaving
  the support is in [24,40]. Births need not equal deaths: their difference
  is the change in state weight.

There is no restriction on the labeled supports, amplitudes, or terminal
state beyond these conditions. Weights and turnover may change arbitrarily
at every update. The frozen rational replay proves

    Pr_setup[E] <= 2^(-106287443107715149/10995116277760)
                 < 2^-9666.

The bound includes the unions over all 251 row counts, their row sets, all
activation times, all allowed histories, and all messages. No selection of
a history independently of setup is assumed.

The statement is over J. As before, an eventual full minimum-distance bound
over J transfers to the extension-field encoder. A particular larger-field
message need not have a J-component with the same restricted support history.

## From one transition to arbitrary histories

Retain the preceding positive local coefficients. For one quartet, a counts
incoming nonzero coordinates, b births, d deaths, and t inactive-row labels.
Its permitted selected-zero polynomial is P(ell). A raw local summand is

    binom(4,a) p^t(1-p)^(4-t) (P(ell)/q^c)^theta,

where q=2^32 and c is covered inactive-row occupancy. The local constraints
and the definition of c are unchanged. Round the fractional power upward
on the same 48-bit dyadic grid. Introduce positive markers H,X,Y for a,b,d.
Let C_(a,b,d) be the resulting coefficient, summed over local supports and
inactive labels and multiplied by H^a X^b Y^d. Put

    Z=sum C_(a,b,d),    rho_(a,b,d)=C_(a,b,d)/Z.

For 128 independent draws from this artificial coefficient law, let
mu(u,b,d) be the probability of the indicated three totals. This is a
generating-function device, not an independence assumption about the actual
state or random partition. The true weighted row sum for a transition from
weight u to weight v with d deaths is bounded by

    k(u,v,d) = Z^128 mu(u,v-u+d,d)
               / (binom(512,u) H^u X^(v-u+d) Y^d).

The binomial denominator averages the actual uniform quartet partition.
For each fixed incoming support, symmetry gives this same row sum. Summing
over a shared intermediate support therefore composes these kernels.

Let A be total state weight and D total deaths after activation. Initially
suppose A+D>=W. A low-weight output needs at least z=A+D-W selected zeros.
The conditional field exponent and selected-zero marker give

    q^[theta(W-K+128r)] ell^(theta W)
      * ell^(-theta A) (q ell)^(-theta D).

Thus the weighted transfer on allowed weights is

    T(u,v)=sum_(d=24)^40 k(u,v,d) q^(-theta d) ell^[-theta(v+d)].

The count can also use this expression when z<0. In that case use the
zero-selection-free conditional count first. Every feasible local P has
constant coefficient one. Since q ell>=1, replacing z=0 by the negative
value only increases this upper bound. It does not request a negative
number of selected zeros.

Choose the positive potential f(w)=X^w ell^(theta w). It cancels all
next-weight markers in T(u,v)f(v)/f(u). Therefore

    sum_v T(u,v) f(v)/f(u)
      <= 17 max_(u,d)
         Z^128 sum_(v=216)^232 mu(u,v-u+d,d)
         / (binom(512,u) H^u ell^(theta u) [XY(q ell)^theta]^d)
      =: beta.

The maximum ranges over u in [216,232] and d in [24,40]. The sum over v
is performed before the remaining death-count maximum. Consequently

    T^h 1 <= beta^h f/min(f)

coordinatewise. This proves coverage of all changing histories in the band;
there is no additional uncounted choice of a weight or turnover sequence.

The initial activation contributes at most

    I = 17 max_(w=216)^232 binom(512,w) ((1-p)X)^w / min(f).

The frozen markers make f decreasing, so min(f)=f(232). The helper checks
this sign using directed logarithm intervals. The free terminal state is
already accounted for by division by min(f), not forced to return to zero.

## Delayed activation and the regional route

If the state remains zero, the recurrence forces the input and output to
be zero too. A nonzero message must activate by step 2048. Otherwise every
outer row has zero evaluations at the first 128 coordinates and hence is
the zero RS word. This argument uses the actual outer and regional route.

Fix an activation step j<=2048. There are h=4096-j+1 active updates and
R=ceil(h/16) touched regions. For an exact inactive-row set of size r, write

    pi_r=binom(8192,r) p^r(1-p)^(8192-r).

Generate independent Bernoulli-p labels and condition each region to have
exactly r inactive labels, as in the preceding proof. A complete zero-prefix
region has no restriction in the relaxed local label count. Its numerator
quota probability is exactly pi_r, cancelling its conditioning factor.
Only the R touched regions need the remaining pi_r^(-R) bound.

This cancellation does not discard the zero-prefix input equations. They
are already present in the shortened outer space used by the conditional
message count. It removes only an unnecessary label-conditioning penalty.

For this j and r, the family contribution is at most

    binom(8192,r) 2^theta pi_r^(-R)
      q^[theta(W-K+128r)] ell^(theta W) I beta^(h-1).

As a function of integer r, its logarithm is an affine function minus
(R-1)log binom(8192,r), hence convex. Both row-count endpoints therefore
bound the entire interval. The helper evaluates all 2048 activation times,
takes the largest directed upper bound, and pays factors 251 and 2048.
It does not infer coverage from sampled activation times or row counts.

The final certificate does not need the preceding Fourier quota cap.
That refinement remains available and is retained in the diagnostic mode.

## Frozen replay and checks

Only proof weights changed. The retained values are

    theta=7/40,
    p=112627/167924,
    H=2042267/1022811,
    X=97905/432928,
    Y=173812/965177,
    ell=13592/664659.

The bounded optimizer proposed these weights but reached its iteration cap.
Convergence is not needed for validity: the frozen mode performs no
optimization and checks the resulting rational bound directly.

The positive integer DP retains a 233-by-57-by-41 coefficient box, using
three uint64 arrays including the multiplication temporary, below 16 MiB.
All omitted degrees exceed queried nonnegative degrees. Coefficients and
each multiplication stage round upward on the 2^-31 grid, with scalar
overflow checks before every stage. Summed coefficient masses use integer
sixteenth powers for an upper logarithm on a 1/16-bit grid. All other logs
use the preceding rational interval implementation. No floating result is
trusted by the replay.

Run from the repository root:

    python -B research/workstreams/one_route_coupled_scan/coded_state/quartet_history_band.py --frozen

It asserts the exact final exponent. Independent checks include 16 small
coefficient cells, 1536 labeled partition cases giving nine varying-weight
transition identities, two complete-path potential inequalities, 406
binomial concavity inequalities, seven integer logarithm checks, five
zero-prefix quota identities, and nine negative-zero-count checks. The prior
224/32 frozen replay still gives its original 46% exponent, -104.10742.

Selected intermediate directed bounds, all for the unchanged construction:

| Covered family | Distance | Log2 upper bound |
|---|---:|---:|
| Weights 216--232, deaths 24--40, r=1850--2100, activation at step 1 | 45% | -33545.74 |
| Same, activation anywhere in steps 1--128 | 45% | -15639.51 |
| Same, all possible activation times, retuned proof markers | 45% | -9666.79 |
| Same first-step band using the original markers | 46% | +72650.25 |

The positive entry is an inconclusive upper bound, not a distance obstruction.
A separate zero-selection-free bound closes late-start bands very strongly;
it is a diagnostic fallback, not needed in the final frozen result.

All runs were serial, file-free, and protected by a 60-second watchdog and
the 20-GiB disk guard. The final replay took about 8.5 seconds. No benchmark,
implementation, library, paper, setup-distribution, commit, or push change.

## Remaining work

This closes one continuous band of histories, including all its activation
times, rather than the whole code. Missing cases include other row counts,
weights and turnover outside the band, paths that leave and return, and
full resets followed by reactivation. The 46% band remains open as well.

Next test a composable extension that includes band exits and returns.
Zero-state stretches must retain their outer shortening cost; a scalar
transfer that admits arbitrary quiet loops while forgetting active-row
requirements can be much too loose. Use cheap envelopes away from critical
weights and exact coefficients near them, rather than storing an unbounded
three-dimensional table. Keep the construction fixed until this test reveals
a genuine structural problem or a clearly identified counting limitation.
