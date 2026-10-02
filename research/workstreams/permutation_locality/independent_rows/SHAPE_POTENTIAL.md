# Shape-preserving continuation bounds

This experiment keeps the four-bit independent-row encoder unchanged.
It tests whether taking separate maxima for individual transition entries
loses useful dependencies. The motivation is the continuation-potential
method in Section 5 and Appendices G--K of the companion lifting paper.
Its odd-characteristic cancellation inequalities are not imported.

The code is `candidates/shape_potential.py`. It rebuilds the existing
outward local bounds, retains packet shapes through eight active windows,
and evaluates selected homogeneous support classes in binary64. The
numerical screen is not an outward certificate or a complete message cover.

## Result (2026-09-28)

The one-step continuation refinement helps, but does not close the
four-bit proof. At K=2^20, two updates, distance 9.25%, union support
u=200, tilt .056 and rho .9, the selected log2 expected-count upper
proposals are:

| Active groups q | Original finite bound | Intersect, then entrywise max | Linear common potential | Whole-shape common potential |
|---:|---:|---:|---:|---:|
| 80 | -288.4113 | -289.0488 | -286.5479 | -317.4912 |
| 96 | +1779.6714 | +1779.0252 | +1781.4630 | +1751.8553 |
| 128 | +8190.6733 | +8190.0227 | +8192.3634 | +8164.7589 |

Lower is better. Positive values are uninformative upper proposals, not
lower bounds on the bad-message count. Even the negative q=80 entry
concerns only the selected homogeneous support class. The original column
reproduces the earlier two-update screen.

Intersection before maximization saves about 0.64--0.65 bits. Whole-shape
continuation saves another 27.60--30.94 bits against the matching linear
common-potential control. After its finite prefactor, the improvement
against the original finite bound is 25.91--29.08 bits. These are reductions
in selected log upper bounds, not new whole-code security margins.

Each comparison uses the original auxiliary support probability: about
.92820728656634, .91447482006172, and .88785445682414, respectively.
The final screen uses 16 potential iterations; the best proposal is the
last one in each case. Coordinate-ratio spreads are below 2.8e-5, but
convergence is not used to justify a witness. The finite terminal factor
is retained (approximately 7761--8494).

An initial screen omitted the shape-specific single-packet activation and
used 12 iterations. It gave -317.4643, +1751.8828 and +8164.7796. Retaining
that exact activation leaves its maximized continuation unchanged, for
the structural reason below. Four more potential iterations improve the
proposals by less than 0.03 bits. Thus this omission does not explain the
remaining gap.

The candidate suite passes 63 tests, including the separate two-epoch
local calculation described below. No production code, parameter choice,
performance number, or completed two-bit certificate changed.

### Next bounded experiment

Test a small history-aware bound that remembers the packet weight that
created a fresh state, or an equivalent two-epoch block bound. Include
both returns to zero and nonzero-state persistence; a scalar return
formula alone does not bound the full continuation. Use the existing
valid bounds for paths not refined. Compare q=96 and q=128 before
attempting another full support cover or a larger census. If this also
gives only a small gain, prioritize inner-parameter retuning rather than
assuming a longer run analysis must close the proof.

## Local continuation inequality

Use the eleven nonnegative envelope coordinates of [BRIDGE.md](BRIDGE.md):
zero, fresh and mature mass, mature density, five uniform expansion classes,
and two mature tail bounds. Density and tail coordinates are auxiliary
bounds, not additional probability mass; their terminal weights are zero.

For an epoch with j active packet windows, a shape h is the sorted tuple
of its j packet weights in {1,2,3,4}. Let A[j,h] be a valid nonnegative
envelope matrix for that shape, including the factor rho for each
weight-four packet. The input-weight tilt is one. Window and lane choices
are averaged as in the existing local inequalities, conditional on the
incoming state and assigned packet weights.

Write x for an incoming row vector of envelope quantities. For each fixed
shape, the next envelope is bounded componentwise by x A[j,h]. For a
nonnegative column v of continuation coefficients, define componentwise

    F_j(v)_i = max_h sum_k A[j,h]_(i,k) v_k.

Then x A[j,h] v <= x F_j(v) for every shape h. This statement allows the
maximizing shape to differ across envelope coordinates; that is an upper
relaxation, not a claim that setup chooses shapes adaptively.

In contrast, the old entrywise envelope first forms
Abar[j]_(i,k) = max_h A[j,h]_(i,k). Positivity gives

    F_j(v) <= Abar[j] v.

The inequality can be strict because one shape must supply the whole
row contribution before the maximum is taken. Each F_j is monotone and
positively homogeneous. It need not be linear.

### Constructing the shape matrices

Start each A[j,h] from the current refined all-shape matrix T_j. Intersect
its entries only with bounds on the same component contribution for h:

- the exact zero-to-fresh activation exp(-tilt*b) rho^(b=4) for one packet;
- output moments from the exact window histogram, through eight windows;
- full feedback counts and translated density, through six windows;
- joint cancellation/output counts, through three windows.

The helper uses the same refinement functions as the existing screen,
with one complete shape record at a time. It retains the lazy/refresh
split when refining a zero-return entry. Density records already include
the lazy probability; the reciprocal all-one weight stays in the outer
measure. No alternative mass-based column decomposition is mixed in.

Every shape must be present: the counts for j=1,...,8 are
4,10,20,35,56,84,120,165. The all-one shape is never removed. Above eight
windows the family is the singleton containing T_j, so all occupancies
remain bounded. The construction checks 0 <= A[j,h] <= T_j entrywise.
It also requires the local-record metadata to match the requested tilt,
penalty and update count.

This produces two separate possible gains. First, intersecting bounds
before maximizing over h may improve the entrywise envelope itself.
Second, applying v before maximizing can retain dependencies between
entries. The experiment reports these gains separately.

There is a deliberate remaining relaxation: the fresh-state coordinate F
does not record the weight of the packet that created it. The local bound
may therefore combine a cheap activation by a weight-one packet with a
later bound derived from a weight-four fresh-state distribution. Retaining
the current packet's shape does not recover that missing history. In
particular, tightening the four individual Z-to-F activation entries
cannot by itself change F_1(v)_Z: that row has only one nonzero target,
and its maximum still occurs at packet weight one. A two-epoch bound or
fresh-state refinement must keep the earlier packet type to address this
specific relaxation. Source-state and longer-run correlations remain open.

## Without-replacement placement

One region has E=64 epochs with W=32 windows each. For a terminal
potential v, let V[e,r] bound the continuation over e epochs containing
r active packets, averaged over their uniform placement into distinct
windows. Initialize V[0,0]=v. The backward recurrence is

    V[e,r] = sum_j H(e,r,j) F_j(V[e-1,r-j]),
    H(e,r,j) = binom(W,j) binom((e-1)W,r-j) / binom(eW,r).

Only feasible j are included. The probability H is the exact
hypergeometric placement probability, not independent-slot sampling.
The local bound is uniform over the packet weights and their assignments
to the epoch, so the continuation may forget their remaining histogram.
This forgetting weakens the result but preserves the bound.

For q active outer groups and an auxiliary support probability p, define

    P_p(v) = sum_r binom(q,r) p^r (1-p)^(q-r) V[64,r].

Here p is a coefficient-extraction witness, not a change to setup.
Auxiliary independent Bernoulli indicators choose whether each group
occurs in each region. Averaging this experiment and then conditioning
each group to occur in exactly u regions costs beta_u(p)^(-q), where
beta_u(p)=binom(256,u)p^u(1-p)^(256-u), as in [BRIDGE.md](BRIDGE.md).
Uniformity over shapes permits the conditional histories needed by this
argument. It does not assume that shapes in the actual code are independent.

## Finite continuation witness

Let tau be the terminal column that selects only mass coordinates. Choose
v>0 with v_0=1 and lambda>0 such that P_p(v)<=lambda v. Set
C=max_i tau_i/v_i. Monotonicity and positive homogeneity imply

    P_p^256(tau) <= C lambda^256 v.

Thus the auxiliary moment from the zero initial state is bounded by
C lambda^256. No convergence claim is needed: any verified positive
super-solution suffices. The numerical implementation iterates positive
potentials to find proposals and retains the best finite bound, including C.
It does not treat P_p as an ordinary matrix.

For the homogeneous class with union size u in each of q active groups,
the existing outer measure a_rho(u) then supplies

    exp(tilt*D) C lambda^256 binom(2048,q)
        * (a_rho(u)/beta_u(p))^q

as a sufficient expected bad-message bound, provided every local inequality
and the finite super-solution are validated. D is the forbidden output
weight. The screen uses the existing valid shell/CDF cap for a_rho(u).

The initial potential is derived from the linear, shape-intersected region
matrix. A linear common-potential control separates the finite C-factor
from the new nonlinear improvement. The existing finite matrix moment is
also retained: a new sufficient witness need not beat that tighter finite
evaluation merely because it beats the linear common-potential control.

## Verification scope and reproduction

Small-geometry tests compare the new recurrence with exact enumeration of
slot subsets and every shape sequence. They check that singleton families
recover the old noncommuting matrix placement, and that the joint bound
does not exceed the entrywise envelope. Additional tests check complete
shape catalogs, terminal factors, and rejection of invalid geometry and
zero potential coordinates. The local construction uses outward Arb
coefficients; conversion and all global recurrence arithmetic are binary64.

```sh
python -B -m unittest discover -s research/workstreams/permutation_locality/independent_rows/candidates -p 'test_*.py'
python -B research/workstreams/permutation_locality/independent_rows/candidates/shape_potential.py --rounds 2 --tilt .056 --penalty .9 --groups 80 96 128 --support 200 --cutoff 193986 --maximum 8 --iterations 16
```

The selected points fix BCH[256,128], K=2^20, N=2^21, IMT(128,19),
two updates, and four-bit independent-row packets. They compare at 9.25%
distance, not against the old 10% cutoff. Each point keeps its original
optimized p for this first comparison; no new p/tilt/penalty grid is claimed.
Raw output stays under ignored `tmp/`, not in Git. Production defaults,
the paper, and the completed two-bit certificate are unchanged.

## Two-epoch return: an exact local building block

The companion paper suggests keeping activation, emitted weight and return
to zero together over a run. The smallest nontrivial SPIN case can be
enumerated directly, without changing the encoder. This calculation is
separate from the global continuation screen above.

Start in state zero and take two consecutive epochs with one nonzero
four-bit packet in each. Let x and x' be their 128-bit input words, of
packet weights a and b. Let E be the expansion map and B the feedback map.
The actual inner recursion is y=x+E(s), s'=T(s)+Bx. After the first
epoch the state is Bx; the second emits x'+E(Bx). With R independent
transvections, alpha=2^(-R), and S=2^19-1 nonzero states, the exact tilted
zero-return moment is

    rho^((a=4)+(b=4)) E[exp(-tilt*(a+|E(Bx)+x'|))
        * (alpha 1_{Bx=Bx'} + (1-alpha)/S)].

Both feedbacks are nonzero. The expectation averages independent uniform
window and lane-mask choices in the two epochs, conditional on a and b.
The refresh term is valid because the IMT refresh is uniform over all S
nonzero states and is independent of the input choices. This does not
treat the emitted output as independent of the feedback.

`candidates/two_epoch_return.py` checks all 480 single-packet inputs of
the actual maps: their feedbacks are nonzero and pairwise distinct. Hence
the lazy branch can return to zero only when x=x', which requires a=b.
Conditional on matching weights, that event has probability
1/(32 binom(4,a)); for unequal weights its probability is zero.
This map property was already used in the fresh-state analysis. The new
local helper retains both epoch weights and their combined emitted weight.

The helper enumerates all 230400 ordered packet pairs with integer output
histograms, then evaluates their moments with 192-bit Arb arithmetic.
At tilt .056, rho .9 and R=2, the largest local moment is approximately
0.0001400885051544294, attained at a=b=4. Separately maximizing the first
activation factor and the conditional second-step return factor gives
approximately 0.0001841286905312212. Their ratio is about 1.3144, or
0.3944 bits per such return. These are local moment factors, not failure
probabilities for the code or a global margin gain.

Thus cross-weight cancellation restrictions are real, but the worst
matching-weight return is not eliminated. This small case alone is not
evidence that a longer-run analysis will close the middle occupancies.
Other outgoing states and occupancies must also be bounded before using
two-epoch blocks in the global placement recurrence. Merely multiplying
the existing envelope matrices is not a new joint-state bound.

Tests compare the histogram expression with a direct enumeration of the
state-transition kernel on a small map, check exact rational return
probabilities at zero tilt, and reject maps without the required feedback
injectivity.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/two_epoch_return.py --tilt .056 --penalty .9 --rounds 2
```
