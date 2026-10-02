# GF-Specific Bounds at Fixed Packet Occupancy

The transferred sparse proof permits adversarial packet weights. GF(16)
randomization supplies a stronger input distribution. This note specifies
new bounds that use that distribution while keeping the same encoder.
These operators are proof tools, not whole-code certificates on their own.

## Conditional Input Distribution

Fix a message before sampling setup. Condition on the active groups in a
region and on the number j of active packets assigned to an inner step.
The j positions form a uniform subset of the W=32 packet positions.
Their values are independent uniform nonzero four-bit vectors. The other
positions are zero. The fresh inner update is independent of this input
and the entering state. Packet values at different steps are independent,
conditional on the placement, because their GF multipliers are independent.

Write X for this step's 128-bit input. For entering state a, output and
next state are Aa+X and T(a)+CX. The output uses the entering state,
before the update. Here T composes r fresh independent IMT maps; each
map samples u uniformly from the nonzero states and v uniformly from
u's orthogonal complement, including zero, then sends a to a+u(v dot a).
These samples are independent of the input, entering state, and earlier
updates. Put L=2^19-1, alpha=2^-r, and beta=1-alpha.

For a fixed nonzero a, one map has distribution one-half point mass at
a plus one-half uniform nonzero state. Indeed, every nonzero b!=a
requires u=a+b and v dot a=1; exactly 2^17 admissible v satisfy this.
The transition probability to b is therefore 1/(2L). The uniform
nonzero distribution is stationary. Composing r such kernels gives
alpha times the point mass at a plus beta times that uniform distribution.
All maps fix zero. Thus r=2 gives (alpha,beta)=(1/4,3/4), whereas r=3
gives (1/8,7/8).

The bounds below retain alpha and beta symbolically. Changing r changes
the encoder and requires regenerating its local operators and full
certificate. The input and feedback counts are unchanged; a two-update
numeric certificate is not reused as a three-update certificate.

For 0<z<1, the analysis weights a transition by z^wt(Aa+X). All moments
below are over the conditional input distribution, with a fixed j.

## Output Moments and Feedback Atoms

For a fixed state a, let h_w count its expansion packets of weight w.
The output moment is

    H_a(j,z) = [x^j] product_{w=0}^4
                 (z^w + x ((1+z)^4-z^w)/15)^h_w / binom(W,j).

Here [x^j] extracts the coefficient of x^j. The inactive factor preserves
the expansion packet; the active factor sums over all 15 nonzero inputs.
The denominator averages the selected packet positions. The implementation
retains each expansion histogram's multiplicity, allowing both a maximum
over nonzero states and exact averages over all states or weight classes.
Replacing z by z^2 supplies the second moments used by Cauchy--Schwarz.

Feedback probabilities use an exact character sum. For a character u,
let c(u) be the number of packet positions on which C^T u is zero.
An active GF packet has character expectation 1 at these positions and
-1/15 at every other position. Thus

    f_u(j) = [x^j] (1+x)^c(u) (1-x/15)^(W-c(u)) / binom(W,j).

The exact zero-feedback probability is p0=2^-19 sum_u f_u(j). For each
nonzero target b, character orthogonality gives both bounds

    Pr[CX=b] <= 2^-19 sum_u |f_u(j)|,
    Pr[CX=b] <= 2^-19 sum_u |f_u(j)-c|,

for every constant c. The implementation chooses a weighted median c
and checks all arithmetic as rational numbers. Denote the resulting
common nonzero-atom upper bound by b_j. It is not a bound asserting that
the feedback is uniform.

### Exact Feedback Counts at Small Occupancies

The preceding absolute-value bound discards cancellation in the character
sum. For j<=8, `feedback_exact.py` instead computes every feedback count
by integer Walsh inversion. Define

    F_j(u) = [x^j] (1+15x)^c(u) (1-x)^(W-c(u)),
    n_j(b) = 2^-s sum_u (-1)^(u dot b) F_j(u).

Here n_j(b) counts the position/value choices giving feedback b. There
are D_j=binom(W,j)15^j equally likely choices. Thus p0=n_j(0)/D_j and
b_j=max_{b!=0} n_j(b)/D_j. The optional `--exact-feedback` flag uses
these exact values through eight packets and the previous bounds thereafter.

For the actual feedback map, the first counts are:

| j | n_j(0) | Largest nonzero count | D_j |
|---|---:|---:|---:|
| 1 | 0 | 1 | 480 |
| 2 | 0 | 5 | 111600 |
| 3 | 27 | 63 | 16740000 |
| 4 | 3637 | 3823 | 1820475000 |

These counts average uniformly over the selected positions. They are
not bounds conditional on one fixed position set. The fixed-occupancy
placement operator uses precisely this averaged distribution.

All Walsh arithmetic is integral. Before each transform, the implementation
bounds every intermediate by the sum of absolute input entries and checks
that the sum fits int64. It then checks divisibility, nonnegative counts,
total mass, and agreement with the independent zero-feedback calculation.

### Cancellation Must Still Emit Weight

Let d_A be the minimum nonzero expansion weight. With j active packets,
wt(X)<=4j. For every nonzero entering state a, the triangle inequality gives

    wt(Aa+X) >= max(0,d_A-4j).

Put e_j=z^max(0,d_A-4j) and nu_j=1-p0. The weighted lazy return to zero
is consequently at most e_j b_j, for every fixed nonzero a. For a uniform
nonzero entering state, its average is at most e_j nu_j/L. Indeed, each
input X has at most one cancelling nonzero state, namely a=CX.

The factor e_j is a pointwise upper bound on the output-weight factor,
including on the cancellation event. Multiplying it by the event
probability does not assume independence. The implementation takes the
minimum of these bounds, the total output moments, and the earlier
Cauchy--Schwarz bounds. For this expansion, d_A=48.

`occupancy_kernel.py` combines these quantities with the three-coordinate
refresh envelope. Its output and feedback are bounded jointly; their
marginal bounds are not multiplied under an independence assumption.

`return_moment.py` additionally counts feedback and wt((AC+I)X) jointly
for j<=4. This gives the exact weighted lazy-return sums before outward
rounding. Its entire feedback marginal is checked against Walsh inversion.
The shared-shuffle sparse verifier uses this census through j=3,
together with the rank/trimmed bounds and distance floor. Dense-cell
witnesses can select the same refinement explicitly.

For a fixed entering state a!=0, a lazy return requires CX=a, and its
emitted word is (AC+I)X. Let C_j(a,w) count the position/value choices
having feedback a and emitted weight w. Their number before fixing a
is D_j=binom(W,j)15^j. The exact weighted lazy-return moment is

    R_j(a,z) = sum_w C_j(a,w) z^w / D_j.

The arbitrary-state coordinate may use max_{a!=0} R_j(a,z); the uniform
density coordinate may use sum_{a!=0} R_j(a,z)/L. After introducing
the birth classes below, a class E_l may use max_{a in E_l} R_j(a,z).
These are bounds on a dominated measure, not assumptions about a
conditioned state's distribution. Multiplying by alpha and adding the
existing refreshed-return bound gives a new candidate for each zero
entry. `fixed_joint_return_probe.py` takes its minimum with the old entry.
All other entries remain unchanged.

The census checks integer mass and every feedback count independently
by Walsh inversion. At j=4 there are 1,820,475,000 choices for the actual
32-packet step; this is below the checked uint32 counter limit. Weighted
evaluation rounds powers upward to dyadic rationals and checks its uint64
accumulation bound. Small tests enumerate all four-packet inputs and
verify the refined zero entries against exact transitions.

## Preserve a Density Cap Across Lazy Steps

`density_kernel.py` adds a cap on the arbitrary nonzero component. Write
its total-mass bound as M and its pointwise bound as D: mu(a)<=D for every
nonzero a. D is not additional probability mass. The terminal coefficient
of D is zero.

The remaining component is bounded pointwise by U/L. Thus the coordinates
are (Z,M,D,U), where Z bounds zero mass. The old three-coordinate bound
instead loses the density information in M after every step.

Let d_A be the minimum nonzero expansion weight. For every nonzero entering
state, z^wt(Aa+X) is at most z^max(0,d_A-wt(X)). Define

    g_j = E[z^max(0,d_A-wt(X))].

The distribution of wt(X) is the j-fold convolution of (4,6,4,1)/15.
For a fixed outgoing state t on the lazy branch, entering state is
uniquely determined by X as a=t+CX. Its contribution is therefore at most
alpha D g_j. This propagates the density cap without replacing it by M.

For a lazy return to zero, a=CX must be nonzero. The multiplier of D is
bounded by alpha times

    k_j = min(g_j, (1-p0) z^max(0,d_A-4j)).

The two alternatives bound the same weighted event. In particular, the
first is not multiplied by 1-p0: the average output moment and feedback
event need not be independent. The uniform component gives the same
lazy-zero and lazy-density bounds divided by L.

The zero-source weighted feedback atom is at most z^j b_j, because every
selected packet has positive weight. Its total moment is
(((1+z)^4-1)/15)^j. Upper and lower bounds on the weighted zero-feedback
event give the separate zero and nonzero mass entries. The refresh
contribution uses the output moment and the pointwise 1/L refresh bound.

## Exact One-Packet Refinement

`single_packet.py` enumerates all 480 position/value choices. For every
outgoing target t it computes an outward dyadic sum of

    E[z^wt(A(t+CX)+X) 1_{t+CX != 0}].

It records the maximum for nonzero t, the maximum including zero, and the
value at zero. These refine the one-packet density and cancellation bounds.
The sum excludes zero entering state; that state has its separate row.
Integer summation is exact and checked to fit uint64. Python-integer
replays check the maximizing target and selected control targets.

For j>=2, condition on the other j-1 packets. The remaining packet can
use W-j+1 positions. Its restricted average is at most W/(W-j+1) times
the unrestricted nonnegative average. Shift t by the other packets'
feedback and bound their output perturbation by z^-wt(X_other). Their
values remain independent uniform nonzero. Consequently, the one-packet
maximum including target zero, multiplied by

    W/(W-j+1) * (((1+1/z)^4-1)/15)^(j-1),

bounds both the lazy-density and cancellation averages. This bound may be
loose; it is combined with the preceding valid bounds by taking a minimum.

## Preserve Expansion Weight During Empty Steps

The `gf-classes` variant replaces U by one component U_v per expansion
weight v. It bounds density by U_v/n_v within the n_v states of that class.
Moments for the class use its actual histogram multiplicities.

An empty lazy step keeps the same state and hence the same class. Its
contribution remains in U_v, multiplied by alpha z^v. A refresh sends
the weighted mass to class w with envelope factor beta n_w/L. This avoids
turning an averaged state into an arbitrary minimum-weight state after
an empty lazy step.

On a nonempty step, the lazy component enters M with its class-average
output moment and D with a bound divided by n_v. A class-average return
to zero has event probability at most min(b_j,(1-p0)/n_v). Both a maximum
output-weight bound and Cauchy--Schwarz with the class second moment apply.
These statements use pointwise domination, not uniformity of a state after
conditioning on earlier outputs.

## Retain the Class of a Newly Activated State

The separate `gf-birth-classes` variant records the expansion weight of
a state created from zero. Unlike `gf-classes`, it makes no density claim
within a weight class. Its coordinates are Z, M, U, and one mass bound
F_l for each set E_l={a!=0:wt(Aa)=l}. For the actual expansion, the five
values of l are 48, 56, 64, 72, and 80. Every coordinate represents mass,
so every terminal coefficient is one.

Condition on exactly j active packet positions, with the distribution
specified above. The weighted mass born into E_l is

    nu_{l,j}(z) = E[z^wt(X) 1_{CX in E_l} | J=j].

To evaluate it, let f_l indicate E_l and let hat(f_l) be its unnormalized
Walsh transform. Put S=2^19. For each character u, let n_r(u) count the
packet positions where C^T u has weight r. Define

    a_r(z) = ((1+z)^(4-r) (1-z)^r - 1) / 15,
    P_u(t,z) = product_{r=0}^4 (1+t a_r(z))^n_r(u).

Character orthogonality gives

    nu_{l,j}(z) = sum_u hat(f_l)(u) [t^j] P_u(t,z)
                  / (S binom(W,j)).

Each coefficient sums over the selected positions; each a_r already
averages the 15 nonzero packet values. The signed census is integral.
Arb encloses the signed polynomial sum before the final upper rounding.
The zero-source row uses these class masses instead of arbitrary
nonzero mass. Its zero-state entry remains unchanged.

For a source in E_l, let H_{l,j} be the maximum output moment over that
class, and let R_{l,j} be its maximum rank/trimmed return bound. Put
h_{l,j}=z^l for j=0 and zero otherwise. The new row has entries

    Z:   min(old arbitrary-state Z entry, H_{l,j},
             alpha R_{l,j} + beta H_{l,j}/L),
    M:   alpha (H_{l,j}-h_{l,j}),
    U:   beta H_{l,j},
    F_l: alpha h_{l,j}.

The empty lazy step preserves the class. A nonempty lazy step can leave
the class and is charged to M. The overlap with the separate return
bound only enlarges the envelope. The existing M and U rows are retained.

`occupancy_birth_classes.py` implements this operator. Exact small tests
check the class masses and mixed sequences of empty and active steps,
including zero and rank-one feedback maps. At 8.7%, q=48, a search over
eleven output tilts covers all supports with 175.1443 bits of margin.
Fresh 256- and 384-bit runs agree to the displayed precision. The original
GF rank operator also closes with this tilt grid, giving 63.0187 bits.
Thus both the grid and the class refinement help; the earlier sparse
failure was not evidence against the encoder. Different stopping points
mean these margins do not measure the refinement's optimal gain.
The first full 8.7% replay passed every sparse occupancy except q=4;
its coarse output-tilt grid missed the assigned stopping budget there.
A fine grid later verifies q=4 at the stronger 9% cutoff, with 89.2548
bits of margin at 384-bit precision. Fresh whole-code assembly at 9%
uses a separate fine grid for q=2--8 and regenerates the rest with the
broader grid. The shared-shuffle whole-code certificates at 5% and 6%
use this variant. Their separate outer counts and complete coverage are
recorded in `SHARED_SHUFFLE.md`.

## Composition and Verification

The existing placement polynomial combines the local operators without
assuming independent slot occupancies. If T_j is the local matrix, a
region containing q active groups is bounded by

    [x^q] (sum_{j=0}^{32} binom(32,j) x^j T_j)^64 / binom(2048,q).

The exact noncommuting-matrix placement tests still apply. For a fixed
outer ensemble, these local refinements leave its union-support CDF
unchanged. The shared-shuffle route must use its own four-row CDF, not
the independent-row CDF. Combine that CDF with a full support cover.
Every occupancy omitted by a dense cover must be verified; selected
support vectors or occupancies do not establish a whole-code certificate.

The small-model tests enumerate every eight-bit input and every entering
three-bit state, for one, two, and three updates. They check exact feedback
probabilities, output-moment polynomials, and eight successive bounded
transitions. Class tests check each class separately. These tests support
the implementation; the conditional-distribution and domination arguments
above justify applying the operators to the actual construction.
