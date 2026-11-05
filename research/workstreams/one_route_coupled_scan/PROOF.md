# Proof interfaces and structural checks

The first task is to establish a useful whole-code bound for the
[specified construction](CONSTRUCTION.md), not to assume refreshed states
are uniform over F. The results below are local lemmas and necessary tests,
with several partial support classes now bounded separately. They do not certify
35% distance for the complete code.

## One common space of message constraints

Fix an outer realization, a route, and a proposed set of active source groups.
Let V be the corresponding F-linear message space. For s groups, dim V=128s.
Each exposed input, state, and output coordinate is a linear form on V.

A branch specifies zero and nonzero coordinates. Let W be the span of imposed
zero forms, and keep every previously required nonzero form in a list B.
Discard a branch if any form in B belongs to W. This is a necessary condition
for a common nonzero witness, not a sufficient condition.

For any tau >= 1, use the potential

    Phi = 1[branch retained] tau^(-dim W).

A nonzero witness implies dim W <= dim V - 1. Thus its existence indicator
is bounded by tau^(dim V - 1) times the final potential, summed over the
branches that cover the event. This counts constraints, not |F|^(dim V)
field-valued messages. Keep the same W across all steps.

## Fixed-mask one-step lemma

Fix the exposed history before D_t, including input and current-state forms.
Fix proposed current-state, output, and next-state support masks H, J, H'.
All current nonzero state forms are required to remain nonzero. Write

    r_j = Q[j],   z_j = (R Q)[j],
    y_j = X[j] + a_j r_j,
    q'_j = z_j + y_j.

First impose coefficient-independent zero equations:

- If r_j is required zero, the output and next-state forms reduce to X[j]
  and z_j+X[j]. Impose whichever are required zero.
- If r_j is required nonzero and both y_j and q'_j are required zero,
  impose z_j=0. If z_j was required nonzero, discard the branch.

Let d be the rank gained in this preliminary step. All these forms are known
before exposing current coefficients. Retain all earlier nonzero requirements.

For each j in H for which at least one of y_j and q'_j is required zero,
there remains exactly one equation

    u_j + a_j r_j = 0.

When both zeros were requested, the equations coincide after z_j=0. Do not
charge two coefficients or two ranks. Let p count these charged coordinates.
Expose their independent coefficients sequentially in the common quotient.

If u_j is outside W+F r_j, the equation adds one rank for every a_j. If
u_j=kappa r_j modulo W, every coefficient except a_j=kappa forces the required
r_j into W and kills the branch. The one potentially surviving value has
probability at most 1/256. An already killed r_j discards the branch.

Each charged coordinate therefore contracts expected potential by at most
1/psi, where psi=min(tau,256). Consequently, conditional on the preceding history,

    E[Phi_after] <= tau^(-d) psi^(-p) Phi_before.

For tau<=256 this is the original tau^(-d-p) bound. The
[separated-weight refinement](RANK_PRICE_SEPARATION.md) permits larger tau
when independent equations and exceptional roots retain their distinct costs.

Reveal any uncharged coefficients without imposing zero equations, and append
the newly required nonzero forms only after their coefficients are exposed.
Doing so cannot increase the potential. Earlier rank gains can alter later
dependence classifications; sequential exposure handles that dependence.

Nonzero requirements restrict the branches being counted. They must not be
used to condition and renormalize the setup distribution. The lemma does not
assert that p equations are jointly independent. It combines rank gains with
exceptional-coefficient probabilities.

The next step is a continuation bound uniform over every incoming quotient.
Support masks alone do not determine the rank d. Any discarded deterministic
rank must be an explicit relaxation, with its loss assessed.

## Exact reset obstruction

Fix every setup choice and restrict the message to s outer groups, of dimension
D=128s over F. Mark every width-w input batch that could contain an output
from those groups. All other complete batches have input identically zero
on this message space, regardless of which message is chosen.

An initial empty prefix emits zero without constraints. For any later empty
interval, imposing Q=0 at its beginning uses at most w homogeneous linear
equations and makes the whole interval emit zero. Choose up to

    r = floor((D-1)/w)

such intervals, or all intervals if fewer exist. A nonzero message remains
because the total rank is below D. Overall injectivity ensures its output is
nonzero. If the chosen intervals have lengths L_i in complete batches and
the initial empty prefix has length L_pre, then

    d_min <= N - w (L_pre + sum_i L_i).

Choosing the longest intervals gives the strongest bound of this form. The
bound is independent of coefficient distribution. Dependence between reset
equations only strengthens the existence argument.

For one group, the guaranteed reset counts are 127, 63, 31, and 15 at widths
one, two, four, and eight. These are dimension counts, not measured distance
ceilings. For w>1, scalar-coordinate gaps are invalid: use whole empty input
batches and account for coordinates that share a batch. Include the terminal
empty suffix as an optional interval; do not require the final state itself
to vanish.

A failed setup found by this upper bound is genuine. Passing it is only a
necessary screen. A random sampled route is not a distribution-wide proof.

## Sparse proof gate and its limitation

One active outer group supplies exactly one allowed input coordinate per
region. Since w divides g, a batch then has at most one nonzero input.

Suppose such a batch extinguishes a nonzero incoming state. If its possible
input coordinate is j, every other coordinate r obeys

    Q[r-1] + a_r Q[r] = 0.

Indices are cyclic. Removing edge j leaves a path connecting all coordinates.
Because every a_r is nonzero, one zero state coordinate would force all state
coordinates to zero. Thus any nonzero state that this batch extinguishes
must have full support.

For a fixed incoming full-support mask and requested zero outgoing state,
all w coefficients are charged by the one-step lemma. The expected rank
potential contracts by at most psi^(-w), even if state forms are dependent.
If the input batch is empty, extinction is impossible.

These local statements do not suffice globally. A cluster of 32 groups
already [rules out width two](WIDTH_TWO_OBSTRUCTION.md) at 35%. Future
continuation bounds must handle group unions from the start, not certify
isolated groups while ignoring this obstruction. Include all possible outer
messages and the terminal state. A fixed-message transition is insufficient.

## Paired inputs and state wiring

Let H be the incoming state support, and I the input support of a reset step.
Since Q'=0 implies X=(R+D)Q, necessarily

    H symmetric_difference R(H) subset I subset H union R(H).

The output support is exactly R(H). For the natural cycle, two adjacent input
coordinates can extinguish a singleton state: if Q[j]=u is the only nonzero
state, choose X[j]=a_j u and X[j+1]=u. The output has just one nonzero record.
Routed physical pairs can supply precisely these adjacent coordinates.

This is a local reset pattern, not a whole-code low-weight witness. Its rank
charge is one fresh coefficient, plus any deterministic equations, not w.

An alternative fixed wiring puts each physical pair at opposite points of
the cycle: cycle orders (0,2,1,3) for w=4 and (0,2,4,6,1,3,5,7) for w=8.
The determinant proof is unchanged. For an input confined to one physical
pair, a proper nonempty reset support must have at least w/2 coordinates:
its cyclic boundary lies in that antipodal pair. This increases the fresh
reset charge to at least w/2. These are proposed wiring variants, not changes
to the initial natural-cycle specification or whole-code certificates.

## What transfers and what remains open

The companion BAA proof supplies coefficient-root counting and common-quotient
rank accounting. Its global permutation normalizations do not describe the
paired regional route. Its prime-field numerical tables are not certificates
for this code. Some prime-specific cancellation refinements need new proofs.

The existing grouped outer remains injective over F, but its GF256 point-density
comparison cannot be relabeled as a uniform-F distribution. A whole-code bound
must use justified shortened dimensions, actual ranks, or an explicit comparison
that preserves the correlations it needs.

The [wider reset audit](WIDTH_RESET_BOUNDS.md) now accounts for batch
collisions and excludes width four at 45%, without excluding 35%. Exact
support-mask enumeration confirms the single-pair wiring benefit, but also
shows that several active pairs can restore cheap local reset patterns.

The [rank continuation](RANK_CONTINUATION.md) supplies a composable finite
operator and the actual regional normalization. Its lookahead lemma charges
a repeated singleton-localization pattern that the first-order operator
misses. The GF256 support-containment bridge supplies a valid way to bound
larger-field shortened dimensions, but its numerical slack is unresolved.

Open obligations are:

1. Extend the partial support-class bounds to the remaining classes.
2. Preserve actual outer support correlations and quantify the bridge's slack.
3. Implement regional bounds with lookahead across boundaries and a free terminal state.
4. Cover partially active states, multiple active pairs, and all group unions.
5. Audit other low-charge histories, including consequences of R^w=I.
6. Evaluate a complete ideal-setup failure bound before any distance claim.

The baseline and its numerical certificates remain unchanged.

## Partial closures and the remaining sparse obstruction

The [dense-support argument](DENSE_SUPPORT_BOUND.md) covers every exact outer
support missing at most 1% of its coordinates. Its failure bound is below
2^-167780 at 35% output distance, uniformly over any fixed injective outer
and route. It uses the ordinary rank transfer and an exact integer comparison.

For sparse inputs, the [empty-span lemma](EMPTY_SPAN_RANK.md) retains the
fact that an empty interval cannot grow the span of the state forms. At
width four, at most three independent rank gains can survive throughout
such an interval. Further cancellation charges must represent exceptional
coefficient values. Joint branch accounting yields a positive envelope
that closes the exact full-support one-group class below 2^-849, summing
all 8192 group choices. The initial lookahead-only numerical screen did not
close that class.

Separating the rank price tau from the actual coefficient count q=256
closes more classes without changing the encoder. At tau=65536, the union
over all 33..64 fully populated columns in one group has failure below
2^-291. The bound retains all group and column choices and a free terminal.

The [partial-column moment argument](PARTIAL_COLUMN_MOMENTS.md) also closes
every support occupying exactly 33 columns in one group, including scalar
holes, below 2^-219. It counts pairs of GF256 codewords under the same
sampled outer, rather than treating their random mixers independently.
These classes overlap; their margins are not whole-code guarantees.

The [all-one-group argument](ALL_ONE_GROUP.md) now covers every scalar
support in one group below 2^-147. It replaces the constant pair rank of
the 33-column argument with polynomial rank profiles. A line-constraint
lemma counts every rank-one evaluation, including those beyond an RS
information set. Four exactly verified envelopes bound the same h-classes;
taking the smallest complete bound per h and then summing closes the class.

For [two adjacent groups](TWO_GROUP_ROUTE.md), the actual route alternates
between a full physical pair and two distinct singleton pairs. The latter
can share a batch with probability 1/4095. Retaining that collision and
the fresh source pivots in the first information-set columns is valid,
but its first broad upper bound remains uninformative.

The [endpoint-tail refinement](ENDPOINT_TAIL.md) further retains the state
support required by a cheap paired reset. It applies a bounded empty-gap
Bellman operator to the reset terminal before collapsing to a common linear
envelope. Positive series products then sum the actual ordered gaps,
including zero gaps, collisions, both pair parities, and the free terminal.
At 35%, the union for two adjacent groups on the same 33..35 fully populated
columns is below 2^-117. Larger shared column counts, different column sets,
partial columns, nonadjacent pairs, and larger group unions remain open.

The [exact-support bridge](EXACT_SUPPORT_BRIDGE.md) gives a second way to
relate large-field shortened dimensions to actual GF256 outer words. Keep
the older containment bound too: with the same point cap, it dominates the
uniform-factor version of the new bridge. Neither bridge transfers the
old inner's complete certificate to this construction.

## Large two-group diagnosis

The small [paired-rank operator](PAIRED_RANK_GATE.md) retains the span of
the old state together with the current inputs. Exact toy-field tests
validate its local cases. At tau=1024, z=.998, and lambda=.993, its broad
adjacent/shared-full-column bound remains about 9.928e207. It is not a
counterexample or a certificate.

The [lazy-outer argument](LAZY_OUTER_COEFFICIENTS.md) uses still-unrevealed
Toeplitz diagonals to charge later input/state matches. A raw-matrix
experiment is conditioned back to the actual invertible distribution once,
after the complete calculation. Its endpoint-aware broad bound is about
5.751e162 at those markers. Removing live-to-quiet transitions diagnostically
only changes the log2 bound from 540.676 to 535.744. Thus reset tails alone
are not the main remaining loss at this point.

The [online-route recurrence](ONLINE_ROUTE.md) supplies a different
conditional-expectation argument. It exposes pair labels and swaps before
each local step, then uses a suffix bound uniform over the remaining labeled
forms and common quotient. It does not exchange a nonlinear operator with
the old fixed-route average. No lookahead constraints are used.

A first combined online/lazy screen retains rank/support within each
region and collapses only at regional boundaries. Its floating-point
log2 upper-bound calculation is about 286.351 for the same broad union,
still uninformative. It changes no code operation. This screen is not a
directed certificate, and the lazy filtration has separate obligations.
Before enlarging the source class, test whether retaining residual rank
between regions removes the large boundary loss. Partial-column, unequal-
column, nonadjacent, and larger-group classes remain unresolved.
