# What the fixed-row cut certificate can cover

The [verified cut](RANDOM_PAIRING_CUT_BOUND.md) bounds every half-full
history for one selected target-row set. This continuation studies the
missing quantifier over row sets. It extends the deterministic statement
to a neighborhood, but rules out two straightforward ways to turn those
neighborhoods into full coverage. These are limitations of the proof
methods, not structural counterexamples to the code.

## A certified neighborhood

Keep the same sampled route and matching geometry. For a half-full
history H and a target set T, write

    g_T(H)=128|T|-C_T(H)-D(H)-83887.

Each outer row has at most256 covered coordinates. If T' removes a rows
from T and adds b rows, then for every H,

    g_T'(H)-g_T(H) <= 128(a+b).

Removing a row can increase the expression by at most256-128; adding a
row can increase it by at most128. The same inequality holds after taking
the maximum over admissible histories.

The existing bound max_H g_T(H)<=-47942 therefore implies

    max_H g_T'(H)<=-70

for every T' differing in at most374 memberships. In particular, it covers
every same-size set obtained by at most187 swaps. This is a deterministic
family statement for the sampled geometry. It is not a setup failure
probability, a field-failure bound, or a minimum-distance certificate.

There are also setup-independent bounds at extreme row counts. Positivity
of C_T and D gives g_T<=128|T|-83887. Conversely, the total visible count
is 2^20+D+J, and the complement of T covers at most256(8192-|T|)
positions. Therefore

    g_T <= 964689-128|T|-2D-J.

All target sets of size at most655 or at least7537 consequently have
g_T<=-47. The middle range still includes the difficult profiles. Absence
of a positive rank gap is only absence of this sufficient obstruction;
it does not establish distance.

## Why neighborhood enumeration is not enough

For a center of size1974, the number of same-size sets within187 swaps is

    V=sum_(j=0)^187 binom(1974,j) binom(6218,j).

There are binom(8192,1974) such sets in total. Exact integer bit-length
bounds show that covering them with neighborhoods of this radius requires
at least2^4425 centers, even before accounting for overlap. Thus exhaustive
coverage by fixed-set cut checks at the current margin is not viable.
This statement does not rule out stronger certificates with much larger
coverage regions or an analytical probability bound.

## A barrier to the obvious joint cut relaxation

One can make row membership variable inside the cut graph. If a_i=1
means row i is outside T, a coverage-OR node z contributes z*(1-a_i).
An edge from z to a_i represents that cost. Row penalties encode
-128|T|. Thus the unconditioned joint objective

    F(H,T)=C_T(H)+D(H)-128|T|

is still representable by a cut. However, the actual problem also requires
every state to be half-full and the target set to have a specified size s.
Enforcing both requirements only by linear Lagrange prices is too weak.

For s<=8192/2, consider a half-and-half mixture of two assignments:

- every state is empty and |T|=2s;
- every state is full and T is empty.

Every layer has mean weight256, the mean target size is s, and the mean
objective is -128s. For s>=8192/2, instead use an empty-state assignment
with all8192 target rows and a full-state assignment with2s-8192 targets.
Its mean objective is -128(8192-s).

The price terms cancel in either mixture. Since a minimum is no greater
than this mixture average, every joint relaxation using only those linear
prices has lower bound at most

    -128 min(s,8192-s).

At s1974 this is -252672, whereas excluding a positive gap would require
a bound of at least-83887. The conclusion is independent of how the
prices are optimized, and persists under a different death price because
both extreme assignments have zero deaths. It also allows different
prices for different time layers.

This is not an admissible half-full history and not a codeword. It is a
convex-mixture obstruction to this particular relaxation. Exact small
examples have hard-constrained minimum zero while the joint dual remains
negative. A joint proof must retain more than the average constraints;
adding row vertices and tuning more prices cannot suffice by itself.

## Screen a flow that is easier to analyze

A possible alternative to selected-set enumeration is an explicit flow
whose lower tail can be analyzed over setup. As a first control, assign
half of each target gate's unit capacity to each incoming coordinate
timeline. Discard its outgoing-coordinate connection. The two assigned
halves use no more than the original gate capacity, so a feasible flow
in these restricted timelines is feasible in the full graph.

On a timeline, give each state vertex supply a/c, each temporal edge
capacity1/2 in each direction, and service (T_i+T_partner)/2 from the next
gate. The exact two-state min-cut recurrence is small. With integer scale
c, supply a, edge capacity c/2, and service v, let d=F1-F0 and keep the
accumulated value min(F0,F1). One update uses

    g0=max(0,-d)+min(0,d+c/2)+a,
    g1=max(0,-d)+min(d,c/2)+v,
    d_next=g1-g0,
    increment=min(g0,g1).

The helper checks this recurrence against486 exhaustive short timelines.
An initial implementation omitted the normalization max(0,-d), producing
impossible negative flow values. That error was corrected before using
the screen; the erroneous outputs are not results about the design.

For an IID diagnostic only, take target probability p=1974/8192 and
service probabilities (1-p)^2, 2p(1-p), p^2 for service0,c/2,c. The final
vertex of each window has no next gate. At supply1/4, the extrapolated
mean dual lower cost is about164116 for32-step windows and177284 for
256-step windows. The corresponding selected-set rank threshold is168785.
The longer-window mean barely exceeds it and is far below the verified
adaptive cut lower bound216727.

These means are not bounds for the actual regional distribution or its
lower tail. Correlations, exact regional counts, the union over target
sets, and the field count are still unresolved. The screen does not
justify a full proof based on the equal-split allocation. It suggests
that the ability to share gate capacity between timelines matters.

## Validation, resources, and next step

Exact checks pass4096 sensitivity cases,360 mixture identities,375 joint
dual-price checks, and486 timeline recurrences. The toy hard-constrained
minima are zero for target sizes1,2,3; they are strictly above the joint
dual barrier.

Local free disk initially measured22.35GiB, then fell to12.37GiB before
the first test. The local guard stopped execution. Used the configured
Peach connection after checking63GiB free in /tmp and no numerical process
under the account. Standard-library sources were streamed over host-key-
pinned stdin, with30s timeout,256MiB address-space cap, Python -B, and a
remote20GiB disk guard. No remote files or packages were created. Runs
were serial; exact scope checks took0.10s and the corrected service screen
0.04s. No cleanup, system changes, or native benchmark was performed.

The next useful proof step is not more neighborhood centers or a joint
linear-price cut. Study a flow allocation that shares service according
to timeline demand, and identify a local interface that could support an
actual setup tail bound. Keep the equal-split scheme as a negative control.
Any new argument must still compose with the field constraints and cover
all relevant histories; the deterministic neighborhood result does not
replace those requirements. The construction remains unchanged.
