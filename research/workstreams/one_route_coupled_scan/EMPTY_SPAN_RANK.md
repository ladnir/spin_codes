# Empty gaps with a finite state-span budget

The sparse gates in [SPARSE_GATE.md](SPARSE_GATE.md) contain long intervals
with zero input. The support-only continuation can repeatedly charge a new
rank without remembering how many state forms remain. This note retains
that budget. It defines an empty-gap operator, not a whole-code certificate.

## 1. The state span cannot grow during an empty gap

Use the construction in [CONSTRUCTION.md](CONSTRUCTION.md). Write
F=GF(2^128), E=GF(256), q=256, and A=theta+E, where theta has degree 16
over E. Each diagonal entry is sampled independently from A. The wiring R
is one fixed w-cycle, with w=4 or w=8. For zero input the recurrence is

    Y_t = D_t Q_t,
    Q_(t+1) = (R + D_t) Q_t.

Every diagonal entry is nonzero. Also R+D_t is invertible for every allowed
diagonal: its determinant is 1+product_j(theta+c[t,j]), a nonzero polynomial
in theta of degree w<16. Thus an empty step preserves the span of the state
coordinate forms, not merely the nonzeroness of a particular state vector.

Fix the outer, the full route, and the inner coefficients preceding the gap.
Work over the common message space V from [PROOF.md](PROOF.md). Let

    U = span_F {Q_start[j] : 0 <= j < w} subset V*.

For every realization of the remaining coefficients, the state coordinates
at each subsequent empty step span exactly U. Keep one common space W of
imposed zero forms and all earlier nonzero requirements. As before, discard
a branch if W contains a required nonzero form. Its global potential is

    Phi = 1[branch retained] tau^(-dim W),       tau >= 1.

The rank weight tau need not equal or be smaller than q. Independent zero
equations cost 1/tau; exceptional coefficients retain probability 1/q.
See [the separated-weight argument](RANK_PRICE_SEPARATION.md).

Define the residual state rank by

    r = dim_F ((U+W)/W).

If the current state support is H, a retained live branch satisfies
1 <= r <= |H|. Every zero equation imposed inside the empty gap belongs to
U. An independent equation increases dim W by one and decreases r by one.
The potential ratio for this branch is therefore 1/tau.

**State-span lemma.** A retained live empty gap entering with residual rank
r_start has at most r_start-1 independent rank gains, regardless of its
length. At every complete step, its residual rank is at most the outgoing
state support size.

Indeed, invertibility preserves U throughout the gap. If residual rank
became zero, W would contain every entrance state form, contradicting an
entrance nonzero requirement. At a complete step, the outgoing coordinates
still span U; only the coordinates required nonzero can survive modulo W.

This lemma does not bound the number of coefficient charges. After the
independent rank budget is exhausted, a cancellation may still survive at
one exceptional coefficient value. Those probability costs must remain.

## 2. Sum the zero and nonzero branches before taking a maximum

Fix a live incoming support H and put K=R(H). Since D_t has nonzero entries,
the output support is exactly H. The next support contains the mandatory
coordinates M=H symmetric_difference K. Its optional coordinates are
O=H intersection K; all coordinates outside H union K remain zero.

Process the coordinates of O in a fixed order. At coordinate j, put

    u = (R Q_t)[j],       v = Q_t[j].

Both forms are required nonzero. The proposed next coordinate is u+a_j v.
Before exposing a_j, classify u and v in the current quotient by W. Let
F_nonzero(r) and F_zero(r) bound the remaining continuation for the two
support choices, relative to the current potential. The three cases give:

1. If u is outside W+F v, the zero branch adds one rank for every a_j.
   The sum of branches is bounded by

       F_nonzero(r) + tau^(-1) F_zero(r-1).

   This case requires r>=2. The two summands enumerate message-support
   branches; they are not a probabilistic split of coefficient values.

2. If u=kappa v modulo W and kappa belongs to A, only a_j=kappa permits
   the zero branch. That branch adds no rank. The bound is

       (1-1/q) F_nonzero(r) + (1/q) F_zero(r).

3. If u=kappa v modulo W and kappa is outside A, only the nonzero branch
   survives. Its bound is F_nonzero(r).

Take the maximum of these three bounds, not their sum. The actual quotient
has one classification. Earlier zero equations can change later
classifications, so recompute the maximum sequentially. If an earlier
equation kills u, v, or another required nonzero form, the actual branch
contributes zero. Ignoring such additional contradictions is conservative.

## 3. A Bellman operator for L empty steps

Fix 0<z<=1. Let v_0(r,H) be a nonnegative terminal bound, uniform over all
retained quotients and forms with residual rank r and state support H.
For a free terminal state, set v_0=1 on valid live states and on (0,empty).
Set it to zero on invalid states.

Define v_L(r,H) to bound the expected sum of terminal contributions,
weighted by z to the total output weight and by Phi_final/Phi_initial,
over L empty steps. This definition applies only to retained initial
branches. The expectation uses the original independent coefficient
distribution. For the quiet state set

    v_L(0,empty) = v_0(0,empty).

It emits zero and remains quiet. A live state cannot reach it during the
empty gap.

For L>=1 and live H, list O={j_1,...,j_m} in the chosen order. Define a
within-step continuation B_i(r,T), where i optional coordinates have been
processed and T contains M together with those selected nonzero. Residual
rank r includes the zero equations already imposed in this step. At i=m,
set

    B_m(r,T) = v_(L-1)(r,T)       if 1 <= r <= |T|,
               0                otherwise.

The rank test is essential: the complete next state spans the same
residual U, and only |T| nonzero coordinates remain. In particular, T=empty
cannot terminate a live gap.

For i<m, write j=j_(i+1) and

    A_r = B_(i+1)(r,T union {j}),
    Z_r = B_(i+1)(r,T).

Set

    B_i(r,T) = max {
        A_r,
        (1-1/q) A_r + (1/q) Z_r,
        A_r + tau^(-1) B_(i+1)(r-1,T)  [only when r>=2]
    }.

The L-step operator is

    v_L(r,H) = z^|H| B_0(r,M).

The output factor uses the incoming support H. Previous nonzero
requirements remain active even when zero equations are added later in
the step. Reveal the remaining, uncharged coefficients before continuing
to the next step; the terminal bound is uniform over their values.

## 4. Exposure and composition

At each optional coordinate, the unexposed a_j remains uniform on A,
conditional on the preceding exposed coefficients. W and the classification
above may depend on that preceding history, but not on a_j. A branch's
nonzero requirements are retained indicators, not conditioning events.
In particular, no survival probability is renormalized away.

For an independent equation, the global potential ratio is exactly 1/tau
before possible branch rejection. For a dependent equation, the unique
root has probability 1/q or is absent from A. These observations prove the
three local bounds. Their sequential maximum is uniform over all incoming
quotients. Induction over optional coordinates and then over empty steps
proves the Bellman bound by the tower property of expectation.

The same W persists outside the gap. Do not replace the global potential
tau^(-dim W) by tau^(-r), and do not reset the rank budget between empty
steps. At the next gap, recompute its entrance state span and residual
rank; preceding rank gains remain in W.

This operator uses ordinary next-state zero equations inside a genuinely
empty interval. A lookahead equation involving the next nonempty input can
lie outside U and must be handled separately at the boundary. Likewise,
active steps may enlarge the state span. Combining this operator with the
sparse gates requires an explicit boundary interface. Section 6 gives one
such interface, without asserting that its numerical bound closes.

## 5. A possible positive envelope for long gaps

Let T map a terminal vector f to the one-step Bellman bound defined above,
restricted to valid live states. T is monotone and positively homogeneous.
If a positive vector v and lambda>0 satisfy

    T(v) <= lambda v,

then every nonnegative terminal vector f satisfies

    T^L(f)(r,H) <= lambda^L v(r,H) max_(s,J) f(s,J)/v(s,J).

The maximum ranges over valid live states. This follows by bounding f by
a multiple of v and applying monotonicity at each step. Keep the quiet
state separate and exact, as before.

A candidate search could use the iteration v_next=1+T(v)/lambda. A finite
fixed point would satisfy a strict envelope inequality. Neither convergence
nor a useful lambda is asserted here. After obtaining a candidate, check
the displayed inequality with exact rational arithmetic or outward-rounded
intervals. The active-step interface must also account for its possible
outgoing residual ranks; each gap boundary may cost a potential ratio.

Unlike the earlier support matrices, T is nonlinear because of its maxima.
Do not insert it into matrix coefficient extraction or interchange it with
an average over routes. In fact T is a maximum of nonnegative linear
continuations and is convex, so T(E[f]) <= E[T(f)]; substituting the average
terminal vector can go in the wrong direction for an upper bound. One safe
option is to average bounds for fixed routes or insertion positions.
Another is a verified common envelope that bounds every relevant boundary
continuation before averaging. Neither option removes the need to retain
the actual paired-route normalization.

## 6. A two-state envelope for the full-support one-group gate

Assume a positive vector v and lambda have passed T(v)<=lambda v for the
chosen tau, z, and wiring. This section gives an active-step interface and
an analytic gate bound. It uses ordinary support transfers throughout,
without lookahead equations. It asserts no numerical value for the bound.

Let Omega consist of all pairs (r,H) with nonempty H and 1<=r<=|H|. Use two
basis functions on retained incoming branches: the quiet indicator and the
function equal to v(r,H) on live branches and zero on quiet branches.
A coefficient vector (c_q,c_l)^T therefore represents the continuation

    f(quiet) = c_q,       f(r,H) = c_l v(r,H).

All coefficients are nonnegative. The empty-step envelope acting on these
coefficients is the linear matrix

    E_env = [[1, 0], [0, lambda]].

### One singleton input can add at most one state dimension

Fix a step whose exact input support is {j}, and write x for its nonzero
input form. Let U be the span of the incoming state forms and let W be
the current common space of zero equations. If the incoming state has
residual rank r, every outgoing state form belongs to U+F x. The new zero
space W' contains W. Consequently its outgoing residual state rank r'
satisfies

    r' <= dim((U+F x+W')/W') <= r+1.

For a live outgoing support H', also 1<=r'<=|H'|. This bound uses the same
global W; it neither resets accumulated rank nor assumes x is independent
of the current state. Earlier nonzero requirements remain active.

Let M_I(z) be the ordinary support matrix from Section 1 of
[RANK_CONTINUATION.md](RANK_CONTINUATION.md). Its entries already include
output weights and the fresh-coefficient charges. Evaluate that matrix at
charge parameter psi=min(tau,q), not tau when tau>q. Define

    a = max_j z v(1,{j}),

    b = max_(j,(r,H) in Omega) M_{ {j} }(z)[H,empty] / v(r,H),

    d = max_(j,(r,H) in Omega) [
          sum_(H' nonempty) M_{ {j} }(z)[H,H']
            max_(1 <= s <= min(r+1,|H'|)) v(s,H')
        ] / v(r,H).

Then a uniform singleton-step envelope is

    A_env = [[0, a], [b, d]].

From a quiet state, the input produces one nonzero output and the same
singleton next state. Its residual rank is one. Thus quiet-to-quiet is
impossible and quiet-to-live is bounded by a. From a live state, M_I bounds
the expected relative global potential for each outgoing support. A quiet
outgoing state has continuation c_q. For a live outgoing state, the inner
maximum in d bounds its continuation for every possible residual rank.
Summing the support branches proves the live row of A_env.

These bounds are uniform over all previous nonzero requirements. If a
new equation contradicts one, the branch contributes zero. Dropping such
additional contradictions in the matrix envelope is conservative. The
maxima over j also make the envelope valid before averaging the singleton's
physical position.

### Average positions only after obtaining the linear envelope

For the width-four one-group gate, a region contains m=2048 steps and
exactly one nonzero input. Its step is uniform among those m positions;
its coordinate is uniform among the four positions. This follows from
the actual uniform pair permutation and independent swap. Regional routes
are independent. The envelopes above already hold for every fixed route
and every realized input form.

Define G_m(lambda)=sum_(ell=0)^(m-1) lambda^ell. Averaging the linear envelope
for each possible insertion position gives

    B_env = (1/m) sum_(ell=0)^(m-1)
                E_env^ell A_env E_env^(m-1-ell)

          = [[0, a G_m(lambda)/m],
             [b G_m(lambda)/m, d lambda^(m-1)]].

No nonlinear Bellman operator is moved through an expectation here. First
bound each fixed input-mask sequence by a product of common linear
envelopes; then average these nonnegative products over the actual route.
Independence of the regional insertion positions gives B_env^256.

Put v_min=min_(r,H in Omega) v(r,H). The terminal vector

    t_free = (1, 1/v_min)^T

dominates the constant-one continuation, so it includes every terminal
state without imposing a reset. Starting from the quiet state, the moment
for one group's exact full source support is at most

    (1,0) B_env^256 t_free.

The shortened message dimension of that class is 128 over F. There are
8192 choices of the group. Let P_full denote the probability that some
nonzero message in this union of exact support classes has weight at most
e. Then

    P_full <= 8192 tau^127 z^(-e) (1,0) B_env^256 t_free,

where e=floor(0.35N), N=2^21. This statement is conditional on any allowed
outer realization and sums all messages in this exact support class.
It does not cover partially supported groups or multiple active groups.
To turn the displayed expression into a finite certificate, verify its
constants and final bound with exact arithmetic or rigorous intervals.

## 7. Verified checkpoint

The [positive-vector checker](theory/check_empty_span_potential.py) verifies
all 32 valid live-state inequalities with exact fractions. Its float loop
only proposes a rational vector; no claim depends on convergence or the
accuracy of the proposal. The one-group checker uses z=499/500, tau=256,
and lambda=993/1000, with an explicit cap of 16384 proposal iterations.

The [full-support gate checker](theory/check_sparse_full_gate.py) also checks
the singleton boundary constants and evaluates the final positive products
with upward rounding at 80 and 120 decimal digits. Both establish
P_full<2^-849 at 35% output distance. The same envelope fails to establish
a useful bound for the 33-full-column class; see [the selected results](SPARSE_GATE.md#6-selected-finite-checks-2026-10-20).

The local Bellman implementation separately passes 840 exact GF16 fixtures
at width two, with one or two empty steps, exact residual rank one or two,
two terminal potentials, and rank weights tau=1,2,8 with q=2. These are bounded implementation tests. The
span lemma and coefficient argument, not the toy tests, justify width four.

The older failed 33-column screen is now superseded by a proof-weight
change. At tau=65536, z=499/500, lambda=993/1000, all fully populated
column choices in one group have combined failure below 2^-291. A different
verified envelope, z=999/1000 and lambda=4981/5000, closes the complete
33-column class including holes below 2^-219 after actual outer averaging.
The coefficient set still has exactly 256 elements in both calculations.
