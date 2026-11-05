# Quiet constraints with a bounded set of live batches

The retained 40% all-quiet receipt does not, by itself, force a low-weight
nonquiet word to have many live batches. A valid extension must account for
both the omitted quiet equations and the difference between input and output.
This note gives that reduction for the exact
[small-state proposal](SMALL_STATE_PROOF_GATE.md). No new probability bound
has been evaluated here, and the frozen quiet receipt is unchanged.

## Which quiet equations can fail

Fix the outer mixers and paired route. Let C_out be their routed outer
code over F=GF(2^32), obtained by scalar extension from E=GF256. Write
b=64, s=16, T=N/b=32768. For a fixed inner setup, the recurrence is

    Y_t=X_t+A Q_t,
    Q_(t+1)=R_t(Q_t+C X_t),       Q_0=0,

where every R_t=D2_t H16 D1_t is invertible. The last state is free.
Define the incoming-live set of a trajectory by

    L={t in {1,...,T-1}: Q_t!=0},       a=|L|.

For t<T-1, if both t and t+1 are outside L, invertibility gives C X_t=0.
Consequently all possibly omitted quiet constraints lie in

    J(L)=(L union (L-1)) intersect {0,...,T-2}.                (1)

Here L-1={t-1:t in L}; it is a set shift, not a cardinality subtraction.
Also Y_t=X_t outside L, so

    wt(X) <= wt(Y)+b a.                                     (2)

These statements hold for every setup and message, including a message
selected after setup.

For a fixed J subset {0,...,T-2}, define the partially quiet outer subcode

    C_J={X in C_out : C X_t=0 for t outside J, t<T-1}.         (3)

All equations defining C_J have coefficients in E. Thus C_J over F is the
scalar extension of the corresponding E-code. Its symbol distances over
the two fields are equal: a nonzero E-component of an F-word has contained
support, and every E-word embeds in F.

For e=734003, (1)--(3) imply the event inclusion

    {some nonzero encoded word has wt(Y)<=e and live set L}
      subset {d_min(C_J(L))<=e+64a}.                         (4)

The right event depends only on outer and route randomness. It is a safe
superset, not an assertion that its witness realizes the proposed live set.
The F-word's E-component may occupy fewer source groups; sum over all
source-group counts rather than assigning it the original word's count.

## Exact geometry of the deleted set

For a>0, let r be the number of runs of consecutive indices in L, and
let epsilon=1 if T-1 is in L and 0 otherwise. Then

    |J(L)|=a+r-epsilon <=2a.                                 (5)

Each run [l,u] contributes [l-1,u], except that T-1 is excluded from J.
Consecutive runs are separated by at least one absent live index, so
these contributed intervals are disjoint, although they can abut.

The number of live sets with fixed a,r,epsilon is exactly

    N(a,r,epsilon)
      = binom(a-1,r-1) binom(T-1-a,r-epsilon),                (6)

with out-of-range binomials interpreted as zero. This follows by composing
the a live indices into r positive runs, then distributing the quiet gaps.
Summing (6) over r and epsilon recovers binom(T-1,a).

The 40% input threshold remains available without increasing its distance
parameter when

    e+64a <= L0=838860,       equivalently a<=1638.            (7)

This range is an eligibility condition for the proposed bound, not a
claim that every value of a in that range closes.

## What the quiet minimum distance gives deterministically

Let C_quiet=C_empty. Suppose d_min(C_quiet)>L0. For a coordinate set S
with |S|<=L0, the map X -> (C X_t)_(t in J) is injective on C_J intersect
F^S: its kernel would be a quiet word supported in S. Therefore

    dim(C_J intersect F^S) <=s|J|.                           (8)

There is a related bound retaining the actual state equations. Fix an
inner setup, an allowed output support W, and a live set L. Enlarge the
trajectory class by requiring Q_t=0 outside L and Y outside W=0, while
allowing additional zeros within L and W. This is a linear message space.
Projection to the a states (Q_t)_(t in L) has as its kernel exactly the
quiet words supported in W. If |W|<=L0, its dimension is consequently at
most s a. These are bounds on one fixed linearized support space. They
do not count possible supports or supply fresh-coefficient probability
charges by themselves.

Neither bound proves that C_J has large distance. In general, a high-distance
code C0 can be extended by one independent unit vector. One added constraint
can recover C0 as its quiet subcode, while deleting that constraint restores
a weight-one word. This is a limitation of a quiet-distance-only argument,
not a construction-specific low-distance witness.

## Transfer through the positive quiet kernel

The old proof uses uniform E-valued reference amplitudes after its proved
outer-density comparison. Put q0=|E|=256. For ell allowed positions in
one 16-coordinate emission block, its exact quiet kernel is

    k_ell(z)=sum_(c=0)^ell binom(ell,c) z^(ell-c)(1-z)^c
                          q0^(-c-min(4,ell-c)),

where 0<z<=1. The unconstrained kernel is

    g_ell(z)=(z+(1-z)/q0)^ell.

Term by term in the nonnegative zero-coordinate expansion,

    k_ell(z) >= q0^(-4) g_ell(z).                            (9)

A full 64-coordinate batch has four such blocks. Removing its quiet
equation therefore increases this reference kernel by at most q0^16.
For every fixed routed support and every fixed deleted set J, deleting
all those equations increases the product of kernels by at most
q0^(16|J|). This is a coefficientwise comparison in the nonnegative
zero-coordinate expansion at the same z. It is not a claim of
coefficientwise domination of the ordinary weight enumerators in z.

Apply this comparison to the positive reference expression for (3), after
the source-density comparison and before its fractional power alpha.
The resulting loss is q0^(16 alpha |J|). It is independent of the location
of J, so the remaining product is exactly the old all-quiet product.
Thus the old paired-route averaging and matching reduction remain available.
The already charged free-terminal loss is unchanged: J never includes T-1.

In particular, let U_(g,f)(L0) be a rigorously evaluated old bound for one
source-count/matching histogram at threshold L0, using the fixed markers
z_g,alpha_g of its source-count band. The same positive witness gives

    U_(g,f)(J,w)
      <= U_(g,f)(L0)
           [q0^(16|J|) z_g^(L0-w)]^alpha_g.                 (10)

The source density, projective factor, histogram count, and route
normalization in U_(g,f) are unchanged. Source-count bands may use
different z_g,alpha_g exactly as allowed by the frozen verifier.

Equation (10) must not be described as multiplying the actual all-quiet
failure probability by a deletion loss. The probability for C_quiet may
be zero while C_J has a low-weight word. The comparison is between the
positive reference bounds which precede the probability estimate.

## Bounded next certificate gate

For a specified a>0, a valid complete bound for the event in (4), summed
over every live set of that size, is

    sum_(r,epsilon) N(a,r,epsilon)
      sum_(g,f) U_(g,f)(L0)
        [256^(16(a+r-epsilon))
          z_g^(838860-734003-64a)]^alpha_g.                  (11)

All live-set multiplicities remain outside alpha. A coarser version
replaces a+r-epsilon by 2a and the outer sum by binom(T-1,a).
The source count g in this sum belongs to the auxiliary E-valued word,
not the actual F-valued encoder message. Keep every nonzero auxiliary
source count. In particular, the direct bounds for actual F messages on
one through six groups do not permit removing g<=6 here: a message on
more groups can have an E-component on six or fewer groups, and that
component need not realize the actual live trajectory. Such an exclusion
requires a separate bound for the corresponding C_J classes or a proved
support-preserving comparison. The present reduction supplies neither.

The finite proposal is to replay the existing positive witnesses with
(10), initially for a small declared live-count range and every auxiliary
source count, and use directed arithmetic
for every new claim. A verified old vertex margin can be reused together
with an outward-rounded transfer factor; stored floating margins alone
are not proof evidence. Do not modify the old source or receipt. A new
checker can read those frozen witnesses and retain its own compact bound.

The deterministic part of a bounded regression should enumerate short
binary live sets to check (1), (5), and (6), then check (9) with exact
rational arithmetic for ell=0,...,16. No new checker or numerical run is
part of this note. If (11) fails, the failure is in this deletion-based
envelope, not evidence of a low-distance codeword. A larger gate would
need sharper partly quiet source ranks, actual paired-route geometry,
or the fresh full-field constraints on nonquiet trajectories.
