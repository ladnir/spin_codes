# Separate the rank weight from the coefficient-set size

The 33-full-column class has shortened message dimension four. Its current
bound uses the same value, 256, for two different purposes: the inverse
probability of an exceptional coefficient and the penalty for one new
independent equation. The joint empty-gap continuation can separate them.
This changes a proof parameter, not the construction or its randomness.

## 1. The global rank weight can be arbitrary

Fix the outer realization, full route, and a nonempty exact source support.
Let V be the corresponding shortened F-linear message space, of dimension
D>=1. Keep the common zero-form space W and all nonzero requirements from
[PROOF.md](PROOF.md). Discard a branch when W contains a required nonzero
form. For any real tau>=1, define

    Phi = 1[branch retained] tau^(-dim W).

If a nonzero message satisfies a branch, dim W<=D-1. Therefore its existence
indicator is at most tau^(D-1) Phi. Summing the covering branches and taking
expectations preserves the bound. No coefficient-distribution assumption
is used in this step, so it does not require tau<=256.

An independent equation multiplies Phi by 1/tau before possible rejection.
This multiplier is an accounting weight. It is not the probability that
the equation holds for a fixed message.

## 2. Exceptional roots retain their actual probability

The construction still samples each fresh coefficient independently from
A=theta+GF(256), a set of size q=256. Consider a required zero equation

    u + a v = 0,

where v is required nonzero and u,v are fixed before exposing a. In the
current quotient by W there are two possibilities.

- If u is outside W+F v, the equation adds one rank for every a. Its relative
  expected potential is at most 1/tau.
- If u=kappa v modulo W, only a=kappa can retain the branch. That value
  occurs with probability 1/q when kappa belongs to A, and never otherwise.
  Its equation adds no rank. Every other value forces v into W and rejects
  the branch.

Earlier exposed coefficients can change this classification. Sequential
exposure in the same W remains valid. Nonzero requirements are retained
indicators; do not condition on them and renormalize the coefficient law.

The [empty-gap Bellman operator](EMPTY_SPAN_RANK.md) already distinguishes
these possibilities. Its optional-coordinate alternatives remain

    F_nonzero(r) + tau^(-1) F_zero(r-1)          [only if r>=2],

    (1-1/q) F_nonzero(r) + (1/q) F_zero(r),

    F_nonzero(r).

Take their maximum. The coefficient q stays equal to 256, independently
of tau. The state-span invariant and the outgoing-rank rejection are
unchanged. Thus this Bellman operator is sound for every tau>=1, not only
tau<=q.

## 3. The coarse active-step matrix uses the smaller parameter

The ordinary fixed-mask matrix merges independent equations and exceptional
roots into one charge. For arbitrary tau>=1, the valid common charge is

    max(tau^(-1),q^(-1)) = psi^(-1),
    psi = min(tau,q).

If a fixed support branch has p fresh-coefficient charges and d known
deterministic rank gains, its expected relative potential is at most

    tau^(-d) psi^(-p) <= psi^(-p).

Accordingly, the ordinary matrix M_I(z) from
[RANK_CONTINUATION.md](RANK_CONTINUATION.md) must be evaluated with charge
parameter psi. In particular, when tau>256, use M_I(z;256), not M_I(z;tau).
Replacing the exceptional-root cost by 1/tau would be invalid.

This gives a composable two-state envelope:

1. Obtain and rigorously verify a positive live vector v satisfying
   T_tau(v)<=lambda v for the empty-gap operator, with root probability 1/q.
2. Compute the singleton constants a,b,d from Section 6 of
   [EMPTY_SPAN_RANK.md](EMPTY_SPAN_RANK.md), using M_I(z;psi) and this v.
3. Use the same linear regional envelopes and free terminal state.
4. Multiply the resulting moment by tau^(D-1) z^(-e).

The source factor uses tau, not psi. Both empty and active steps still act
on the same global potential tau^(-dim W). The coarser active matrix merely
overestimates some independent rank costs. Fixed-route linear domination
still precedes route averaging.

## 4. Why this may help the four-dimensional class

For the 33-full-column class, D=4 for every allowed outer realization.
The source factor is therefore tau^3. Increasing tau suppresses independent
zero branches in the empty-gap operator without weakening exceptional-root
costs. This may reduce the positive vector's ratios across residual ranks,
and therefore reduce the active continuation constant d.

For example, replacing tau=2^8 by tau=2^16 costs 24 bits in the source factor.
It reduces each independent rank multiplier by a factor 256. These are
accounting identities, not a predicted final margin. A net improvement
requires a new verified potential and the complete finite gate calculation.

The refinement does not impose the exact global limit of three surviving
independent rank gains. Tracking that remaining budget explicitly remains
a possible further improvement. Nor can tau be increased without cost:
the source factor grows, and exceptional roots always retain probability
1/256. A failed bound after this change is not a low-distance witness.

## 5. The GF256 source-counting bridge is unchanged

The existing bridge uses

    tau^(D-1) <= (q^D-1)/(q-1),       1<=tau<=q.

That stated range remains essential. The bridge cannot be reused unchanged
when tau>q. For the present fixed-D class, the direct tau^3 source factor
needs no enumerator comparison, so this restriction is irrelevant.

Extending the refinement to all source supports requires their actual
shortened dimensions, a valid dimension-dependent bound, or a new bridge.
It does not follow from closing one fixed-D support class. The original
coefficient distribution, rate, target distance, and retained certificates
are unchanged.

## 6. Verified use of the refinement

With tau=65536, z=499/500, and lambda=993/1000, the exact positive-envelope
check passes. The union over all 33..64 fully populated columns in one
active group has failure below 2^-291, including all 8192 group choices.
See [the finite calculation](SPARSE_GATE.md#7-separate-rank-and-probability-costs-wider-sparse-closures).

The [shared-outer pair bridge](PARTIAL_COLUMN_MOMENTS.md) supplies a valid
source factor at tau=q^2 for the minimum 33-column class, including holes.
With its separately verified marker, this larger class closes below 2^-219.
Neither result uses the unimplemented global-rank refinement or certifies
all message supports.
