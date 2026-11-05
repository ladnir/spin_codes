# A direct full-field rank gate for seven groups

The [retained small-group gate](SMALL_STATE_ROUTE_GATE.md) bounds the new
GF32 refresh using the older, weaker coefficient probabilities. This note
instead uses the actual GF32 distribution. It gives a positive operator
for all messages on seven actual outer groups, without replacing their
source space by an independent uniform-input density.

The local rank bounds below are analytic. They do not assert that the
resulting finite operator is below a target failure probability. The new
[floating diagnostic](theory/screen_small_state_direct_rank.py) implements
this interface; a passing directed replay would still be required.

## Fixed construction and source set

Use the construction in [the small-state gate](SMALL_STATE_PROOF_GATE.md).
Let E=GF256 and F=GF(2^32), with q=2^32. There are 8192 outer groups.
Each group contains four RS[64,32] rows and independent invertible
E-valued column mixers. The retained paired route has 256 regions and
alternating adjacent-pair matchings. The inner emits 64 symbols per batch:

    Y_t = X_t + A Q_t,
    Q_(t+1) = D2_t H16 D1_t (Q_t + C X_t),     Q_0 = 0.

Here C=A^T, CA=0, and H16 is superregular. Both diagonals have independent
entries uniform in F minus zero. Setup is sampled once and then used for
every message. There are n=32768 emitted batches; the final state is free.

Fix a set S of seven actual source groups. Initially fix the outer mixers
and route too. Its ambient source space is V=F^896. It is safe to count
all nonzero vectors in V, including those using fewer than seven groups.
The target event is existence of such a vector with output weight at most
e=734003. No step below unions over q^896 independently chosen messages.

For a batch let I be its allowed input coordinates. Let ell_1,...,ell_4
be their occupancies in the four RS[16,4] emission blocks. Define

    c = rank_F C_I = sum_b min(4,ell_b),
    o = #{b : ell_b>0}.

Thus c<=7. These quantities concern allowed coordinates, not the support
of a particular message after setup.

## Separate input equations from fresh diagonal pivots

Fix every D1 matrix and a proposed sequence of exact state supports H_t.
The remaining randomness consists of all D2 entries. For a batch, write
H=H_t, H'=H_(t+1), h=|H|, j=|H'|, and m=16-j. Put R=H16 D1_t and

    M = R[H'^c,H],             N = R[H'^c,:] C_I.

The outgoing zero coordinates require

    M Q_t[H] + N X_t[I] = 0.                                  (1)

Superregularity gives rank M=min(m,h). Eliminate these incoming-state
columns first. The remaining equations are pure input equations, of rank

    d_actual = rank[M,N] - rank M.                             (2)

Collect all such equations on the original message coordinates into one
matrix B. Impose B before exposing any D2 entry. Although B may now be
F-valued, it is independent of every D2. No E-valued scalar-shadow argument
is required.

Use the common-quotient potential

    P_tau(d) = (tau^d-1)/(tau-1),       1<tau<=q-1.

Retain all required nonzero state forms in this quotient. If any becomes
identically zero, the branch has potential zero. Otherwise a surviving
nonzero source witness has potential at least one. An independent pivot
reduces dimension and contracts the potential by at least tau. A dependent
pivot permits at most one value of its fresh nonzero D2 coefficient, so
its survival probability is at most 1/(q-1)<=1/tau.

Consequently, for the fixed D1 matrices and support path, the inherited
common-quotient argument bounds its existence probability by

    P_tau(896-rank B) * tau^(-sum_t min(16-|H_(t+1)|,|H_t|)).   (3)

All support requirements remain intersections inside the probability
space. One must not condition and renormalize on the future support path.
The averaging over D1 comes only after (3), jointly with its input-rank
tags. This ordering prevents double charging the same constraint.

## Actual GF32 rank probabilities

Write P_q(k)=(q^k-1)/(q-1), with P_q(0)=0. Let W be a fixed k-dimensional
subspace and let an a-row MDS check act on diag(D)W, where D has independent
uniform nonzero F entries. For each projective nonzero vector in W, its
image vanishes with probability at most (q-1)^(-a). If its support has
size at most a, that event is impossible. Otherwise choose a pivot
coordinates and expose the remaining coordinates first.

If the resulting map has rank d<k, its kernel contains P_q(k-d) lines.
Expected kernel-line counting therefore gives the rank-class bound

    cap_q(k,a,d) = min(1, P_q(k)(q-1)^(-a)/P_q(k-d)),  d<k,
    cap_q(k,a,k) = 1.                                          (4)

These are separate bounds for disjoint actual rank classes. Their sum
may exceed one; that is a safe positive overcount. They are not multiplied
by another marginal activation probability.

### Activation from a quiet state

For h=0, apply (4) to W=im C_I, k=c, and a=m. Tag the actual rank d of
R[H'^c,:] restricted to W. Activation requires a nonzero kernel vector,
so d<c. Its safe range is

    max(0,c-j,rho(ell,m)) <= d <= min(m,c-1),                    (5)

where

    rho(ell,m) = max sum_b max(0,min(4,ell_b)-4+u_b),
                  0<=u_b<=4, sum_b u_b<=m.

To prove this floor, intersect each block image with a fixed u_b-coordinate
subspace. The total intersection is supported on at most m coordinates,
where the MDS check is injective. Also, a nonzero activation requires
j>=max(1,17-4o). Since X_t is nonzero, its emitted weight is at least one.

Quiet-to-quiet transitions instead impose C_I X_t[I]=0. Their input-rank
tag is c and their probability factor is one. For seven groups these
equations do not imply X_t=0.

### A live incoming state

Suppose h>0 and m>h. The column span of R[H'^c,H] does not depend on
D1_t[H], because those entries only rescale columns. Choose a fixed
left-kernel matrix L of H16[H'^c,H]. Then

    L H16[H'^c,H^c]

is an (m-h)-row MDS check on 16-h coordinates. Any m-h additional columns,
together with H, form an invertible m-by-m minor of H16.

Before exposing D1_t[H^c], choose any fixed
k=max(0,c-h)-dimensional subspace of im C[H^c,I]. Row deletion guarantees
that this choice exists. Apply (4) with

    a=m-h,        n'=16-h,
    max(0,a+k-n') <= d <= min(a,k).                             (6)

The selected rank d is a lower bound on (2), so it can serve as the input
rank tag. Its subspace is fixed before the diagonal entries used by (4).
There is no assumption that an adaptively chosen subspace remains random.

If k=0 or h>=m, tag zero with mass one. For extinction j=0, the projected
check is square and invertible; tag k with mass one. This tag may be
smaller than the actual input rank, which only weakens (3).

## The local positive polynomial

The operator has 17 cardinality states h,j in {0,...,16}. A monomial also
records a nonnegative input-rank tag d. Fix 0<z<1 as an output-weight
marker. Sort the four occupancies decreasingly when computing

    b_h=ceil(h/4),
    w_min(h,ell)=max(0,13 b_h-sum_(b=1)^b_h ell_b).               (7)

Every live Q_t occupies at least b_h of its four coefficient blocks.
Their RS outputs each have weight at least 13; adding X can erase at most
the occupied input positions. This proves (7) without treating A as a
superregular coefficient matrix.

The coefficient associated with rank tag d is as follows:

| Transition | Coefficient before its rank tag |
|---|---|
| h=j=0 | 1, with tag c |
| h=0,j>0 | binom(16,j) z cap_q(c,16-j,d), for (5) |
| h>0 | binom(16,j) z^w_min tau^(-min(16-j,h)) times the rank-class cap from (6) |

For a live transition retain only these necessary support conditions:
extinction requires c>0 and h<=4o; a nonzero outgoing state requires
h+j>=17-4o. Summing binom(16,j) possible outgoing supports is valid because
the local bounds are uniform over their identities. The initial state is
h=0. All-quiet histories are included.

Do not replace the explicit rank classes by the older collapsed estimate
E[tau^(-d)1[d<c]]<=(q-1)^(-m). A sufficient condition previously used for
that shortcut was tau>=q+1, incompatible with the present tau<=q-1.
The larger field permits a larger tau in (3); it does not require one.
Its source potential also grows with tau, so the best marker is not
determined by the local root probability alone.

## Exact paired route and the source cutoff

The route uses the retained striped layout. Within a half-batch of 16
pairs, s singleton pairs and f full pairs have block occupancies

    (f+u, f+s-u),  0<=u<=s,

with multiplicity multinomial(16;16-s-f,s,f) binom(s,u). Convolve the two
halves to obtain four block occupancies. Each pair straddles two different
RS emission blocks. Replacing this law by independent eight-pair blocks
would change the route being analyzed.

For a region, extract singleton/full-pair counts (s,f) after 128 ordered
batches. Divide by

    binom(4096,s) binom(4096-s,f) 2^s.                          (8)

Retain each local input-rank tag during coefficient extraction. For the
fixed source set, s=7-2f. Its two alternating matchings have fixed counts
(f_even,f_odd), each in {0,1,2,3}. They do not get resampled at each column.
An outer column combines four regional operators in even,odd,even,odd
order. There are 64 outer columns.

Within a column, the pure input constraints act on disjoint routed
coordinates, so their ranks add. Its invertible mixers preserve rank.
Any 32 RS columns jointly realize all their 4*7 evaluations. Therefore,
if d_l is the sum of the rank tags in column l,

    rank B >= top32(d_0,...,d_63),       0<=d_l<=28.              (9)

This argument holds for F-valued constraints, not just E-valued ones.
The right side is a lower bound even when some live tags are only
guaranteed ranks.

The strict cutoff top32<896 excludes zero-dimensional source quotients.
Use the identity

    top32(d) = min_v [32v + sum_l (d_l-v)_+].                   (10)

A path below the cutoff has a minimizing integer v in {0,...,27}.
For each threshold put D_v=896-32v and retain
k_v=sum_l(d_l-v)_+<D_v. Its source potential is exactly

    tau^D_v/(tau-1) * tau^(-k_v) * (1-tau^(-(D_v-k_v))).        (11)

Summing the 28 threshold operators is a safe overcount. A separate z
marker can be fixed for each threshold. Include z^(-734003) and
binom(8192,7). Take the maximum over the 16 fixed matching profiles,
not an additional union over profiles: each source set has one profile.
All these operations preserve nonnegativity.

## The discarded terminal state

The operator may use all n transitions and sum every final syndrome
support H_n. This does not impose a final reset. If the implementation
omits D1_(n-1), augment the probability space with an independent copy.
It only partitions the discarded syndrome and cannot change the emitted
word or the bad-output event.

In each fixed H_n branch, its pure input equations may enter B. The final
pivot charge concerns the incoming D2_(n-2), used to form Q_(n-1).
No charge uses D2_(n-1). Summing every H_n is exhaustive, including zero
and nonzero discarded states. Thus the final multiplication by the
all-ones column is a free-terminal bound, not a termination assumption.

## Status and next bounded check

The rank-class formulas and their exposure order have been independently
reviewed. A read-only review of `local_terms` and `rank_cap` in the new
diagnostic found that they implement this interface. No numerical job was
run for this note.

The bounded comparison is now complete. At tau=65280 the inherited-marker
floating bound has log2 score +1473.40; a complete two-pass coordinate
search over all 28 thresholds and 16 profiles improves it to +815.65.
The untuned tau=2^32-1 endpoint is worse. These are failed upper bounds,
not low-weight witnesses. No directed certificate is asserted for this
rank-potential operator.

The separate [exact-route whole-state gate](SMALL_STATE_EXACT_ROUTE_GATE.md)
closes seven through sixteen groups with directed arithmetic. Further
seven-group tuning here is unnecessary. Preserve the source-rank interface
and the following support-sum refinement for later source classes.

## Sum block-support patterns before applying the weight marker

There is a stronger replacement for (7) that preserves the same 17-state
operator. Equation (7) places the incoming support in as few emission
blocks as possible, and places those blocks at the largest input
occupancies. Neither worst case needs to hold for every support in the
union over paths.

For an exact incoming support H, let H_b be its restriction to the four
state coordinates of emission block b. Each nonempty H_b produces a
nonzero RS[16,4] word. Adding X in that block can erase at most ell_b
positions. Therefore

    wt(Y_t) >= sum_(b:H_b nonempty) max(0,13-ell_b).             (12)

Define the positive degree-16 polynomial

    F_ell(x,z) = product_(b=1)^4
          (1 + z^max(0,13-ell_b) ((1+x)^4-1)),
    E_h(ell,z) = [x^h] F_ell(x,z) / binom(16,h).               (13)

Within a block, choosing a support of size j has multiplicity binom(4,j).
Thus the coefficient in (13) is exactly

    sum_(H subset of {1,...,16}, |H|=h)
        z^(sum_(b:H_b nonempty) max(0,13-ell_b)).               (14)

For live incoming states, replace z^w_min in the local table by E_h.
Leave quiet transitions, activation weight z, rank tags, and all other
coefficients unchanged. Formula (13) is symmetric in the four ell_b,
so the existing sorted-occupancy geometry remains sufficient.

This replacement is valid even though the actual state support is not
uniform. To see why, first fix the route, the sequence of support sizes,
and every local rank tag. Apply the preceding rank and probability caps
to each exact support path. Those caps depend on h,j,ell and the tag,
but not on the identities of incoming or outgoing support coordinates.
The top32 source factor also depends only on the tags. The only remaining
dependence on an incoming support H_t is its emitted-weight factor (12).

Now sum the positive bound over all exact supports in that fixed size
path. For t=1,...,n-1, its sum is (14). The final support H_n has no
emission after it, so its sum contributes only binom(16,|H_n|). Equivalently,
each local transition contributes binom(16,j) as before, while its live
incoming emission factor is divided by binom(16,h). These neighboring
factors cancel exactly along the path. The initial support is empty and
has multiplicity one. This proves the replacement by (13), including
the exhaustive free-terminal convention.

One may also retain 35 states (h,b), where b is the number of nonempty
state blocks. Their outgoing multiplicity is

    n(h,b)=binom(4,b)[x^h]((1+x)^4-1)^b.

Because every other local cap ignores b, that operator collapses to a
17-state average. Formula (13) is stronger: it retains which emission
blocks are nonempty during the averaging, not only their number. It also
uses the sum of the four separate nonnegative block floors.

All facts needed for this collapse are present in the current operator:
activation floors depend only on the input occupancy; projected live
rank caps use h rather than the support's block pattern; and extinction
feasibility uses h and the number of occupied input blocks. A later cap
that depends on support identities cannot be combined with (13) by
multiplying separate averages. It would need a joint support sum.

This refinement is an analytic upper bound only. Its first discriminating
check is the existing seven-group diagnostic with (13) substituted for
the live emitted-weight factor, while keeping the other markers fixed.
No new rank-tail estimate or assumption of generic matrix rank is needed.
