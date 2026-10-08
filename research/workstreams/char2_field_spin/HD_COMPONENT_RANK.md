# Component classes: a certified sparse term and a failed common mark

The positive block-rank expansion from `HD_AFFINE_PROOF.md` admits an exact
small transfer after component contraction. A shared Chernoff mark cannot
certify the construction, but marks chosen separately for different component
classes avoid that obstruction.

The main finite result is positive: at K20, the entire permitted-singleton
term of this expansion is below 2^-84.93786249. The result is uniform over
all fixed regionwise row orders. It covers both individual-row affine
placement and four-row-group placement, through the same post-inner control.

The subsequent [pair analysis](HD_COMPONENT_PAIRS.md) bounds the entire
two-permitted-row term below 2^-41.11393101. This is not a minimum-distance
certificate: the terms for three or more permitted rows remain unbounded
in aggregate. A permitted-row term is not the probability of failure from
codewords with that many active outer rows.

Subsequent work bounds the [dense tail](HD_COMPONENT_DENSE.md), but also
proves that rare actual-code kernels defeat the uncapped first moment in
the middle range. Separate class marks do not remove that obstruction.
Use this component transfer within the new [routing-conditioned approach](HD_ROUTING_CONDITION.md),
not as a standalone path to a complete certificate.

## Splitting the positive expansion

Retain R=32768, b=64, rho=32, q=65536, N=bR, and e=838860. Write
Q=q-1 and z=Q^(-rho). The prefix carries globally, as do both post states.
The route consists of independent uniform permutations within the b regions.

Let J be a set of permitted source rows. Require all source differences
outside J to vanish. A row inside J may also vanish; J is not an exact active
support. Contracting the corresponding prefix equalities gives b|J| free
component values. Let C_J map those values through the actual regional route
and post-inner to the output.

For 0<v<=1, define the expected weight enumerator

    W_J(v) = E_(route,post) sum_{t in F^(b|J|)} v^wt(C_J t).

The zero assignment contributes one. The prefix equalities are constraints
on one global path, so free components may cross region boundaries.

The positive row expansion gives the following bound, with an independently
chosen mark v_a for each class size a:

    Pr[distance <= e]
      <= 1/Q sum_{a=1}^R z^a (1-z)^(R-a)
           v_a^(-e) sum_{|J|=a} (W_J(v_a)-1).

Every term is nonnegative. Thus a mark can be chosen after fixing the class
size without changing the bound. An even finer choice per J would also be
valid. The factor 1/Q counts scalar lines rather than their nonzero multiples.

This remains the block envelope: a nonzero row is assigned weight z, including
row supports that the actual MDS code forbids. Exact MDS support feasibility
can strengthen the bound, but the singleton result below does not need it.

## Why one mark for every class fails

The common-mark inequality from `HD_AFFINE_PROOF.md` cannot succeed. This is
an obstruction to that inequality, not to every class-specific envelope.

For any fixed setup, take a source difference supported at the final physical
prefix position. There is one active source row. Its prefix output, routed
input, and post output are confined to the last region, so the final weight
is at most R. The block envelope therefore contributes at least

    z v^(R-e)

to its common-mark Chernoff bound. For v<=1/2, this exceeds 2^805580 at K20.
The actual MDS membership probability of that source difference is zero.

For v>=1/2, count all output vectors with every coordinate nonzero. There are
Q^(N-1) scalar lines. Each has at most R active inverse source rows, so its
block weight is at least z^R. Their common-mark contribution is at least

    Q^(K-1) v^(N-e) >= 2^14470333.

The two inequalities cover every positive mark. They do not survive the
class split: the small permitted-row classes need a mark close to one,
whereas dense classes can use a much smaller mark. The calculation below
confirms that this distinction is operational, not merely formal.

## An exact regional transfer for components

Fix J with |J|=a>0. In region c, list its cut positions as

    0 <= j_1 < ... < j_a < R.

The component preceding the first cut is carried from the preceding region.
Its live flag is sigma. The component after the final cut has flag sigma'.
The a-1 interior components have lengths j_(i+1)-j_i. Their weighted support
polynomial is

    P_c(x) = product_{i=1}^{a-1} (1+Q x^(j_(i+1)-j_i)).

Every live free component has Q possible field values. Components are distinct
variables even when two nonzero assignments happen to have the same value.

Let F_h(v) be the four-state post kernel averaged over all input masks of
weight h in a region. Its normalization is C(R,h), not C(N,h). The post state
is not reset between regions. The transition from sigma to sigma' is

    Q^sigma' sum_l [x^l]P_c(x)
       F_(j_1 sigma + (R-j_a) sigma' + l)(v).

There are eight frontier states: the carried component flag and the two
post-state live flags. Start with all three flags zero, multiply the regional
transfers, and sum all final states. The result is W_J(v).

This transfer is exact for the conditioned four-product post ensemble.
For fixed component values, its local support law depends only on whether
the input and state components are zero. Fresh post coefficients give that
law even when an input value repeats at many positions. Conditional on a
component assignment, the independent regional permutations produce uniform
support masks of their respective occupancies. These facts justify summing
component values through their live flags here; no new independent values
are assigned to their occurrences.

For a singleton, P_c=1 and the transition simplifies to

    Q^sigma' F_(j_c sigma + (R-j_c) sigma')(v).

For two permitted rows, the sole interior factor is 1+Q x^(j_2-j_1).
Thus a specified pair remains inexpensive once the kernels F_h are available.

## A shared two-state post control

For scalable bounds, put p=1/Q and beta=v+(1-v)p. The positive matrices

    T0 = [[1, 0], [0, beta]],
    T1 = [[0, v], [v p^2, beta]]

index the zero state vector and a nonzero state vector. They dominate the
conditioned four-product law after aggregating its three nonzero states.

For an old state with exactly one live component, the output-zero probability
on live input is p, and return to the zero vector is impossible. For two live
components, return probability is p^2 and that return has nonzero output.
Its live-to-live marked mass is beta-p^2. These cases establish the stated
entrywise control. The small tests verify it directly.

The same matrices dominate the distinct three-product ensemble in
`theory/hd_three_product.py`, which sets c=d1 and needs no conditioning.
Consequently every bound obtained from this common control applies to both
post ensembles. This does not transfer the earlier fresh-coefficient rank
pivot proof to the three-product map; its output and updated first state
are correlated.

Replacing the exact four-state kernels by this control reduces the component
frontier to four states. The replacement is an upper bound, not an equality.

## Computing a controlled kernel without a degree-R polynomial

Fix h live input positions. Their h+1 zero gaps are a uniform weak composition
of R-h. If G is the sum of L selected gaps, its positive factorial expansion is

    E[beta^(-G)]
      = sum_{ell=0}^{R-h} C(R-h,ell)
          (L)_(ell)/(h+1)_(ell) (1/beta-1)^ell,

where the parenthesized quantities are rising factorials. This also defines
the cases L=0 and L=h+1, which give 1 and beta^(-(R-h)), respectively.

Let M_L denote this moment. The two-state controlled kernel has the run sums

    F01 = sum_{t=0}^{floor((h-1)/2)}
            C(h-t-1,t) v^(2t+1) p^(2t) beta^(R-2t-1) M_(t+1),
    F10 = p^2 F01,
    F11 = sum_{t=0}^{floor(h/2)}
            C(h-t,t) v^(2t) p^(2t) beta^(R-2t) M_t,
    F00 = sum_{t=1}^{floor(h/2)}
            C(h-t-1,t-1) v^(2t) p^(2t) beta^(R-2t) M_(t+1).

These formulas count live-state excursions and the zero gaps between them.
For h=0, use diag(1,beta^R). Exact rational tests compare every coefficient
with direct polynomial multiplication for all h at region sizes through eight.

The directed K20 evaluation keeps t=0 through 3 in F01 and F11, and t=1
through 4 in F00. Put lambda=h p^2. A valid omitted-tail bound for the first
two entries is

    lambda^4 / (4! (1-lambda/5)).

For F00, multiply this bound by p^2. For F10, multiply the F01 bound by p^2.
The combinatorial coefficient is at most h^t/t!, while the remaining marked
path factor is at most one. This proves the tail bound when lambda/5<1.
At K20 that condition holds with ample margin.

The gap-moment series is also positive. After term ell, every subsequent
term ratio is at most (R-h)(1/beta-1)/(ell+1). Once this number is below one,
a geometric tail bound supplies a directed upper. No floating saddle-point
approximation or global-interleaver normalization is used.

## A uniform K20 bound for the entire permitted-singleton term

It is unnecessary to inspect all 32768 row profiles. The controlled kernels
have bounds uniform over every positive occupancy h. Define

    eta = (1+p^2)^R - 1,
    B = v(1-beta^R)/(R(1-beta)),
    D = p^2(1+p^2)^(R-1),
    G = [[D, B+eta], [p^2(B+eta), beta^R+eta]].

Then F_h is entrywise at most G for 1<=h<=R. To see this, separate paths with
no return to zero from paths with a return. Their total additional mass is
at most eta. Without returns, the birth term is largest at h=1, where the
first live position is uniform and gives B. The live-to-live no-return term
is beta^R. A path ending at zero must return on its last live input; its
mass is at most D. Finally, the run sum gives F10=p^2 F01.

Let E=F_0, let L be the controlled F_R, and let G0 be G with its upper-left
entry replaced by one. The following block matrix dominates every singleton
cut position j, with block indices sigma,sigma':

    K = [[E, Q G], [G0, Q L]].

The reason is precise: the two diagonal occupancies are always zero and R;
the 0-to-1 occupancy R-j is positive; and the 1-to-0 occupancy j may be zero.
Keeping the full-occupancy block L is important for the finite bound.

Let H=e_(component0,post0)^T K^64 1 - 1. Summing the R possible permitted
singletons and dropping the factor (1-z)^(R-1)<=1 gives the upper bound

    R * Q^(-33) * v^(-838860) * H.

A 128-bit directed-Arb computation at the rational mark v=2949/2950 proves

    log2(singleton-term upper) < -84.93786249.

The per-profile uniform upper is below 2^-99.93786249. These statements hold
for any fixed row permutation within every coordinate region, including the
two currently considered placements. They concern the a=1 permitted/free-row
term of the positive expansion, not an event consisting of actual one-active-
outer-row codewords. The common two-state control makes them valid for both
post ensembles described above.

For comparison, the sharper run kernels were evaluated at six frozen affine
profiles: rows 0,1,2,8192,16384,32767. Their individual upper logarithms range
from -112.494604 to -112.491151. These spot checks are not the justification
for multiplying by R; the uniform block matrix is.

## A bounded two-row feasibility check and the remaining work

At v=1474/1475, the specified affine pairs (0,1), (0,16384), and (1,32767)
have individual permitted-pair contributions below 2^-214.79295. A fourth
profile uses cuts (0,1) in every region, as in an aligned within-group pair.
Its individual contribution is below 2^-201.74253.

These are only four profiles. They do not establish a maximum over pairs,
and no factor C(R,2) is applied to them. They show that the next sparse term
is numerically viable at this resolution, including a clustered pair profile.

The [subsequent all-pairs bound](HD_COMPONENT_PAIRS.md) completes this step
using uniform occupancy bounds, rather than extrapolating these four profiles.
It preserves the exact all-live occupancy R block and the interior factor
1+Q x^length. Sharper first-return sums and minimum-two occupancy bounds
control the losses after multiplication by Q^2.

A complete 40-bit result still requires the sum of every a>=3 term, or a
stronger envelope that bounds the omitted classes. The next step should retain
component lengths and their joint constraints when bounding these terms.
Neither the common-mark obstruction nor the first two classes settles them.

Ten tests pass. They cover both local post laws, the run sums and their
directed uppers, exact component flags and value counts, the positive
permitted-row expansion, uniform occupancy bounds, the common-mark
obstruction, the K20 singleton receipt, and the scope of the pair checks.

```text
python -B research/workstreams/char2_field_spin/theory/hd_component_rank_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_component_rank.py
```

Only this new note and the two new `hd_component_rank` helpers were written.
No native build, remote call, or performance benchmark was run.
