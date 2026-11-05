# Count live quartet selectors before state amplitudes

This note gives a quotient enumeration for the unchanged quartet
construction. It groups actual state histories by their live quartet
sets and forced births, then counts amplitudes through allowed coordinate
runs. The construction, free terminal state, and 45%/40-bit target are
unchanged. This is an audited counting interface, not a full distance
certificate; the remaining global geometry sum has not been bounded.

Use the once-sampled setup and projective field argument in
[the rank reset](QUARTET_GLOBAL_RANK_RESET.md). Write q=2^32, n=512,
T=4096, and N=nT. The recurrence and output are

    Q_0=0,   Q_t=A_t Q_(t-1)+X_t,
    Y_t=Q_t+delta Q_(t-1),

where delta is nonzero. Every live quartet has an independent raw uniform
GF(q^4) multiplier. The accepted setup has density loss at most
gamma<=1024/989. Fix the routing and all quartet partitions before the
field count below.

## Selectors, births, and allowed states

For an actual state history, let U_t be the union of the quartets at
update t meeting supp(Q_(t-1)). Let B_t=supp(Q_t)\U_t. Set

    U_1=empty,          U_(T+1)={1,...,n},
    B_t subset U_(t+1)\U_t,
    A_t=(U_t intersect U_(t+1)) union B_t,
    V_t=U_t union B_t.

The formal final set U_(T+1) imposes no terminal restriction and has no
random quartet partition or additional update. Every actual history has
supp(Q_t) subset A_t. Conversely, arbitrary states supported in A need
not realize the selected live sets U; this distinction is essential.

Fix a candidate pair (U,B), with each U_t a union of current quartets.
Write f=sum_t|B_t|. For a nonzero actual history f>0: its first nonzero
state occurs in a quartet whose incoming state is zero. We require every
declared B position to have a nonzero state value. These positions are
the forced births outside U, not all coordinate activations.

The inclusion A_(t-1) subset U_t follows from the definition of A and
the restriction B_(t-1) subset U_t. Thus every unselected quartet has
zero incoming state for every allowed amplitude assignment. Outside U,
the input equals the state on B and is zero elsewhere.

## Apply the field rank before enlarging the amplitudes

For outer row i, let ell_i and b_i be its numbers of positions in U and
B. Put

    R_U=sum_i min(128,ell_i),
    R_B=sum_i (min(128,ell_i+b_i)-min(128,ell_i)),
    R=R_U+R_B.

Fix a candidate output Y and recover Q from
Q_t=Y_t+delta Q_(t-1), independently of the random quartet maps. At this
step retain all the exact selector conditions: every selected quartet
has a nonzero incoming state, and Q_B is nonzero in every declared B
position. On each selected quartet, multiplication of that fixed nonzero
incoming vector by a raw uniform GF(q^4) element is uniform in GF(q)^4.
The input coordinates on U are therefore jointly uniform.

Let b denote the vector Q_B. Eliminating the input coordinates on U from
the outer parity checks leaves R_B independent homogeneous equations on
b. If b satisfies these residual equations, the raw-field probability
of the remaining checks is q^(-R_U); otherwise it is zero. This is the
same RS rank decomposition as in the rank reset. In particular, there
are at most q^(f-R_B) compatible vectors b with every entry nonzero.
Some footprints have none, which only strengthens the bound.

Only after obtaining this field probability may the amplitude count be
enlarged to all Q supported in A with the declared nonzero values on B.
The enlarged family can have zero incoming state in a selected quartet.
The factor q^(-R_U) is not asserted for those additional assignments;
they merely upper-bound the number of assignments in the exact family
to which that factor was already applied.

## Exact output weights of allowed coordinate runs

For each coordinate, decompose its allowed times A into maximal runs of
consecutive updates. A run of length ell ending before T is closed;
one ending at T is open. A run is forced if its first position belongs
to B and unforced otherwise.

Every B position starts a run. Indeed, B_t is disjoint from U_t, whereas
A_(t-1) is contained in U_t. Consequently a run contains at most one B
position, at its start. Fixing the nonzero value at a forced start does
not change that run's weighted count: multiplying all states in the run
by a nonzero scalar transports any such value to any other and preserves
which outputs are zero. Thus the count is the same for every compatible
nonzero vector b, even though the residual RS equations can couple its
coordinates.

For 0<x<=1 define

    a=1+(q-1)x,    c=1-x.

The weighted enumerators, where each nonzero output contributes x, are

| Allowed run | Weighted amplitude count |
| --- | --- |
| Unforced, closed, length ell | (a^(ell+1)+(q-1)c^(ell+1))/q |
| Forced start fixed nonzero, closed, length ell | x(a^ell-c^ell)/q |
| Unforced, open, length ell | a^ell |
| Forced start fixed nonzero, open, length ell | x a^(ell-1) |

To derive the closed formulas, rescale states and outputs by powers of
delta. The ell allowed states then give ell+1 output symbols with one
zero-sum equation. The q-ary single-parity-check weight enumerator is
(a^(ell+1)+(q-1)c^(ell+1))/q. Fixing its first output to a nonzero value
contributes x; the remaining ell outputs must sum to a prescribed
nonzero value, whose enumerator is (a^ell-c^ell)/q. For an open run the
map from its states to its ell output symbols is bijective, so there is
no final equation. A forced first output leaves ell-1 free symbols.

A closed run includes its closing output at the following forbidden
state position. Distinct maximal runs have disjoint output intervals,
so their enumerators multiply. There is no closing output after update
T. Let E_(U,B)(x) be this product, with all forced starts fixed to any
chosen nonzero values. It is independent of those choices.

As a check, E_(U,B)(1)=q^(sum_t|A_t|-f): each forced start fixes one
otherwise free state value, for open and closed runs alike.

## Projective normalization and the resulting bound

The exact candidate family and residual equations are homogeneous.
Because f>0, every candidate is nonzero and its scalar orbit has q-1
elements. All elements of an orbit have the same output support and
the same field-membership event. Count one representative per orbit.

For fixed geometry and (U,B), the probability under the accepted field
setup that a nonzero output of weight at most W realizes this pair is
bounded by

    min(1,
        gamma/(q-1) * q^(f-R) * x^(-W) * E_(U,B)(x)).   (1)

Indeed, x^(-W) times the weighted amplitude sum bounds the number of
outputs of weight at most W. There are at most q^(f-R_B) compatible
nonzero birth vectors, each with enlarged weighted count E_(U,B)(x).
Multiply by q^(-R_U), divide by q-1, and pay the accepted-setup density
loss gamma once. The normalization here is gamma/(q-1), not gamma/q.
At x=1, (1) reduces to gamma*q^(sum|A_t|-R)/(q-1), before the cap at one.

Every actual history has one pair (U,B), so a union bound over these
pairs covers all histories without separately enumerating their state
support masks. The live-set choices are at most 2^(128(T-1)), since the
first layer is empty. At fixed f, a further crude birth-mask count is
binom(N,f), although the inclusions B_t subset U_(t+1)\U_t permit a
sharper count. These are interfaces, not a claim that their crude product
with (1) closes the target.

The remaining task is to control the global sum of (1) over selectors
and births under the random partitions and routing. The allowed-run
geometry retains temporal intersections of consecutive live sets; it
must be kept when estimating that sum. No additional field-rank factor
may be attached to the enlarged amplitude family without a separate
proof. Actual state volume Ls and death count D are not fixed by (U,B),
and need not equal the corresponding allowed-set statistics. Prior
classes involving Ls or D must be imposed before enlargement, or their
indicators must be retained in a refined amplitude count; their gains
cannot be inferred from the allowed-run lengths alone. The full-code
certificate remains open.

## Retain actual state volume before a fractional power

The same quotient can retain the actual state volume Ls without a
separate union over individual support masks. Define

    A0=sum_(t=1)^(T-1) |U_t intersect U_(t+1)| + |U_T|,
    m=Ls-f,           G=sum_t |U_t|=F-f.

Here A0 counts allowed positions outside B, so sum_t|A_t|=A0+f.
For each fixed compatible nonzero birth vector, the number of allowed
amplitude assignments with exactly Ls nonzero states is

    binom(A0,m) (q-1)^m,

with value zero unless 0<=m<=A0. This counts every choice of the m
nonzero positions and their values once. It ignores the output-weight
restriction, so it also bounds the desired low-weight class. Retain
the exact-live filter through the field-rank step, then enlarge this
amplitude count as before. The resulting bound is

    P_field(U,B;Ls)
      <= min(1, gamma/(q-1) q^(Ls-R) binom(A0,m)).       (2)

The discarded factor is (1-1/q)^m<=1. Formula (2) restores actual Ls;
it does not identify Ls with the allowed volume A0+f.

Choose 0<=theta<=1, y>0 and b>=1. Apply theta to the completed
fixed-(U,B) bound, after summing its amplitudes. On the feasible range,
the binomial coefficient and the capacity indicator give

    1_(A0>=m) binom(A0,m)^theta
      <= b^(-m) y^(-theta*m) beta^A0,
    beta=b(1+y)^theta >=1.                             (3)

For theta=0 the left side means the capacity indicator. Thus one may
use beta=b and the correction beta^(-m). For theta>0, the specialization
b=1 uses y=beta^(1/theta)-1, requiring beta>1. Applying the fractional
power separately to individual actual support masks would lose the
grouping achieved by (2).

## Separate the fixed-profile route factor

Fix an actual profile (h,R,Ls,f,F), with G=F-f, and put

    C_H=F-R+128h,      C_out=R-128h.

For any fixed U,B after the partitions are exposed, the necessary
route event has the uniform bound

    Q_route=min(1,
      binom(F-f,C_H-f) binom(N-F,256h-C_H)
                         / binom(8192,h)^255).          (4)

Every actual forced birth is in an eligible row. For one candidate H,
count the eligible-slot assignments that contain B and have exactly
C_H covered positions. Dropping the individual regional quotas leaves
the two numerator binomials. The true uniform regional law has
denominator binom(8192,h)^256; union over candidate H contributes one
binom(8192,h). Thus (4) already includes the H-set union.

The field bound in (2) is uniform over routes in this profile, and A0
does not depend on the route. Apply (4) while the fixed G,R,F restriction
is still present. Only afterward enlarge the selector sum by positive
markers. There is no second H union or separate route penalty.

For fixed G, births lie in the complement of U, which has N-G positions.
Their count is at most binom(N-G,f). The additional condition
B_t subset U_(t+1) can be dropped for this upper count. Neither this
birth bound nor (4) asserts that a relaxed selector pair is realizable.

## The exact 129-state selector overlap moment

Let i,j denote numbers of selected quartets in consecutive layers, and
let z>0 mark the volume G. For independent uniform coordinate subsets
of sizes 4i and 4j, define

    H_(i,j)(beta)=E beta^|U intersect V|
      =sum_l binom(4i,l) binom(4j,l)
                        (beta-1)^l / binom(512,l),
    K_(i,j)=H_(i,j)(beta) binom(128,j) z^(4j),
    t_i=beta^(4i).

Summing over chosen quartet subsets and averaging their independent
partitions gives exactly

    M_sel(z,beta)=e_0 K^(T-1) t
                =sum_selectors E_partitions z^G beta^A0. (5)

There are T-1 random selector layers because U_1 is fixed empty. The
terminal t contributes |U_T| to A0. It contributes no new partition,
quartet choice, or volume marker. For beta>=1, replacing t by beta^512
times the all-ones vector is a valid looser endpoint bound.

Combining (2)--(5) yields the fixed-profile probability upper bound

    Q_route binom(N-G,f)
      * (gamma/(q-1))^theta q^(theta*(Ls-R))
      * b^(-m) y^(-theta*m) z^(-G) M_sel(z,beta).       (6)

This is a union over selectors and births after grouping amplitudes.
It changes no part of the construction. It retains actual Ls but
continues to ignore output weight and actual death count in this
amplitude bound.

## An exact rational specialization

Take theta=1, b=1, y=3, beta=4, and volume marker z=1/4.
Use the positive potential v_i=16^i. Given a fixed old selected set,
a new quartet containing c old positions contributes

    1+z^4*16*beta^c = 1+4^(c-2).

For c=0,...,4 the exact inequality

    1+4^(c-2) <= (17/16) 2^c                           (7)

holds, with equality at the two endpoints. Multiplying over all 128
new quartets, whose old occupancies sum to 4i, proves

    K v <= (17/16)^128 v.

This holds for every partition before averaging. The exact endpoint
satisfies t_i/v_i=16^i<=16^128=2^512, and v_0=1. Therefore

    M_sel(1/4,4) <= 2^512 (17/16)^(128(T-1)).           (8)

No floating matrix or numerical witness is needed. Since
gamma/(q-1)<2^-30, (6) gives the explicit profile exponent

    log2 binom(N-G,f) + log2 Q_route -30 +32(Ls-R)
      +2G -(Ls-f)log2(3)
      +128(T-1)log2(17/16) +512.                       (9)

The coarse constant -30 is safe for the accepted density and the
q-1 projective denominator. Directed bounds for the binomial terms and
the two rational logarithms turn (9) into an integer certificate.

[quartet_selector_rank_certificate.py](quartet_selector_rank_certificate.py)
implements this bound without a production matrix or floating witness.
It also checks a geometry-only control with theta=0, beta=2, volume
z=1/2, and potential v_i=4^i. The local comparison is

    (1+2^(c-2))^2 <= (25/16)2^c,

giving a moment bound 2^256(5/4)^(128(T-1)) and the capacity correction
2^(-(Ls-f)). The serial exact replay passed 10 local comparisons,
10 small overlap rows, 16 free-endpoint comparisons, and 2142 directed
entropy comparisons. Independent source audit passed both branches.

For the actual profile

    (h,R,Ls,f,F)=(1536,753664,730726,4,851968),

the directed integer exponents are

| Bound | Integer log2 probability upper bound |
| --- | ---: |
| Grouped actual-L rank bound, theta=1 | -305581 |
| Geometry-only capacity control, theta=0 | +126588 |

The first value is an exact fixed-profile upper bound. The positive
control is a limitation of that control, not a codeword witness.

Equations (6) and (9) are fixed-profile bounds. They still require a
union over the actual profile coordinates. The previously established
strict 2^96 bound on the number of (h,R,Ls,f,F) profiles applies; an
integer exponent at most -172 therefore gives a qualifying profile
class with total failure below 2^-76. A finite tested grid does not
establish that every remaining profile qualifies. The full 45%/40-bit
certificate remains open.

## Closed allowed runs and actual deaths

There is an exact accounting relation, but selector departures alone do
not lower-bound the number of closed allowed runs. Define

    C_U=sum_t |U_t\U_(t+1)|=G-A0,
    H_A=A0+f-Ls,
    C_A=number of closed coordinate runs in A.

Here H_A counts allowed positions at which the actual state is zero.
The state-zero count at live coordinates satisfies

    E=F-Ls=C_U+H_A.                                    (10)

Every actual coordinate death occurs either at an allowed zero position,
or at the closing edge of an allowed run. Thus

    D<=H_A+C_A,
    C_A>=max(0,D-H_A).                                 (11)

In particular, actual D cannot be replaced by C_A without accounting
for allowed holes.

For a sharper description, let I count closed coordinate runs of U of
length one that have no forced birth immediately before them. The
formal final U_(T+1) is full, so C_U counts closed U runs exactly. Each
closed U run gives one closed A run unless it is one of these I
singletons. Therefore

    C_A=C_U-I.                                        (12)

The singleton term can grow with T. As a selector/support example,
start four core coordinates with four forced births. At every later
update select four quartets, each containing one core coordinate and
three fringe coordinates. Alternate two disjoint twelve-coordinate
fringe sets and keep the actual state on the core. Every selected
quartet is live, all allowed runs are open, and C_A=0, whereas
C_U=12(T-2). This rules out a general bound C_A>=C_U minus only a bounded
endpoint allowance. The example does not impose RS compatibility and
is not a codeword witness.

Finally, wt(Y)=Ls+D-Z_survivor, where Z_survivor counts zero outputs
between two nonzero states. For D<=W-Ls, the output-weight condition
therefore places no additional restriction on the amplitudes. An
output-aware quotient can still help where larger deaths dominate its
count, but that gain cannot be assumed for low-D histories.

## Exact verification record

The parent-run standalone
[quotient checks](quartet_live_selector_quotient_checks.py) passed 140
run-enumerator comparisons, 240 run-envelope inequalities, 138 fixed-birth
fibres, 63 coupled residual linear systems, and 84,752 support-containment
cases. The amplitude checks use GF(2) and GF(4), including nontrivial
delta, and verify projective normalization and the free terminal.
A zero-incoming negative control confirms why the field rank must be
applied before enlarging the amplitude family. These checks validate the
counting identities; they do not bound the remaining global geometry sum.
