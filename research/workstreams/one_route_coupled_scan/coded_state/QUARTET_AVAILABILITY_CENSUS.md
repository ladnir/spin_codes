# Recover the footprint from availability

The fixed quartet construction and full 45%/40-bit target are unchanged.
This checkpoint strengthens the global rank argument; local refinements
remain parked. The whole low-rank remainder is not yet certified.
Notation and the all-run field bound are from
[the closure checkpoint](QUARTET_ALL_RUN_CLOSURES.md).

## A canonical recovery identity

Fix the route and quartet partitions. Let S_t be an actual state support,
and write Gamma_t(U) for the union of update-t quartets meeting U.
Its covered footprint is V_t=Gamma_t(S_(t-1)) union S_t. Availability is
A_t=V_t intersect the union of next-step quartets wholly covered by
V_(t+1), with A_0 empty and A_T=V_T. Every actual state satisfies S_t<=A_t.

For every realizable footprint,

    V_t = Gamma_t(A_(t-1)) union A_t.                 (1)

The inclusion from left to right follows from S<=A. Conversely, A_t<=V_t,
and every quartet meeting A_(t-1) is wholly contained in V_t by the
definition of availability. Thus an availability history determines at
most one realizable footprint for fixed partitions. This is not the
claim V=A; nor does it apply to arbitrary relaxed footprints.

Two independent audits agree. Exhaustive small-mask checks verify (1)
for 200704 support histories across the width-4 and width-8 toy geometries.

## Count histories by their runs

Write L=sum|A_t|, D=sum_(t<T)|A_t\A_(t+1)|, and E=|A_T|. A coordinate
starts unavailable and has a free terminal state. Its exact generating
function is

    G_T(s,z,u) = (1,0) [[1,s],[z,s]]^T (1,u)^tr.

The coefficient of s^L z^D u^E in G_T^512 counts all availability masks
with those statistics. For one nonempty coordinate history, the exact
coefficient is binom(L-1,D+E-1) binom(T-L,D), E in {0,1}.

Append a fixed zero after each coordinate history and concatenate. This
injects the histories into length-(N+512) strings ending in zero, with L
ones and D+E runs. Therefore, for L>0,

    count(L,D) <= B(L,D)
      := min{binom(N,L),
          sum_(E=0)^512 binom(L-1,D+E-1)binom(N-L+512,D+E)}.   (2)

Invalid binomial indices contribute zero. There is no binom(512,E)
factor: prescribed sentinel positions identify the terminal coordinates.
The count summed over L also satisfies

    count(D) <= min{2^N, binom(N+512,2D+512)}.          (3)

One way to prove (3) is to use the one-coordinate coefficient
binom(T+1,2D+1), then apply Vandermonde across coordinates. These counts
replace a footprint union by (1); they are not extra multiplicative
factors on top of the old 2^N union.

For nonempty profiles B(L,D)>0, put h(L,D)=ceil(log2 B(L,D)). The class

    32R-16L+15D >=16W+h(L,D)+135                     (4)

has total failure below 2^-76. For each fixed profile, the field bound
and census give less than 2^-119. There are fewer than 2^43 profiles.
This conservative allowance avoids needing a sharper profile-count
constant. Condition on arbitrary geometry first, then average.

## Linear exclusions from a weighted census

A weighted version gives linear remaining-event constraints. For integer
j>=0, the matrix [[1,1],[2^(-2j),1]] has positive eigenvector
(1,2^-j)^tr with eigenvalue 1+2^-j. The free terminal vector (1,1)^tr is
bounded entrywise by 2^j times this eigenvector. Consequently,

    sum_A 2^(-2j D) <= 2^(512j)(1+2^-j)^N <= 2^c_j. (5)

Use c_0=N. For j>=1 a valid integer choice is
c_j=512j+ceil(3N/2^(j+1)), since ln(2)>2/3 and ln(1+x)<x.
The first inequality follows by integrating 1/t>4/(1+t)^2 on (1,2).
At j=3 the sharper exact comparison 9^64<2^203 gives c_3=361984.

For each fixed j, BOTH following classes have failure below 2^-76:

    32(R-L)-2jD >= c_j+45,                          (6)
    32R-16L+(15-2j)D >=16W+c_j+92.                 (7)

For (6), the x=1 field bound is less than 2^-31 q^(L-R), because
gamma(1+2/q)^N <=1048576/1011747<2. On the class, its summand is bounded
by 2^(-c_j-76)2^(-2jD); apply (5). For (7), use the previously audited
bound 2^[16+16W-32R+16L-15D] and the same weighted census.

These are whole-history class theorems. They do not assume a particular
eligible-row set or add a union beyond the availability census. The
free terminal is retained in (5). They concern J-field histories; their
individual class labels need not survive extension-field basis projection.

Examples of exact constants are:

| j | c_j | rank-class threshold c_j+45 | closure-class threshold 16W+c_j+92 |
|---:|---:|---:|---:|
| 3 | 361984 | 362029 | 15461564 |
| 6 | 52224 | 52269 | 15151804 |
| 8 | 16384 | 16429 | 15115964 |
| 12 | 6912 | 6957 | 15106492 |

Combining sixteen j values in both families costs at most
32*2^-76=2^-71. A choice of j depending on the footprint cannot silently
reuse one 2^-76 allowance. Preserve the earlier certificates as well.

The helper `quartet_availability_census.py` passes 211 exact profile-count
checks, 20 weighted-census comparisons, the recovery checks above, and
exact rational constant checks. No floating timing or diagnostic is used
to justify (6) or (7).

## Remaining-event substitutions

The integer complement of (7) has cap C=16W+c_j+91. A log2 marker tau>=0
multiplies its canonical count by 2^[tau(C-32R+16L-(15-2j)D)]. It scales
the rank weight a by 2^(32tau), m and b by 2^(16tau), and closure weight
alpha by 2^[-(15-2j)tau], with global factor 2^[tau(C-4096h)].

For (6), the cap is c_j+44. Its marker scales a,m,b by 2^(32tau) and
alpha by 2^(2j tau), with global factor 2^[tau(c_j+44-4096h)]. Apply all
markers before positive relaxation. Survival-tag methods require the
resulting alpha<=1; the implementation rejects a violation.

Selected five-tilt tagged probes still fail: the closure j=3 marker at
tau=1/64 gives +206414, and rank j=8 markers at tau=1/128,1/64 give
+153687,+207149, versus +134694 without those markers. These FLOAT
values do not refute the class theorems. They show that the tested
positive relaxation does not yet exploit the smaller remainder well.

## Direct availability transfer

Identity (1) also supplies a different global counting method. Fix a
candidate availability history BEFORE averaging the independent quartet
partitions. Reconstruct V by (1); apply the grouped field bound on the
event that it is a genuine fixed point of the availability definition,
then drop that indicator positively. Do not condition a fresh partition
on the random realized availability mask from the previous step.

The canonical forced births are exactly A_t\Gamma_t(A_(t-1)). For incoming
availability weight u and outgoing weight v, put a_birth=mp chi and r0=r chi.
If k quartets meet the old mask, their distribution is

    P_(u,k) = binom(128,k)[z^u]((1+z)^4-1)^k/binom(512,u).

The positive 513-state kernel is

    K_(u,v)=sum_k P_(u,k) r0^(4k)
      [y^v](alpha+b y)^u(1+b y)^(4k-u)(1+a_birth y)^(512-4k).  (8)

It retains the closure weight alpha exactly; no auxiliary tag rounding
is needed. Old-hit quartets are wholly covered; new positions outside
them are forced births and receive the eligible-row factor p. The free
terminal must satisfy Gamma_T(A_(T-1))<=A_T=V_T, giving

    t_u=sum_k P_(u,k)(r0 b)^(4k)(1+a_birth)^(512-4k).  (9)

The complete count is e_0 K^4095 t. The actual-coordinate-count basis is
common across quota tilts; regional conditioning uses the same eligible
label coefficient argument as before. Exact empty transitions and all
nonempty endpoints must be retained. Unlike tags, alpha in this direct
kernel is an ordinary positive closure marker and may exceed one.

`quartet_direct_availability.py` implements (8) using positive Taylor and
homogenization convolutions. Exact small-mask checks pass 530 matrix
entries, 14 terminals, 14 partition laws, two zero operators, eight free
terminal compositions, and the regional recurrence. The initial h=2048
nine-tilt FLOAT result is +748551, substantially worse than the tagged
control. Its relaxed midpoint has about 396 available positions per
step, with an almost full covered footprint. Such masks can violate the
availability fixed-point condition that was dropped. The census itself
is still valid; the direct transfer is a looser numerical control.

## Retain global consistency before adding more states

Canonically V_(t-1) intersect Gamma_t(A_(t-1))=A_(t-1). Any covered old
position in a quartet hit by old availability must itself be available,
because that quartet is wholly covered next. Therefore the unavailable
old covered set E_(t-1)=V_(t-1)\A_(t-1) is disjoint from Gamma_t(A_(t-1)).
For each t>=2, their total size is at most 512. Summing, and using
E_T=empty and sum|Gamma_t(A_(t-1))|=F-f, proves

    2F-L-f <= N-512.                                (10)

The canonical forced-birth count here is the same f used by the original
footprint count. A log2 marker tau>=0 applies (10) by scaling coverage
chi by 2^(-2tau), m by 2^(2tau), b by 2^tau and the global weight by
2^[(N-512)tau]. Apply it canonically before enlarging either transfer.
The positive direct transform remains valid: mp chi stays unchanged
while b increases. The weaker condition 5F-L<=4N follows by charging
every unavailable covered position to one of at most four positions
in a next-step quartet containing a coverage hole.

Both identities pass independent audit and the exact recovery checks.
For the direct count, tau=1/4,1/2,1 gives +641425,+959653,+1994404 at
the same selected point. Thus the scalar marker improves one result but
does not make the direct count competitive. An ordinary coverage-cap
marker of -1 also fails (+1636806). These are failed FLOAT bounds, not
counterexamples to the code or the run-census theorems.

The next useful refinement must retain the unavailable covered positions
jointly with old availability, so their exclusion from old-hit quartets
cannot be relaxed independently at every update. It should preserve
the global rank grouping, exact quotas and free endpoint. The tagged
bound remains the better tested numerical control. No full low-rank
certificate is claimed; all prior checkpoints remain available.

The same scalar exclusion at tau=1/8 in the broad six-level tagged count
also fails (+265482 versus +105850 without it). Final checks of the census,
direct transfer and tagged transfer all pass. Every numerical process
has exited; no remote experiment or local artifact remains running.

## A joint availability/covered-hole interface

There is an audited way to retain the missing consistency exactly before
any computational compression. Track u=|A_prev| and e=|E_prev|, with the
two sets disjoint. Use x,z to mark old A,E and y,w to mark new A,E. Put
a_birth=mp chi and r0=r chi. The local positive polynomial is

    P = r0^4 { [w+b y+x(alpha w+b y)]^4 - [w+b y]^4 }
        + (1+z)^4 { (1+a_birth y)^4-(a_birth y)^4 }
        + (a_birth y)^4.                            (11)

Each displayed difference removes a full coefficient slice, leaving
nonnegative coefficients. The first term describes a quartet hit by
old availability: old E is forbidden, every new position is covered,
and old A moving to new E pays alpha. The second term describes an
unhit quartet with partial new birth coverage; old E is allowed. The
last term permits fresh full coverage only when both old sets are empty.
Thus old A lies in full new quartets, old E lies in nonfull new quartets,
and E_new=Gamma(A_prev)\A_new.

The coefficient of x^u z^e y^v w^f in P^128, divided by
multinom(512;u,e,512-u-e), is the transition in the actual two-count basis.
Start at (0,0) and restrict the free terminal to E_T=empty, summing every
A_T. This enforces the availability fixed point, but does not assert
realizability by actual nonzero field states; that remains an upper count.
Fix candidate masks before partition averaging as in (8).

All m,b weights are charged on new A, so the full expression is
e_(0,0) K^4096 t with t_(u,e)=1[e=0]. The zero-coverage operator changes:
Z_(u,e),(0,0)=1[u=0], including incoming states with e>0. In vector form
Z=v e_(0,0)^tr, v_(u,e)=1[u=0], and Z^2=Z. Do not replace this by the
single-state zero operator used for (8); that would drop valid histories
whose old unavailable coverage ends before an empty region. Regional
composition and occupied-region counting must use this general rank-one Z.

The uncompressed basis has 131841 states, so a dense implementation is
not appropriate. The next task is a positive compression or matrix-free
evaluation of this GLOBAL interface, with exact small checks. Merely
independently multiplying a separation probability into (8) is invalid:
the same partition controls both separation and the outgoing weights.
