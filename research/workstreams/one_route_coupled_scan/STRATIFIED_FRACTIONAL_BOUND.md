# A finite fractional-moment bound for stratified routing

This note gives a finite bound for fixed MDS outer blocks, one stratified
permutation, and the uniform scalar recursive inner. It avoids charging
the full multiplicity of every rare route. The
[finite checkpoint](STRATIFIED_CLOSED_CERTIFICATE.md) evaluates this bound
and its batched extension. This note itself derives the bound; it contains
no implementation timing.

## 1. Construction and shortening

Fix a field F_q, an even block length B, k=B/2, M blocks, and N=MB.
Each outer block is a fixed injective [B,k] MDS code over F_q. The blocks
may share the same encoder. A Reed--Solomon encoder over a subfield is
permitted if it remains MDS over F_q.

There are B consecutive regions, each of length M. Region j receives
coordinate j from every outer block, in an independently sampled uniform
permutation of the M block labels. This defines one permutation of the
N outer coordinates. Write C for the routed outer code.

Independently sample every scalar feedback coefficient uniformly from F_q.
For memory L>=1, encode a routed word x by

    y_t = x_t + sum_(j=1)^L alpha[t,j] y_(t-j),
    y_t=0 for t<=0.

All coefficients and permutations are fixed once for every message. The
inner is invertible and has a free terminal state. Distance counts nonzero
F_q symbols. Fix a failure threshold e in {0,...,N}.

For a fixed coordinate set U, let u=|U| and let t_i count its coordinates
from outer block i. The MDS property gives the exact shortened dimension

    D(U) := dim(C intersect F_q^U) = sum_i (t_i-k)_+.       (1)

Here F_q^U denotes vectors supported inside U, and z_+=max(z,0).
In particular, D(U)=0 excludes every nonzero input supported in U.

## 2. Output intervals and their counts

For a nonzero output y, let A be the positions whose preceding L output
symbols are not all zero. Put U=A union supp(y). Cluster successive nonzero
outputs whenever their gaps are at most L. Each cluster gives an interval
from its first nonzero position to L positions after its last, truncated
at N. The c intervals are disjoint but may be adjacent. Their first
positions are the activation positions, so

    |A|=u-c.                                              (2)

For fixed y and route, its inverse input is uniform on an affine cube of
dimension u-c. Therefore

    Pr_inner[I^(-1)(y) in C] <= q^(D(U)-u+c).              (3)

The exact affine-intersection dimension is D(A); replacing it by D(U)
only enlarges this bound. If D(U)=0, the probability is zero: the first
activation input equals the nonzero first output.

Separate two classes according to the last nonzero output position.

| Class | Forced zero coordinates | Free coordinates v | Representation count n(u,c) |
|---|---:|---:|---|
| Complete: last nonzero at most N-L | Lc | u-Lc | binom(v-1,c-1) binom(N-u+c,c) |
| Truncated: last nonzero greater than N-L | L(c-1) | u-L(c-1) | binom(v-1,c-1) binom(N-u+c-1,c-1) |

In the complete class, each interval has length at least L+1. Subtract
its L forced zeros and count compositions of v into c positive parts.
The c+1 outside/inter-interval gaps are nonnegative and total N-u.
In the truncated class, the last interval ends at N, so its final gap is
fixed to zero. The same argument leaves only c nonnegative gaps. The
counts may include assignments that do not realize the proposed clusters;
this only enlarges the bound. A binomial coefficient with invalid integer
arguments, and any class with v<c, contributes zero.

Define the q-ary volume

    V_q(v,e) = sum_(j=0)^min(v,e) binom(v,j) (q-1)^j,

with value zero if v<0 or e<0. Each representation has at most

    W(u,c) := (q-1)^c V_q(v-c,e-c)                        (4)

eligible output assignments: its c activation symbols must be nonzero.
The weaker V_q(v,e) is also a valid choice. Keep q in this volume even
when introducing a fractional moment below.

## 3. An elementary route moment bound

Fix U before sampling the route. Suppose U contains f whole regions and
h partially occupied regions. In partial region j, let p_j be its selected
fraction of positions. For a fixed block i,

    t_i = f + sum_(j=1)^h Z[i,j],
    Z[i,j] independently distributed as Bernoulli(p_j).

Indicators in one region are not independent across blocks: they form a
uniform subset of prescribed cardinality. The following comparison keeps
that dependence on the lower side of the bound.

Let theta>0 and g(s)=exp(theta*(s-k)_+). Condition on all regions except j.
Each block factor has the form

    g(v_i+Z[i,j]) = a_i+b_i Z[i,j],       a_i,b_i>=0.

Expand the product over i. For any d distinct block labels, their joint
inclusion probability in a uniform subset of size s_j is

    (s_j)_d/(M)_d <= (s_j/M)^d,

where (s)_d is the falling factorial. Thus replacing this region by
independent Bernoulli(p_j) indicators cannot decrease the expectation.
Repeat for every partial region. All block factors then become independent,
giving

    E_route exp(theta D(U))
       <= [E exp(theta*(f+sum_j Z[j]-k)_+)]^M.             (5)

No independence between the original block occupancies is assumed.

## 4. Reducing partial-region probabilities to one binomial

Put mu=u/M. Then sum_j p_j=mu-f. For any convex function g, averaging two
Bernoulli probabilities at fixed sum increases E g of their sum. Indeed,
the coefficient of their product is the nonnegative second difference
of g, and their product is maximized when the probabilities are equal.
Repeated averaging, followed by continuity, gives

    E exp(theta*(f+sum_j Z[j]-k)_+)
       <= E exp(theta*(f+Bin(h,(mu-f)/h)-k)_+).            (6)

At most two regions per interval are partial. Consequently use

    r = min(2c,B)             in the complete class,
    r = min(2c-1,B)           in the truncated class.       (7)

The second bound uses the last interval's endpoint N. It is not necessary
to maximize over every h<=r. Add deterministic Bernoulli variables until
there are r: an empty region contributes 0, while a whole region
contributes 1 and decreases f by one. This leaves the occupancy distribution
unchanged. Apply (6) to the resulting r probabilities.

For 0<=s<=B, define

    H_r(s,theta)
       = max_f E exp(theta*(f+Bin(r,(s-f)/r)-k)_+),        (8)

where f ranges over integers satisfying

    0<=f<=B-r,       0<=s-f<=r.

The expectation in (8) is an explicit sum of r+1 nonnegative binomial
terms. Equations (5)--(8) imply

    E_route exp(theta D(U)) <= H_r(mu,theta)^M.           (9)

There is also a deterministic zero test. If any region is partial, then
f<=ceil(mu)-1; otherwise t_i=mu for every block. In either case,

    ceil(mu)+r-1 <= k  implies  D(U)=0.                   (10)

For example, a single terminal interval with u<=kM cannot support a
nonzero outer word. Apply (10) before evaluating a moment bound.

## 5. Preserve the fixed-size identity through complementation

The route always satisfies sum_i t_i=u. At rate one half, (1) gives

    D(U) = u-N/2+D(U^c).                                 (11)

The complement has exactly the same partial regions, although it may
have more interval components. Therefore its bound in (8) uses the same
r from (7), not a larger component-based value. Applying (9) to U^c gives

    E_route exp(theta D(U))
       <= exp(theta*(u-N/2)) H_r((N-u)/M,theta)^M.         (12)

Take the minimum of (9) and (12), and denote it by K(theta;u,c).
Both are upper bounds on the same moment, so K>=1 mathematically.
Identity (11) must be applied before the independent-Bernoulli comparison;
that comparison does not preserve the total occupancy.

If u>N/2, the positive shift in (11) already makes D(U)>0. Outcomes with
D(U^c)=0 still contribute. Do not discard those complement outcomes.

## 6. Fractional moment and the complete finite sum

Fix one interval representation and fix its routed outer. Summing (3)
over its eligible outputs gives a conditional failure bound

    1[D(U)>=1] min(1, W(u,c) q^(D(U)-u+c)).              (13)

For 0<nu<=1, min(1,z)<=z^nu. Set theta=nu ln q. Averaging (13) over
the actual route therefore gives

    exp(theta*(c-u+log_q W))
          * E_route[exp(theta D(U)) 1[D(U)>=1]].          (14)

This is a change of proof weight, not a change of field or coefficient
distribution. The bound applies to the existence of a bad output, not
to one message chosen before setup.

For every nonnegative integer D,

    exp(theta D) 1[D>=1]
       <= (exp(theta D)-1)/(1-exp(-theta)).

Consequently define

    J(theta;u,c)
       = min(K(theta;u,c),
             (K(theta;u,c)-1)/(1-exp(-theta))).           (15)

This bounds the expectation in (14). If u>N/2, use J=K: subtracting
one from the unshifted complement moment would incorrectly discard
positive-D(U) outcomes. Equivalently, (15) already selects K because
D(U)>=u-N/2>=1.

For each class and valid pair (u,c), set its contribution to zero if
W=0 or (10) holds. Otherwise let

    F(u,c) = min(1,
        inf_(0<theta<=ln q)
          exp(theta*(c-u+log_q W(u,c))) J(theta;u,c)).      (16)

Then the full setup-failure probability satisfies the finite bound

    Pr[d_min<=e]
       <= min(1, sum_classes sum_(u,c) n(u,c) F(u,c)).     (17)

All sums are finite: 1<=u<=N and 1<=c<=min(e,N). Invalid counts vanish.
The parameter theta may be selected separately for each class and (u,c).
Every selected theta yields a valid upper bound; a finite tested set of
values does not need to attain the infimum in (16).

## 7. What a numerical certificate must cover

A complete certificate must upper-bound the entire sum (17), including
every valid u and c. A grid of selected activity sizes is only a diagnostic
unless interpolation or a monotonic envelope covers the omitted values.
The maxima in (8) must cover every feasible integer f. The field volume,
binomial terms, exponentials, and final sum require certified upper bounds.
Near K=1, evaluate the difference in (15) without uncontrolled cancellation.

No small-coefficient inner follows from this argument. Its inverse step
uses independent uniform symbols over the full field. The linked finite
checkpoint records audited arithmetic and complete coverage for selected
parameters; changing the coefficient distribution requires a new argument.
