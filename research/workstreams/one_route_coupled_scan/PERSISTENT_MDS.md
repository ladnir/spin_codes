# Persistent MDS state with uniform feedback

Status: an analytic construction with complete finite distance certificates
listed below. Subsequent [native implementations and timings](NATIVE_PERSISTENT.md)
cover its GF32 profiles and the additive-coordinate implementation.
The [window-state certificates](STRATIFIED_CLOSED_CERTIFICATE.md) and their
[native measurements](NATIVE_BASELINE.md) remain separate checkpoints.

The proposed change replaces a window of previous outputs by a smaller
persistent state. Resetting that state imposes MDS constraints on the output
symbols since its last reset. This supplies a different counting argument;
it does not inherit the window state's forced zero tails.

## 1. Field, outer, route, and setup

Let F=GF(2^d), where d is a positive multiple of 32, and write q=|F|.
Identify its subfields E=GF(256) and J=GF(2^32). Fix a primitive element
alpha of J and set

    R = (2^32-1)/(2^8-1) = 16843009.

Fix an even outer length B<=256, k=B/2, and N=MB. The outer is a direct
sum of M copies of one fixed [B,k] MDS code over E, extended linearly to F.
Thus K=N/2. An RS code with distinct subfield evaluation points suffices.

Use the stratified route from the finite checkpoint: coordinate j from each
outer block enters region j, and each region independently permutes its M
block labels uniformly. Concatenating the B regions gives one permutation
of the N outer coordinates. There is no second global routing stage.

Choose a batch size b dividing N, and set the state dimension s=b. Require
b<=255 and n:=N/b<=R. Fix distinct nonzero beta_0,...,beta_(b-1) in E.
Define the deterministic b-by-b matrices

    P = diag(1,alpha,...,alpha^(b-1)),
    C[i,j] = beta_j^i,                 0<=i,j<b.

Both matrices are invertible. For each batch t, independently sample an
F-linear map T_t:F^b->F^b with the following property:

    For every fixed S!=0, T_t S is uniform on F^b.             (1)

These maps are independent of the route. All setup randomness is sampled
once and fixed for every message. Property (1) concerns its setup
distribution, not fresh randomness during encoding.

Start with S_0=0. On input batch X_t, emit and update

    Y_t = X_t + T_t S_t,
    S_(t+1) = P S_t + C Y_t,              0<=t<n.             (2)

Emit all N symbols and discard the terminal state. No return to zero is
required at the end. The state update uses only fixed local arithmetic.

### Two exact realizations of the feedback

A uniform dense b-by-b matrix satisfies (1): its rows are independent,
and each row's dot product with a fixed nonzero S is uniform on F.

For 2b<=256, there is also a spectral realization with 2b variable field
products per batch. Fix 2b distinct points Omega in E and a degree-graded
polynomial basis P_0,...,P_(2b-1), with deg P_i=i. Interpret S as the first
b coefficients of a polynomial p_S in this basis. Let G evaluate p_S at
Omega. Let H interpolate a degree<2b polynomial and retain its highest b
coefficients in the same basis. Independently sample r_a from F for each
a in Omega, including zero, and define

    T = H diag(r) G.                                            (3)

Every b columns of H are independent. In monomial coordinates, those
columns are row-equivalent to a Vandermonde matrix with nonzero weights
1/L'(a), where L(X)=product_(a in Omega)(X-a). Indeed, the highest b
coefficients of L(X)/(X-a) are unit-triangular combinations of
1,a,...,a^(b-1). A degree-graded basis changes these b rows invertibly.

For S!=0, p_S has degree at most b-1 and is nonzero at at least b+1 points
of Omega. Consequently H diag(GS) has rank b. Applying it to the independent
uniform r_a proves joint uniformity in (1), not only uniform marginals.
Different batches use independent diagonals. This is the same rank
principle as the [spectral window map](SPECTRAL_BATCH.md), with a different
state size.

## 2. Inverse and global reset constraints

Given Y, compute S by the second equation of (2), then recover

    X_t = Y_t + T_t S_t.

Thus the inner is F-linear and invertible for every setup. The complete
code has rate exactly one half. For fixed Y, its entire state sequence is
deterministic and independent of every T_t. This is the property needed
by the inverse-output proof below.

Suppose S_a=0 and inspect the state after batches a,...,z. Unrolling gives

    S_(z+1) = sum_(h=a)^z P^(z-h) C Y_h.                       (4)

The column for output coordinate (h,j) in this equation is

    (1,xi_(h,j),...,xi_(h,j)^(b-1))^T,
    xi_(h,j) = alpha^(z-h) beta_j.

All these points are distinct. A coincidence for h!=h' would put
alpha^(h-h') in E*. The order of alpha modulo E* is R, whereas
0<|h-h'|<n<=R. At the same time h, distinct beta_j give distinct points.

It follows that any w selected output coordinates in (4) have column
rank min(b,w). In particular:

- A nonzero output segment starting and ending at zero state has weight
  at least b+1.
- On a fixed support of size w>=b+1, its reset equations have rank b.
- Disjoint closed segments impose these constraints on disjoint output
  variables, so their ranks add.

Zero output cannot extinguish a live state because P is invertible.
Zero input does not have this deterministic guarantee: P+C T_t may be
singular. The proofs below account for such resets through their setup
probabilities instead of excluding them.

## 3. Excursions and their representations

Fix a nonzero output Y. Ignore batches with S_t=Y_t=0. An excursion starts
at a batch with S_t=0 and Y_t!=0. Since C is invertible, its next state is
nonzero. Include subsequent batches through the first return to zero,
or through batch n-1 if the state never returns.

These intervals are disjoint and ordered; consecutive intervals may abut.
Every interval except possibly the last is closed. Each closed interval
has at least two batches and at least b+1 nonzero output symbols. The last
interval has at least one batch and at least one nonzero output symbol.
An open last interval actually reaches the end of the stream. Allowing
other endpoints only enlarges the representation count below.

For a representation with c excursions and v total batches, let U contain
all scalar coordinates in those intervals. Delete the first batch of each
interval to obtain A. Write r=c-1. Then

    u:=|U|=bv,       a:=|A|=b(v-c),       r=c-1.                (5)

The entering state is nonzero exactly at batches in A. The input is zero
outside U. If the output weight is at most e, the feasible indices satisfy

    1<=c,       2c-1<=v<=n,       (b+1)(c-1)+1<=e.            (6)

For fixed (v,c), the number of possible interval representations is at most

    J(v,c) = binom(v-c,c-1) binom(n-v+c,c).                    (7)

To see this, assign minimum lengths two to the first c-1 intervals and
one to the last. Composing the remaining length gives the first factor.
The c+1 nonnegative gaps, including prefix and suffix, give the second.
This count permits abutting intervals, a closed singleton last interval,
and open intervals ending early. Those impossible cases are harmless
overcounts; every actual nonzero output is included.

For a fixed route, write C_out for the routed outer code and define

    D(S) = dim_F(C_out intersect F^S).

If D(U)=0, no nonzero output with this representation is possible, since
its inverse input belongs to C_out intersect F^U.

## 4. Inverse-output envelope

Fix a route, a representation, and an output y realizing it. At every
batch in A, S_t(y)!=0. Property (1) and independence across batches imply
that the inverse input is uniform on an affine translate of F^A. The
translate lies in F^U. Its intersection with C_out has either zero or
q^D(A) points. Therefore

    Pr_inner[I^-1(y) in C_out] <= q^(D(A)-a).                  (8)

D(A)=0 does not imply an empty affine intersection. That outcome cannot
be removed from the moment used in this envelope.

Let V_q(u,e)=sum_(w=0)^min(u,e) binom(u,w)(q-1)^w. A uniform upper bound
on the number W of eligible outputs for the representation is

    W <= (q-1)^(-br) V_q(u,e).                                (9)

For each exact support of weight w that can occur, all r guaranteed closed
intervals have at least b+1 support coordinates. Choose b pivot coordinates
in each of them. The remaining w-br coordinates must be nonzero; the reset
equations determine the pivots uniquely. Thus at most (q-1)^(w-br)
assignments realize that support. Summing all supports, including invalid
ones, proves (9). The exact support restrictions or the exact MDS weight
enumerators can replace this relaxation later.

For numerical evaluation, the same argument permits q choices for each
free support coordinate. Together with the dimension of the full reset
space, this gives the two alternative bounds

    W <= q^(-br) sum_(w=0)^min(u,e) binom(u,w)q^w,
    W <= q^(u-br).

The first counts valid supports using q^(w-br) assignments, then includes
invalid supports as an overcount. It does not replace q-1 by q in the
negative exponent of (9). Write

    Ent(u,e) = u ln u-e ln e-(u-e)ln(u-e),

with 0 ln 0=0. For 0<e<u and q>=u, the Chernoff choice
z=e/[q(u-e)]<=1, the bound sum_w binom(u,w)<=2^u, and the full-space bound
give

    log_q W <= min(u, e+u/d, e+Ent(u,e)/(d ln 2))-br.          (9a)

For e>=u use u-br. All feasible classes in (6) have e>=1. Formula (9a)
is a separate valid upper bound and is suitable for directed evaluation.

Equations (8)--(9) give, conditional on the route,

    p_vol <= 1[D(U)>0] min(1, W q^(D(A)-a)).                 (10)

For any 0<nu<=1 and tau=q^nu, dropping the indicator and using
min(1,x)<=x^nu yields

    E_route p_vol <= W^nu tau^(-a) E_route tau^D(A).          (11)

One may also divide the actual output count by q-1 before this step,
because all nonzero scalar multiples have the same representation and
membership event. Formula (11) does not rely on that optional improvement.

## 5. Forward envelope from reset probabilities

This second bound keeps source dimensions that the output-volume count
can discard. Fix the route, a representation, and a nonzero input
x in C_out intersect F^U before exposing the inner maps.

At the last batch of each of the first r excursions, the entering state
must be nonzero. Condition on all earlier maps. If the required earlier
conditions have failed, this branch contributes zero. Otherwise, (1)
makes Y_t uniform on F^b, even though X_t is fixed. Since C is invertible,
there is exactly one value

    Y_t = C^-1 P S_t

that resets the state. It has probability q^-b. Sequential exposure at
the r distinct terminal batches therefore bounds the probability of all
required resets by q^(-br). Additional state-pattern and output-weight
requirements can only reduce that probability.

All events here are invariant under nonzero F-scaling of x. There are

    L_D = (q^D-1)/(q-1) = sum_(j=0)^(D-1) q^j,
    D=D(U),

projective input lines, with L_0=0. A union bound gives

    p_reset <= min(1, L_D q^(-br)).                            (12)

This is a whole-code existence bound: its representatives are fixed after
the route and before the independent inner setup. It does not restrict
the eventual bad message to have been chosen before setup.

For tau=q^nu, subadditivity of x^nu on nonnegative terms gives

    L_D^nu <= sum_(j=0)^(D-1) tau^j = (tau^D-1)/(tau-1).

Hence

    E_route p_reset
      <= tau^(-br) (E_route tau^D(U)-1)/(tau-1).             (13)

Unlike (11), this envelope retains an exact zero-dimensional correction.
The reset probability only needs C invertible and uniform feedback. The
global MDS property supplies the stronger weight restriction in (6) and
the output count in (9).

## 6. Reuse of the regional route moments

For the actual fixed MDS outer and stratified route,

    D(S) = sum_i (t_i(S)-k)_+,
    D(S) = |S|-N/2+D(S^c),                                  (14)

where t_i counts selected coordinates from source block i. These identities
are unchanged. Both U and A are unions of at most c intervals, so each meets
at most min(2c,B) partial regions.

The [regional comparison](STRATIFIED_FRACTIONAL_BOUND.md) supplies the
following moment envelope. Put rho=min(2c,B), and define

    H_rho(mu,tau)
      = max_f E tau^((f+Bin(rho,(mu-f)/rho)-k)_+),

where the maximum covers integer f satisfying
0<=f<=B-rho and 0<=mu-f<=rho. For a union S of at most c intervals,
complementing before comparison gives

    E_route tau^D(S) <= Moment_c(|S|,tau),

    Moment_c(t,tau)
      := tau^((t-N/2)_+)
         H_rho(min(t,N-t)/M,tau)^M.                         (15)

Padding partial-region indicators justifies the same rho for fewer actual
partial regions. Fixed-cardinality sampling and the complement identity
must be retained; replacing the route by independent input positions is
not part of this construction.

For a fixed (v,c), equations (11), (13), and (15) give two uniform bounds:

    U_vol(v,c;tau) = W^(log_q tau) tau^(-a) Moment_c(a,tau),

    U_reset(v,c;tau)
      = tau^(-b(c-1)) (Moment_c(u,tau)-1)/(tau-1).           (16)

The first can equivalently use the factor tau^(-min(a,N/2)) and the
complemented H moment. The subtraction in the second must be evaluated
with a valid upper bound on Moment_c-1; cancelling a near-one quantity in
uncontrolled floating point is not a certificate.

If the geometry of U deterministically gives D(U)=0, set the entire term
to zero. For example, u/M+2c<=k is a sufficient coarse test. No such zero
test is justified merely from D(A)=0.

For e>=1, the finite whole-code failure bound is

    sum_(v,c satisfying (6)) J(v,c)
      * min(1, inf_(1<tau<=q) U_vol(v,c;tau),
                 inf_(1<tau<=q) U_reset(v,c;tau)),            (17)

with the deterministic-zero exclusions above. For e=0, injectivity makes
the failure probability zero. A finite set of verified markers may replace
either infimum. Different markers and the two envelopes may be selected
separately for each class: each bounds the same class before the final
union. No independence between the two envelopes is asserted or needed.

The complete directed evaluations below use both envelopes from (17).
The earlier window-state certificates remain separate: they do not certify
this construction.

### Complete finite checkpoints

All rows have K=2^20, N=2^21, rate 1/2, and a free terminal state. The
probability is over ideal setup, fixed once for every message. The margin
is a bound on the probability that any nonzero message has output weight
at most the stated threshold.

| Field | B | Batch=state | Bad-weight threshold | Distance guarantee | Margin at least | Cover rectangles |
|---|---:|---:|---:|---:|---:|---:|
| GF(2^32) |64|64|734003|>35%|46 bits|2846|
| GF(2^32) |128|32|734003|>35%|44 bits|2939|
| GF(2^128) |256|128|943718|>45%|48 bits|710|

The helper [persistent_mds_certificate.py](theory/persistent_mds_certificate.py)
first proposes a rectangular cover in floating point. Each rectangle fixes
one of the two valid envelopes and a rational marker. The independent
replay evaluates the frozen exact rational route moments and outward-rounded
Decimal arithmetic, at both 80 and 120 digits. Every pair (v,c) in (6) is
covered, not just a grid of selected supports.

Each individual (v,c) contribution is at most 2^-75. There are fewer than
2^29, 2^31, and 2^27 possible pairs respectively. This gives the conservative
margins in the table. The worst directed rectangle exponents are
-76.73874519393029, -84.48709388494824, and -75.75965463754208 respectively.
These decimal excerpts are diagnostic; the helper verifies the inequalities
using the full directed values.

Reproduce from this workstream's theory directory:

```text
python -B persistent_mds_certificate.py --batch 64 --field-bits 32 --term-bits 75 --seconds 90
python -B persistent_mds_certificate.py --B 128 --batch 32 --field-bits 32 --term-bits 75 --seconds 90
python -B persistent_mds_certificate.py --B 256 --batch 128 --field-bits 128 --distance 9/20 --term-bits 75 --seconds 90
```

All runs were serial. Floating coverage took about 7, 13, 8 seconds; each
directed replay took under 3 seconds. These are proof-checking times, not
encoder timings. The helper caps coverage at 50000 checks and stores no
generated files. Exact route caching is bounded by the frozen 4096-entry
cache. Independent analytic review checked the new counting and replay.

The separate [exact toy tests](theory/test_persistent_mds.py) pass over
GF(16), using two GF(4) points and the analogous multiplicative-coset
construction. They check 340 reset-column subsets, 510 nonzero-state
spectral ranks, 352 recurrence/inverse fixtures, 920 unrolled prefixes,
and all 247 interval representations through seven batches. Negative
controls verify the coset-order limit and the quiet-kernel trap when s<b.
These checks validate local proof interfaces, not the large-field native
kernel. Run `python -B test_persistent_mds.py` from the same directory.

Source SHA-256 receipts:

```text
persistent_mds_certificate.py 62f1905de9d3f582be8e97a09a5a21773681ea7938c54a7f4b4a73402ac173ee
test_persistent_mds.py       fe875a535d8fc5237e3290d7620385036eeabe1448058013602b47114dd455b6
```

The following smaller-state probes failed the same upper-bound test:

| B | Batch=state | Field bits | Distance | First failed (v,c), log2 envelope |
|---|---:|---:|---:|---|
|64|16|128|35%|(47077,1190), -71.7233|
|64|32|128|35%|(23048,106), -72.8165|
|256|64|128|45%|(15072,384), -70.3228|

These are failures to establish the uniform 2^-75 per-class bound. They
are neither measured low-weight words nor upper bounds on actual distance.
The first failed point does not locate the worst class. No stronger claim
about these profiles is justified by these probes.

## 7. Cost and two unproved continuations

The spectral realization uses 2b full-F products per b output symbols,
plus two size-2b transforms. The persistent update adds b-1 products by
fixed elements of J, a fixed subfield Vandermonde moment map, and XORs.
For power-of-two b<=128, the beta_j can be an affine subspace of E avoiding
zero, allowing a transform-based moment map. Its scaling and basis
conversion must be implemented and checked, not assumed free.

These are algebraic operation counts. They are not cycle counts or a
measured speedup. In particular, the route and the outer still perform
the same work as the current full-field baseline.

### Smaller state than the batch

The recurrence can be defined with s<b by taking the first s rows of C
and an s-dimensional P. The Vandermonde reset lemma then uses rank s.
However, C has a nontrivial kernel. A nonzero Y in ker C can start and
finish at zero state in one batch, without exposing any feedback map.
Such quiet singleton excursions invalidate the q^-s reset charge for
every closed interval and the length restrictions used above.

The code may still have good distance, but it needs a new bound on the
outer/route intersection with these deterministic local kernels. Neither
(12)--(17) nor the window-state certificate is a proof for that variant.

### One observed scalar instead of uniform feedback from every live state

A cheaper proposed feedback is

    Y_t = X_t + r_t ell(S_t),       ell(S)=sum_i S_i,

with independent uniform r_t in F^b and the same persistent update. It
uses one full-F product per output and is jointly uniform when ell(S_t)!=0.
For a nonzero state S followed by zero outputs,

    ell(P^h S) = sum_i S_i (alpha^h)^i.

This is a nonzero polynomial of degree below s evaluated at distinct
points. It vanishes at at most s-1 times in that zero-output run. This
local observation bound does not control complete codewords: nonzero
outputs alter the state, and unobserved live states create holes in the
inverse random-coordinate set. A new rank/observability argument is
required. No distance or implementation claim is made for this follow-up.

The next discriminating step is a matched native implementation of the
two GF32 checkpoints, with direct scalar oracles for feedback and C Y.
Four SIMD field lanes can share one sampled GF32 code in each 128-bit
record. A nonzero record word contains a nonzero lane, so the same
symbol-distance bound gives the record-distance bound without a lane
union. This is not the GF128 construction or a GF256-linear embedding.
Preserve the window construction as the exact-uniform control while
investigating the two weaker interfaces separately.
