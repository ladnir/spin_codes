# Additive-coset coordinates for the persistent MDS inner

This construction folds the persistent state update into the transforms
used by the feedback. It preserves the two properties needed by the
[persistent-state proof](PERSISTENT_MDS.md): exact joint-uniform feedback
on every nonzero state, and MDS constraints at every state reset.
Consequently the same finite certificate bounds apply. This note records
the algebra and its operation count. Independent native checks and timings
are recorded in [Native persistent-state gate](NATIVE_PERSISTENT.md).

The first implementation target is F=GF(2^32), with four payload lanes
per 128-bit record. Uniform GF(2^128) feedback is a separate extension
specified below. The outer and the single stratified route are unchanged.

## 1. Fixed points and the Gray-code schedule

Write E=GF(256) and J=GF(2^32). The construction works over F=J or an
extension F containing J; set q=|F|. Use the established AES-byte embedding
of E in J. Let gamma_0,...,gamma_7 be the fixed Cantor basis of E, with

    gamma_0=1,       gamma_i^2+gamma_i=gamma_(i-1).

Choose b=2^h with 1<=h<=7. For 0<=j<b, define

    beta_j = sum_(i=0)^(h-1) bit_i(j) gamma_i,
    V_h = {beta_j : 0<=j<b}.

In particular beta_0=0. Zero is a valid evaluation point; the multiplicative
construction's requirement that the beta_j be nonzero is not needed here.

Fix an F2-basis of J whose first eight vectors are gamma_0,...,gamma_7.
In the degree-four byte-tower representation, the remaining vectors can
be the raw coordinate bits 8,...,31. Form the quotient basis

    eta_0,...,eta_(31-h)
      = gamma_h,...,gamma_7, raw_bit_8,...,raw_bit_31.

The phrase raw_bit refers to this fixed tower representation. It is not
an assumption about an arbitrary application's polynomial basis for J.
These vectors, together with gamma_0,...,gamma_(h-1), form a basis of J.

Let n=N/b and require n<2^(32-h). For the reflected binary Gray code
g(m)=m xor (m>>1), define

    a_(-1)=0,
    a_t=sum_l bit_l(g(t+1)) eta_l,       0<=t<n,
    Omega_t=a_t+V_h.

The n output cosets and the predecessor coset Omega_(-1)=V_h are distinct.
Every scalar output position (t,j) therefore receives a distinct point
a_t+beta_j in J. Define

    delta_t=a_t+a_(t-1)=eta_(v_2(t+1)),                        (1)

where v_2(m) is the exponent of two dividing the positive integer m.
Every delta_t lies outside V_h, including delta_0. The initial batch does
not require a degenerate identity transition or an exceptional proof case.

## 2. The local interpolation matrix

For delta outside V_h, let D_delta send the values of a degree<b polynomial
on V_h to its values on delta+V_h, preserving the beta_j order. For any
fixed degree-graded polynomial basis, this is

    D_delta = Eval_delta Interp_0.                            (2)

These are b-point transforms. With the existing novel basis, no conversion
to or from monomial coefficients is required between them.

Let s_h(X)=product_(v in V_h)(X-v), and put kappa=s_h'(X).
The derivative is a nonzero constant. For the chosen Cantor basis it is
one. Lagrange interpolation gives the entrywise formula

    (D_delta)[i,j]
      = s_h(delta)/(kappa*(delta+beta_i+beta_j)).             (3)

Indeed s_h(beta_i+delta)=s_h(delta), and the Lagrange denominator at each
beta_j is kappa. All denominators in (3) are nonzero since delta is outside
V_h. Equation (3) proves that D_delta is symmetric.

Evaluating at delta+V_h is translation of the underlying polynomial by
delta. Translation preserves degree<b and is its own inverse in
characteristic two. Thus

    D_delta^T=D_delta,       D_delta^2=I.                     (4)

The 2b points in V_h and delta+V_h are distinct. Consequently

    G_delta = [ I ; D_delta ]

generates an MDS [2b,b] evaluation code: any b rows are independent, and
G_delta Z has weight at least b+1 for Z!=0. Equivalently, every square
minor of D_delta is nonzero. This implies that any b columns of

    H_delta = [ D_delta , I ]

are independent. To verify the last implication, select identity columns
first; the remaining determinant is a square minor of D_delta.

## 3. State coordinates and the transpose issue

For each output coset, define its b-by-b moment matrix

    C_t[i,j]=(a_t+beta_j)^i,          0<=i,j<b,

with row zero identically one. Define C_(-1) similarly on V_h.
These matrices are invertible. Consider a conceptual accumulated syndrome

    S_0=0,       S_(t+1)=S_t+C_t Y_t.                        (5)

The implementation will not form these global moment matrices. It stores
the state in the preceding coset's coordinates:

    Z_t=C_(t-1)^(-1) S_t.

The resulting transition matrix is C_t^(-1)C_(t-1). Its equality with
the evaluation transform in (2) needs justification: moments and
evaluations normally differ by a transpose.

Let M(a) be the lower-triangular translation matrix

    M(a)[i,k]=binom(i,k) a^(i-k),       k<=i,

with binomial coefficients interpreted in characteristic two. Then
C_t=M(a_t)C_(-1), and M(a)M(a')=M(a+a'). Hence

    C_t^(-1)C_(t-1)=C_(-1)^(-1) M(delta_t) C_(-1).

In monomial coordinates, the evaluation transfer is the transpose of
this expression. Equation (3) makes that transfer symmetric, so

    C_t^(-1)C_(t-1)=D_(delta_t).                             (6)

Equation (6) is independent of the polynomial basis used to implement the
two transforms. It is not a generic rule that a moment map and evaluation
map coincide; the common additive-subspace derivative is essential.

There is another way to check the basis convention. Let E_a evaluate a
fixed degree-graded basis on a+V_h, and let J_rev reverse coefficient order.
The rows of J_rev E_a^(-1) evaluate fixed polynomials of successive degrees
0,...,b-1. Those polynomials are independent of a: changing the coset
changes only the constant term of its vanishing polynomial, whereas each
coefficient of its Lagrange quotient uses only nonconstant terms. This
gives an invertible common change of moment rows. In these coordinates,
the transition is directly E_(a_t) E_(a_(t-1))^(-1), again (2).

## 4. Feedback and the actual recurrence

Independently for each batch, sample two vectors r_(t,0),r_(t,1) from F^b.
All coordinates are independent and uniform, including zero. These coins
are independent of the outer route. Sample them once and reuse them for
every encoded message. The Gray schedule and all D matrices are fixed.

Start with Z_0=0. Set D=D_(delta_t) and compute

    V_t       = D Z_t,
    Y_t       = X_t + D(r_(t,0) .* Z_t) + r_(t,1) .* V_t,
    Z_(t+1)   = V_t + Y_t,                                  (7)

where .* denotes coordinatewise multiplication in F. The last equation
is (5) expressed through (6). Emit all N symbols and leave the terminal
state unconstrained.

For fixed Z!=0, the random feedback in (7) is

    H_delta diag(r_0,r_1) G_delta Z.                         (8)

At least b+1 coordinates of G_delta Z are nonzero. Their independently
scaled values are uniform field elements, and any b corresponding columns
of H_delta have rank b. Therefore (8) is jointly uniform on F^b.
It is independent of earlier batches' keys, conditional on the state
before exposing the current keys.

The resulting random b-by-b matrix is not sampled with independent uniform
entries. That stronger property is unnecessary. Equation (8) proves the
exact uniform-vector property used by the persistent-state proof.

For a fixed output Y, the inverse reconstructs Z deterministically using
Z_(t+1)=D_(delta_t)Z_t+Y_t, then recovers X from the middle line of (7).
Thus the inner is F-linear and invertible for every setup. Its state for
a fixed output does not depend on any random feedback key.

## 5. Global MDS resets and certificate transfer

The coordinate changes in (6) preserve whether the state is zero. A closed
excursion from batch a through z satisfies, in the conceptual coordinates,

    sum_(t=a)^z C_t Y_t=0.                                  (9)

Every column in this equation is the length-b Vandermonde column at its
output position's point a_t+beta_j. Those points are globally distinct.
Consequently a support of size w has rank min(b,w), and every nonzero
closed excursion has at least b+1 nonzero output symbols.

The remaining interfaces of PERSISTENT_MDS.md also hold without change:

- A nonzero output batch at zero state activates the state, because
  Z_(t+1)=Y_t there. A nonzero quiet singleton cannot occur.
- Zero output cannot erase a live state, because D_delta is invertible.
- For a fixed output, all inverse batches with nonzero entering state are
  independently uniform on F^b. This gives the same affine-fiber bound.
- For a fixed input and exposed past, a required reset at a live terminal
  batch has probability q^(-b): Y_t is uniform and must equal D_delta Z_t.
- The terminal state remains free. The same excursion representations,
  first-batch activation costs, and output-volume counts apply.

Therefore both the inverse-volume and projective-reset envelopes, their
minimum, the regional moment bounds, and the complete finite sum are
identical functions of (N,B,b,q,e). Every completed certificate for those
parameters in PERSISTENT_MDS.md transfers to (7).

This transfer is a proof of the same failure bound for a new setup
distribution. It is not equality of the sampled matrix distributions,
nor a claim that corresponding deterministic benchmark seeds give the
same code. It introduces no additional global permutation.

## 6. Transform schedule and cost

For each batch, equation (7) applies D_delta twice: once to Z and once to
r_0 .* Z. Each application is one b-point interpolation on V_h followed
by one b-point evaluation on delta+V_h. The total is

    four b-point transforms,
    2b variable field products,
    local XORs and the retained b-symbol state.

There is no separate C Y moment transform, no monomial-basis conversion,
and no separate multiplication by a diagonal persistent-state map. The
computed V_t is used both in the feedback and in the state update.

The inverse transform on V_h always uses GF256 coefficients. The forward
transform uses coefficients in GF256 whenever delta lies in GF256; other
shifts require GF32 coefficients. These coefficients are deterministic
and can be prepared for the finite set of possible Gray-code increments.
Their arithmetic and preparation are not assumed free.

The quotient-basis order in Section 1 puts its 8-h GF256 increments first.
By (1), a non-GF256 increment occurs exactly when t+1 is divisible by
2^(8-h). Among n transitions, their number is floor(n/2^(8-h)). When n
is divisible by this period, the fractions are

| Batch b | GF256 shifts | Other GF32 shifts |
|---|---:|---:|
|32|7/8|1/8|
|64|3/4|1/4|
|128|1/2|1/2|

All shifts remain outside V_h, including the cheaper GF256 shifts.
Changing coefficient arithmetic by this deterministic distinction does
not change either the code distribution or the proof.

These counts identify work removed from the earlier persistent kernel.
They do not establish latency or a speedup. In particular, the cost of
the less frequent GF32-coefficient transforms must be measured together
with the outer, route, key traffic, and state updates.

## 7. Payload lanes and a later GF128 instantiation

For the first target F=J, a 128-bit record can carry four J-symbol lanes
encoded by the same sampled J-linear code. If that code has distance d_min,
every nonzero record-valued word has at least d_min nonzero records: choose
a nonzero lane, whose support is contained in the record support. There
is no four-lane failure union when all lanes share the same sampled code.

One may also extend that fixed J-linear code to GF(2^128) without changing
its distance. This extension does not turn J-valued keys into independent
uniform GF(2^128) keys or obtain a new larger-alphabet distance bound.

A genuine F=GF(2^128) instantiation instead samples all r coordinates
uniformly from F. The schedule and D coefficients remain in J. Equations
(8)--(9) then hold over F and transfer the corresponding q=2^128 certificate.
Its full-field product cost must be measured separately; the GF32 SIMD
timing cannot be assigned to it.

The GF32 implementation now passes direct Cauchy evaluation, complete
recurrence, and conceptual-syndrome checks. Matched scalar-route timings
are about23--24 ms for B64/b64, versus34--35 ms for the multiplicative
control. The scalar-observation and s<b variants remain separate, unproved
changes. Packet routing has its own proof and certificates in PACKET_ROUTE.md.

## 8. Full Cantor schedule: cheaper deterministic coefficients

The optional `SPIN_OR_FULL_CANTOR` mode replaces the raw coordinate bits
in Section1 by an extension gamma_8,...,gamma_31 of the existing Cantor
chain. It changes the deterministic output points, not the route or the
distribution of random feedback diagonals. The original schedule remains
the default control.

Such an extension exists in J. The F2-linear map Delta(x)=x^2+x satisfies
Delta^32=0 and has one-dimensional kernel, so it has one Jordan chain of
length32. Any initial chain beginning at1 extends to that length. In
particular, for i<=30, gamma_i has absolute trace Delta^31(gamma_i)=0,
so Delta(x)=gamma_i is solvable. Either root extends the chain and preserves
binary independence. Setup chooses one root by bounded binary elimination,
then verifies all relations and the full rank32 basis.

Use quotient directions gamma_h,...,gamma_31 in the same Gray schedule.
They still enumerate distinct cosets of V_h. The proof of exact uniformity,
global MDS resets, and each finite certificate therefore remains valid.
No additional permutation or per-message randomness is introduced.

This schedule makes Delta^i(gamma_k)=gamma_(k-i) a byte-field element as
soon as i>=k-7. In a transform, only its bottom k-7 levels then require
GF32 products; upper levels use GF256 products. The raw-bit schedule had
this shortcut at k8 but generally not at later directions. The reduction
is exact, not an approximation of a uniform coefficient. Native timings
must still establish whether the saved arithmetic improves total latency.
