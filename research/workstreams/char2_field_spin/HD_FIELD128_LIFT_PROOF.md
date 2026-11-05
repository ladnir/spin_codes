# Proof obligations for the literal GF(2^128) lift

This note audits the outer counting measure and inner transition for the
literal lift of [the certified paired construction](HD_SHIFTED_CLOSED_CERTIFICATE.md).
The formulas below are valid at field order q=2^128. They do not supply a
new numerical distance certificate or a runtime measurement.

## Field and setup

Let F0=GF(256), let F=GF(2^128), and fix a field embedding F0 into F.
Embed every coefficient of the four RS[64,32] constituents, the four
RS[16,4] emission blocks A, C=A^T, and the interpolation matrix H16.
All dimensions count F symbols: K=2^20, N=2^21, and the state has 16 symbols.

Independently for every outer group and column, sample a uniform invertible
4-by-4 Toeplitz matrix over F. Independently for every inner step, sample
both diagonals D1,D2 with entries in F minus zero. All these choices and
the 256 regional permutations and fair pair swaps are mutually independent.
The initial state is zero, the final state is unrestricted, and state
persists across regions. The recurrence is

    Y = X + A Q,
    Q_next = D2 H16 D1 (Q + C X).

This is an F-linear random code. Its randomness is not restricted to F0.
The old componentwise distance argument therefore does not certify this
ensemble. A concrete implementation must use a field embedding; interpreting
an AES-field byte as an arbitrary polynomial-basis field element is insufficient.

Every minor of a fixed F0 matrix has the same zero/nonzero status in F.
Consequently all constituent MDS properties, total nonsingularity of H16,
the identity CA=0, and all principal nullities of I+AC survive the lift.
Since (AC)^2=0, the map I+AC remains its own inverse.

## Outer counting measure

Write q=|F|. First sample an unconditioned Toeplitz matrix T using seven
independent elements of F. For fixed x != 0 in F^4, the map from these
seven elements to Tx has rank four. Indeed, fix the last nonzero coordinate
of x and expose the other diagonals; the remaining four diagonal variables
give a triangular map with nonzero diagonal. Thus Pr[Tx=y]=q^-4.

The number of invertible 4-by-4 Toeplitz matrices is (q-1)q^6.
This finite-field identity is stated and proved by
[Kaltofen and Lobo, ISSAC 1996](https://kaltofen.math.ncsu.edu/bibliography/96/KaLo96_issac.pdf).
Rejection sampling therefore gives, for every fixed x != 0 and y,

    Pr[Tx=y | T invertible] <= q^-4/(1-1/q)
                            = 1/((q-1)q^3) = rho_T.

For y=0 the conditional probability is zero. For nonzero y this is an
atom bound, not a uniform-image claim. The seven entries cease to be
independent after conditioning; the bound does not assume otherwise.
Separately sampled accepted matrices remain independent.

Fix a nonzero complete target group y, and let h count its nonzero
four-symbol columns. Invertibility forces the source columns to vanish
exactly where y vanishes. If h<33, no nonzero RS[64,32] constituent can
have this support. Otherwise the constituent polynomials must vanish on
64-h specified nodes. Each shortened constituent has dimension h-32.

Choose any h-32 nonzero target columns and reveal only their accepted
Toeplitz matrices. Their inverse images, together with the 64-h zero
columns, specify 32 evaluations of each constituent. They determine at
most one source group. The other 32 accepted matrices are still independent
and have their original distributions. Conditional on the revealed matrices,
each unrevealed column either has zero chance or costs at most rho_T.
Hence the expected counting measure at y is at most rho_T^32.

For a fixed set U of r active groups, the expected outer counting measure
is therefore dominated by beta_T(q)^r times uniform measure on F^(256r), where

    beta_T(q) = q^256 rho_T^32
              = q^160/(q-1)^32
              = q^128 (1-1/q)^-32.

The outer is injective, so membership and counting at a fixed target agree.
The comparison measure may also contain zero values in potential positions.
It does not assert that the actual outer output coordinates are independent.
The inherited bound beta_old(q)=q^256/(q-1)^128 is also valid because
rho_T <= (q-1)^-4. Use either bound consistently in each new receipt.

## Local transition

Fix the potential input set for one packet before sampling its refresh.
Under the outer comparison measure, X is uniform on this coordinate subspace
and is independent of the incoming state and the fresh refresh. Let
l=(l_0,l_1,l_2,l_3) be its four block occupancies, with 0<=l_i<=16.
Fix 0<z<1, put t=z+(1-z)/q, and define

    k_l = sum_(c=0)^l binom(l,c) z^(l-c)(1-z)^c q^(-c-min(4,l-c)),
    f_n = sum_(c=0)^n binom(n,c) z^(n-c)(1-z)^c q^-min(4,c),
    K = product_i k_(l_i),             T = t^(sum_i l_i),
    d = q^16-1,
    V = T (q^16 product_i f_(16-l_i)-1)/d.

Here K is the weighted probability of CX=0 from incoming state zero.
V is the output moment averaged over a uniform nonzero incoming state.
The MDS ranks used in both expressions are unchanged over F.

For a fixed nonzero refresh input of support size a, consider a target
of support size b. Set R=D2 H16 D1. The first diagonal makes the a active
input values uniform and nonzero. The s-b output-zero constraints, with
s=16, have full possible rank by total nonsingularity of H16. If a+b<=s,
the event is impossible. Otherwise at most (q-1)^(a+b-s) input assignments
satisfy them. The second diagonal hits the b nonzero target values with
probability (q-1)^-b. Therefore

    Pr[R x = y] <= (q-1)^-16,       x != 0.

Let B=d/(q-1)^16. An incoming weighted measure with zero mass at most Z
and nonzero atoms at most W/d maps to a measure bounded by the row-vector
product (Z,W)M, where

    M = [[K,       B(T-K)],
         [Jbar/d,  B V]].

To define Jbar, let P=I+AC on one 16-coordinate block and let nu_c be the
maximum nullity of P[D,D] over all c-element coordinate subsets D. For
0<=l<=16 and rho=min(4,l), define

    R_l = sum_(c=0)^l sum_(e=0)^min(3,16-l)
          binom(l,c) binom(16-l,e) z^(16-c-e)(1-z)^(c+e)
          * min(1-q^-rho,
                q^-max(0,c-min(nu_c,max(0,4-l+c))),
                q^(max(0,l-c-4)-l) (q^min(rho,4-e)-1)),
    j_l = min(t^l-k_l, (1-q^-rho) z^max(0,13-l), R_l),
    Jbar = min(product_i(k_(l_i)+j_(l_i))-K, d V).

These are the formulas in `hd_common_refresh_packets.py`. The return event
requires Q=CX != 0 and has output (I+AC)X. The three terms in R_l bound
its probability, the rank of active output-zero constraints, and the number
of nonzero syndromes consistent with inactive output-zero constraints.
All rank statements use fixed matrices over F0. In particular, the exact
principal-nullity table in `hd_common_refresh_principal.py` applies unchanged.
Its implementation performs F0 arithmetic despite calling a GF65536 helper.

The lower-right entry uses B V, without subtracting Jbar/d. Subtracting an
upper bound on the returning mass would reverse the required inequality.
The density B applies only to transitions into a nonzero state.

For stable numerical evaluation, expand the local difference positively:

    t^l-k_l = sum_(c=0)^(l-1) binom(l,c) z^(l-c)(1-z)^c q^-c
              * (1-q^-min(4,l-c)).

Then form T-K by a positive telescoping product. Directly subtracting
rounded logarithms can erase this positive mass when z is much smaller
than 1/q. The new FLOAT screen uses the positive expressions; a directed
replay still needs rigorous rounding of every operation.

For an empty packet, retain the existing improvement

    M = diag(1, min(B V, L0(q,z))).

Here L0(q,z)=max_(y != 0) E_R[z^wt(A R^-1 y)]. It is bounded by
`hd_grouped_quiet.quiet_cap(z,q=q,state=16,width=64,h=4)`.
The inverse of H16 is totally nonsingular by complementary minors.
For inverse-refresh input weight a and output weight b, put k=a+b-16.
The exact transition probability is zero if k<=0, and otherwise equals

    binom(16,b) (q-1)/(q-1)^a
        * sum_(j=0)^(k-1) (-1)^j binom(a+b-1,j) q^(k-1-j).

Conditional on b, the support is uniform among the b-element subsets.
Each active four-coordinate state block emits at least 13 symbols.
Averaging z^(13 times the number of active blocks) gives the implemented
L0 bound. This proof uses exact MDS support counts, not uniform refresh.

## Reusing the route proof

The cyclic support counts, feasibility tests, regional multinomial
denominators, swap probabilities, coefficient degrees, and support boxes
are independent of q. Keep the existing alternating route and its exact
independent swaps. For 0<alpha<=1, fix route choices, take entrywise
fractional powers of the local nonnegative operators, then average over
the swaps. Moving the averaging inside the fractional power is a different
bound and must not be imported silently.

For a support class C with active-group weight r and forbidden output
weight L, the same proof has the form

    Pr[exists a bad message with support in C]
      <= (beta_T(q)^r / ((q-1) z^L))^alpha
         * sum_(U in C) e0 E_route[product_t M_(U,t)^circ(alpha)] 1.

The expectation includes the independent pair permutations and swaps.
The projective divisor is q-1 and is inside the fractional power.
Every nonzero F-scalar multiple has identical support and output weight.
Support multiplicities stay outside the power. Adjacent matching statistics
are not independent; retain the existing joint counts or common-potential
argument. State remains continuous across the 256 regional operators.

For 35%, set L=734003. A concrete route to a new certificate is to retain
the old disjoint support cover, recompute its local operators and all
coefficient bounds at q=2^128, and aggregate the new ceilings. Reoptimize
marks, fractional exponents, and state potentials wherever replay fails.
The cover geometry can be reused; its GF256 numerical ceilings cannot.
The required stopping condition is complete coverage and a directed union
ceiling at most 2^-40. The same procedure can then raise L along the
40%, 45%, 47%, 48%, and 49% ladder.

Three implementation issues remain material. First, field order is
hardcoded in existing packet wrappers, density functions, prefixes, caches,
construction identifiers, and receipt validators. Second, q^16=2^2048
overflows binary64; positive log-domain screens and exact/directly rounded
replays must avoid large intermediate powers and cancellation. Third, a
larger q does not prove that all support classes pass the same witnesses.
No monotonicity theorem for the complete optimized comparator has been
established here. New complete coverage is still required.

The intended guarantee concerns ideal independent setup over F and symbol
Hamming distance. It does not certify a concrete pseudorandom seed,
bit-Hamming distance, or an implementation that shares random coefficients.
