# A complete one-route certificate over GF(16)

The one-route construction now has a complete finite distance certificate.
At K=65536 and N=131072, the probability that its minimum symbol distance
is below 13108 is less than 2^-42. The bound concerns the setup ensemble,
not an empirical search for low-weight words. It covers every nonzero
source word and every active-outer-row count.

The decisive change was to use GF(16) coefficients and four state symbols.
The earlier GF(256) candidate used two state symbols. Both have 65536
possible state vectors, but the GF(16) outer has far fewer base-field
messages. The smaller coefficient field removes the specific large-subspace
obstruction to the previous word first moment. The certificate below,
not that heuristic comparison, establishes the new distance bound.

## Construction and probability space

Let F=GF(16), L=8192, and t=16. Independently encode each of L message
rows with the same GRS[16,8] code over F. Every row contributes one symbol
to each of sixteen regions. Independently sample a uniform permutation
of the L symbols in each region. No permutation preserves a packet of
outer rows. Independently multiply every routed symbol by a uniform
element of F\{0}.

Process the region-major stream in t-symbol steps, without resetting
the state at region boundaries. For j in F, define the row

    A_j=(1,j,j^2,j^3),       C=A^T.

The ordering of F determines the sixteen lane labels. Initialize Q=0
in F^4. At each step, independently sample M uniformly from GL(4,F)
and compute

    R=M Q,       Y=X+A R,       Q'=R+C X.

Discard the final state. Setup samples all permutations, diagonals, and
matrices once; the resulting encoder is linear and fixed for every
message. Any implementation that changes these distributions needs a
separate argument.

All generator entries lie in F. Extending scalars to GF(2^m), for any m
divisible by four, preserves minimum symbol distance. To see this,
expand an extension-field source word in an F-basis. A nonzero component
produces a base-field codeword whose support is contained in the complete
output support. Conversely, base-field codewords are still available
after extension. In particular, the result applies to m=8,32,64,128
with an explicit compatible subfield representation. It does not treat
an arbitrary extension-field state as a uniform element of F^4.

## A local transfer for one fixed source word

Fix a base-field source word before exposing setup. Conditional on its
input support in a step, the fresh diagonal makes its nonzero input
symbols independent uniform elements of F\{0}. Conditional on a nonzero
incoming state, the fresh M makes R uniform over F^4\{0}, independently
of those input symbols. We therefore retain only the two state categories
zero and nonzero for this fixed-word calculation.

Fix z in (0,1], and mark an output of weight v by z^v. Write q=16,
s=4, and

    f(z)=z+(1-z)/(q-1).

For input support size w, the probability that C X=0 is

    p0(w) = (q-1)^(-w)
            * sum_(i=0)^w (-1)^i binom(w,i) q^max(0,w-s-i).

This follows by inclusion-exclusion from the MDS column ranks of C;
the formula includes p0(0)=1. Let n_v be the number of nonzero words
of a [t,s] MDS code with one prescribed support of size v. Inclusion-
exclusion gives, for v>=1,

    n_v = sum_(i=0)^v (-1)^i binom(v,i) q^max(0,v-(t-s)-i).

Set n_0=0 because this enumerator excludes the zero word.

A nonzero A R has at most s-1 zero coordinates. Its unconditional
marked output moment, after adding X, is

    G_w(z) = 1/(q^s-1)
      * sum_(j=0)^(s-1) n_(t-j)
      * sum_a binom(w,a) binom(t-w,j-a)
               z^(t-w-j+2a) f(z)^(w-a),

where max(0,j-(t-w))<=a<=min(w,j). Here j counts zero coordinates of
A R, and a counts those lying in the active input support.

Let J_S(z) be the joint moment of the output and the event Q'=0 for
an input support S. Its unmarked return probability is

    p_ret=(1-p0(w))/(q^s-1).

The exact two-state transfer has rows and columns ordered zero, nonzero:

    [[z^w p0(w),       z^w (1-p0(w))],
     [J_S(z),          G_w(z)-J_S(z)]].

The certificate uses an entrywise upper envelope. It puts an upper
bound on J_S in the return entry and subtracts a lower bound from G_w
in the nonreturn entry. This does not assume that return and output
weight are independent.

### Bounds on the joint return moment

The MDS weight of A R gives

    p_ret z^t <= J_S(z) <= p_ret z^max(0,t-s+1-w).

For a singleton input, averaging over its uniformly selected lane gives
the sharper exact value

    J_1(z)=(13 z^15+z^16+2 z^14)/(16*(16^4-1)).

This identity uses the stated monomial A and C=A^T. It is not asserted
for an arbitrary generator basis of the same MDS code.

Dense inputs need a second upper bound. Since the lane set is all of F,
its power sums through degree six vanish. Hence C A=0 and
B=I+A C satisfies B^2=I. On a return, R=C X and Y=B X. For any output
zero set Z and input support S, invertibility gives

    rank(B_(Z,S)) >= max(0, |Z|-(t-|S|)).

Each independent pivot determines at most one nonzero input value.
Consequently, for w=|S|,

    Pr[(B X)_Z=0] <= (q-1)^(-max(0,|Z|-t+w)).

Expand the output marker with nonnegative coefficients:

    z^wt(Y) = sum_Z z^(t-|Z|) (1-z)^|Z| 1[Y_Z=0].

Dropping the condition C X!=0 yields the valid return upper bound

    J_S(z) <= 1/(q^s-1)
      * sum_(h=0)^t binom(t,h) z^(t-h) (1-z)^h
                          (q-1)^(-max(0,h-t+w)).

The certificate takes the smaller of this bound and the preceding MDS
bound. At w=t, the new bound is f(z)^t/(q^s-1). This avoids the constant
upper bound that made the earlier dense-step estimate ineffective.

## Regional normalization and complete source coverage

Let T_w(z) denote the local upper-envelope matrix averaged over uniform
input masks of weight w. The singleton averaging above is exact; all
other joint-return bounds hold for every mask individually. Define

    F(u,z)=sum_(w=0)^t binom(t,w) T_w(z) u^w.

For a region with h active positions, the entrywise transfer bound is

    R_h(z)=[u^h] F(u,z)^(L/t) / binom(L,h).

The denominator is regional, not global. It follows from the uniform
permutation of individual symbols in each region. Multiplying region
matrices retains the state across their boundaries.

Fix a set of a active outer rows. For a prescribed support of size
w>=9 in one GRS[16,8] row, its exact multiplicity is the MDS formula
above with t=16,s=8. For a=1,2,3,4, the certificate sums all binary
a-by-16 support patterns and uses these exact multiplicities. A
positive dynamic program retains each row's support size.

For larger a, use the upper multiplicity q^(w-8). Assigning this same
positive envelope to impossible row supports only increases the bound.
Define

    M_a(z)=sum_(h=0)^a binom(a,h) q^h R_h(z).

The complete a-class failure contribution is at most

    binom(L,a)/(q-1) * q^(-8a) * z^(-13107)
          * [1,0] M_a(z)^16 [1,1]^T.

The one factor q-1 removes scalar multiples of a complete source word.
It is not applied separately to outer rows. A separate rational z may
be chosen for each class or for each interval of classes.

### Large-middle classes without large coefficient arrays

For positive x, define

    H(x,z)=sum_(w=0)^t binom(t,w) T_w(z) (q x)^w (1+x)^(t-w).

The identity

    binom(a,h)/binom(L,h)=binom(L-h,a-h)/binom(L,a)

implies

    M_a(z)=[x^a] H(x,z)^(L/t) / binom(L,a).

All matrix-polynomial coefficients are nonnegative. Evaluating at a
positive rational x therefore gives the class bound

    binom(L,a)^(-15) * q^(-8a) / (q-1)
      * z^(-13107) * x^(-16a)
      * [1,0] H(x,z)^8192 [1,1]^T.

For fixed x,z, its a-dependent part is log-convex: the binomial
coefficient is log-concave, and the remaining factors are geometric.
Thus the larger endpoint bound covers every integer in an interval.
Eight fixed rational witness pairs cover a=128 through 5997.

### Dense tail

Condition on any incoming refreshed state. On every active input lane,
the independent diagonal gives an output moment at most f(z). Inactive
lanes contribute at most one. Thus a fixed source word of outer weight
W has total output moment at most f(z)^W, regardless of routing.

Let P(v) be the exact nonzero-word enumerator of the GRS[16,8] outer.
The a-class bound is therefore

    binom(L,a) P(f(z))^a z^(-13107)/(q-1).

Successive terms have decreasing ratio
(L-a)P(f(z))/(a+1). At z=1/52, an exact integer comparison bounds the
entire tail a=5998 through 8192 by 2^-64.

## Finite result and reproduction

All certification operations use rational witnesses, exact integer
counts, and 80-digit directed Decimal arithmetic. Positive polynomial
and matrix products round upward. Positive denominators round downward.
Floating-point searches selected the witnesses but do not participate
in their replay.

| Active outer rows | Method | Certified contribution, rounded upward here |
|---:|---|---:|
| 1 through 127 | Direct regional coefficients; exact outer counts for 1 through 4 | 1.763460e-13 |
| 128 through 5997 | Eight coefficient-Chernoff intervals | 3.740909e-22 |
| 5998 through 8192 | Exact geometric-tail comparison | 2^-64 |
| All nonzero words | Sum of all three blocks | 1.763460e-13 < 2^-42 |

The unrounded combined upper bound is approximately
1.7634597907074477e-13, corresponding to about 42.37 bits. The verified
integer claim is 42 bits. The final encoder is injective: given Y and
the previous state, recover X=Y+A M Q and then update Q. This recovers
every input step, and the outer is injective. Thus the certificate
concerns ordinary minimum distance, not only nonzero outputs after
quotienting a kernel.

From the repository root:

```sh
python -B research/workstreams/char2_field_spin/theory/gf16_certificate.py --certify-full
python -B research/workstreams/char2_field_spin/theory/gf16_certificate.py --verify PATH_TO_RECEIPT
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_gf16_certificate.py -v
```

The verifier checks the construction, all coverage boundaries, every
rational sparse witness, both endpoints of every middle interval, and
the exact dense-tail comparison. It rejects understated bounds and
overstated margins. Independent algebra and rounding checks are in
`theory/test_state_screen.py` and `theory/test_coupled_chernoff.py`.

The earlier GF(256) rank-based attempt is preserved separately in
[ONE_ROUTE_PROOF.md](ONE_ROUTE_PROOF.md). It is not used in this certificate.
