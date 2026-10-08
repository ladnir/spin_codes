# A complete K18 certificate with six GF(16) state symbols

The regional construction has a complete certificate at K=262144,
N=524288, rate 1/2, and minimum symbol distance at least 52429.
Its setup-failure probability is at most

    8.272540769181098e-14 < 2^-43.

This rounded display is an upper bound. The certificate uses the full
unrounded value in its receipt. Its approximate margin is 43.4587 bits.
The probability is over the ideal independent setup distribution below;
it is not a per-seed certificate or an empirical distance measurement.

This profile changes the state dimension from four to six. The earlier
K16/state4 certificate remains separate and unchanged. Both profiles use
the same outer code and the same regional routing rule.

## Construction

Work over F=GF(16). Let L=32768. Encode each of L message rows with
GRS[16,8]. Place coordinate j of every outer row in region j. Independently
sample a uniform permutation of the L individual symbols in each region.
Independently multiply every routed symbol by an element sampled uniformly
from F minus zero. There are no intact packets of outer rows.

Process each region in sixteen-symbol steps. The state continues across
all region boundaries. Label the sixteen lanes by j in F and define

    A_j=(1,j,j^2,j^3,j^4,j^5),       C=A^T.

Initialize Q=0 in F^6. At every step, sample an independent uniform
M in GL(6,F), then compute

    R=M Q,       Y=X+A R,       Q'=R+C X.

Discard the final state. Setup samples every permutation, scalar, and
matrix once; the encoder then uses this fixed linear map for every message.
The certificate requires these independent distributions.

All coefficients lie in F. Scalar extension preserves minimum symbol
distance, so the same conclusion holds over GF(2^m) whenever four divides m.
This includes m=8,32,64,128 with a compatible subfield representation.
Expand an extension-field word in an F-basis: every component output has
support contained in the complete output. A nonzero component therefore
gives a base-field witness of no larger weight. The reverse inequality
follows because base-field words remain available after extension.

## Local transfer and its state dependence

The fixed-word calculation follows [GF16_PROOF.md](GF16_PROOF.md).
Conditional on the active input support, the independent diagonal makes
its nonzero symbols independent uniform elements of F minus zero.
Conditional on a nonzero incoming state, the fresh M makes R uniform over
F^6 minus zero. The transfer therefore retains the two categories zero
and nonzero, not every state value.

The parameterized local calculation substitutes s=6 in the MDS kernel
counts and the [16,s] output spectrum. Every count is an integer.
The columns of C have the MDS rank property. Also C A=0: every relevant
power sum over F has degree at most ten and vanishes. Consequently
B=I+A C satisfies B^2=I. The invertible-B bound on the joint return
moment remains valid.

The singleton return moment depends on the actual generator basis.
For a singleton at lane j, a return gives R=C_j X_j and output B X.
The implementation enumerates the sixteen possible lanes exactly and
counts the nonzero entries of B's column j. It does not reuse the
state4 singleton formula. All sixteen lane counts are independently
checked against explicit matrix multiplication.

    J_1(z)=(4 z^12+8 z^13+3 z^14+z^15)/(16*(16^6-1)).

Let T_w(z) be the resulting entrywise upper transfer matrix for a
uniform input support of weight w. Here 0<z<1 marks output weight.
The return entry uses an upper bound on the joint return moment; the
nonreturn entry subtracts a lower bound from the complete output moment.
No independence between return and output weight is assumed.

## Retaining more of the outer spectrum

The original positive source envelope loses useful constants over GF(16).
For one prescribed nonzero outer support of size w, write n_w for its
exact GRS[16,8] multiplicity. For 9<=w<=16,

    n_w=sum_(i=0)^w (-1)^i binom(w,i) 16^max(0,w-8-i).

Set n_w=0 for w<9. For any rational rho>0, define the exact rational

    eta(rho)=max_(9<=w<=16) n_w / (16^(w-8) rho^w).

Thus every actual nonzero row support has multiplicity at most
eta(rho) rho^w 16^(w-8). Assigning this positive envelope to impossible
support sizes only increases the resulting sum. A separate rho may be
chosen for each active-row class or interval.

For a region with h active positions, its transfer is bounded by

    R_h(z)=[u^h] (sum_(w=0)^16 binom(16,w) T_w(z) u^w)^(L/16)
                    / binom(L,h).

The normalization follows from the uniform regional permutation of
individual symbols. Fix a set of a active outer rows and define

    M_a(z,rho)=sum_(h=0)^a binom(a,h) (16 rho)^h R_h(z).

Summing all row-support patterns gives the class contribution

    B_a <= binom(L,a)/15 * eta(rho)^a * 16^(-8a) * z^(-52428)
                  * [1,0] M_a(z,rho)^16 [1,1]^T.

The factor 15 removes scalar multiples of the complete nonzero source
word. It is not applied independently to outer rows. For a=1,2,3,4,
the certificate instead retains every row's exact support size and n_w.

## Covering all active-row classes

For a=1 through 223, the certificate computes the regional coefficients
directly. Classes 1 through 4 use exact outer counts. The remaining
classes use the positive rho/eta envelope above.

For the large middle, define

    H(x,z,rho)=sum_(w=0)^16 binom(16,w) T_w(z)
                                      (16 rho x)^w (1+x)^(16-w).

For positive rational x, the coefficient identity and positivity give

    B_a <= binom(L,a)^(-15) * eta(rho)^a * 16^(-8a) /15
               * z^(-52428) * x^(-16a)
               * [1,0] H(x,z,rho)^L [1,1]^T.

For fixed x,z,rho, the right side is log-convex in a. The binomial
coefficient is log-concave; every other a-dependent factor is geometric.
The larger endpoint bound therefore covers an entire integer interval.
Eight rational witnesses cover a=224 through 23991.

The dense tail uses only the fresh diagonal. Put f(z)=z+(1-z)/15 and
let P be the exact nonzero-word enumerator of GRS[16,8]. The complete
a-class contribution is at most

    binom(L,a) P(f(z))^a z^(-52428)/15.

At z=1/52, successive terms have decreasing ratio below one. An exact
integer comparison bounds their sum for a=23992 through 32768 by 2^-80.

| Active outer rows | Certified contribution, rounded upward here |
|---:|---:|
| 1 through 223 | 8.272540709781916e-14 |
| 224 through 23991 | 5.931646348493063e-22 |
| 23992 through 32768 | 2^-80 |
| Complete union | 8.272540769181098e-14 < 2^-43 |

The causal inner is injective: from Y and the previous Q, recover
X=Y+A M Q and then update Q. The outer is also injective. Thus the
bound concerns ordinary minimum distance of the complete encoder.

## Reproduction and remaining scope

The fixed rational witnesses are in `theory/gf16_size_witnesses.py`.
Replay uses exact integer counts and 80-digit directed Decimal arithmetic.
Floating searches selected witnesses but do not enter certification.
From the repository root:

```sh
python -B research/workstreams/char2_field_spin/theory/gf16_sizes.py --certify
python -B research/workstreams/char2_field_spin/theory/gf16_sizes.py --verify PATH_TO_RECEIPT
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_gf16_sizes.py -v
```

The verifier regenerates every class and interval endpoint. It checks
continuous coverage through all L active-row counts, the exact dense
comparison, displayed endpoint bounds, all contribution sums, and the
claimed margin. The checkpoint receipt is
`tmp/char2-field-spin/gf16-s6-full-k18-d10.json`; it is reproducible output,
not a dependency of the proof program.

The same envelope did not close K18 with four or five state symbols.
For K20 it remained negative in middle classes even with eight state
symbols. Those are unsuccessful upper-bound screens, not lower bounds
on code-failure probability. They do not show that those codes have
poor distance. A subsequent exact return-spectrum screen through input
weight four did not close K18 with four states; the a=64 class remains
negative. See [the progress record](ONE_ROUTE_PROGRESS.md) for this result
and the next outer/inner tradeoffs. This new screen does not change the
six-state certificate above.
