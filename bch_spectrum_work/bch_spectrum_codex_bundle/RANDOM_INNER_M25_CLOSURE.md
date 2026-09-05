# A deterministic distance certificate with 25 state bits

Updated: 2026-09-04.

The same BCH-derived [256,128] outer meets the 40-bit distance target under
ideal RandomStepConv with **25 state bits**. This result uses deterministic
spectrum bounds, not the three statistical shell assumptions.
The complete outward bound has approximately 41.90289 bits.

This changes the inner memory from 22 to 25. The original M22 result remains
conditional. This note does not substitute M25 for that target without an
explicit parameter choice, and it does not claim a result for RM2Sub.

## Encoder and theorem

Let P and Q be the even extensions of the primitive narrow-sense binary BCH
codes of length 255 and designed distances 37 and 39. Work over GF(256) with
modulus 0x14d. Index the extended coordinates by field elements, and let
p_37(c) be the sum of the 37th powers of the support coordinates of c.
Fix the five-dimensional binary subspace S={0,...,31} and set

\[
C:=\{c\in P:p_{37}(c)\in S\}.
\]

The map p_37 on P has kernel Q and an eight-dimensional image. Thus C has
dimension 128. It is even, contains the all-one word, and has minimum distance
at least 38. The nonzero quotient cosets have equal spectra, so the bounds
below also apply to any five-dimensional choice of S.

Fix a linear isomorphism G from F_2^128 to C. Apply G to each of 8192 message
rows. Setup independently samples a uniform coordinate permutation for every
row and a uniform permutation for every one of the 256 regions. The route
permutes each row, transposes the array, permutes each region, and serializes
the result in region-major order. The serialized length is N=2^21.

The inner starts in the zero state of F_2^25. At each position i, setup samples
an independent uniform binary 26-by-26 matrix. It maps the current state and
one input bit to the next 25-bit state and one output bit. Discard the final
state. All matrices and permutations are independent and are fixed once for
every message. Write E_theta for the resulting linear encoder.

**Theorem (computer-assisted, ideal setup).** For this fixed C and this setup,

\[
\Pr_\theta\bigl[\exists x\in\mathbb F_2^{2^{20}}\setminus\{0\}:
  \operatorname{wt}(E_\theta(x))\le209716\bigr]
\le U_{25}<2^{-41}<2^{-40}.
\]

The exact rational U_25 is stored in
`generated/bch256_m25_unconditional_distance.json`.
Outside the bad event, E_theta is injective and has minimum distance at least
209717, hence relative distance greater than 0.1 and rate 1/2.
Injectivity is part of this probabilistic conclusion, not a property of every
possible inner setup. No variance estimate is needed.

The theorem uses the published exact endpoint spectra described below.
It makes no unproved assumption about the three low BCH shells. As usual for
this computer-assisted argument, correctness of the arithmetic implementation
and the explicit mathematical reductions is part of the verification boundary.

## Deterministic BCH bounds without search-based lower bounds

Let q_w=A_w(Q), and let h_w count weight-w words in a fixed nonzero Q-coset
inside P. Coordinate scaling acts transitively on nonzero quotient syndromes,
so

\[
A_w(P)=q_w+255h_w,\qquad A_w(C)=q_w+31h_w.
\]

Both spectra are complement-symmetric. Half-spectra therefore suffice.
The code L of dimension 71 and designed distance 59 is contained in Q.
The code U of dimension 187 and designed distance 19 contains P.
Direct generator divisibility checks establish L <= Q <= P <= U.

Fujiwara and Kusaka's IEICE 2021 paper supplies the complete spectrum of L
in Table 7 and of U's dual in Table 10. Every retained entry was matched
against the source PDF. Exact MacWilliams transforms recover the other two
spectra. The dual of L has minimum distance 16, so Q and every Q-coset are
orthogonal arrays of strength at least 15.

These facts give the following linear constraints on q_w and h_w:

- nonnegativity, the unique zero word, and the minimum-distance support;
- degree-at-most-15 moment identities for Q and its cosets;
- coefficientwise containment bounds from L <= Q <= P <= U;
- the reversed dual containments, expressed by the MacWilliams transform.

The new checker independently reconstructs the Krawtchouk coefficients by a
three-term recurrence and validates all 396 retained LP rows. It reuses exact
dual witnesses but replaces every orbit-search lower bound by zero. This
weakening preserves dual feasibility and avoids any need to trust a search
enumeration. It also omits the earlier lattice-rounding improvement at weight 38.

The resulting integer caps are:

| Weight | Deterministic cap |
|---:|---:|
| 38 | 773397230534890 |
| 40 | 4210950731605482 |
| 42 | 38281539412412374 |
| 44 | 1144487705125679705 |
| 46 | 7169795974927066729 |
| 48 | 56162197984147928451 |
| 50 | 645473488311357207401 |

These caps are much larger than the statistically tested caps. M25 tolerates
them. At every listed weight, the new cap is at most r=1+2^-20 times the
corresponding cap in the earlier occupation certificates.

For other even weights w<=128, use the existing Johnson certificates at
52 and 54 and the bound

\[
A_w(C)\le\min\left\{2^{128},
 \left\lfloor\binom{256}{w-18}/\binom{w}{w-18}\right\rfloor\right\}.
\]

The packing bound follows because distinct weight-w words cannot share
w-18 support coordinates: such a pair would have distance at most 36.
Complement symmetry supplies weights above 128, and the all-one shell has
one word. No statistical test outcome is an input to these inequalities.

## Occupation one

For a positive rational s, put b=(1+exp(-s))/2 and epsilon=2^-25. Inactive and
active states give the two moment matrices

\[
Z=\begin{pmatrix}1&0\\\epsilon b&(1-\epsilon)b\end{pmatrix},
\qquad
T=\begin{pmatrix}\epsilon b&(1-\epsilon)b\\
                  \epsilon b&(1-\epsilon)b\end{pmatrix}.
\]

Let R_0=Z^8192 and let R_1 be the average length-8192 product with one T and
8191 copies of Z. For outer weight w, let M_w be the average product of 256
region matrices with w copies of R_1 and 256-w copies of R_0. The initial
row vector is (1,0); the terminal column is (1,1)^T. Chernoff's inequality gives

\[
F_1\le\sum_{w>0} A_w(C)\,8192\,e^{209716s_w}
       (1,0)M_w(s_w)(1,1)^\mathsf T.
\]

The witness s_w may vary with w. Binary64 screening chooses witnesses, which
the certificate then treats as exact dyadic rationals. The primary calculation
uses 192-bit Arb, unnormalized support counts, and exact binomial denominators.
It includes all 92 possible nonzero shells, including weight 256.

The independent calculation uses 256-bit Arb, a normalized support recurrence,
and binary powering of a degree-one matrix polynomial for the region law.
All 92 coefficients lie below the stored outward bounds. Their aggregation
uses exact rational numbers.

## Reusing all higher occupations

For a fixed input sequence, represent the marginal inner law using independent
fair bits B_i and uniform real variables V_i. Starting from a_0=0, put

\[
u_i=a_i\lor x_i,\quad y_i=u_iB_i,\quad
a_{i+1}=u_i\,\mathbf1\{V_i\ge2^{-m}\}.
\]

A uniform linear map sends each fixed nonzero state-input vector to a uniform
state-output vector, so this is the correct law for memory m. Coupling two
memories with the same B_i and V_i preserves a_i and y_i in increasing order
as m increases. Thus each fixed message's lower-tail probability and its moment
at exp(-s) decrease as memory increases. This is a marginal coupling, not
pointwise monotonicity of one encoder on all messages.

Consequently the old M22 transfer bounds remain upper bounds at M25. There is
one further allowance: we weakened some deterministic spectrum caps above.
The occupation-two and occupation-three bounds are positive homogeneous
expressions of degree q in their row multiplicity bounds, so enlarging each
row bound by at most r costs at most r^q.

The same factor works for the occupation tail. Every old cumulative cap was
the minimum of a cumulative shell bound and independent polynomial or total-
mass bounds, possibly propagated from a later prefix. The true cumulative
count is at most r times each such retained cap. Hence the old reference factor
Gamma can be replaced by r Gamma, costing r^q at occupation q. No reference
row is conditioned to be nonzero.

For q<=8192, an elementary binomial estimate gives

\[
r^q\le(1+2^{-20})^{8192}
\le\sum_{j\ge0}(8192\cdot2^{-20})^j=128/127.
\]

Write V_2,V_3,V_tail for the earlier certified M22 bounds computed from the
old deterministic comparison envelopes. Although those comparison envelopes
included the search lower bounds, their numerical values can be reused
without assuming those lower bounds: the factor above accounts for the
weakened, independently justified envelopes. The complete bound is

\[
U_{25}=U_{1,25}+\frac{128}{127}(V_2+V_3+V_{\rm tail})<2^{-41}.
\]

This covers every occupation from 1 through 8192. Markov's inequality applied
to the count of bad nonzero messages proves the theorem. Different messages
need not have independent outputs.

## Relation to the current SPIN draft

The audited draft is
`C:/Users/peter/.codex/worktrees/8c0d/permute_conv/output/pdf/spin_codes_draft.pdf`,
SHA-256 `86493c53e04119af1ebfd7252f5c3bfffa5d786c1926d7f7fc810bc02a98211b`.

The distance event matches the general first-moment framework in Theorem 3.1:
noninjectivity or an image word of weight at most the cutoff. The fixed outer
also satisfies the independence requirement. The route is exactly the
block-permute/transpose/region-permute route of Section 6.4.

However, this is a separate finite ideal-inner specialization, not a proof of
the draft's Random SPIN theorem or Structured SPIN theorem. The draft's random
convolution, equation (19), stores trailing output bits and is invertible for
every setup. RandomStepConv refreshes the entire state with a random linear
map and need not be invertible. Their transfer laws differ. The current
theorem targets 10% distance; the draft's main asymptotic theorems target 11%
or 11.002%. Those claims must not be conflated.

For the general framework, classes may retain the row-weight profile. Our
calculation sums these classes and then groups them by occupation. It does
not assume that occupation alone determines the conditional failure law.
The full result establishes the finite distance property, not an asymptotic
family theorem, a benchmark, or a proof-system security theorem.

## Verification and next decision

Run these commands from the bundle directory:

    python -B code/audit_bch_closure_envelope.py --verify
    python -B code/certify_bch_m25_closure.py --verify
    python -B code/verify_bch_m25_closure.py

The first re-derives the BCH LP constraints and exact dual bounds. The second
recomputes all 92 Q1 transfers independently. The third checks dependencies,
the higher-occupation comparisons, and the complete rational sum. The existing
certificates retain independent arithmetic checks for all higher occupations.

The next choice is whether to adopt M25 as the random-model endpoint or retain
M22 as a separate tightening goal. No performance cost has been measured for
this change. No PRG replacement or transfer to the draft's other inners follows
from this ideal-setup theorem.

External spectrum source: Toru Fujiwara and Takuya Kusaka, *The Weight
Distributions of the (256,k) Extended Binary Primitive BCH Codes with k<=71
and k>=187*, IEICE 2021, DOI 10.1587/transfun.2020EAP1119; local file
`paper/e104-a_9_1321.pdf`, Tables 7 and 10, printed pages 1326-1327.
