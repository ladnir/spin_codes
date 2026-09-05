# The original M22 distance target is closed

Updated: 2026-09-05.

The fixed BCH-derived [256,128] outer meets the finite SPIN distance target
under ideal RandomStepConv-M22. The computer-assisted bound covers 2^20
message bits, rate 1/2, and relative distance greater than 0.1. Its exact
failure upper bound U satisfies

\[
U < \frac{61}{100}\,2^{-40}<2^{-40}.
\]

The diagnostic value of -log2(U) is 40.72497637925699. This is an upper
bound on setup failure, not an estimate of the actual failure probability.
No statistical shell assumption or variance estimate enters this result.
The inner memory remains 22; the earlier M25 fallback is unnecessary here.

The final receipt is
`generated/bch256_m22_unconditional_closure_audit.json`.
The exact joint LP witness is in `generated/shift_rank_oa29_joint`.
Older notes and frozen receipts retain their historical, weaker conclusions.

## Fixed outer and random setup

Work over F256 with modulus 0x14D and primitive element alpha=2. Let P and Q
be the even extensions of the binary primitive narrow-sense BCH codes of
length 255 and designed distances 37 and 39. Their dimensions are 131 and
123, and Q is contained in P. Index extended coordinates by field elements.
For a binary word c, define

\[
p_{37}(c):=\sum_{a:c_a=1} a^{37}.
\]

Fix the five-dimensional binary subspace S={0,...,31}, using the polynomial
basis representation, and define

\[
C:=\{c\in P:p_{37}(c)\in S\}.
\]

On P, this syndrome map is surjective onto F256 and has kernel Q. Hence C
has dimension 128. It is even, contains the all-one word, and has minimum
distance at least 38. Coordinate scaling acts transitively on nonzero
syndromes. Thus, if q_w counts weight-w words of Q and h_w counts weight-w
words in any fixed nonzero Q-coset inside P, then

\[
A_w(P)=q_w+255h_w,\qquad A_w(C)=q_w+31h_w.
\tag{1}
\]

In particular, the argument also applies to every other five-dimensional
choice of S. No complete spectrum for C is assumed or computed.

Fix a linear isomorphism G:F2^128 -> C. Split a message of K=2^20 bits into
8192 rows and apply G to each row. Setup samples the following independently:

1. A uniform permutation of the 256 coordinates of each encoded row.
2. A uniform permutation within each of the 256 regions obtained by transposing
   the array. Each region has length 8192.
3. An independent uniform binary 23-by-23 matrix M_i at each of N=2^21 positions.

Serialize the permuted regions in region-major order. Start the inner encoder
with state a_1=0 in F2^22. At position i, apply M_i to the state and input bit;
interpret the result as the next 22-bit state and one output bit. Discard
the final state. Write E_theta for this linear encoder, where theta denotes
all setup randomness. Every message uses the same theta.

## Finite distance theorem

The theorem uses two published exact endpoint spectra: Table 7's k=71 column
and Table 10's k=187 column in Fujiwara and Kusaka, IEICE 2021,
DOI 10.1587/transfun.2020EAP1119. They describe the dimension-71 extended BCH
code and the dual of the dimension-187 extended BCH code. These are explicit
mathematical inputs. The code constructions, containments, derived spectra,
rank witnesses, and numerical inequalities are checked by the retained programs.

**Theorem (computer-assisted, ideal RandomStepConv-M22).** For the fixed C
and random setup above, with D=209716,

\[
\Pr_\theta\left[
  \exists x\in\mathbb F_2^{K}\setminus\{0\}:
  \operatorname{wt}(E_\theta(x))\le D
\right]
\le U < \frac{61}{100}\,2^{-40}.
\tag{2}
\]

The exact rational U is the `full_first_moment_upper` field of the final
receipt. Outside the event in (2), E_theta is injective and its image has
minimum distance at least 209717. Its rate is 1/2 and relative distance
exceeds 0.1. Injectivity is part of this probabilistic conclusion; it is
not asserted for every possible setup.

## Proof: only the required weighted spectrum functional

Let Z_theta count nonzero messages whose encoded weight is at most D. For
q=1,...,8192, let F_q be the expected count with exactly q nonzero message
rows. The expectation is over theta. The union bound gives

\[
\Pr[Z_\theta\ge1]\le\mathbb E[Z_\theta]
=\sum_{q=1}^{8192}F_q.
\tag{3}
\]

No independence between different messages is needed. In particular, neither
a variance estimate nor a randomized certification guarantee is an input.

For one nonzero row, let p_w be the failure probability for a fixed outer
word of weight w, averaged over its route and inner setup. The route makes
its occupied regions a uniform w-subset. Its one position within each
occupied region is independently uniform. The transfer certificates give
nonnegative rational numbers c_w with 8192 p_w <= c_w. They cover every
possible nonzero weight: even weights 38 through 218, and weight 256.
Therefore

\[
F_1\le\sum_w c_w A_w(C).
\tag{4}
\]

For weights 38, 40, and 42, the sharper coefficients come from the cellwise
active-time comparison in [M22_ONE_SHELL_REDUCTION.md](M22_ONE_SHELL_REDUCTION.md).
All other coefficients come from the full-Arb Chernoff transfer certificate.
Both use the same M22 marginal law. For a fixed message, an active state-input
vector maps to an independent fair output bit and a uniform next state.
The cell calculation brackets active time, then bounds the conditional
Binomial(active time,1/2) tail. Its positive recurrence uses directed rounding;
Arb encloses the region coefficients and analytic tail factors.

The remaining task is to upper-bound (4) for the actual fixed outer. Let
T={38,40,42} and define b_w=c_w+c_(256-w) for w in T. Complement symmetry and
(1) make the contribution of these six shells

\[
\sum_{w\in T}b_w(q_w+31h_w).
\tag{5}
\]

We bound (5) jointly. Maximizing three separate shell counts would permit
incompatible worst cases and lose the required margin.

### Why the spectrum satisfies the joint LP

The final model has 196 nonnegative variables and 1163 normalized rows.
Its four half-spectra are q, h, the spectrum r of the hull P intersect P-perp,
and the common nonzero coset spectrum s of the nested hull quotient.
The following facts justify its constraints:

- Exact generator algebra gives the BCH containments and quotient identities.
  Published endpoint spectra and exact MacWilliams transforms bound both
  primal and dual spectra coefficientwise.
- Binary elimination computes the nested hulls, their dimensions 85 and 93,
  their doubly-even property, and the required intersections. Exact quadratic
  forms give the signed weight sums. Their Fourier support and sign counts
  give further linear inequalities. See
  [BCH_HULL_SPECTRUM_CONSTRAINTS.md](BCH_HULL_SPECTRUM_CONSTRAINTS.md).
- The hull of the dimension-71 anchor has dimension 69. Its spectrum is
  precisely the weight-0-mod-4 part of that anchor. Its dual contains the
  larger hull's dual. Complete affine orbit certificates prove that the
  latter has 0, 16592, and 0 words at weights 14, 16, and 18. Completeness
  follows by matching the containing code's exact shell totals; a partial
  search is not used as a complete enumeration.
- A finite Fourier case cover proves d(Q-perp)>=30 directly. Every case
  has a checked rank or BCH root-run witness. One case splits exhaustively
  according to whether F(23)=0. Parity and verified affine invariance transfer
  the bound to the extension. Thus Q and every coset are orthogonal arrays
  of strength at least 29. See
  [BCH_SHIFT_RANK_CERTIFICATES.md](BCH_SHIFT_RANK_CERTIFICATES.md).

The last item also justifies all inherited strength-21 equations. Neither
the formerly inaccessible SchaubPlus table nor an external dual-distance
table remains a required mathematical input. The published endpoint spectra
are still used. Search-based shell lower bounds and statistical shell caps
are not assumed.

### Exact joint witness and the rest of the sum

Round each b_w upward to a dyadic number bar_b_w with 192 fractional bits.
The checker verifies exactly that

\[
b_w\le\bar b_w<b_w+2^{-192}.
\]

After exact variable scaling and objective normalization, the retained dual
prices have the correct row signs and dominate every objective coefficient.
Their exact rational objective value therefore upper-bounds (5). Let B denote
that value after undoing normalization. The checker also verifies primal
feasibility and primal-dual equality, although dual feasibility alone suffices
for this upper bound. Solver status and binary64 tolerances are not trusted.

For every other shell, use the deterministic baseline cap a_w. These caps
come from the exact BCH LP, constant-weight packing, or exact Johnson-scheme
LP witnesses. Set

\[
R_1:=\sum_{w\notin T\cup(256-T)}c_w a_w.
\]

This includes the unique all-one word. Consequently F_1 <= B+R_1.

Let V_2, V_3, and V_tail be the retained M22 bounds for occupations two,
three, and four through 8192. Each deterministic baseline cap is at most
r=1+2^-20 times its old comparison cap. The same comparison holds for
cumulative caps used by the polynomial-square tail envelopes. Positivity
and degree-q homogeneity therefore cost at most r^q. Since q<=8192,

\[
r^q\le\sum_{j\ge0}(8192\cdot2^{-20})^j=128/127.
\]

Hence

\[
\sum_{q=2}^{8192}F_q
\le R_{\ge2}:=\frac{128}{127}(V_2+V_3+V_{\rm tail}).
\]

The prefix-domination, marginal monotonicity, and positive transfer arguments
are given in [RANDOM_MODEL_FULL_DISTANCE_BOUND.md](RANDOM_MODEL_FULL_DISTANCE_BOUND.md).
The envelope-inflation argument is given in
[RANDOM_INNER_M25_CLOSURE.md](RANDOM_INNER_M25_CLOSURE.md), under “Reusing all
higher occupations.” Here every transfer already uses memory 22, so no
change-of-memory argument is required. The older statistical assumptions
in the conditional note are not imported into this theorem.

Combining these inequalities proves (2), with U=B+R_1+R_ge2. For orientation,
the contributions as fractions of the 2^-40 budget are:

| Contribution | Fraction of target, diagnostic |
|---|---:|
| Joint paired weights 38, 40, and 42: B | 0.5154835704434867 |
| All remaining shells and higher occupations | 0.08952337965156525 |
| Full upper bound U | 0.6050069500950520 |

The comparison U < 61/(100*2^40) uses exact rational arithmetic, not these
displayed decimals. This proves the needed weighted inequality without
proving the earlier sufficient condition A38<=10^13 or finding the spectrum.

## Verification record and reproduction

The final audit ran successfully and then reproduced its receipt unchanged.
It independently recomputed all 92 original M22 transfer coefficients using
256-bit Arb, a normalized support recurrence, and binary polynomial powering.
It reconstructed the joint model, checked its rational witness, verified
78 JSON evidence files and their recorded hashes, and reaggregated every
occupation exactly. The separate cell-DP replay reproduced every saved array
bit-for-bit and passed the 6400-path toy comparison.

The retained higher-occupation evidence includes 83 independent full-Arb
region recurrences, 8106 direct full-length matrix powers, and 5397 exact
prefix comparisons. The final audit checks these receipts and their hashes;
it does not claim to rerun those long computations. It does recompute the
Q2/Q3 aggregates and the sum of all 8189 tail bounds.

The inherited OA21 and signed-weight model audits also replayed successfully.
The endpoint source pages were visually rechecked, and all retained entries
were matched programmatically. These are checks against the published
spectra, not new enumeration of the dimension-71 code.

Run the following sequentially from the bundle directory. Do not use Python
optimization flags, which disable assertions. No command below runs a solver
or overwrites a frozen receipt.

```text
python -B code/certify_bch_shift_rank.py --verify
python -B code/certify_bch_shift_rank_split.py --verify
python -B code/certify_bch_oa21_probe.py --verify
python -B code/audit_bch_mod4_cap.py --verify
python -B code/check_bch_hull_fourier.py --verify
python -B code/certify_bch_q1_activity_cells.py --verify
python -B code/verify_bch_m22_closure.py
```

The final command itself reconstructs the nested hull and low-shell models
and verifies the joint LP. Its `source_sha256` manifest preserves the evidence
boundary. This is a computer-assisted mathematical result, not a proof
formalized in a proof assistant. Correctness of the published spectra,
explicit reductions, checker implementations, and arithmetic libraries remains
the ordinary trust boundary. There is no statistical false-accept allocation
to add to (2).

## Application scope and next work

The finite event and route match the general first-moment framework of
Theorem 3.1 and the block-permute/transpose/region-permute route in Section 6.4
of the audited SPIN draft. Its SHA-256 remains
`86493c53e04119af1ebfd7252f5c3bfffa5d786c1926d7f7fc810bc02a98211b`.

This theorem concerns the ideal refreshed-state inner defined above. It does
not prove the draft's equation (19) random convolution, RM2Sub, a PRG
replacement, 11% distance, an asymptotic family theorem, or an end-to-end
proof-system security theorem. No claim of measured performance follows.
Those are separate questions, not outstanding assumptions within (2).

The next recommended task is independent review and paper integration of
this finite theorem. Transfer to the chosen practical inner should remain a
separate, explicit obligation. The original ideal-M22 closure goal is complete.
