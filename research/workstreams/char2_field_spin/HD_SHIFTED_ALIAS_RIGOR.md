# Directed replay of the shifted-pair alias bound

The positive alias bound has a directed replay for the class with 512
full pairs and 128 runs under the even matching. Its contribution to the
established support union is at most

\[
  2^{-756.596785624611}.
\]

A second replay covers the general profile \((512,512,400,400)\), with
singletons under both matchings, below \(2^{-208.246106966782}\). The final
section specifies that profile and its witness. Each receipt covers one
support class; neither establishes a whole-code distance claim.

The implementation is
[`theory/hd_shifted_alias_certificate.py`](theory/hd_shifted_alias_certificate.py).
The numerical proof begins with exact rational packet caps. It does not
treat the floating coefficients from the screen as exact inputs.

## Exact inputs and regional operators

Fix a support that consists of \(a\) full even pairs in the cycle of 4096
pairs. Let \(c\) be its number of contracted runs. Under the odd matching,
the four pair types have counts

\[
  (n_{00},n_{10},n_{01},n_{11})=(4096-a-c,c,c,a-c).
\]

The exact number of indexed supports in this class is

\[
  C_{a,c}=\frac{4096}{c}
  \binom{a-1}{c-1}\binom{4096-a-1}{c-1}.
\]

Let \(z\in(0,1)\) and \(\alpha\in(0,1]\) be rational witnesses. For each
sorted occupancy \(\ell\in\{0,\ldots,16\}^4\), let \(E_\ell(z)\) be the
nonnegative common-refresh matrix returned by the exact packet formulas.
The implementation includes the GF256 diagonal density factor
\((256^{16}-1)/255^{16}\) in its birth and live entries. For zero occupancy,
the live entry also uses the exact inverse quiet cap. These are the same
rational ingredients used by the grouped-subfield certificate.

Each entry is rounded upward to a dyadic rational with 60 mantissa bits.
Write \(\widehat E_\ell\) for the resulting matrix. Entrywise powers of
\(\widehat E_\ell\) and ordinary products of these matrices remain sound
majorants because all entries are nonnegative.

For an odd half-step, write \(h=(h_{10},h_{01},h_{11})\), where each component
is nonnegative and \(|h|\le16\). Define

\[
  w(h)=\binom{16}{16-|h|,h_{10},h_{01},h_{11}},
  \qquad
  L(h)=(h_{10}+h_{11},h_{01}+h_{11}).
\]

Concatenating \(L(h)\) and \(L(k)\), then sorting their four coordinates,
defines \(\ell(h,k)\). The odd-step polynomial is the matrix polynomial

\[
  P(x,y,u)=\sum_{h,k}w(h)w(k)
  \widehat E_{\ell(h,k)}^{\circ\alpha}
  x^{h_{10}+k_{10}}y^{h_{01}+k_{01}}u^{h_{11}+k_{11}}.
\]

Here \(\circ\alpha\) denotes entrywise exponentiation. Arb encloses each
rational power before coefficient summation. The two singleton orientations
remain separate polynomial coordinates. Exchanging the halves preserves
the matrix, so the implementation sums unordered half pairs with the exact
multiplicity factor two when the halves differ.

The ideal odd regional operator is

\[
  R_o=\frac{[x^cy^cu^{a-c}]P^{128}}
  {\binom{4096}{4096-a-c,c,c,a-c}}.
\]

The even operator retains two-step composition before the entrywise alpha
power. For a full-pair half profile \(p\le q\), use occupancy \((p,p,q,q)\),
degree \(p+q\), and multiplicity
\(\binom{16}{p}\binom{16}{q}(2-\mathbf1_{p=q})\).
For every ordered pair of these profiles, multiply their matrices first.
Raise each product entry to \(\alpha\), attach the product multiplicity,
and sum into the corresponding coefficient of a polynomial \(H(x)\).
Then

\[
  R_e=\frac{[x^a]H^{64}}{\binom{4096}{a}}.
\]

Arb polynomial multiplication computes this coefficient. Truncating above
degree \(a\) preserves it because all polynomial degrees are nonnegative.
The rational tilt \(a/(4096-a)\) cancels through the coefficient normalizer.
The implementation takes upper endpoints to obtain entrywise majorants
before composing the regional operators.

## Positive alias extraction

Choose positive rational markers \(s,f\). Let \(D=\operatorname{diag}(v)\)
for an arbitrary positive dyadic vector \(v\). Choose a positive dyadic
\(\lambda\) at least every row sum of \(D^{-1}P(s,s,f)D\), using Arb upper
endpoints. Define

\[
  Q(x,y,u)=\lambda^{-1}D^{-1}P(sx,sy,fu)D.
\]

A floating Perron calculation proposes \(v\). No correctness claim depends
on that proposal being an exact eigenvector. Directed arithmetic determines
\(\lambda\) and all normalization factors. The normalized polynomial has
nonnegative coefficients and row sums at most one at \((1,1,1)\).

Fix a power-of-two grid size \(g\ge64\). Evaluate \(Q\) at triples of
\(g\)-th roots of unity, take the ordinary matrix power 128 at each triple,
and apply the normalized inverse DFT. At target \(t=(c,c,a-c)\), the result is

\[
  A_g(t)=\sum_{j\in\mathbb Z_{\ge0}^3:\ j\equiv t\pmod g}
       [x^{j_1}y^{j_2}u^{j_3}]Q^{128}.
\]

Every summand is entrywise nonnegative. Therefore \(A_g(t)\) majorizes the
target coefficient. The directed replay encloses \(A_g(t)\), takes the
upper endpoints, undoes \(D\) and \(\lambda\), and divides by
\(s^{2c}f^{a-c}\) and the exact multinomial normalizer. This gives a directed
entrywise majorant for \(R_o\). No bound on the mass of the extra aliases is
needed.

The numerical enclosure uses these operations:

1. Arb encloses the exact rational powers, sums, and normalizations.
2. Each real input and root component is converted to binary64 endpoints.
   Arb verifies containment after conversion.
3. Explicit radix-2 butterflies propagate rectangular complex intervals.
   Each separate real addition, subtraction, or multiplication expands its
   result outward with `nextafter`.
4. Seven interval matrix squarings enclose the 128th power at each grid point.
5. A directed sum extracts the chosen inverse-DFT coefficient. The generic
   path also includes directed phase factors for nonzero target residues.
6. Arb carries the final scaling, normalizers, regional product, and logarithm.

These intervals require IEEE-754 binary64 arithmetic with gradual underflow.
The implementation checks the array type and subnormal arithmetic, and
rejects nonfinite results. NumPy supplies the individual arithmetic ufuncs;
its FFT routines are not part of the proof. There is no numerical cushion.

For this witness, \(t=(128,128,384)\equiv(0,0,0)\pmod{64}\), so inverse
extraction is a directed average. The imaginary intervals contain zero.
That check diagnoses implementation failures; real upper endpoints provide
the mathematical bound.

## Targeted receipt

The exact witness is

```text
active_pairs = 512
contracted_runs = 128
z = 6022500333374007/12500000000000000
alpha = 4734086816355307/50000000000000000
singleton_marker = 49459/1000000
full_marker = 4669/25000
grid = 64
packet_mantissa_bits = 60
even_composition_steps = 2
field_order = 256
forbidden_weight = 734003
```

Let \(T=734003\), and use the established outer comparison measure and
regional composition from `HD_GROUPED_PACKET_PROOF.md`. The class bound is

\[
 C_{a,c}\,2^{\alpha(2a(2048-128\log_2 255)-\log_2 255-T\log_2 z)}
 \left((R_eR_o)^{128}_{00}+(R_eR_o)^{128}_{01}\right).
\]

The order \(R_eR_o\) follows the alternating schedule. The scalar prefactor
retains the outer density factor and projective normalization.

Local directed replays at 256 and 384 Arb bits both return `accepted=true`
for threshold \(-80\). Their log2 upper majorants begin

```text
-756.59678562461149706511757632065141863007952384322891983379735620333655
```

At 384 bits the reported final Arb radius is below \(4.32\cdot10^{-109}\).
This is the enclosure of the computed upper majorant, after the interval
FFT has already supplied its outward upper endpoints. It is not a claim
that the FFT error or aliasing loss has that size.

Both precisions give the following binary64 intervals for the normalized
real alias entries, in row-major order:

| Entry | Lower endpoint | Upper endpoint |
|---|---:|---:|
| 00 | 1.0607676633319353e-05 | 1.0607676633338726e-05 |
| 01 | 6.2025771508733215e-06 | 6.202577150883904e-06 |
| 10 | 1.0681670048671724e-05 | 1.0681670048689903e-05 |
| 11 | 6.15910921586584e-06 | 6.159109215876494e-06 |

The JSON receipt also includes an exact rational ceiling on its log2 bound:

```json
{
  "log2_upper_rational": "-756596785624611/1000000000000",
  "accepted": true,
  "directed": true,
  "full_under_even_matching_only": true,
  "whole_code_certificate": false
}
```

To replay from the repository root:

```powershell
python research/workstreams/char2_field_spin/theory/hd_shifted_alias_certificate.py
python research/workstreams/char2_field_spin/theory/hd_shifted_alias_certificate.py --precision 384
python research/workstreams/char2_field_spin/theory/hd_shifted_alias_certificate_tests.py
```

Eight tests check exact rational arithmetic enclosures, Arb input conversion,
FFT inversion, cyclic matrix polynomial aliases, inverse extraction,
orientation-preserving half multiplicities, rejection paths, and upward
packet rounding. The cyclic reference explicitly enumerates integer matrix
polynomial products on a small grid, including noncommuting coefficient
matrices and nonzero target residues. An independent read-only audit also
checked the coefficient construction, two-step normalization, and interval
scaling against this reference.

The remaining proof task is directed coverage across all feasible support
classes, followed by aggregation. The targeted receipt resolves the
numerical alias obligation for the displayed class; it does not resolve
those coverage obligations.

## Arbitrary oriented profiles

The separate implementation
[`theory/hd_shifted_general_certificate.py`](theory/hd_shifted_general_certificate.py)
extends the directed replay to supports that have singleton pairs under
both matchings. It retains the earlier full-even certificate and its
two-step even-region bound.

Write \(e,o\) for the active even and odd positions, and \(a,b\) for the
full-pair counts under the even and odd matchings. Put \(r=e+o\). In the
ordered coordinates \((00,10,01,11)\), the regional histograms are

\[
 H_0=(4096-r+a,e-a,o-a,a),\qquad
 H_1=(4096-r+b,o-b,e-b,b).
\]

The second histogram exchanges the parity labels because an odd pair
starts at an odd position. The replay uses three independent positive
rational markers for each region. Its coefficient array weights degree
\((i,j,k)\) by \(m_{10}^i m_{01}^j m_{11}^k\), with no equality assumption
between the singleton markers.

Both regional operators use the one-step polynomial \(P\) defined above,
raised to the ordinary matrix power 128. For histogram \(H\), extract its
three nonempty counts and divide by \(\binom{4096}{H}\). The implementation
also divides by the marker monomial when undoing the numerical tilt.
Monomials exceeding a target coordinate can be removed before the power:
their degrees cannot contribute to that target coefficient.

Let \(C(e,o,a,b)\) be the exact joint support count from
`hd_shifted_general_count.exact_joint_count`. With the same outer comparison
measure, the general class bound is

\[
 C(e,o,a,b)\,
 2^{\alpha(r(2048-128\log_2 255)-\log_2 255-T\log_2 z)}
 \left((R_0R_1)^{128}_{00}+(R_0R_1)^{128}_{01}\right).
\]

The implementation also accepts the exact-integer upper counts
`joint_count_upper` and `joint_run_count_upper`; the receipt names the
selected method. The regional product always retains the schedule order.
When both histograms and marker triples agree, one deterministic numerical
replay supplies both enclosures. This reuse does not couple the independently
sampled regional permutations or change the comparison measure.

For the profile \((e,o,a,b)=(512,512,400,400)\), both histograms equal
\((3472,112,112,400)\). The rational witness is

```text
z = 603269865746597/1250000000000000
alpha = 64134464781859/625000000000000
markers for each matching =
  (41044534/928270931, 12934738/292533347, 195385295/968824422)
support_count_method = exact
packet_mantissa_bits = 60
precision = 256
```

The exact count has log2 approximately 2194.415539248103. At grid 64, the
directed log2 upper bound is approximately \(-39.82744412205585\), and the
receipt correctly rejects threshold \(-40\).

At grid 128, the directed replay returns

```text
log2_upper_enclosure =
  [-208.2461069667823305256190214054107692055743028497024283233379505110439709
   +/- 3.78e-71]
```

Its exact rational ceiling and scope flags are

```json
{
  "profile": [512, 512, 400, 400],
  "grid": 128,
  "precision": 256,
  "log2_upper_rational": "-104123053483391/500000000000",
  "accepted": true,
  "threshold": -40,
  "full_under_either_matching": false,
  "directed": true,
  "whole_code_certificate": false
}
```

The normalized real alias intervals at grid 128 are

| Entry | Lower endpoint | Upper endpoint |
|---|---:|---:|
| 00 | 7.3491344233922064e-06 | 7.349134423434767e-06 |
| 01 | 4.377879091536198e-06 | 4.377879091560127e-06 |
| 10 | 7.519190872703577e-06 | 7.519190872744705e-06 |
| 11 | 4.296293272396932e-06 | 4.296293272420811e-06 |

All imaginary intervals contain zero. Their largest absolute endpoint is
below \(2.10\cdot10^{-17}\). As in the first receipt, the final Arb radius
encloses the computed majorant after interval endpoints have been taken.
It does not measure the loss from positive aliases.

The general API is `regional(z, alpha, histogram, markers, ...)` for one
operator and `receipt(profile, z, alpha, markers, ...)` for an ordered
pair of operators. The latter receives two marker triples.

```powershell
python research/workstreams/char2_field_spin/theory/hd_shifted_general_certificate.py --grid 128
python research/workstreams/char2_field_spin/theory/hd_shifted_general_certificate_tests.py
```

Five additional tests cover three distinct markers, compatibility with the
original equal-marker interface, target truncation, two different regional
histograms, and reuse of equal regional enclosures. An exact integer test
also checks the schedule order with noncommuting matrices. An independent
read-only audit checked the general orientation, normalization, truncation,
and reuse arguments and found no defect.

## Selective paired composition and certified neighborhoods

The further helper `hd_shifted_good_pair_certificate.py` composes two
low-impurity steps before the entrywise fractional power. Here a step has
low impurity when it contains at most one singleton pair. Write G for the
one-step polynomial of these patterns and D for the other patterns, after
the prior entrywise power. Write C for the good-good polynomial obtained
by composing the original matrices before taking that power. The paired
polynomial is C+GD+DG+DD. It has nonnegative coefficients and is bounded
coefficientwise by (G+D)^2. The count and decomposition are detailed in
[the structural progress note](HD_STRUCTURAL_ITERATION.md).

The directed helper builds C from 425^2 ordered exact-cap geometry pairs.
It retains both singleton orientations by expanding `(x10+x01)^k` with
exact binomial multiplicities. For a positive dyadic vector v, put
V=diag(v). Choose a positive binary64 dyadic rho whose square dominates
every row sum of `V^-1 (C+GD+DG+DD)(1,1,1) V`; Arb verifies this inequality.
Transform G and D with similarity V and scale rho, and C with similarity
V and scale rho^2. The interval expression C+GD+DG+DD then encloses the
paired polynomial with the common similarity and scale rho^2. Its 64th
matrix power requires the final factor rho^128. Degree 64 is folded into
residue zero at grid 64, never discarded.

For `(E,O,a,b)=(512,480,450,420)`, the directed point replay gives

    log2 contribution <= -394206727488049/1000000000000.

This replaces a positive one-step bound without changing the encoder.
The fixed rational weight witnesses are

    z = 4446635480191217/9007199254740992
    alpha = 1765869599810503/18014398509481984
    even markers = (19486566/827921929, 4949092/433906991, 151838032/701946107)
    odd markers = (16539401/719028543, 31644277/898709963, 23695721/116405276).

`hd_shifted_box_certificate.py` extends the replay to the inclusive box

    E in [510,514], O in [478,482], a in [448,452], b in [418,422].

All 625 profiles are feasible. The support-count envelope uses the exact
rational transfer-matrix markers

    (78107/24499979, 59183/54662624, 79848470422/92864511, 17740896315/66165169).

The scalar normalizer is bounded uniformly by convex log-factorial chords
in `hd_shifted_box_normalizer.py`. For each region, the full interval inverse
transform supplies coefficient upper endpoints for all residues that its
histogram box can request. Their entrywise maxima bound every profile's
unnormalized regional matrix. The product retains the order R0 R1. Multiplying
by the total support-count upper bound and the uniform scalar envelope
bounds the whole box union; the maximizing profiles need not coincide.

Both 256- and 384-bit replays give

    log2 box contribution <= -12092295652059/100000000000.

The point and box receipts set `directed=true` and
`whole_code_certificate=false`. Exact-reference tests cover the good-good
coefficient identity at alpha=1, separate marker powers, degree-64 folding,
noncommuting products, inverse-transform normalization, and wrapped target
intervals. Independent read-only audits found no gap in either replay or
in the combination of the count, scalar, and matrix envelopes.

Run from the theory directory:

```text
python -B hd_shifted_good_pair_certificate.py
python -B hd_shifted_box_certificate.py --radius 2
python -B hd_shifted_box_certificate.py --radius 2 --precision 384
```

The remaining obligation is a disjoint cover of every nonempty support
profile and a directed sum of all accepted box contributions. These local
receipts alone do not establish the requested whole-code margin.
