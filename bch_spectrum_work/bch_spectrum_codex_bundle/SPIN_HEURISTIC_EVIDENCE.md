# Heuristic evidence for the BCH-256 SPIN constituent

Application correction, 2026-09-04: the weight-38 target is one allocation in
a joint low-shell budget. Bounding that shell alone does not close the 40-bit
Q1 screen with the current other-shell LP caps. The directed budget audit
gives 36.8688 bits with the target cap and 36.9955 bits with its contribution
removed. This supersedes the later historical claims of weight-38-only closure.
See `generated/endpoint_application_budget_audit.json` and
`ENDPOINT_SAMPLING_CERTIFICATE.md` for the corrected boundary.

## Question and scope

The finite SPIN implementation needs a fast rate-half outer constituent.  The
candidate studied here is the extended BCH-derived code

\[
  C=[256,128,38].
\]

The proposed implementation reuses one fixed copy of \(C\) in every outer
row.  Independent permutations randomize the coordinates within each row and
the positions within each transposed region.  This note uses RandomStepConv-M22
as a temporary transfer model.  A separate computation will replace that model
with the final RM2Sub transfer.

This note does not claim a distance theorem for that implementation.  The
complete weight enumerator of \(C\) is unknown.  The objective is narrower: to
assess whether the authenticated spectra of smaller BCH constituents support
using \(C\) as a practical heuristic choice.

The Expanded Accumulant construction remains the proof-backed alternative.
Its known spectrum and variance support a rigorous selection argument.  The
reported implementation is approximately four milliseconds, or a factor of
1.4, slower than the BCH path.  BA-3 is not treated as a proof fallback here.
Its sampled constituent would require a variance or concentration calculation
that is not currently available.

## Comparable curve experiment

The diagnostic experiment fixes

\[
  M(k)=\lceil\log_2 k\rceil+2,
  \qquad N=2k,
  \qquad D=\lceil0.10N\rceil.
\]

A rate-half constituent of length \(B\) is repeated in

\[
  L=\frac{2k}{B}
\]

outer rows.  Each authenticated BCH spectrum is compared with the expected
spectrum of a full-rank random \([B,B/2]\) map.  The same map is reused in all
rows.  The reported quantity is the binary64 occupation-one margin
\(-\log_2 Q_1^{\mathrm{ub}}\).

Holding \(L\) fixed is important.  It compares block sizes without changing
the number of active-row locations or crossing the one-row routing boundary.

| rows \(L\) | block \(B\) | \(\log_2 k\) | BCH margin | matched-random margin | BCH advantage |
|---:|---:|---:|---:|---:|---:|
| 2 | 8 | 3 | -1.12 | -1.99 | +0.87 |
| 2 | 32 | 5 | 1.32 | 0.81 | +0.51 |
| 2 | 64 | 6 | 7.62 | 7.16 | +0.46 |
| 2 | 128 | 7 | 20.09 | 19.71 | +0.38 |
| 2 | 256 | 8 | unknown | 46.51 | unknown |
| 16 | 8 | 6 | 0.43 | -3.63 | +4.07 |
| 16 | 32 | 8 | 4.11 | 0.87 | +3.24 |
| 16 | 64 | 9 | 11.92 | 8.46 | +3.46 |
| 16 | 128 | 10 | 29.69 | 24.77 | +4.93 |
| 16 | 256 | 11 | unknown | 58.74 | unknown |
| 8192 | 8 | 15 | -7.81 | -12.48 | +4.68 |
| 8192 | 32 | 17 | -3.78 | -7.76 | +3.98 |
| 8192 | 64 | 18 | 4.43 | 0.02 | +4.41 |
| 8192 | 128 | 19 | 23.10 | 16.70 | +6.40 |
| 8192 | 256 | 20 | unknown | 51.57 | unknown |

The extension to \(\log_2 k<8\) was recomputed from the existing evaluator.
The source spectra and transfer calculation are unchanged.

Three observations are stable across these slices.

1. No authenticated BCH constituent is worse than its matched random control.
2. The advantage is several bits once the construction has many rows.
3. At the target \(L=8192\), the advantage rises from 4.41 bits at \(B=64\)
   to 6.40 bits at \(B=128\).

The last row supplies the most direct extrapolation.  The random \(B=256\)
control clears the 40-bit screen by 11.57 bits before assigning any positive
BCH advantage.  Extending the observed 4--6 bit advantage would predict a
margin near 56--58 bits.  This range is a trend estimate, not an interval or a
probabilistic guarantee.

## Low-weight sensitivity

At \(B=256\), \(L=8192\), and the neutral RandomStepConv-M22 schedule, the
random control is dominated by weight 29.  Its expected multiplicity at that
weight is only 0.4413.  The BCH candidate contains no nonzero word below
weight 38 and contains only even-weight words.

A separate shell calculation assigns 82.24 bits of margin to one weight-38
word.  If weight 38 were the only nonzero shell, that shell could contain about
\(2^{42.24}\) words before its margin fell to 40 bits.  The exact relation

\[
  A_{38}(C)=3968m
\]

would reach that shell-only threshold only near \(m=1.3\cdot10^9\).  Higher
shells also contribute to the union bound, so this calculation is a
sensitivity scale and not a sufficient condition.

The application-matched RM2Sub-S19 diagnostic is consistent with the neutral
curve.  At \(N=2^{21}\) and relative distance 0.11, a random-even modeled
\([256,128,\ge38]\) spectrum has 61.47 occupation-one bits.  Weight 38 is the
dominant shell.  The authenticated \([128,64,22]\) spectrum has 26.70 bits in
the corresponding direct RM2Sub calculation.

These results support the following working hypothesis:

> The missing low-weight shells of \(C\) are unlikely to erase the complete
> 11.57-bit random-control cushion, and the known BCH rungs suggest a positive
> rather than negative correction to the random prediction.

The hypothesis would fail if \(C\) placed an exceptional fraction of its
\(2^{128}\) words in the first few admissible shells.  Minimum distance alone
does not exclude that failure mode.

## Quantified low-shell stress test

The minimum distance and evenness of \(C\) permit a more relevant null model
than an unrestricted random linear code.  The model assigns binomial mass to
the even weights from 38 through 218.  It reserves one word at weights 0 and
256.  This model has a RandomStepConv occupation-one margin of 57.55 bits.

The first two admissible shells dominate the result:

| weight | null-model \(\log_2 A_w\) | shell contribution \(\log_2 Q_1\) | shell-only \(\log_2 A_w\) limit at 40 bits |
|---:|---:|---:|---:|
| 38 | 24.28 | -57.96 | 42.24 |
| 40 | 29.20 | -59.79 | 48.99 |
| 42 | 33.95 | -62.48 | 56.43 |
| 44 | 38.54 | -65.29 | 63.83 |

If all other modeled shells remain fixed, the combined multiplicities at
weights 38 and 40 must increase by a factor of approximately
\(2^{17.60}=1.99\cdot10^5\) to erase the 40-bit margin.  The following stress
points show the remaining margin.

| common multiplier for \(A_{38},A_{40}\) | margin |
|---:|---:|
| 2 | 56.58 bits |
| 16 | 53.60 bits |
| 256 | 49.60 bits |
| 4,096 | 45.60 bits |
| 65,536 | 41.60 bits |

Exact smaller BCH spectra provide a scale for this multiplier.  For each code,
the comparison below conditions the even binomial model on the same minimum
distance and total code size.

| code | exact minimum-shell count | conditioned-binomial count | ratio |
|---|---:|---:|---:|
| \([8,4,4]\) | 14 | 14.00 | 1.00 |
| \([32,16,8]\) | 620 | 321.27 | 1.93 |
| \([128,64,22]\) | 243,840 | 324,498.61 | 0.75 |

The largest observed inflation is less than a factor of two.  The target would
need an inflation near two hundred thousand across both dominant shells.  This
comparison is heuristic because three smaller codes do not define a
probability distribution for the missing spectrum.

## What the current linear constraints can prove

The exact properties of \(C\) include complement symmetry, minimum distance
38, and orthogonal-array strength at least 15.  A degree-14 Krawtchouk
majorant uses these facts without estimating individual shells.  This
relaxation gives only -13.44 bits of worst-case margin.  Thus OA-15 alone does
not certify the target.

A normalized linear program also includes the exact BCH sandwich, the dual
sandwich, the dual-code inclusion, and the current affine-orbit shell floors.  Its
binary64 worst-case margin is 6.29 bits.  This result also fails to certify 40
bits.  The maximizing relaxation places about \(2^{97.5}\) words of \(Q\) at
weight 44.  The conditioned-binomial model assigns about \(2^{38.5}\) words to
the corresponding shell of \(C\).  The 59-bit discrepancy identifies a loose
direction in the available constraints; it is not evidence that the BCH code
has that spectrum.

These two failures clarify the evidence boundary.  The current result is a
strong heuristic stress test, but not a proof from the known BCH constraints.

## A proof-sized low-weight target

The RandomStepConv coefficients decay fast enough that a certificate does not
need the complete spectrum.  Consider one shell of weight \(w\).  Two distinct
supports in that shell cannot share a subset of size \(w-18\).  Otherwise, the
corresponding codewords would have distance at most 36.  Therefore

\[
 A_w(C)\binom{w}{w-18}\leq\binom{256}{w-18}.
\]

Applying this packing bound separately through weight 128 controls every
remaining shell.  The trivial \(2^{128}\) bound is smaller from weight 70
onward.

The weight-52 packing term was the dominant limitation.  Johnson-scheme
Delsarte programs now strengthen the terms at weights 52 and 54.  The programs
have 34 and 36 distance variables, respectively.  QSopt_ex solved both over
the exact rationals.  Their exact optima give

\[
\begin{aligned}
 A_{52}(C)&\leq 2979928035058718446462766169<2^{91.268},\\
 A_{54}(C)&\leq 47955805775385501045909339217<2^{95.276}.
\end{aligned}
\]

The Johnson bounds improve the corresponding packing caps by 4.35 and 4.44
bits.  Combining them with the shellwise packing bounds gives 45.80 bits of
RandomStepConv margin for every shell of weight at least 52.  The combinatorial
bounds and Johnson programs are exact.  The reported margin still uses
binary64 transfer coefficients.

The same exact program at weight 38 gives only

\[
 A_{38}(C)\leq971390652389758748<2^{59.753}.
\]

The application cap below is \(2^{41.80}\).  Thus distance-only Delsarte
constraints miss the required minimum-shell bound by about 18 bits.  This
negative result rules out applying the same generic program mechanically to
the seven low shells.  A useful upper bound there must exploit the linear BCH
structure.

The affine invariance of \(P\) gives one intermediate refinement.  Its
weight-38 shell is a 2-design.  Fixing a coordinate pair contained in a
weight-38 word and deleting that pair leaves a constant-weight
\((254,36,38)\) code.  An exact Johnson-scheme program bounds the pair
incidence by

\[
 \lambda_{38}(P)\leq22789092036927452.
\]

Combining this bound with the 2-design lattice gives

\[
 A_{38}(C)\leq128630323189940096<2^{56.837}.
\]

This gains another 2.92 bits over the direct distance-only program, but it
still misses the application cap by 15.04 bits.

### Exact BCH-sandwich bounds

The full BCH sandwich gives substantially stronger bounds.  The exact model
couples the spectra of \(Q\), every nonzero \(Q\)-coset, \(P\), and \(C\).  It
uses the endpoint spectra, MacWilliams identities, dual-code inclusions,
OA-15 equalities, support restrictions, and the rigorous orbit-search floors.
Direct exact solution was ill-conditioned.  An exactly equivalent export
replaces each shell variable by its conditioned-binomial scale and divides
each row by a positive integer.  QSopt_ex then solved each program over the
rationals.  An independent checker reconstructs all 396 rows and verifies
primal feasibility, dual signs, all dual column inequalities, and exact
primal--dual objective equality.

The resulting rigorous caps are:

| weight | exact BCH-sandwich cap on \(A_w(C)\) | \(\log_2\) cap | excess over sufficient application cap |
|---:|---:|---:|---:|
| 38 | 773397229452928 | 49.458 | 7.659 bits |
| 40 | 4210950726984874 | 51.903 | 5.181 bits |
| 42 | 38281539335776556 | 55.087 | 3.612 bits |
| 44 | 1144487702470434098 | 59.989 | 3.924 bits |
| 46 | 7169795968863348791 | 62.637 | 2.137 bits |
| 48 | 56162197953802356662 | 65.606 | 0.825 bits |
| 50 | 645473488233155997799 | 69.129 | 0.212 bits |

At weight 38, \(q_{38}=0\), \(A_{38}(C)=31h_{38}\), and the 2-design lattice
forces \(h_{38}\) to be a multiple of 128.  Applying that lattice after the
continuous optimum gives the displayed integer cap.  The weight-38 result is
10.29 bits stronger than the direct Johnson bound and 7.38 bits stronger than
the residual-design bound.

These caps do not yet prove the target.  Substituting all seven caps into the
binary64 RandomStepConv transfer, together with the exact high-shell bounds,
gives only 32.706 bits of margin.  The weight-38 term alone gives 32.782 bits;
the other six shells barely change the total.  The computation therefore
isolates one concrete remaining problem: improve the weight-38 multiplicity
bound by 7.66 bits.  We no longer need comparable progress on every low
shell.

It remains to bound only

\[
  A_{38}(C),A_{40}(C),\ldots,A_{50}(C).
\]

Even after charging the generic worst-case bound for weights 52 and above,
the seven low shells can jointly exceed their conditioned-binomial values by
about \(2^{17.52}=1.88\cdot10^5\) before the total margin reaches 40 bits.
One simple sufficient condition applies this common multiplier to every low
shell:

| weight | sufficient upper cap on \(\log_2 A_w(C)\) |
|---:|---:|
| 38 | 41.80 |
| 40 | 46.72 |
| 42 | 51.48 |
| 44 | 56.07 |
| 46 | 60.50 |
| 48 | 64.78 |
| 50 | 68.92 |

These caps are not necessary conditions.  They allocate the remaining failure
budget in proportion to the conditioned-binomial shell contributions.

This reduction suggests two proof-oriented computations.

First, a PAC model counter can estimate the seven shell sizes with an explicit
multiplicative error and failure probability.  Directly counting \(C\) misses
useful symmetry.  Instead, count the corresponding shells of \(P\) and \(Q\).
Both codes are affine-invariant, so each shell is a 2-design.  If
\(\lambda_w(P)\) counts weight-\(w\) words containing one fixed coordinate
pair, then

\[
 A_w(P)=\lambda_w(P)\frac{\binom{256}{2}}{\binom{w}{2}},
\]

and the same identity holds for \(Q\).  The exact quotient relation then gives

\[
 A_w(C)=\frac{31A_w(P)+224A_w(Q)}{255}.
\]

ApproxMC provides a factor-\((1+\varepsilon)\) estimate except with probability
\(\delta\).  A union bound can allocate one total failure probability across
the fourteen counts.  A preliminary ordinary-CNF weight-38 count did not
finish in ten minutes.  A native-XOR satisfiability probe with the fixed-pair
reduction did not finish in two minutes.  The method is mathematically sound,
but the current encoding is not yet an interactive computation.

ApproxMCPB was also built locally to preserve the cardinality constraint as a
pseudo-Boolean equality and the parity checks as native XORs.  The fixed-pair
encoding exactly recovers the known incidence counts 3 and 35 for the
length-8 and length-32 codes.  The corresponding known length-128 count is
6930, but that calibration did not finish after several solver configurations.
The failure already at length 128 makes this route unsuitable for a credible
length-256 estimate without a materially better encoding or counter.

Second, a deterministic search can enumerate affine-orbit representatives in
the same seven shells.  A complete orbit list plus independently checkable
membership, stabilizer, and exhaustion certificates would give exact counts.
This route should start from the known Wambach representatives rather than ask
a generic solver to rediscover minimum words.

### Calibrated systematic-basis search

The deterministic search has now been implemented for the minimum shells.  It
constructs every inequivalent primitive-element representation, puts its
generator matrix in systematic form, and enumerates sums of at most (r)
rows.  The hot loop uses four fixed 64-bit limbs.  Each discovered word is
mapped back to the standard coordinates and checked for code membership.
Affine orbits are identified exactly by mapping each ordered pair of support
points to ((0,1)).  Distinct orbit closures therefore give rigorous lower
bounds, although the search radius does not prove completeness.

The method has a useful exact calibration.  At radius (r=4), its affine
closure recovers the complete published minimum shell at every smaller
rate-half BCH rung:

| code | exact minimum shell | recovered shell | recall |
|---|---:|---:|---:|
| extended BCH ([8,4,4]) | (A_4=14) | 14 | 100% |
| extended BCH ([32,16,8]) | (A_8=620) | 620 | 100% |
| extended BCH ([128,64,22]) | (A_{22}=243840) | 243840 | 100% |

For the length-256 sandwich codes, radius four finds 14 affine orbits in the
weight-38 shell of (P), and 25 affine orbits in the weight-40 shell of (Q).
Radius five is not saturated: it finds 115 weight-38 orbits in (P) and 222
weight-40 orbits in (Q).  The resulting rigorous bounds are

\[
 A_{38}(P)\ge7507200,
 \qquad
 A_{38}(C)\ge\frac{31}{255}7507200=912640,
\]

and

\[
 A_{40}(Q)\ge12990720.
\]

The radius-five \(P\) search also proves \(41583360\) weight-40 words outside
\(Q\).  Combining their affine closures with the independent \(Q\) closure
gives the rigorous union bound

\[
 A_{40}(C)\ge12990720+\frac{31}{255}41583360=18045952.
\]

These are lower bounds, not substitutes for the upper bounds needed in a
certificate.  They nevertheless provide calibrated evidence about scale.  A
Chao1 extrapolation from the number of singleton and doubleton orbit hits gives
the heuristic values

\[
 A_{38}(C)\mathrel{\approx}9.69\cdot10^5,
 \qquad
 A_{40}(C)\mathrel{\approx}1.81\cdot10^7.
\]

The row-search observations are not independent samples, so this extrapolation
has no PAC guarantee.  Its value is diagnostic.  The search was then extended
through weight 44.  The application-scale results are summarized below.

| weight | rigorous \(A_w(C)\) floor | Chao1 heuristic | heuristic / random | application cap / heuristic |
|---:|---:|---:|---:|---:|
| 38 | 912640 | \(9.69\cdot10^5\) | 0.0953 | \(3.95\cdot10^6\) |
| 40 | 18045952 | \(1.81\cdot10^7\) | 0.0586 | \(6.42\cdot10^6\) |
| 42 | 42289664 | \(4.73\cdot10^7\) | 0.00569 | \(6.61\cdot10^7\) |
| 44 | 267281408 | \(3.04\cdot10^8\) | 0.00152 | \(2.48\cdot10^8\) |
| 46 | 1247821056 | \(1.45\cdot10^9\) | 0.000334 | \(1.13\cdot10^9\) |
| 48 | 5525363968 | \(6.40\cdot10^9\) | 0.0000760 | \(4.96\cdot10^9\) |
| 50 | 1574477312 | \(1.84\cdot10^9\) | 0.00000124 | \(3.03\cdot10^{11}\) |

Thus every probed shell is far below both its matched random expectation and
the sufficient application cap.  In particular, the weight-44 search gives no
support for the normalized LP's pathological \(2^{97.5}\)-word spike.  The
weight-50 row uses radius four; every earlier row uses radius five.

Charging all presently proved words through weight 50 to the transfer, and no
other nontrivial words, leaves 62.19 bits of margin.  Replacing those floors by
the calibrated Chao1 estimates leaves 62.12 bits for the seven low shells.

The proof-shaped heuristic combines those seven estimates with the exact
combinatorial bounds for every weight at least 52.  Under the binary64
RandomStepConv transfer, the combined margin is 45.799588 bits.  The low-shell
estimates can increase jointly by a factor of \(4.47\cdot10^6\) before that
margin reaches 40 bits.

This calculation is not a certificate.  Chao1 does not upper-bound affine
orbits that the systematic-basis search never observes.  The reported transfer
also uses nearest binary64 arithmetic.  The high-weight bounds still set the
45.80-bit result almost entirely, but they now leave a 5.80-bit cushion.  The
result therefore gives strong evidence for the low shells.  Turning the seven
low-shell estimates into certified upper bounds, and enclosing the transfer
arithmetic, are the remaining proof obligations.

The exact Johnson receipts are the `generated/johnson_n256_w38_d38.*`,
`generated/johnson_n256_w52_d38.*`, and
`generated/johnson_n256_w54_d38.*` files.  The reproducible search summary is
`generated/wambach_search_evidence.json`.  The compact JSON and TSV orbit
receipts retain one canonical representative, orbit size, stabilizer size, and
search-hit multiplicity for every discovered orbit.

## What the literature supports

The published boundary is consistent across the available sources.  Complete
weight distributions are known for the length-128 extended primitive BCH
family.  At length 256, published complete distributions cover dimensions at
most 71 and at least 187, not the middle dimensions near 128.  The Okayama
[weight-distribution index](https://isec.ec.okayama-u.ac.jp/home/kusaka/wd/)
and the [Fujiwara--Kusaka paper](https://globals.ieice.org/en_transactions/fundamentals/10.1587/transfun.2020EAP1119/_p)
record this boundary.

Sala and Tamponi developed linear programs for unknown BCH distributions using
known subcodes, supercodes, MacWilliams identities, and Pless identities.  Of
particular relevance, their Method D maximizes the downstream error functional
directly instead of reconstructing every shell.  This is the same strategy
needed for the SPIN transfer.  See [A Linear Programming Estimate of the
Weight Distribution of BCH(255,k)](https://eprints.soton.ac.uk/254459/1/46it06-shrt2.pdf).

Classical results show that broad portions of many primitive BCH spectra are
close to binomial.  Krasikov and Litsyn give a quantitative error estimate in
[On Spectra of BCH Codes](https://bura.brunel.ac.uk/bitstream/2438/3321/1/On%20spectra%20of%20BCH%20codes.pdf).
Their displayed parameter regime does not cover the low tail of the present
\(t=18\), length-255 code.  The theorem therefore supports central-spectrum
plausibility but does not justify the weight-37 or weight-38 multiplicity.

Augot, Charpin, and Sendrier analyze minimum BCH words through their locator
polynomials and Newton identities in
[Studying the Locator Polynomials of Minimum Weight Codewords of BCH Codes](https://www-rocq.inria.fr/secret/Pascale.Charpin/AugChaSen-i3e.pdf).
Their results establish the true minimum distances at length 255 but do not
enumerate the weight-37 shell needed here.  Their method nevertheless exposes
the right next reduction.  For a weight-37 support in \(B(255,37)\), the first
36 power sums vanish.  Newton identities force every odd locator coefficient
below degree 37 to vanish.  Because \(\gcd(37,255)=1\), a unique cyclic shift
normalizes the leading locator coefficient.  After applying the inverse
Frobenius map, each normalized word is therefore equivalent to a polynomial
\(g\) of degree at most 18, with fixed constant term, that agrees with a fixed
monomial on 37 distinct points of \(\mathbb F_{256}^{*}\).  This
Reed--Solomon-coset formulation is a deduction from the locator equations, not
a counting theorem claimed by the paper.

The list-decoding problem has a useful endpoint property.  For every candidate
\(g\), the substitution \(z=u^{128}\) identifies its agreement points with the
roots of \(z^{37}+g(z^2)\).  This is a nonzero polynomial of degree 37.
Consequently, no candidate has more than 37 agreements, and the desired list
is exactly the endpoint \(N(g)=37\), rather than the tail \(N(g)\geq37\).

Let \(L_{37}\) denote the number of normalized polynomials.  The affine
2-design identity and the cyclic normalization give

\[
 L_{37}=\frac{A_{37}(P_{\mathrm{punct}})}{255}
       =\frac{19}{128}h_{38},
 \qquad
 A_{38}(C)=\frac{3968}{19}L_{37}.
\]

Consequently, the occupation-one target follows from
\(L_{37}<2^{34.093}\).  This target is 7.66 bits below the list-size cap
induced by the exact BCH-sandwich program.
The executable check `code/verify_locator_reduction.py` verifies all 36
syndromes, the unique cyclic normalization, the locator form, and exactly 37
agreements on Wambach's published minimum word.  Its receipt is
`generated/locator_reduction_wambach.json`.

This formulation also identifies a proof-sized incidence target.  Each listed
polynomial contributes \(\binom{37}{24}\) pairs \((g,T)\), where \(T\) is a
24-point agreement set.  On a fixed \(T\), the 18 free coefficients of \(g\)
face six more field equations than unknowns.  The random-rank baseline for the
number of admissible \(T\)'s is therefore \(\binom{255}{24}/256^6\).  If the
actual count is at most 6.069 times this baseline after charging the
derivative-zero family below, double counting gives
\(L_{37}<2^{34.093}\).  Thus the next theorem needs a constant-factor point
bound for one six-equation incidence variety.  It does not require the full
BCH spectrum.

The support is a super-Vandermonde set in the terminology of Sziklai and
Takáts.  Their
[classification paper](https://doi.org/10.1016/j.ffa.2008.06.004) proves that
the vanishing power sums are equivalent to a fully split polynomial of the
form \(Y^{37}+g(Y)^2\).  The available classification covers sizes below the
characteristic and above \(q/2\).  It does not cover size 37 in
\(\mathbb F_{256}\).  Blokhuis, Marino, Mazzocca, and Polverino later
[classified the adjacent sizes](https://doi.org/10.1007/s10623-016-0254-z)
\(p+1\) and \(q/p-1\), which specialize here to 3 and 127.  Their result also
leaves size 37 open.

The six-equation incidence model nevertheless yields two rigorous structural
facts.  Let \(R_T\) be the remainder of \(u^{146}\) modulo
\(\prod_{t\in T}(u+t)\).  The incidence equations set the constant coefficient
of \(R_T\) to one and its coefficients of degrees 19 through 23 to zero.

First, suppose that the corresponding polynomial \(g\) satisfies
\(g'\ne0\).  Differentiating the interpolation equations shows that each
Jacobian column is a scaled Lagrange polynomial.  After invertible row and
column operations, any six columns with \(g'(t)\ne0\) form the Vandermonde
matrix on \(1,t,\ldots,t^5\).  Since \(\deg g'\le16\), at least eight columns
survive.  The Jacobian therefore has rank six at every nondegenerate
admissible 24-set.

Second, suppose that \(g'=0\).  Then \(g=h^2\) for a polynomial \(h\) of degree
at most nine, and \(h(u)=u^{73}\) on the 37 agreement points.  Any nine
nonzero points determine \(h\).  Double counting the nine-point subsets gives
the exact bound

\[
 L_{37}^{\mathrm{deg}}
 \leq
 \left\lfloor\frac{\binom{255}{9}}{\binom{37}{9}}\right\rfloor
 =87550900<2^{26.384}.
\]

This degenerate family is more than seven bits below the locator target.
After charging it in full, the smooth incidence count may be at most 6.069
times the random-rank reference.  The remaining proof obligation is thus a
point bound for the smooth, codimension-six locus.

As a consistency check, the 138 certified affine orbits generate 5244
distinct normalized locators.  Every locator has exactly 37 agreements.
All 20976 sampled admissible 24-sets have Jacobian rank six.  These observations
follow the theorem above but do not bound undiscovered points.  The receipt is
`generated/locator_incidence_jacobian.json`.

A separate importance-sampling experiment probes the point count itself.
Choose 18 nonzero points \(B\) uniformly and interpolate the unique candidate
\(g_B\) through them and \(g_B(0)=1\).  If \(E_B=N(g_B)-18\), then

\[
 \frac{\mathbb E_B\binom{E_B}{6}}
      {\binom{237}{6}/256^6}
 =
 \frac{\#\{\text{admissible 24-sets}\}}
      {\binom{255}{24}/256^6}.
\]

Five million independent samples give a ratio of 1.032; the normal 95-percent
upper endpoint computed from the sample variance is 1.136.  The application
allows 6.069 after charging the degenerate family.  This is direct evidence
for the exact statistic that remains to be bounded, but it is not a rigorous
confidence bound for unseen algebraic components.  The source is
`code/sample_locator_extensions.cpp`, and the receipt is
`generated/locator_extension_monte_carlo.json`.

The same statistic now has a size ladder. Exhaustive enumeration gives ratios
1.829 at \(q=8\) and 1.994 at \(q=32\), at the highest available moment order.
Five-million-sample estimates of the sixth-order ratio give 0.784 at \(q=128\)
and 1.032 at \(q=256\). In contrast, the exact fully split endpoint at
\(q=128\) is 96.94 times its random reference. The broad random-splitting
model is therefore false, but the narrower incidence model remains consistent
with every available size. The formal target lemma, heuristic assumption,
and proof routes are separated in `LOCATOR_INCIDENCE_DUAL_TRACK.md`.

The application side of that target is now rigorous. Arb enclosures bound the
transcendental terms, positive binary64 transfer recurrences are rounded
outward, and the shell aggregation uses exact rational arithmetic. The result
is the sufficient integer cap

\[
 A_{38}(C)\leq3{,}827{,}351{,}840{,}403.
\]

In incidence normalization, the certified lower endpoint for the permitted
factor is \(6.098404530135757\). Thus \(R_6\leq6\) suffices with 0.02347 bits
of margin. The application threshold is no longer conditional on nearest
binary64 arithmetic.

The first character-sum refinement has also been tested exactly. The
exceptional ordinary Walsh ridge applies directly to affine perturbations
\(u^{146}+cu+1\). Its signed contribution, summed over all \(c\), is \(-64\).
Moreover, exhaustive enumeration shows that perturbations of degree at most
one, two, and three have at most 4, 7, and 10 roots, respectively. All three
slices therefore contribute zero to the 18-base and 24-point locator
incidences.

The 5,244 normalized locators already obtained from the certified affine
orbits are concentrated at the other endpoint: 5,206 have degree 18, 31 have
degree 17, and seven have degree 16. This is consistent with the low-degree
exclusion, although it remains evidence only for the discovered components.

The uniform-base sampler now records matching strata. Degree-17 interpolants
contribute 0.00463 to the total sixth-moment ratio, while degree-18
interpolants contribute 1.02708. Conditioning on whether the linear
coefficient lies in the affine ridge image changes the estimated conditional
ratio from 1.024 to 1.063; the normal 95-percent upper endpoint in the ridge
image is 1.295. Neither stratum exhibits consequential concentration.

An exact exceptional-family calculation reaches the same conclusion without
sampling. Exhausting all \(2^{18}\) interpolants whose coefficients lie in
\(\mathbb F_2\) finds 10 endpoint locators, but the entire family's
contribution is only \(3.33\cdot10^{-9}\) on the \(R_6\) scale. On the set
side, 112 of the 5,365 Frobenius-invariant 24-sets are admissible, despite an
independent-rank expectation near \(1.9\cdot10^{-11}\) sets. The enrichment is
real and extreme, but its absolute \(R_6\) contribution is only
\(1.05\cdot10^{-17}\).

This enumeration also finds four new affine weight-38 orbits, each of size
65,280. The rigorous floors improve to
\(A_{38}(P)\geq7{,}768{,}320\) and \(A_{38}(C)\geq944{,}384\).

The next nested subfield has now been exhausted rather than sampled. The map
\(x\mapsto x^4\) partitions \(\mathbb F_{256}^{*}\) into 3 singleton, 6 pair,
and 60 four-element orbits. Of all 267,516,561 invariant 24-sets, exactly
67,900 are admissible. They recover all 40,712 \(\mathbb F_4\)-coefficient
interpolants with at least 24 agreements. The resulting contribution is
\(R_6^{(\mathbb F_4)}=1.134\cdot10^{-8}\); 99.786 percent is supplied by the
34 endpoint locators. Thus the larger algebraic family remains negligible on
the application scale, although endpoint domination is pronounced.

Those 34 locators occupy 15 affine orbits. Ten are new and each has size
65,280. The rigorous floors therefore improve again to
\(A_{38}(P)\geq8{,}421{,}120\) and \(A_{38}(C)\geq1{,}023{,}744\). The exact
receipt is generated/frobenius4_invariant_incidence.json.

The \(\mathbb F_{16}\) stratum requires sampling because it contains
\(3.488\cdot10^{17}\) invariant 24-sets. An invariant 18-point base determines
a unique \(\mathbb F_{16}\)-coefficient interpolant. Weighting each agreeing
six-point extension by the reciprocal number of bases in the resulting
24-set gives an unbiased estimator for this stratum.

Ten million samples give ratio 1.00624 with standard error 0.00412. The normal
95-percent upper endpoint is 1.01432. This estimate is evidence rather than a
rigorous upper bound. The ten sampled endpoint polynomials are exact
witnesses. They occupy nine new affine orbits and certify
\(A_{38}(P)\geq9{,}008{,}640\) and
\(A_{38}(C)\geq1{,}095{,}168\).

The nine affine closures contain 64 new \(\mathbb F_{16}\)-coefficient
locators. In total, 190 certified normalized locators in 30 affine orbits have
coefficients in \(\mathbb F_{16}\). This inventory is a lower bound; the
GF(16) sample is not exhaustive.

This exact result narrows, but does not finish, the proof. After polynomial
coefficient orthogonality, the full 24-point character sum is indexed by a
six-dimensional Vandermonde kernel. Its multipliers contain
\(P_T'(x_i)^{-1}\), which couples the 24 points. The ordinary one-variable
Walsh spectrum therefore does not factor the full sum. A successful
character-sum proof must exploit this kernel form or control generalized
higher-degree perturbations.

The eighteen universal Reed--Solomon factorial moments have also been
optimized exactly. They permit \(R_6=3.3072\cdot10^7\), so generic MDS moment
identities are insufficient. The useful monomial-specific replacement is the
symmetric system

\[
 e_{24}h_{122}=1,\qquad h_{123}=\cdots=h_{127}=0.
\]

The characteristic-two identity \(H=EH^2\) halves the indices in this system
by Frobenius squaring. This is now the preferred basis for a
computer-assisted proof.

## Strongest next evidence

The next computations should answer the application functional rather than
attempt the complete spectrum.

1. Recompute the exact-spectrum \(B=32,64,128\) ladder with the final
   RM2Sub transfer at fixed row counts.
2. Recompute the shell stress test with the final RM2Sub coefficients.
3. Attack the normalized weight-37 locator-polynomial list through the
   six-dimensional Vandermonde-kernel character sum. A bound
   \(L_{37}<2^{34.093}\), about 7.66 bits below the present exact LP cap, is
   sufficient; exact enumeration is unnecessary. The clean target is
   \(R_6\leq6\).
4. Add any locator-derived inequalities to the exact scaled sandwich LP and
   regenerate its rational primal--dual certificate.
5. Evaluate occupations beyond one with the same fixed-code spectrum model.
   For a fixed constituent, row-weight profile counts factor through the
   ordinary enumerator; a biweight enumerator is not required after the
   independent row permutations.

The present evidence rates the RandomStepConv occupation-one hypothesis as
strong.
It does not yet rate the all-occupation claim.  Expanded Accumulant should
remain the formal fallback until the fixed-BCH transfer survives the
high-occupation diagnostics.
