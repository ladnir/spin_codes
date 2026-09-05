# A direct statistical test for the BCH weight-38 cap

Updated: 2026-09-04

The fixed BCH-derived code C has parameters [256,128,38]. The current
RandomStepConv low-shell allocation assigns weight 38 the target

\[
 A_{38}(C)\le U,\qquad U:=3{,}827{,}351{,}840{,}403.
\]

This is an allocated shell target, not a sufficient condition for the entire
occupation-one functional. The 2026-09-04 application audit found that the
other current LP caps still prevent closure. With this weight-38 target,
their combined directed upper bound has a margin of 36.8688 bits, below 40.
Even removing weight 38 entirely leaves only 36.9955 certified bits.
See `generated/endpoint_application_budget_audit.json`.

The sixth incidence moment is one sufficient route to this shell cap. A direct
endpoint test gives another route with an explicit statistical error bound.
It does not assume a random-like BCH spectrum or estimate a variance.
The fixed experiment is retained under `generated/endpoint_certificate_20260904`.
It completed with **conditional statistical acceptance**: 150,361,035 samples,
zero endpoint hits, and a complete matching replay using the second algorithm.
The Python checker independently verified all 407 retained audit vectors.
That directory's `certificate.json` records the source hashes, byte-stream
hash, exact probability inequality, and randomness qualification.

The rational upper bound for the log2 false-accept probability is approximately
-40.00000016075008 under ideal independent bytes. Actual execution used Windows
CNG cryptographic pseudorandomness, so this is not an unconditional algebraic
proof. The archive contains 2,812,280,832 bytes, of which 2,812,026,464 were
consumed. The primary pass took 243.613 seconds and the sequential full replay
took 1,053.709 seconds; these are run durations, not controlled benchmarks.

The primary histogram also gives the diagnostic estimate R_6=0.9753118802.
This estimate is distinct from the accepted endpoint test and is not a
certified upper bound on R_6.

## The endpoint probability

Let D be the 255 nonzero elements of GF(256). For each 18-subset B of D,
let g_B be the unique polynomial of degree at most 18 satisfying

\[
 g_B(0)=1,\qquad g_B(u)=u^{146}\quad(u\in B).
\]

Define N_B as its number of agreements on D. The polynomial
z^37+g_B(z^2) has degree 37, so N_B is at most 37.
Let L_37 count the normalized polynomials with exactly 37 agreements.
The existing BCH locator reduction gives

\[
 A_{38}(C)=\frac{3968}{19}L_{37}.
\]

For a uniform 18-subset B, define p:=Pr[N_B=37]. Then

\[
 \boxed{p=L_{37}\frac{\binom{37}{18}}{\binom{255}{18}}.}
 \tag{1}
\]

To prove (1), fix a normalized endpoint polynomial g. Exactly C(37,18)
bases lie within its agreement set, and every such base interpolates to g.
Two different polynomials cannot share a base: the 18 agreement conditions
and the condition at zero determine a polynomial of degree at most 18 uniquely.
Thus endpoint bases form disjoint families of equal size. Dividing their
total size by C(255,18) proves the identity.

This counts all endpoint polynomials, regardless of coefficient subfield,
derivative, or previously discovered affine orbit. GF(16)-invariant bases
have a different law and cannot be substituted into (1).

## A fixed-sample acceptance test

Since A_38(C) and L_37 are integers and gcd(3968,19)=1, L_37 is a multiple
of 19 and A_38(C) is a multiple of 3968. Therefore violation of the cap forces

\[
 L_{37}\ge L_{\rm bad}:=
 19\left(\left\lfloor\frac{U}{3968}\right\rfloor+1\right)
 =18{,}326{,}533{,}524.
\]

Define the exact rational probability

\[
 p_{\rm bad}:=L_{\rm bad}\frac{\binom{37}{18}}{\binom{255}{18}}
 \approx1.84395410537\cdot10^{-7}.
\]

Fix a sample count n before drawing fresh independent uniform bases
B_1,...,B_n. Count an endpoint whenever N_{B_i}=37. Accept the shell cap
if and only if the endpoint count is zero; otherwise return inconclusive.

For any fixed code violating the cap, the false-accept probability satisfies

\[
 \Pr[\text{accept}] = (1-p)^n \le (1-p_{\rm bad})^n.
 \tag{2}
\]

Consequently, the following sample counts suffice when no endpoints are seen.
They are the smallest integers satisfying the right-hand bound in (2).

| Maximum false-accept probability | Fresh independent samples |
|---|---:|
| 1/20 | 16,246,240 |
| 1/100 | 24,974,428 |
| 2^-40 | 150,361,035 |
| 2^-80 | 300,722,069 |

These are error guarantees for a randomized decision procedure. They are not
posterior probabilities that this fixed code is good. Acceptance does not
prove R_6<=6; it statistically certifies the shell condition directly.
The error probability is separate from the SPIN setup failure probability.
Closing this shell also does not cover occupations beyond one.

A nonzero endpoint count does not refute the shell cap. A future test could
use the full binomial lower tail to allow hits. Such a test needs its own
predeclared sample count, acceptance rule, and error calculation.

## What the old sample establishes

The retained five-million-sample histogram has no endpoints: its largest
extra agreement count is nine, whereas an endpoint has nineteen.
Even under the ideal IID model, a shell at the first prohibited lattice
value would produce zero hits with probability about 0.397732.
Thus this sample does not exclude a violating endpoint family at 95% confidence.

A conservative IID 95% zero-hit confidence rule uses
p<=log(20)/n. Combined with (1) and the shell lattice, its hypothetical upper
cap for the old run is A_38(C)<=12,436,016,063,872. This is still above U.

The existing sampler uses std::mt19937_64 with the fixed seed 188203230.
Its partial Fisher-Yates procedure generates a uniform subset if its integer
draws are independent and uniform. Each partial shuffle is undone before
the next sample. However, a fixed deterministic seed supplies no statistical
probability space, and independent IID draws do not follow merely from
seeding MT19937 randomly. The old run remains diagnostic evidence.

A cryptographic random-bit generator would instead introduce a computational
randomness assumption. Under an ideal random-bit model, unbiased subset
sampling gives the information-theoretic guarantee (2). A report must state
which randomness model its execution uses.

## Arithmetic and checks

Run:

    python code/analyze_endpoint_sampling.py

The script writes generated/endpoint_sampling_analysis.json. It hashes the
directed application threshold, original histogram, smaller-field receipt,
and its own source. Every sample-count decision uses exact rational arithmetic.
Decimal values are display-only approximations.

For 1<=x<=2, put z=(x-1)/(x+1). The script encloses log(x) using the first
m terms of 2 sum_{j>=0} z^(2j+1)/(2j+1). Its positive remainder is at most

\[
 \frac{2z^{2m+1}}{(2m+1)(1-z^2)}.
\]

Range reduction handles log(20), log(100), and integer powers of two.
For -log(1-p), truncating sum_{j>=1} p^j/j after m terms leaves at most
p^(m+1)/((m+1)(1-p)). Rational enclosures verify both
n[-log(1-p_bad)]>=log(1/alpha) and failure at n-1.

The checker also verifies the endpoint identity against the stored complete
q=8 and q=32 base enumerations. Exact small Bernoulli examples independently
check the sample-count logic against direct rational powers.

## Next experiment and continuing proof work

Freeze the target, source hashes, sample count, randomness model, and
acceptance rule before fresh sampling. Audit exact subset selection and
finite-field evaluation. Count every endpoint, including known ones.
Retain the input provenance, total draw count, hits, and run parameters.
Do not pool invariant-base experiments with the uniform-base experiment.
Do not repeat the test until it accepts without budgeting the combined error.
Schedule the run without another concurrent benchmark.

The implementation is `code/test_bch_endpoints.cpp`, built with MSVC C++20,
optimization, and AVX2. `code/run_endpoint_certificate.py` writes an exclusive
run manifest before drawing any experimental bytes. It freezes the sample
count, shell cap, source hashes, executable hash, and acceptance rule.

The sampler calls Windows `BCryptGenRandom` in one-MiB blocks and saves every
block. Within each base, it rejects zero bytes and duplicate nonzero bytes
until it has 18 distinct points. Under independent uniform bytes, the next
accepted point is uniform among the remaining points, so the resulting base
is uniform. The raw stream permits complete deterministic replay.
Microsoft documents this interface and its cryptographic generator in the
[BCryptGenRandom reference](https://learn.microsoft.com/en-us/windows/win32/api/bcrypt/nf-bcrypt-bcryptgenrandom).
This computational generator does not remove the randomness qualification above.

The primary algorithm uses Newton divided differences and AVX2 evaluation
through precomputed monomial contributions. The replay uses incremental
interpolation with a vanishing polynomial and scalar Horner evaluation with
eight independent chains. It also uses a separate duplicate-detection layout.
The driver compares full histograms, endpoint counts, byte consumption, and
rejected-byte counts. Every recorded audit vector is checked a third time by
Python Gaussian elimination on its Vandermonde matrix. Source tables are
checked against an independent carryless-product field implementation.
The replay checks the same random samples; it does not increase the
statistical sample count.

The end-to-end preflight is stored in `generated/endpoint_preflight_20260904`.
Its 4,096 samples are explicitly diagnostic and excluded from the statistical
run. It passed full replay and 257 independent Python audit-vector checks.

To launch one new fixed experiment into a new directory:

    python code/run_endpoint_certificate.py execute generated/endpoint_certificate_20260904

The driver samples and replays sequentially. It refuses an existing run
directory, and does not retry an inconclusive result. To recheck the retained
manifest, full-replay comparison, audit vectors, hashes, and rational inequality:

    python code/run_endpoint_certificate.py verify generated/endpoint_certificate_20260904

The algebraic track remains open. It seeks a deterministic bound on L_37
or on the broader incidence count. The direct statistical test provides a
finite, quantitative intermediate target while that work continues.
