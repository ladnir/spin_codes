# Finite comparison at message dimension \(2^{20}\)

**Status:** resumed on 2026-08-31 after the implementation comparison. The
independent-row B=240 transposed map is the first finite proof target. The
exact setup law is frozen in `FINITE_K20_INDEPENDENT_SETUP.md`.

## Result

The current finite analysis does not yet certify 40 bits of setup failure
for the complete Golay--BA-3/RM2Sub-S19 ensemble. It does identify a useful
concrete design point.

Set

\[
  k=2^{20},\qquad B=240,\qquad L=8832,\qquad N_+=BL=2{,}119{,}680.
\]

The parent code has dimension \(N_+/2=1{,}059{,}840\). Shorten 11,264
fixed input coordinates to obtain dimension \(k\). The resulting rate is
\(0.4946859903\), and the target bad event is

\[
  \mathcal F:=\{d_{\min}\le 233{,}164\}.
\]

Thus \(\mathcal F^c\) gives minimum distance at least 233,165, which is
slightly greater than \(0.11N_+\).

For independently sampled BA rows under the revised \([23,217]\) conditioning
window, the binary64 one-active diagnostic gives

\[
  \log_2 \mathbb E[Z_{233164,1}]\le -46.4966.
\]

The union over occupations 2 through 64 gives a further diagnostic bound of
\(-84.6173\) in base-two logarithms. The one-active class therefore remains
the bottleneck in the computed range.

The occupation-one calculation has since been translated to one-sided
binary64 arithmetic. Its proved upper bound is

\[
  \mathbb E[Z_{233164,1}]
  \le \mathtt{0x1.6f3c66666f370p-47}<2^{-46}.
\]

Occupations 2 through 64 have also been translated outward. Their aggregate
is at most

\[
  \mathtt{0x1.4dc4b8b530941p-1}\cdot2^{-84}<2^{-84}.
\]

The union through occupation 64 remains below \(2^{-46}\). Occupations above
64 remain incomplete.

## Probability spaces

One BA draw consists of two independent uniform permutations in the
Golay--BA-3 encoder. Define \(\mathcal G_B\) as the event that the sampled BA
code has no nonzero word outside the weight interval

\[
  [\lceil0.104B\rceil,\lfloor0.896B\rfloor].
\]

The finite comparison uses three distinct ensembles.

1. **Reused conditioned outer.** Sample one BA code conditioned on
   \(\mathcal G_B\), then reuse it in all \(L\) rows. Sample every route
   permutation and RM2Sub multiplier independently of that code.
2. **Independent conditioned outers.** For each row, sample an independent BA
   code conditioned on \(\mathcal G_B\). Sample the route and inner randomness
   as above.
3. **Ideal random comparator.** Sample a uniform random binary linear map from
   \(k\) input bits to \(N_+\) output bits. This is a distance benchmark, not
   a linear-time construction.

The one-active expectation applies to either conditioned BA ensemble. The
multi-active products apply directly only to the independent conditioned
ensemble. For the reused ensemble, one needs an authenticated spectrum bound
for the selected BA code or higher moments of its sampled spectrum.

For each block size below, the expected BA spectrum and Markov's inequality
give a lower bound on \(\Pr[\mathcal G_B]\). Rejection sampling has the stated
expected number of trials only if setup can decide \(\mathcal G_B\).

## Admissible lengths and shortening

This finite family requires

\[
  24\mid B,\qquad 128\mid L,\qquad N_+=BL,
  \qquad k\le N_+/2.
\]

For a requested \(k\), choose an admissible \(B\), then set

\[
  L:=128\left\lceil\frac{2k}{128B}\right\rceil.
\]

Shorten \(N_+/2-k\) fixed input coordinates. Shortening restricts the parent
message space, so it cannot add a bad codeword. All calculations below use
the full parent message space except the ideal-random comparator. They are
therefore conservative for the shortened structured code.

Puncturing to exactly \(2k\) output coordinates is unattractive here. It
would lose up to \(N_+-2k\) units of distance. For \(B=240\), an 11-percent
guarantee after puncturing would require parent distance about 11.95 percent.

## Block-size comparison

The table reports expected first-moment margins. A margin \(m\) means an
expected bad-codeword count of at most \(2^{-m}\) in the stated diagnostic.

| \(B\) | \(L\) | \(N_+\) | Rate after shortening | Length overhead over \(2k\) | \(\Pr[\mathcal G_B]\) lower bound | Trials upper bound | Unconditioned \(Q=1\) margin | Conditioned \(Q=1\) margin |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 216 | 9,728 | 2,101,248 | 0.499026 | 0.195% | 0.902235 | 1.10836 | 4.123 | 41.780 |
| 240 | 8,832 | 2,119,680 | 0.494686 | 1.074% | 0.993661 | 1.00638 | 5.034 | 46.497 |
| 264 | 8,064 | 2,128,896 | 0.492544 | 1.514% | 0.919553 | 1.08749 | 5.867 | 54.795 |

The nonmonotone good-event bound comes from integer endpoints in the window
and the finite BA spectrum. It is not a monotonicity failure in the code.

The ideal random comparator has margins 2,235.40, 11,454.52, and 16,062.58
bits at the three shortened rates. Most of that increase comes from
shortening. At the unshortened parent rate of one half, the corresponding
ideal margins are 187.40, 190.52, and 190.58 bits.

The random-code benchmark therefore has much more distance slack. A dense
random generator would store about 256.5 GiB at the \(B=216\) point and would
not provide the linear-time encoder sought here.

## Occupation diagnostics

At \(B=216\), the finite two-state RM2Sub transfer gives the following
conditioned diagnostics.

| Active rows \(Q\) | Coverage | Aggregate or pointwise margin |
|---:|---|---:|
| 1 | exact sum over outer weights | 41.780 bits |
| 2--64 | every integer \(Q\) | 68.786 bits aggregate |
| 65--512 | \(65,96,128,192,256,384,512\) | at least 2,520.44 bits |
| 768--9,728 | ten sampled values | unresolved by the current envelope |

At \(B=240\), the revised-window conditioned one-active margin is 46.479 bits
outward. The union over every \(Q\) from 2 through 64 has an outward margin of
84.617 bits. The earlier 48.019- and 89.360-bit values used the narrower
\([25,215]\) law and are retained only as historical comparisons.

The dense-range failure is a limitation of the comparison. A single
worst-weight BA envelope assigns the likelihood of weight 23 to every active
row. A full-spectrum Holder calculation reduces this loss but remains
vacuous at some large occupations. In particular, the endpoint moments
dominate its large Holder powers. These vacuous bounds are not evidence that
the code has low distance.

## Running time and setup size

For \(B=240\), one online encoding performs the following structured work:

- 88,320 constant-size Golay encodings;
- 4,239,360 accumulator bit updates;
- a constant number of permutation passes over 2,119,680 bits; and
- 16,560 fixed-width RM2Sub epochs.

Both the reused and independent-outer variants therefore have \(O(N_+)\)
ordinary and transposed encoding work. Independent BA rows do not add
arithmetic work. They replace two reused BA permutation tables by two tables
per row, which worsens memory traffic.

With 32-bit permutation entries, the reused variant stores approximately
16.17 MiB of route and BA indices. The independent variant stores about
32.34 MiB. Packed indices can reduce both values. A pseudorandom permutation
generator can trade stored tables for setup and evaluation work.

At \(B=240\), the revised lower bound gives at most 1.00638 BA trials per row
under rejection sampling. Across 8,832 rows, this is at most 9,257 expected
BA draws. This statement does not include the cost of deciding
\(\mathcal G_{240}^{23}\). No linear-time verifier for that event is supplied.

## Proof status

The following facts are proved in exact finite form.

- The BA expected spectrum has the combinatorial formula used by the
  evaluator.
- Markov's inequality gives the displayed lower bounds on
  \(\Pr[\mathcal G_B]\).
- Shortening cannot increase the bad-codeword count.
- The fixed RM2Sub transfer matrices are valid entrywise envelopes under the
  assumptions recorded in the source receipts.
- Independent conditioned BA draws permit products of their conditional
  expected spectra.

The numerical values in this note remain diagnostics because the evaluator
uses nearest-binary64 arithmetic. A complete 40-bit claim also requires all
of the following steps.

1. Outward-round the conditioned one-active calculation. The \(B=240\) point
   has 8.019 diagnostic bits of slack.
2. Replace the dense worst-shell comparison by a finite weight-coupled
   transfer, then cover every \(65\le Q\le L\).
3. Give an efficient setup method for conditioning on \(\mathcal G_B\), or
   authenticate one selected BA spectrum and reuse that code.
4. Bind the final receipts and all frozen RM2Sub inputs in a certificate
   manifest.

The first-moment margin concerns the probability, over code setup, of
selecting a code with a low-weight nonzero word. It is not a decoding attack
cost or a cryptographic security level against an adversary.

## Active handoff

The implementation measured both outer layouts. The independent-row
transposed map took 10.407942 ms, 2.699% less than the same-binary B=256 BCH
baseline. It was only 1.411% slower than reuse. This result fixes independent
rows as the first proof interface because their spectrum expectations
factorize.

Occupations 1 through 64 are now outward-certified. The principal
mathematical task is a weight-coupled cover of every occupation from 65
through 8,832. Conditional setup and ordinary-encoder performance remain
open.
