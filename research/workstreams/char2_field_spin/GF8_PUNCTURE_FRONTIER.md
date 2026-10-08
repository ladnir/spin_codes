# Local puncturing near the GF256 rate frontier

The useful construction change is to reallocate redundancy, not enlarge the
state. Keeping six outputs of the existing width-eight/state-eight packet makes
outer GRS[128,80] a strong nominal-rate-5/6 candidate. A short public puncturing
period then gives finer rate increments without new field arithmetic or routing.

The isolated [frontier helper](theory/gf8_puncture_frontier.py) contains floating
screens only. It imports the existing local transfer unchanged and multiplies
different packet matrices in chronological order. The exact dense refinement
has its own [lemma and receipt](GF8_EXACT_SPECTRUM_DENSE.md). This note does not
replace the frozen [128,84]/keep-seven certificate.

## Six-output construction and proof audit

Use the same full-Cauchy maps F,A,C, eta, input diagonals, and state diagonals as
the certified seven-output construction. Retain the first six output rows instead
of seven. Initialize state to zero and discard terminal state. The route still
permutes each of the 128 outer-coordinate regions independently, exactly once.

For requested K=2^20 and outer dimension 80, the row count is 13,120. Thus
Nin=1,679,360, Nout=1,259,520, and the forbidden output threshold is 125,951.
The padded rate is 5/6; the shortened requested rate is smaller.

The selected Cauchy-rank arguments hold for any fixed subset of six top rows.
The birth and return minimum becomes max(0,7-w) for input support size w.
Both w=7 and w=8 can therefore emit zero. The transfer retains the z=6 return
term in both cases and the six-root live continuation at zero input.
It does not remove an all-zero retained packet merely because state remains live.

The source minimum row weight is 49. Sparse coverage through A=1028 joins middle
coverage from W=ceil(.03 Nin)=50,381 because 49*1029=50,421. Full row rank of the
retained inner follows from its causal six-by-six packet pivot minors, but does
not by itself establish injectivity on the outer code. The distance proof must
continue to include output zero.

The independent floating point screens give:

| Active outer rows A | Sparse margin, bits |
| ---: | ---: |
| 1 | 264.56 |
| 8 | 1732.59 |
| 16 | 1452.82 |
| 32 | 1581.36 |
| 64 | 2008.19 |
| 128 | 3131.12 |
| 256 | 5845.65 |
| 512 | 12512.55 |

All sampled middle densities .03, .05, .10, .14, .30, .50, .60, .80, .90,
and .98 are positive. The scores at .03 and .98 are 20,891 and 110,709 bits.
The old dense point bound starts passing near .792595. These point checks motivated
the separate directed sparse/middle/dense replays; they are not such replays.

## Public periodic puncturing

Let m_0,...,m_(r-1) be a public period, where packet t retains its first
m_(t mod r) outputs. Every packet still consumes eight symbols and updates all
eight state coordinates. The output length is

\[
N_{\rm out}=\frac{N_{\rm in}}{8r}\sum_{a=0}^{r-1}m_a.
\]

Require the regional packet count to be divisible by r. This aligns the period
in every region without resetting state. The natural row multiple 64 permits
periods 1, 2, 4, and 8 at width eight.

For marks x,y,tau, let M_m be the existing two-state transfer for m retained
rows. The correct period transfer is the ordered product

\[
M_{m_0}M_{m_1}\cdots M_{m_{r-1}}.
\]

In the sparse calculation, take the same ordered product of matrix polynomials
before raising it to the regional period count. Do not average packet matrices
or reset their live/zero state between phases. Each local bound has the same
conditional interface and continuation cone, so the positive comparison composes.
The source, route conditioning, and discarded terminal state are unchanged.

Rotating which rows are retained at a fixed count does not improve this proof:
the selected Cauchy ranks depend only on the count. Any benefit from such rotation
would require a new argument retaining more information than the present transfer.

Two narrow candidates emerged:

| Outer / puncture period | Nominal rate | Evidence from this exploration |
| --- | ---: | --- |
| [128,84] / (6,6,6,7) | 84.000% | A=1,64,128,256 and sampled middle gates pass; old dense cutoff about .838261 |
| [128,74] / (5,6) | 84.091% | Exact-spectrum dense component closes from .9 with 184 bits; separate agents own the sparse/middle replays |

For the first period, sparse margins at A=1,64,128,256 are 232.65, 1220.09,
1633.75, and 2905.45 bits. Middle margins at .03, .05, .10, .50, and .84 are
11,566, 26,954, 85,779, 734,854, and 464,771 bits. These are floating results.
For the second period, this helper independently found middle scores of 37,591
and 2500 bits at .98 and .99. The exact dense component avoids needing either
density as a join.

## The apparent plateau is not a counterexample

Put rho=Nout/Nin and R=Kpad/Nout. At full inverse-input density, the old dense
margin per raw symbol, before the finite conditioning factor, is exactly

\[
\rho\bigl(8(1-R)-B\bigr)-\log_2(256/255),\qquad B=1.27.
\]

The extra logarithm came from an upper bound on nonzero MDS support counts.
It is a proof loss, not an obstruction imposed by nonzero input diagonals.
The exact-spectrum positive-excess lemma removes it to negligible accuracy in
these full-support examples. It adds no construction randomness or work.

The remaining first-moment scale is approximately

\[
N_{\rm out}\bigl(8(1-R)-B\bigr)-128\log_2(L+1).
\]

For [128,81]/keep-six, R=84.375%, even the improved full-density point has
-26,694 bits with B=1.27. This is a failed proof bound, not a bad-code witness.
It is above the usual GF256 random-linear-code entropy benchmark

\[
1-\frac{h_2(.1)+.1\log_2(255)}8\approx .84144613.
\]

That benchmark is not a universal upper bound on code rate. It describes the
standard random-code first moment. Approaching it with the present proof also
requires accounting for finite conditioning and output-ball slack.

The exact rational Chernoff comparison accepts B=1.2685; it rejects B=1.268.
Using the former would recover about 1875 bits at the (5,6) K20 geometry.
The current complete dense receipt deliberately retains the old B=1.27.
At increasing N, its ball exponent alone approaches a rate ceiling of .84125
for this bound. Sharpening the exponent can narrow that small remaining gap.

## Performance consequences and next gate

All three candidates retain the checked factored packet arithmetic. The adjoint
sets omitted input coordinates to zero, evaluates the same maps, and carries
the full state in reverse chronology. No second global route is needed.
Arithmetic work per raw packet does not fall in proportion to retained outputs.

At fixed requested source size, work scales primarily with the outer row count.
The six-output/dimension-80 candidate uses about 4.59% more raw rows than the
frozen seven-output/dimension-84 encoder. Period (6,6,6,7) at dimension 84 uses
the same raw row count. Period (5,6) at dimension 74 uses about 13.27% more.
These are operation/traffic comparisons, not timing predictions.

The mixed layouts need explicit packed serialization. Four (5,6) packets emit
22 records, so successive four-packet groups alternate 64-byte alignment.
Eight packets emit 44 records, an integral eleven ZMM-sized groups. Packing
those eight packets together can preserve aligned local loads/stores.
For (6,6,6,7), four packets emit 25 records; sixteen packets emit 100 records.
The extra retained coordinate can be packed locally across those four groups.
Either choice needs independent chronology, zero-fill, alignment, and alias tests.

Larger packet/state geometries are not justified by the present evidence: the
existing 8/8 map already supports proof gates near the random-code benchmark.
Do not pay for a larger state merely to repair the old dense counting loss.
The next implementation gate should compare the certified fixed-six profile
against a certified periodic winner, using whole calls and the same addresses.
No native implementation or benchmark was performed in this exploration.
