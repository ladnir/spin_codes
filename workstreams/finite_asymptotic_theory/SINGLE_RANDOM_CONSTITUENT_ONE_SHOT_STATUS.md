# One-shot random constituent: finite status

## Closure update

The one-shot (B=512) route is now an outward finite certificate. A single
uniform (256)-by-(512) generator is sampled and repeated in all 4,096
outer rows. Setup performs no acceptance test. At RandomStepConv-M22 and
(D=\lceil0.109N\rceil=228{,}590), the combined verifier proves

\[
 \Pr[d_{\min}<D]<2^{-42.2593129905}.
\]

The spectrum-event complement contributes at most (2^{-42.2614729204}).
Conditional on that event, all occupations contribute at most
(2^{-51.6439589890}). The theorem is stated in
`FINITE_K20_ONE_SHOT_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md`.

The remainder of this note records the search path. Statements below that
describe the shared-category calculation as open are historical.

## Compact pairwise-systematic replacement

The uniform 256-by-512 matrix is no longer needed for the outer comparator.
Identify each message row with \(u\in\mathbb F_{2^{256}}\), sample one
independent pair \(a,b\gets\mathbb F_{2^{256}}\), and encode

\[
 u\longmapsto (u,au+bu^2).
\]

For distinct nonzero \(u,v\), the determinant of the coefficient matrix for
their parity halves is \(uv(u+v)\ne0\). Their parity halves are therefore
independent and uniform. This gives
\(\operatorname{Var}(A_w)\le\mathbb E[A_w]\). The systematic half makes
every sampled constituent injective.

The original integer shell caps remain valid. Reusing the outward sparse and
dense cap-conditional transfers proves the same 42.2593129905-bit overall
margin. The new outer sampler consumes 512 unbiased bits rather than
131,072. Evaluation costs two \(\mathbb F_{2^{256}}\) multiplications per
row. See
`FINITE_K20_PAIRWISE_SYSTEMATIC_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md`.

## Purpose

This note asks whether setup can sample one random rate-half constituent and
reuse it in every outer row, without enumerating its codewords or resampling
until a spectrum test passes. The target is

\[
 N=2^{21},\qquad D=\lceil0.109N\rceil=228590,
\]

with the region-permuted bit transpose and RandomStepConv-M22. The
calculations below began as nearest-binary64 diagnostics. The closure update
identifies their outward replacements.

## Setup distribution and spectrum event

Let (B\in\{512,1024\}), let (K=B/2), and let (L=N/B). Setup samples a
uniform binary (K\)-by-(B) matrix (G). The same map (u\mapsto uG) is
used in every outer row. Define

\[
 A_w(G)=\left|\{u\in\mathbb F_2^K\setminus\{0\}:
                    \operatorname{wt}(uG)=w\}\right|,
 \qquad
 \mu_w=(2^K-1)\binom Bw2^{-B}.
\]

For each shell, choose an integer cap (T_w). If (mu_w\le\delta), set
(T_w=0) and use Markov's inequality. Otherwise choose the integer cap used
by the receipt so that Cantelli's inequality gives a shell-failure upper
bound at most (delta). Distinct nonzero binary messages are linearly
independent, so their images under an unconditioned uniform generator are
pairwise independent. Hence

\[
 \operatorname{Var}(A_w)\le\mu_w.
\]

The event (mathcal E_B) requires full row rank and (A_w(G)\le T_w) for
every (w). This event is used only in the proof. Setup samples one uniform
matrix and does not test the event, so its complement is charged directly to
the theorem's failure probability.

For (B=512), (delta=2^{-51}) gives

\[
 \Pr[\neg\mathcal E_{512}]
 \le1.8968325343082488\mathbin\cdot10^{-13}<2^{-42.26147}.
\]

The nonzero caps have support (42\le w\le470). For (B=1024),
(delta=2^{-52}) gives

\[
 \Pr[\neg\mathcal E_{1024}]
 \le1.8397365882919309\mathbin\cdot10^{-13}<2^{-42.3050},
\]

and the nonzero caps have support (98\le w\le926). Thus an ordinary
uniform generator already satisfies the stated spectrum event with
probability greater than (1-2^{-40}). No biased sampler or rejection test
is needed for this event. The compact pairwise-systematic replacement above
provides the same guarantee with a 512-bit outer sampler.

## Conditional distance decomposition

For fixed (G\in\mathcal E_B), let (Z_Q) count nonzero messages that use
exactly (Q) outer rows and whose encoded word has weight below (D). The
desired conditional statement is

\[
 \sum_{Q=1}^{L}\mathbb E[Z_Q\mid G]<2^{-40.34}. \tag{1}
\]

The expectation in (1) is over the independent row-coordinate
permutations, region permutations, and RandomStepConv maps. Combining (1)
with the event-complement probability would give an overall setup-or-distance
failure probability below (2^{-40}).

The cap event makes the proof independent of joint weight-spectrum facts.
For (Q=1), the number of weight-(w) local choices is at most (T_w). For
(Q=2), the number of ordered image pairs of weights (a,b) is at most
(T_aT_b). For larger (Q), products of the caps give a valid fixed-code
counting measure.

## Results at (B=512)

All entries in this table use (M=22).

| Occupations | Method | Aggregate margin |
|---|---|---:|
| (Q=1) | exact shell transfer | 75.3684 bits |
| (Q=2) | exact shell-pair transfer | 113.0370 bits |
| (Q=3,\ldots,8) | exact three-band region recurrence | 54.6796 bits |
| (Q=9,\ldots,16) | exact three-band region recurrence | 175.1893 bits |

The three bands are

\[
 [42,79],\qquad[80,432],\qquad[433,470].
\]

Their Bernoulli majorants have parameters approximately

\[
 (p,\log_2c)=(0.1542968751,61.8430),
 (1/2,256.0675),(0.8457031249,61.8430).
\]

Every exact sparse slice is dominated by the pure low-band composition. Its
margin increases from 54.6796 bits at (Q=3) to 175.1893 bits at (Q=9).

At (Q=4096), the sampled three-band composition cover has 5035.5203 bits
of margin. Sampled compositions at (Q=128,256,512,1024,2048) also pass.
These sampled slices are evidence, not a cover of all integer compositions.

## Failed dense relaxations

The caps imply a global Bernoulli-half envelope only with factor
(F=2^{50.4033453111}). Charging (F^Q) fails even at (Q=L). This is a
weakness of the envelope, not evidence of a low-distance code.

The categorical three-band reference independently conditions on the same
row-category composition in each of the (B) regions. It therefore pays the
conditioning loss (B) times. This relaxation fails at intermediate
occupations. At (Q=32), its worst listed composition has 2647.4490
negative bits of margin. Increasing the inner memory from 22 to 30 does not
remove the loss because the reference bad-event bound has already saturated
for the offending mixtures. Rebalancing the bands also does not fix it.

## Next-power-of-two check

At (B=1024,K=512,L=2048), the one-shot method gives 210.7876 bits for
(Q=1) and 96.6856 aggregate bits for the coarse three-band (Q=2,3)
calculation. The bands are

\[
 [98,159],\qquad[160,864],\qquad[865,926].
\]

The categorical (Q=32) relaxation still fails on a narrow mixed-tail
strip. Its worst listed composition is ((27,4,1)). Thus increasing the
constituent from 512 to 1024 bits improves the true sparse margins but does
not cure the proof artifact.

## Resolved shared-category obligation

The former obligation was a shared-type mixed-band transfer bound for the
middle occupations. For one composition ((q_-,q_0,q_+)), the row category
is sampled once and is shared across all (B) regions. A sufficient lemma
must retain that shared category while averaging the independent region
permutations. It must not replace the shared category by (B) independent
categorical draws.

The closing proof merges the symmetric tails by a monotonicity lemma for the
RandomStepConv moment. An exact two-category recurrence covers
(3\le Q\le159). An outward convex cell cover handles
(160\le Q\le4096). Exact shell transfers retain (Q=1,2).

## Evidence

The cap construction is implemented by
`evaluate_single_random_constituent_highprob_renyi.py`. The exact sparse
programs are `evaluate_single_random_constituent_highprob_q1.py`,
`evaluate_single_random_constituent_highprob_q2.py`, and
`evaluate_single_random_constituent_highprob_sparse_bands.py`.
`evaluate_single_random_constituent_highprob_bands.py` implements the
categorical diagnostic, and
`evaluate_single_random_constituent_highprob_uniform.py` records the rejected
global-envelope test.

The principal (B=512) receipts are
`single_random_constituent_B512_highprob_q1_s22.json`,
`single_random_constituent_B512_highprob_q2_s22.json`,
`single_random_constituent_B512_highprob_sparse_q3_8_s22.json`,
`single_random_constituent_B512_highprob_sparse_q9_16_s22.json`, and
`single_random_constituent_B512_highprob_band_q4096_s22.json`.

The (B=1024) fallback receipts are
`single_random_constituent_B1024_highprob_q1_s22.json`,
`single_random_constituent_B1024_highprob_sparse_q2_3_s22.json`, and
`single_random_constituent_B1024_highprob_band_q32_s22.json`.
