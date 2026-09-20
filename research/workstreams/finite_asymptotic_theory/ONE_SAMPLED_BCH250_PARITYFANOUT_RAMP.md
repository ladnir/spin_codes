# One-sampled parity fanout for the shortened BCH250 outer

## Question

The shortened BCH code has exact parameters ([250,125,\ge38]), but its
weight distribution is unknown. Empirical evidence from the smaller
extended BCH code suggests a random-like spectrum. This note asks whether a
small random linear wrapper can make that comparison provable without an
exact enumerator.

The wrapper must be sampled once and reused for every outer row. Sampling a
fresh transformed code per row is outside the selected construction.

## Setup experiment

Fix the shortened BCH code (C\subseteq\mathbb F_2^{250}). For disjoint
sets (S,T\subseteq[250]), define

\[
 F_{S,T}(x)
 :=x+\left(\sum_{i\in S}x_i\right)1_T.
 \tag{1}
\]

The current diagnostic uses ((|S|,|T|)=(31,33)). Because (S\cap T) is
empty, (F_{S,T}) is an involution. It is therefore an invertible linear
map.

For a layer count (r\), sample independent pairs

\[
 (S_j,T_j),\qquad j=1,\ldots,r,
\]

from the disjoint-set law. Define

\[
 F:=F_{S_r,T_r}\circ\cdots\circ F_{S_1,T_1},
 \qquad C_F:=F(C).
 \tag{2}
\]

The setup publishes one realization of (F). Every outer row uses the same
code (C_F). Independent coordinate permutations may follow the fixed
encoder; those permutations do not change its spectrum.

## Exact one-word transition

Fix a word of weight (w). Let (a) be its overlap with (S). If (a) is
even, one fanout layer leaves the word unchanged. If (a) is odd, let (b)
be the overlap of its remaining support with (T). The output weight is

\[
 h=w+33-2b.
 \tag{3}
\]

Both overlap laws are hypergeometric. Thus (3) defines an exact stochastic
matrix (P) on the weights (0,\ldots,250). The expected spectrum after
(r) independent layers is

\[
 \mathbb E_F[A(C_F)]=A(C)P^r.
 \tag{4}
\]

Equation (4) is an expectation over the one sampled composition. It is not a
concentration statement.

## Source spectrum polytope

The source spectrum satisfies four exact constraints.

1. (A_0=1).
2. Every nonzero source weight is even and lies in ([38,218]).
3. The nonzero mass is

   \[
     \sum_{w>0}A_w=2^{125}-1.
     \tag{5}
   \]
4. Each shell satisfies the constant-weight packing bound in
   `FINITE_BCH256_SHORTENING_OPTIONS.md`.

For each output weight (h), maximizing

\[
 \sum_w A_w(P^r)_{w,h}
 \tag{6}
\]

over these constraints is a continuous linear program. The optimizer fills
the source shells in decreasing order of ((P^r)_{w,h}). This greedy rule is
exact for the relaxation.

`analyze_one_sampled_bch250_parityfanout.py` implements (3)--(6). It uses
nearest-binary64 transition probabilities. Different output shells use
different maximizers, so the resulting output envelope need not be jointly
realizable.

## Numerical results

At (L=8448), (N=2{,}112{,}000), and relative distance 11%, the exact
dense RM2Sub reference has 177.405 diagnostic bits of all-active margin. A
uniform per-row comparison may therefore spend at most

\[
 \frac{177.405}{8448}=0.02100
 \tag{7}
\]

bits per row before consuming the complete all-active margin.

One ParityFanout-31x33 layer does not erase the source uncertainty. The worst
pointwise output gap remains 48.032 bits. Its one-active RM2Sub calculation
still closes with 11.239 bits because occupation one can afford that loss.

Layering changes the conclusion. With 64 independent layers, the
mass-coupled expected envelope gives:

- 52.698 bits of occupation-one margin;
- 52.713 bits under the random-even source model;
- 0.012699 central excess bits per row on weights 41 through 209;
- 0.009629 central excess bits per row on weights 61 through 189; and
- 0.003587 central excess bits per row on weights 101 through 149.

Thus the expected envelope is nearly random-like in the central weights. On
weights 61 through 189, its uniform charge is approximately 81.35 bits at
full occupation. The ideal dense reference retains about 96.06 bits before
the two tails and finite prefactors are added.

The scalar XOR proxy for (r) ParityFanout-31x33 layers is

\[
 r(31+33-1).
 \tag{8}
\]

For (r=64), (8) equals 4,032 XORs per scalar row. A coarse sweep over eight
fanout shapes found no materially cheaper mixing point. Smaller target sets
need more layers; larger sets need fewer but more expensive layers. This
proxy does not include circuit fusion with the BCH encoder.

### Failure of the first band relaxation

The first partition was

\[
 [1,60],\qquad[61,189],\qquad[190,249].
 \tag{9}
\]

The central vertex closes with 96.057 diagnostic bits. The two tail vertices
do not: their margins are approximately -145,183 and -142,922 bits. Narrowing
the endpoint bands to weights 1 through 36 and 214 through 249 makes all
three vertices close, and the central vertex retains 64.866 bits. This
vertex check is insufficient. With two low-tail rows and 8,446 central rows,
the same relaxation has margin -228.861 bits. The symmetric high-tail case
has margin -228.813 bits. The loss is the 250-fold categorical-conditioning
charge on a thin mixed face.

The edge failure persists for every tested symmetric cutoff from 12 through
60. Splitting weights 17 through 36 into singleton bands does not remove it.
For example, two weight-36 rows and 8,446 central rows still have margin
-116.335 bits after 128 layers. The obstruction is therefore the transfer,
not the width of the outer bands.

### The 128-layer random-model target

After 128 layers, the expected envelope is essentially random in the bulk.
For the endpoint set

\[
 E:=\{1,\ldots,16\}\cup\{234,\ldots,249\},
 \tag{10}
\]

the source-polytope linear program gives

\[
 \log_2\max_A\mathbb E_F[A_E(C_F)]
 \approx-41.40363381.
 \tag{11}
\]

Thus Markov's inequality gives the same binary64 exponent for the event that
the sampled code contains any word in (E). This is a statement about setup
failure, not a deterministic minimum-distance audit.

On the complementary central band ([17,233]), the expected pointwise
half-Bernoulli envelope costs about (1.87\cdot10^{-6}) bits per row. Its pure
all-active diagnostic retains 177.390 bits. Combining only that pure
conditional calculation with (11) would leave about 41.404 bits. That
combination is not yet valid: a sampled code can contain a small number of
words of weights 17 through 36 or 214 through 233, and those thin mixed
faces are exactly where the current conditioning transfer fails.

## Why the result is not yet a one-sampled theorem

For occupation one, the bad-word count is linear in the sampled spectrum.
Equation (4) therefore suffices after outward verification.

For (Q\ge2), the same code (C_F) appears in every active row. A
pointwise first moment does not control products such as

\[
 A_{w_1}(C_F)\cdots A_{w_Q}(C_F).
 \tag{12}
\]

The expectation of (9) depends on high moments of the sampled spectrum. In
equivalent message-space terms, it depends on the rank and joint support type
of the selected codewords. Replacing (9) by a product of expected spectra
would silently change the setup to independent fanout draws per row.

Consequently, the present receipts prove no reusable-code distance claim.
They establish that the fanout walk has the correct one-word mixing target
and that the required arithmetic margin is plausible.

## Why a scalar random-spectrum slack is insufficient

Suppose a fixed constituent satisfies the pointwise comparison

\[
 A_w\le 2^{125+\varepsilon}\binom{250}{w}2^{-250}
 \tag{13}
\]

on every retained shell. Applying (13) independently in all 8,448 rows
costs (8448\varepsilon) bits. The ideal all-active calculation has only
177.405 bits of margin. A 40-bit target therefore requires

\[
 \varepsilon < \frac{177.405-40}{8448}\approx0.01626
 \tag{14}
\]

before finite composition and occupation costs. This explains why a code
can look very close to random and still defeat a uniform comparison. For the
known EBCH128 spectrum, the best body comparison costs 0.067435 bits per
row; repeated 16,384 times, that is about 1,105 bits.

A random model with slack remains viable, but the slack must be
weight-dependent and probabilistic. Rare shells need discrete count or
setup-failure bounds. Bulk shells need a simultaneous concentration theorem
with a relative slack below the budget in (14). The setup randomness is paid
once. It must not be charged as an independent spectrum draw in every row.

## Closure routes

Two routes preserve one sampled code.

### Fanout-walk concentration

Prove that the random composition in (2) maps every admissible source
spectrum into a good three-band profile with high probability. The theorem
must control a transfer-weighted functional or the relevant band masses. A
union bound over independent shell expectations is too expensive.

The proof may instead analyze the fanout walk on subspaces of
\(\mathbb F_2^{250}\). A sufficiently mixed image of the fixed
125-dimensional code behaves as one sampled random linear code. The current
weight-chain calculation is necessary evidence for this approach, but it
does not prove subspace mixing.

### Statistical setup audit

Sample one fanout composition and then sample uniform messages from its fixed
code. Use an exact binomial or multinomial acceptance test for a small set of
band probabilities. The false-accept probability becomes part of setup
failure. Conditioned on acceptance, reuse the same code for every row.

This route requires a dense certificate stated only in terms of the audited
bands. It must not require the complete shell spectrum. The audit cost and
sample count remain open.

## Next theorem target

Use 128 layers and separate three logically different tasks.

1. Treat (E) from (10) as setup failure and outward-certify (11).
2. Treat weights 17 through 36 and 214 through 233 as defect shells. For a
   fixed number of defect rows, preserve their exact uniform-shell law across
   all 250 regions. Do not replace a rare defect row by a categorical
   Bernoulli type independently in every region.
3. Prove a simultaneous bulk comparison of the form (13) on weights 37
   through 213 with an explicit setup-failure probability and
   (\varepsilon\le0.01) bits per row.

The fixed-code distance theorem must then cover every defect count and every
occupation. The one-sampled setup theorem must bound the probability that
the fixed transformed BCH code violates either the defect-count bounds or
the bulk comparison. A sampled BA constituent can use the same interface:
only the proof of the spectrum event changes.

## Reproducible artifacts

- `analyze_one_sampled_bch250_parityfanout.py` constructs the exact-form
  one-word transition and the mass-coupled expected shell envelope.
- `sweep_bch250_parityfanout_mixing.py` produces
  `bch250_parityfanout_mixing_sweep.json`, including the 64- and 128-layer
  endpoint calculations.
- `diagnose_bch250_parityfanout_three_band_qL.py` implements the band and
  edge diagnostics.
- `bch250_parityfanout31x33_l64_three_band_qL_d11_diagnostic.json` records
  the broad-tail failure.
- `bch250_parityfanout31x33_l64_cutoff36_edges_qL_d11_diagnostic.json`
  records the thin-face failure after narrowing the tails.
- `bch250_parityfanout31x33_l128_cutoff16_pure_qL_d11_diagnostic.json` and
  `bch250_parityfanout31x33_l128_cutoff16_edges_qL_d11_diagnostic.json`
  separate the pure-bulk promise from the one-defect obstruction.
- `bch250_parityfanout31x33_l128_sparse17_36_central_spokes_qL_d11_diagnostic.json`
  records the failure of singleton tail splitting under the same transfer.

All numerical receipts use nearest-binary64 optimization. The final transfer
evaluation uses 100 decimal digits where the diagnostic script says so, but
it is not outward rounded.

## Status

- **Exact algebra:** invertibility of each fanout, the hypergeometric
  transition, source mass, and the continuous shell optimizer.
- **Rigorous combinatorics with binary64 display:** the source packing caps.
- **Diagnostic failure:** the broad three-band transfer and all tested
  symmetric-cutoff variants on thin mixed faces.
- **Diagnostic target:** the 128-layer endpoint exponent, pure bulk margin,
  and proposed defect/bulk split.
- **Open:** defect-row transfer, one-sample bulk concentration, all
  occupations, outward arithmetic, implementation, and performance.
