# Finite \(k=2^{20}\) proof obligations

## Workstream status

**Status:** active on 2026-08-31. The implementation comparison selected the
independent-row B=240 interface for the first proof attempt.

The power-of-two EBCH32--ParityFanout--BA candidate now has a complete outward
first-moment certificate at relative distance 10%.  For
\(k=2^{20}\), \(N=2^{21}\), and \(d=209715\), it proves
\(\Pr[d_{\min}\le d]<2^{-51}\).  The dense range uses a 19-interval one-band
cover with 468.225918 bits of outward margin.  See
`FINITE_K20_D10_CERTIFICATE.md`. A six-attempt bounded rejection algorithm,
using exact exhaustive row tests, proves the unconditional event "setup aborts
or the accepted code has distance at most the target" with 45.323934 bits of
margin. The 11% target and practical setup remain open.

The current preferred point is

\[
  k=2^{20},\quad B=240,\quad L=8832,\quad
  N_+=2{,}119{,}680,\quad d=233{,}164.
\]

The revised-window outward occupation-one certificate proves more than 46 bits of
first-moment margin for that class. No complete finite 40-bit claim has been
proved because occupations 2 through 8,832 and conditional setup remain open.

The measured independent-row transposed map is 2.699% faster than the B=256
BCH baseline and 1.411% slower than reuse. Reuse remains a later
memory-optimized target; it is not the probability space of the first finite
certificate.

## Target statement

After the implementation choice is fixed, the finite theorem should state one
fully specified setup algorithm \(\mathsf{Setup}_{240}\). The algorithm must
output an encoder

\[
  E:\mathbb F_2^{2^{20}}\longrightarrow
  \mathbb F_2^{2{,}119{,}680}.
\]

Define

\[
  Z_d(E):=
  \bigl|\{x\in\mathbb F_2^{2^{20}}\setminus\{0\}:
  \operatorname{wt}(E(x))\le d\}\bigr|.
\]

The desired finite claim is

\[
  \Pr_{E\gets\mathsf{Setup}_{240}}[Z_{233164}(E)>0]
  \le 2^{-40}.
  \tag{1}
\]

Equation (1) implies \(d_{\min}(E)\ge233{,}165>0.11N_+\), except with the
stated setup probability. The 40-bit value is a setup-failure margin. It is
not a decoding attack cost.

The theorem must also state ordinary and transposed encoding costs for the
same setup distribution. An asymptotic \(O(N_+)\) statement is insufficient
for the concrete performance comparison.

## Obligation ledger

An item is closed only when a proof and a reproducible verifier or exact
derivation cover its stated probability space.

## Repeated BA with random Toeplitz convolution

**Finite distance status: closed for one fixed setup.**

`FINITE_K20_BA240_REPEATED_TOEPLITZ_CERTIFICATE.md` and
`ba3_B240_repeated_toeplitz_prefix_outward.json` freeze one repeated
Golay--BA-3 code and all routing orders. Exact prefix ranks and a 256-bit Arb
calculation prove

\[
 \Pr[d_{\min}\le233164]<2^{-180}
\]

for the full parent code over one shared random lower-triangular Toeplitz
kernel. The same statement therefore holds for the specified
\(2^{20}\)-dimensional zero-shortened subcode.

The remaining obligations are implementation obligations, not gaps in this
finite probability bound:

1. select an algorithm and complexity model for full-length binary Toeplitz
   multiplication;
2. implement and benchmark its ordinary and transposed maps; and
3. if strict linear time is required, prove that a bounded-state replacement
   retains a sufficient all-message prefix transfer.

No theorem here assigns a sampling probability to the fixed outer/routing
setup. The setup is explicit and authenticated instead.

The first RM2Sub bridge audit sharpens item 3. Let \(B\) be the selected
19-by-128 syndrome map. The subcode whose 16,560 epoch blocks all lie in
\(\ker B\) has shortened dimension at least 733,936. RM2Sub is the identity
on this subcode for every multiplier schedule. Thus an RM2Sub certificate
must establish both of the following properties:

1. the silent subcode has minimum distance above 233,164; and
2. the suffix-silent subcodes satisfy a weight-sensitive transfer bound that
   can be joined to the existing four-state RM2Sub recurrence.

Raw coordinate-prefix ranks do not imply either property. See
`FINITE_K20_RM2SUB_FIXED_OUTER_BRIDGE.md`.

### O1. Freeze the concrete construction

**Status: closed as a mathematical interface.**

`FINITE_K20_INDEPENDENT_SETUP.md` fixes independent conditioned BA rows, the
final 11,264 shortened parent-message coordinates, every route and region
permutation law, every nonzero RM2Sub multiplier, the zero initial state, the
unconstrained terminal state, and the no-failure conceptual rejection
sampler. Efficient implementation of that sampler remains O2.

Choose exactly one of the following outer interfaces.

**Independent-row interface.** For each \(j\in[L]\), sample a Golay--BA-3
code independently, condition that draw on \(\mathcal G_{240}^{23}\), and sample
an independent coordinate permutation.

**Reused-outer interface.** Fix one authenticated Golay--BA-3 code satisfying
a stated spectrum bound. Reuse that code in every row, but sample the row
coordinate permutations independently.

For either interface, the construction specification must also fix:

- the 11,264 shortened input coordinates;
- the distributions of the \(B\) region permutations;
- the distributions of all RM2Sub multipliers;
- the RM2Sub initial and terminal state conventions;
- the representation or generation method for every permutation; and
- the failure behavior of setup.

Changing any item above changes the probability space in (1).

### O2. Close the outer-code selection step

For the independent-row interface, prove or implement an exact sampler for
the conditional BA distribution. Rejection sampling is sufficient only if
setup can decide \(\mathcal G_{240}^{23}\). The proof must establish:

- correctness of the tail-free test;
- termination and expected setup cost;
- the conditional law of every accepted BA draw; and
- independence across accepted rows.

The current Markov calculation proves

\[
  \Pr[\mathcal G_{240}^{23}]\ge0.9936608994.
\]

It does not supply an efficient test for \(\mathcal G_{240}^{23}\).

For the reused-outer interface, produce one selected BA code and authenticate
a pointwise spectrum or transfer-weighted spectrum upper bound. A minimum-
distance statement alone is insufficient for the multi-active union bound.
If the proof instead averages over one sampled reused BA code, it must control
the required higher spectrum moments. Products of expected multiplicities do
not apply to a shared random code.

### O3. Outward-certify occupation one

**Status: closed, conditional on O6's source-interface audit.**

`certify_golay_ba_rm2sub_finite_one_active.py` and
`golay_ba3_rm2sub_finite_B240_q1_outward_k20_d11.json` cover every weight from
23 through 217 and prove

\[
  \mathbb E[Z_{233164,1}]
  \le \mathtt{0x1.6f3c66666f370p-47}<2^{-46}.
\]

The verifier uses positive one-sided binary64 arithmetic. Its diagnostic
input chooses legal Chernoff tilts only; no diagnostic bound is trusted.

Replace every nearest-binary64 step in the \(Q=1\) evaluator by exact or
outward-rounded arithmetic. The certificate must cover:

- the exact finite BA expected spectrum;
- the conditioning factor for \(\mathcal G_{240}^{23}\);
- every permitted outer weight from 23 through 217;
- the finite RM2Sub entrywise transfer;
- the selected Chernoff tilts;
- the sum over 8,832 possible active rows; and
- all integer rounding in the threshold \(d=233{,}164\).

A fixed rational tilt for each weight is sufficient. The proof need not show
that the recorded tilt is optimal. The final outward upper bound must be at
most \(2^{-40}\) after later occupation classes are added.

### O4. Outward-certify occupations 2 through 64

**Status: closed, conditional on O6's source-interface audit.**

`certify_golay_ba_rm2sub_finite_q2_64.py` and
`golay_ba3_rm2sub_finite_B240_q2_64_outward_k20_d11.json` cover every integer
occupation from 2 through 64. Their aggregate bound is

\[
  \mathbb E\!\left[\sum_{Q=2}^{64}Z_{233164,Q}\right]
  \le \mathtt{0x1.4dc4b8b530941p-1}\cdot2^{-84}<2^{-84}.
\]

The verifier retains the occupation, Bernoulli-envelope probability, and
Chernoff tilt. It uses a scaled mantissa/exponent representation for the
240-region products, preventing proof-path underflow and overflow.

The present binary64 receipt checks every integer \(Q\in[2,64]\) and reports
84.617 bits of aggregate margin under the revised window. Convert its witnesses to an outward
certificate. The verifier must retain the dependence on \(Q\), the candidate
Bernoulli probability, and the Chernoff tilt. It must sum all 63 occupations,
not only the displayed dominant rows.

### O5. Cover every occupation from 65 through \(L\)

**Status: open; exact reduction formulated.**

`FINITE_K20_TWO_TRACK_DENSE_PLAN.md` now records the active split. The exact
11% common-norm column-Hölder shortcut is numerically vacuous and has been
rejected. At 10.9%, all 116 sampled five-band compositions for \(Q=L\) close,
with at least 6,016 diagnostic bits. A complete barycentric simplex cover and
the occupations below \(L\) remain open.

`FINITE_K20_CONDITIONING_WINDOW.md` records the revised \([23,217]\) setup
law and the finite endpoint obstruction. `FINITE_K20_DENSE_TRANSFER_TARGET.md`
gives a finite weight-band mixture
identity that retains mixed BA row weights and reduces the remaining work to
a multitype RM2Sub transfer plus a positive composition-sum compression
lemma. Neither proposed compression route is proved yet.

This is the principal mathematical gap. Derive a finite weight-coupled
transfer that uses the complete conditioned BA spectrum. The proof must cover

\[
  65\le Q\le8832
\]

without sampled gaps. The current worst-shell envelope and the corrected
Holder diagnostic become vacuous at some dense occupations. They cannot be
used as certificate components.

Acceptable routes include:

- a finite coefficient transfer with one fugacity per BA weight band;
- a finite analogue of the positive-density joint exponent with explicit
  subexponential factors; or
- an exact selected-code spectrum combined with a sharper change of measure.

The proof must include the binomial choice of active rows and the full union
over \(Q\). It must not charge one endpoint likelihood to every active row.

### O6. Audit the finite RM2Sub interface

Bind the transfer calculation to the implemented RM2Sub-S19 recurrence. The
audit must verify:

- the selected degree-7 constituent and its authenticated spectrum;
- the zero-state activation bounds;
- the pre-add-multiply update order;
- the law and independence of field multipliers;
- the behavior of zero and live states;
- the first and last epoch conventions; and
- the correspondence between 16,560 epochs and \(N_+\) output bits.

Every imported receipt must have a verified hash. A modeled constituent or a
nearest-binary64 source receipt cannot support the final theorem silently.

### O7. Complete the first-moment argument

After O3--O5, sum the outward bounds for every \(Q\in[1,L]\). Prove

\[
  \mathbb E[Z_{233164}]\le2^{-40}.
\]

Apply Markov's inequality in the setup probability space fixed by O1. State
all conditioning explicitly. For the shortened code, include the deterministic
inclusion of its message space in the parent message space.

If setup can fail or restart, define whether (1) is conditioned on successful
setup or includes setup failure. Account for that choice in the final bound.

### O8. Produce the certificate package

The final package must contain:

- an outward verifier for each occupation regime;
- machine-readable receipts with rational or directed-rounding witnesses;
- hashes for every local and frozen dependency;
- a clean-checkout command sequence that reproduces every accepted receipt;
- an exact statement of interpreter and arithmetic-library requirements; and
- a manifest whose claim matches equation (1).

The existing binary64 receipts remain exploration inputs. Their status labels
must remain diagnostic.

### O9. Validate concrete performance

**Status: transposed online comparison closed; ordinary and setup costs
open.**

The implementation task measured 10.407942 ms for independent rows and
10.696591 ms for the same-binary B=256 baseline on one pinned Zen 4 core.
The external receipt SHA-256 is
`4b612e86c267a914b6cb5f6b231579f2ceb8dcfd7aa5981ee0b5807ce4fa01ee`.
It excludes setup, conditioning, allocation, destruction, and ordinary
encoding.

The completed comparison covers transposed online throughput, setup memory,
the baseline, hardware, compiler, and thread count. A complete performance
claim still requires:

- ordinary encoding throughput;
- setup time, including any tail-free test;
- memory bandwidth and cache behavior;
- the cost of stored tables versus generated permutations; and
- confirmation that a conditioned setup uses the same online representation.

Independent BA rows and a reused BA code have similar arithmetic counts but
different table sizes and access patterns. A benchmark of the reused variant
does not settle the independent-row proof model.

## Current continuation rule

O1, O3, and O4 are closed at the mathematical-transfer level. Continue with
O5 while preserving the independent-row conditional law. Do not promote a
40-bit theorem until O2, O5, O6, O7, and O8 are closed. Reconsider \(B=264\)
only if the dense cover consumes the seven-bit integer slack remaining after
the proved occupation-one bound.

## Historical 11% power-of-two fanout checkpoint (2026-09-01)

This checkpoint supersedes the B=240 route as the active finite candidate,
but it does not upgrade the earlier B=240 statements. The candidate has

\[
 (B,L,N,k,d)=(256,8192,2^{21},2^{20},230686)
\]

and uses eight genuine \([32,16,8]\) extended-BCH constituents, an independent
ParityFanout-31x33 draw per row, and two independently interleaved terminated
accumulators. `FINITE_K20_EBCH32_PARITYFANOUT_SETUP.md` defines the exact
probability space and accepted-row event \(G_{256}\), with permitted nonzero
row weights \([24,232]\).

The following statements are outward-certified:

- `ebch32_parityfanout31x33_ba3_B256_setup_outward.json` gives rejected
  expected-word count at most `0x1.362062013db23p-10`, acceptance probability
  at least `0x1.ff64efceff611p-1`, and expected trials per row at most
  `0x1.004d9f9ac6153p+0`.
- `ebch32_parityfanout31x33_ba3_rm2sub_B256_q1_outward_d11.json` covers every
  permitted row weight at \(Q=1\). Its expected bad-word count is at most
  `0x1.1ce93915da2eap-52`, for 51 certified integer bits of margin.
- `ebch32_parityfanout31x33_ba3_rm2sub_B256_q2_64_outward_d11.json` covers
  every integer \(Q=2,\ldots,64\). Its aggregate scaled bound has mantissa
  `0x1.585706232b360p-1` and exponent \(-95\), for 95 certified integer bits
  of margin.

The dense range is not certified. The resumable diagnostic receipt
`ebch32_parityfanout31x33_ba3_B256_three_band_cover_all_q_d11.json` has
processed 8,060 tetrahedra, accepted 3,957, left 149 pending, and rejected
none. The exact cell-local box rule is stated in
`FINITE_K20_DENSE_CELL_COVER.md`. An exact tree audit proves that the accepted
and pending leaves partition the explicit three-cell root triangulation.
The Arb verifier covers 3,716 nondelegated nonempty accepted cells and gives
an outward aggregate margin of 43.525570 bits; 234 accepted cells have empty
admissible lattice boxes.

Seven accepted geometric leaves delegate the external type
\((3344,1,21)\) to a local four-band refinement. A separate 17-interval Arb
certificate covers all 3,345 internal split counts with aggregate margin
198.029788 bits. The partition verifier checks that each delegated leaf
contains only that external type. These outward partial sums still omit the
149 pending cells and do not constitute a full dense bound.

The remaining ordered obligations for this route are:

1. finish the gap-free tetrahedral cover of all \(Q=65,\ldots,8192\);
2. re-evaluate every accepted cell with exact or outward interval arithmetic,
   including the Renyi band moments, RM2Sub matrix power, convexity premises,
   composition-count union bound, and the delegated split-low receipt;
3. sum the \(Q=1\), \(Q=2\ldots64\), and dense outward bounds and compare the
   sum, not the individual terms, with \(2^{-40}\);
4. audit the RM2Sub-S19 receipt-to-frozen-source interface and the region
   permutation law for this geometry;
5. implement and authenticate an efficient exact or one-sided-safe test for
   \(G_{256}\); and
6. implement the EBCH32--ParityFanout--BA outer in the frozen-style ordinary
   and transposed encoders and benchmark online encoding, setup, and memory.

Items 1--3 are the remaining distance-certificate gate. Items 4--6 are
independent theorem-to-implementation and linear-time gates. In particular,
the setup acceptance probability does not by itself provide an efficient
setup algorithm, and the frozen Structured SPIN benchmark measures a
different outer pipeline.

## Power-of-two 10% bounded-setup closure (2026-09-01)

This section supersedes the preceding ordered distance obligations at target
10%; it does not close the 11% route. Set

\[
 (B,L,N,k,d)=(256,8192,2^{21},2^{20},209715).
\]

The outward occupation receipts cover every \(Q=1,\ldots,8192\) and prove

\[
 \Pr[d_{\min}\le d\mid\text{all rows accepted}]<2^{-51}.
\]

For each row, the bounded setup samples at most six independent BA candidates.
It enumerates all \(2^{128}\) messages of each candidate and accepts exactly
when every nonzero output weight lies in \([24,232]\). It returns failure if
any row exhausts six candidates. The setup receipt and Markov's inequality
bound rejection of one candidate by

\[
 q\le\mathtt{0x1.362062013db23p-10}.
\]

First-success selection, conditioned on success, has exactly the independent
conditional row law used by the distance receipts. The exact-rational bounded
setup verifier proves

\[
 8192q^6+\Pr[d_{\min}\le d\mid\text{all rows accepted}]<2^{-45}.
\]

Thus the finite mathematical setup-and-distance claim is closed with more
than 40 bits of margin. The exhaustive test requires at most
\(6\cdot8192\cdot2^{128}\) row-message evaluations, so this result does not
close practical setup. The remaining production obligations are:

1. replace exhaustive enumeration with an efficient exact or one-sided-safe
   \(G_{256}\) test, or authenticate a precomputed accepted-row set;
2. audit the proof's RM2Sub-S19 and permutation interfaces against the frozen
   baseline in
   `constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/`;
3. implement the exact EBCH32--ParityFanout--BA row law in the optimized
   ordinary and transposed encoders; and
4. measure setup, memory, and both online encoder directions.

### Practical-setup diagnostics

Removing row conditioning does not retain the requested margin. At bad-weight
threshold \(209715\), the nearest-binary64 occupation-one calculation gives
only 11.089936 bits of margin. Lowering the relative threshold from 10% to 9%
raises that diagnostic margin only to 11.506321 bits. These values are not
certificates, but they rule out an outward translation with 40 bits of slack
under the same occupation-one bound.

Two exact-feasibility backends were tested on one deterministically generated
row. Z3 did not resolve the query for a nonzero word of weight at most 16 in
60 seconds. Gurobi did not resolve the query for a nonzero word of weight at
most 23 in 60 seconds. Both results are `unknown`; neither is evidence that a
forbidden word exists. They show only that the current generic formulations
are not yet suitable for an 8,192-row setup test.

The practical setup work should next exploit the BA trellis and accumulator
structure, produce checkable selected-row certificates offline, or replace the
conditional row law with a structured outer whose admissibility is easy to
verify. Any replacement changes the probability space and requires the dense
occupation receipts to be rechecked.

## Repeated EBCH128 at 11%

For the candidate that repeats one fixed \([128,64,22]\) extended BCH code
over \(L=16384\) outer positions, the following items remain open.

1. Prove a high-occupation transfer for every \(101\le Q\le16384\). It must
   retain all 128 transposed regions and preserve the even-row dependence.
2. Cover mixtures of non-all-one BCH words and the unique all-one BCH word.
3. Outward-certify the existing occupation-one and occupation-2-through-100
   witnesses and sum them with the high-occupation contribution.
4. Bind the published BCH spectrum and the frozen RM2Sub-S19 receipts by
   hash in a final manifest.

The current low-occupation margin is 26.335 bits in nearest-binary64
arithmetic. It is not a complete certificate. The 127-region parity-drop
diagnostic fails at the all-active endpoint, and the exact-marginal
column-Hölder alternative is too loose. See
`FINITE_K20_REPEATED_EBCH128_D11_STATUS.md`.

The RandomStepConv-M30 comparison is closed at
\((L,N,D)=(16560,2119680,233165)\). The outward verifier proves a complete
failure bound below \(2^{-26.1921836158}\). It disperses the uniform-even
parity pivots, retains their multinomial load law, and uses two positive
coefficient bounds. The exact spectrum covers occupation one. An exact
aligned-pivot recurrence covers occupations 2 through 99. The dispersed
transfer covers occupations 100 through 16384.

No finite-distance obligation remains for this RandomStepConv model. The
remaining design obligation is an efficient structured inner with a proved
comparison to this transfer. See
`FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_CERTIFICATE.md`.

## Shortened BCH250 at 11%

The fixed repeated outer is the explicit six-coordinate shortening of the
extended ([256,131,\ge38]) BCH parent. It has exact parameters
([250,125,\ge38]). Use (L=8448), (N=2{,}112{,}000), and
(d=232{,}320). The full-rate analysis covers (1{,}056{,}000) message bits;
zero padding embeds exactly (2^{20}) input bits into that instance.

The exact parameter audit and constant-weight packing envelope are closed.
The nearest-binary64 occupation-one and occupation-2-through-100 ledgers have
13.777 and 19.044 bits of margin, respectively. The ordered remaining
obligations are:

1. prove a high-occupation coefficient transfer for every
   (101\le Q\le8448), preserving the row-weight mixture across all 250
   regions;
2. use both the exact total row mass (2^{125}-1) and the shell packing
   bounds, rather than charging one worst shell independently to every row;
3. outward-certify the low-occupation and high-occupation ledgers and sum
   them;
4. audit the permutation and RM2Sub probability space against the selected
   implementation; and
5. implement and benchmark the ordinary and transposed encoders for the exact
   shortened constituent.

The current pointwise dense envelope fails by 446,724 bits at (Q=L). An
exact-mass correction alone leaves about 32,500 bits of failure. This is a
proof-architecture gap, not a rounding item. See
`FINITE_BCH256_SHORTENING_OPTIONS.md`.

## One-sampled parity fanout on BCH250

For the variant (C_F=F(C_{250})), sample one composition of fanout maps and
reuse (C_F) in every outer row. The one-word expected transition and the
mass-coupled shell relaxation are available. They do not control products of
the sampled spectrum.

The original 64-layer three-band target is rejected. Its broad tail vertices
fail, and narrowing the bands leaves thin mixed faces with negative margins.
The worst tested 64-layer edge has two endpoint rows and 8,446 central rows;
its diagnostic margin is -228.861 bits. This is not a rounding issue.

The revised target uses 128 layers. Its expected endpoint mass on weights
1 through 16 and 234 through 249 is at most approximately
(2^{-41.40363381}). Its expected pure-bulk all-active calculation retains
177.390 bits. These two numbers do not form a certificate because the
weights adjacent to the endpoint create thin mixed faces and because an
expected spectrum is not a spectrum bound for one reused code.

The ordered obligations are:

1. outward-certify the endpoint expectation and its Markov setup-failure
   consequence;
2. build a defect-row transfer for weights 17 through 36 and 214 through
   233 that preserves their uniform-shell law across all 250 regions;
3. prove a simultaneous bulk-spectrum comparison on weights 37 through 213
   with at most 0.01 bits of pointwise half-Bernoulli slack per row and an
   explicit setup-failure probability;
4. cover every defect count and every occupation for the resulting fixed
   code, then sum all occupation and setup-failure terms;
5. convert every transition, spectrum, and RM2Sub calculation to outward
   arithmetic;
6. audit the construction probability space against the frozen baseline;
   and
7. fuse the sampled composition with the BCH encoder and benchmark the
   ordinary and transposed circuits.

The 128-layer scalar proxy is 8,064 XORs per row before circuit fusion. The
current calculation is a proof-model diagnostic, not a performance claim.

## BCH250-124 with independent row-local fanout

The current strongest finite 11% route uses one fixed
$[250,124,\ge 38]$ subcode of the audited BCH250 constituent in all 8,576
rows. It samples an independent 56-layer ParityFanout-31x33 wrapper in each
row. The parent has dimension 1,063,424 and length 2,144,000. Fixing 14,848
parent input coordinates gives exactly $2^{20}$ message bits.

The optimized proof-model certificate is complete. A shell-sensitive Arb
receipt bounds the probability of any wrapped-row word outside weights 13
through 237 by $2^{-45}$. On the central shells, it proves a pointwise
density excess below $389/1250$ bits. The outward occupation bounds are
$2^{-53}$ for $Q=1$, $2^{-52}$ for $2\le Q\le31$, and $2^{-41}$ for
$32\le Q\le8576$ on the central event. Their total is below $2^{-40}$.

The remaining obligations concern proof-compatible setup and equivalence:

1. replace the benchmark schedule generator with a proof-compatible sampler
   or authenticate a table produced by such a sampler;
2. implement the ordinary forward encoder;
3. prove ordinary, transposed, and proof-model equivalence; and
4. repeat the integrated benchmark on the target implementation host.

The unfused fanout proxy is 3,528 scalar XORs per row and 30,256,128 over
all rows. A 52-layer candidate fails the present tail-plus-central proof.
This failure does not rule out a stronger multi-band argument.

The packed 56-layer action is implemented in
`benchmark_bch250_rowlocal_fanout56.cpp`. A pinned i7-13700H measurement gives
2.177353 ms for 56 layers and 0.006059 ms for the post-BCH zero-layer copy.
The isolated increment is 2.171294 ms. See
`FINITE_K20_BCH250_FANOUT56_IMPLEMENTATION_BENCHMARK.md`. This wrapper-only
measurement does not include BCH encoding, coordinate permutations,
transpose, region permutations, or RM2Sub-S19.

The replacement proxy modifies the frozen 256-by-8192 packed, tiled, fused
path. On a pinned Ryzen 9 7950X processor, zero layers take 10.650813 ms and
56 layers take 20.454284 ms. The difference of medians is 9.803471 ms, and
the ratio is 1.920443. Staged zero-layer and 56-layer checks pass.

The replacement remains a BCH256 performance proxy. It does not discharge
the exact BCH250 implementation, forward-encoder, setup-sampling, or
equivalence obligations. See
`FINITE_K20_BCH256_FANOUT56_PROXY_BENCHMARK.md`.

See
`FINITE_K20_BCH250_124_ROWLOCAL_FANOUT_STATUS.md`.
