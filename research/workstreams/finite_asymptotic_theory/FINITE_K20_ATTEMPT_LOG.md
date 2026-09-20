# Finite \(k=2^{20}\) attempt and decision log

## Purpose

This document records the finite-distance routes investigated before any
further construction change. It is a decision log, not a theorem statement.
It separates certified results, conditional arguments, numerical diagnostics,
performance measurements, and failed proof techniques.

The present design goal is one linear-time encoder with all of the following
properties.

1. The message dimension is \(k=2^{20}\).
2. The minimum-distance target is strictly greater than \(0.11N\).
3. The construction-level failure probability is at most \(2^{-40}\).
4. One fixed outer constituent is selected once and repeated at every outer
   position. Row-dependent outer codes or row-dependent wrappers are not the
   desired final interface.
5. Setup is practical. In particular, setup may not enumerate a
   \(2^{128}\)-element message space.
6. Ordinary and transposed encoding have linear work. The current practical
   target for the online transposed map is approximately 11 ms on the Peach
   Ryzen 9 7950X reference system.

Some investigations used a 10%, 10.9%, or asymptotic target. Such a result is
recorded at its actual target and is not promoted to the current 11% target.

## Status vocabulary

- **Certified:** an outward-rounded verifier and its stated probability space
  prove the displayed bound.
- **Proved asymptotically:** a theorem holds along the stated parameter
  schedule as \(N\) tends to infinity. This does not imply a finite
  \(k=2^{20}\) certificate.
- **Conditional:** the implication is proved, but an explicit hypothesis such
  as a constituent spectrum bound remains unauthenticated.
- **Diagnostic:** nearest-binary64 or sampled numerical evidence. It guides
  the search but proves no probability bound.
- **Measured:** a benchmark result. It proves no distance claim.
- **Failed argument:** the tested upper bound does not close. This does not
  prove that the underlying code has small distance.
- **Rejected interface:** the mathematical construction does not satisfy the
  desired single-repeated-constituent or practical-setup requirement.

## Executive ledger

| Route | Strongest result | What failed or remains open | Present disposition |
|---|---|---|---|
| Random outer and random convolution | Clean random-ensemble baseline | Not an explicit structured instance | Reference model only |
| Random outer plus bitwise RandomStepConv | Exact two-state all-message bound; finite 11% diagnostic closes at \(M=10\) with 1,422.7 bits | Binary64 rather than outward arithmetic; full random outer | Strong bounded-memory reference model |
| One repeated random \([512,256]\) outer plus RandomStepConv-M30 | **Certified** literal 11% distance; conditional failure below \(2^{-108}\); 28-attempt bounded setup closes 40 bits | Exact spectrum test and random step maps are impractical | Complete mathematical comparator; redesign target |
| One repeated Golay--BA-3 plus bitwise RandomStepConv | Conditioned \(Q=1\) diagnostic closes with 41.6332 bits at \(M=23\) | Occupations \(Q\geq2\), setup selection, and outward arithmetic remain open | Active finite baseline |
| Random outer and fixed RM2Sub-S19 | Proved asymptotic 11% distance with \(B=9\log_2N+O(1)\) | Finite structured outer not supplied | Reference theorem |
| Random block plus one accumulator | Proved asymptotic linear distance for logarithmic blocks | Uniform interleaver is not the desired structured permutation | Reference theorem |
| Bit-transpose accumulator variants | Exact obstruction and region-shuffle transfer formulas | Pure bit transpose retains damaging adjacency; no finite target certificate | Do not use pure transpose |
| One repeated Golay--BA-3 plus RM2Sub-S19 | Proved asymptotic 11%; finite \(B=240\) low occupations close | Finite dense occupations did not close at 11% or 10% with the attempted bounds | Preserve; revisit only with a new mechanism or proof idea |
| One repeated Golay--BA-3 plus random Toeplitz convolution | Frozen \(B=240\) setup; complete outward all-message 11% certificate with 180 claimed bits | Toeplitz complexity and performance remain open | Finite proof model closed; not a linear-time construction |
| Independent-row Golay--BA-3 | Finite \(Q=1\) and \(2\le Q\le64\) close; measured 10.407942 ms | Dense range and practical conditional setup open; wrong final interface | Diagnostic implementation only |
| EBCH32--ParityFanout--BA | Complete finite 10% first-moment certificate | Bounded setup requires an impractical exact test; 11% dense cover incomplete; independent rows | Mathematical checkpoint, not a deployable construction |
| Repeated exact EBCH128 plus RM2Sub-S19 | Exact low-occupation diagnostics: about 26.335 bits combined | Dense all-one/mixed-tail transfer open; no 40-bit or complete 11% certificate | Stopped |
| Repeated exact EBCH128 plus RandomStepConv-M30 | **Certified** literal 11% distance with failure below \(2^{-26.192}\) | RandomStepConv representation and online work are impractical | Fixed-outer proof comparator closed |
| Unpadded repeated EBCH128 plus RandomStepConv-M20 | **Certified** 10.9% distance at exact \((k,N)=(2^{20},2^{21})\), with failure below \(2^{-7.584}\) | Margin is below 40 bits; random step maps remain impractical | Power-of-two proof comparator closed |
| One random rate-half constituent reused in every row | At B=256 and M22: 46.59 bits for Q=1 and 64.16 bits for Q=2 at 10.9% | Occupations three and above remain open; arithmetic is diagnostic | B=256 is the first plausible BCH--accumulator target |
| Shortened or subcoded BCH256 near rate one half | Exact parameters; useful packing envelopes and random-like diagnostics | Required full spectrum unavailable; current dense envelope fails badly | Stopped pending authenticated spectrum or stronger invariant |
| One sampled BCH250 with one repeated fanout wrapper | Central spectrum becomes nearly random in the model | Categorical transfer loses on thin mixed faces; setup concentration open | Failed argument, not distance refutation |
| BCH250-124 with independent row-local Fanout-56 | **Certified** finite 11% and failure below \(2^{-40}\) | Effective row code varies; 56-layer online cost is far above 11 ms | Proof benchmark only |
| BCH256 Fanout-56 production-structure proxy | Measured 10.650813 ms at zero layers and 20.454284 ms at 56 layers | Performance exceeds budget; proxy is not the exact proof-model code | Reject 56-layer implementation route |
| Frozen Structured SPIN instance | About 10.823872 ms; modeled-spectrum proof has strong displayed margin | Actual outer spectrum hypothesis remains unproved | Conditional, not a distance certificate |
| One repeated degree-14 block EC constituent | Exact-size \([512,256]\) first moment: 57.7008 kernel bits; dense \(z=0.5\) moment is within 0.0016 bits of random | Pair moment and high-probability shell caps remain open | Primary new constituent route |
| Block Expand--\(t\): degree-14 expander plus \(t\) accumulators | \(t=5\), conditional on \(\operatorname{Var}(A_w)\le2\mathbb E[A_w]\), has a complete outward 42.6254-bit result | Length-512 factor-two variance lemma remains open | Primary expander-constituent target |
| Global two-sided regular Expand--Convolute | **Certified** exact \([2^{21},2^{20}]\), 10.9% distance, and 50.2054 failure bits | Different construction; no decoder theorem or implementation benchmark | Strongest current finite linear-time alternative |

## Detailed record

### A. Random-ensemble baselines

The random-outer investigations established the correct probability-space and
first-moment framework. A random rate-half outer paired with the fixed
RM2Sub-S19 inner has a proved asymptotic 11% distance theorem with
\(B=9\log_2N+O(1)\). The fully random outer/random-inner model is an even
cleaner reference. These results explain how much margin an ideal spectrum
provides and identify the occupation ranges that a structured outer must
control.

These routes worked as theory baselines. They did not produce the desired
explicit, single-repeated-constituent finite encoder. See
`RANDOM_OUTER_RM2SUB_CERTIFICATE.md`, `RANDOM_OUTER_RM2SUB_RAMP.md`, and
`RANDOM_SPIN_PROOF_AUDIT.md`.

### B. One accumulator and structured permutations

For a random block outer followed by one accumulator, logarithmic outer blocks
give asymptotic linear distance with explicit constants. Replacing the uniform
interleaver by a pure bit transpose is not benign: adjacent output pairs retain
structured correlations. Adding independent permutations within transpose
regions gives an exact two-state transfer, but no finite \(k=2^{20}\), 11%
certificate was obtained for the desired structured construction.

The proved asymptotic facts remain useful. The failed pure-transpose step is an
actual structural obstruction, whereas the absence of a finite certificate for
the region-shuffled version is only an open proof problem. See
`ACCUMULATOR_SPIN_ASYMPTOTIC.md`, `BIT_TRANSPOSE_ACCUMULATOR_OBSTRUCTION.md`,
and `REGION_SHUFFLED_TRANSPOSE_ACCUMULATOR.md`.

### C. One repeated Golay--BA-3 constituent

The strongest scalable structured theorem selects one Golay--BA-3 code and
reuses it at every outer position. It proves asymptotic 11% distance with

\[
 B=(39/4)\log_2N+O(1)
\]

and linear ordinary and transposed work. This is the correct single-code
probability space and is recorded in `SINGLE_SAMPLED_BA_RM2SUB_D11.md`.

At finite \(k=2^{20}\), the main candidate used

\[
 B=240,\qquad L=8832,\qquad N_+=2{,}119{,}680.
\]

Under the revised conditioning window, the outward \(Q=1\) margin is 46.479
bits and the outward aggregate margin for \(2\le Q\le64\) is 84.617 bits.
The attempted dense proof did not close. In particular, the common-norm
Hölder bound at 11% misses by more than 1.7 million bits. At 10.9%, all 116
sampled five-band compositions at \(Q=L\) close, with at least about 6016
diagnostic bits, but the barycentric cover and lower occupations remain open.
The search also did not obtain a complete 10% certificate for the unchanged
Golay--BA construction.

The finite failures above are failures of the attempted relaxations, not
counterexamples to distance. Nevertheless, after several hours they did not
produce a usable finite certificate. The measured repeated-code implementation
is about 10.26 ms, leaving approximately 0.74 ms under the 11 ms target. This
route should be revisited only if the construction or proof mechanism changes
materially. See `FINITE_K20_COMPARISON.md`,
`FINITE_K20_TWO_TRACK_DENSE_PLAN.md`, and
`FINITE_K20_DENSE_TRANSFER_TARGET.md`.

### D. Repeated Golay--BA-3 with random Toeplitz convolution

The proof-of-concept enlarged the repeated constituent to

\[
 B=720,\qquad L=2944,\qquad N=2{,}119{,}680
\]

and replaced RM2Sub-S19 by \(g=1\) random lower-triangular Toeplitz
convolution. The diagonal coefficient is one. For each fixed nonzero
difference, the activation output is one and the later outputs are
independent fair bits. The parent dimension is 1,059,840; zero-shortening
11,264 input coordinates gives \(k=2^{20}\).

The exact-in-form binary64 \(Q=1\) calculation has 20.245 bits under the
naive row union and 24.305 bits after grouping all row placements of the same
reused outer word. Outer weights 2 through 5 dominate. The same reuse-aware
calculation reaches 40.58 bits only after deleting every outer shell below
weight 13.

The later prefix-rank identity covers every parent message. Independently
sampled BA rows have 20.244610 bits at \(B=720\), with occupations \(Q\ge2\)
contributing only \(2^{-41.491}\). Four fixed repeated-BA setups give
190.521692 bits at \(B=720\). Four more give the same margin at the current
\(B=240\) design point. The rare BA low shells dominate the ensemble average;
they do not describe the tested fixed setups.

The follow-up froze one \(B=240\) setup with a stable
SplitMix64/Fisher--Yates specification and outward-certified the complete
prefix sum. Its computed margin is 190.521691555 bits, and its claimed
integer margin is 180 bits. The first prefix-rank deficit is at 10,081, and
the prefix kernel becomes zero at 1,163,547. See
`FINITE_K20_BA240_REPEATED_TOEPLITZ_CERTIFICATE.md`,
`FINITE_K20_BA240_REPEATED_TOEPLITZ_MANIFEST.json`,
`ba3_B240_repeated_toeplitz_prefix_outward.json`,
FINITE_K20_BA_IDEAL_CAUSAL_INNER_POC.md,
FINITE_K20_RANDOM_CONV_PREFIX_ANALYSIS.md,
evaluate_ba_ideal_causal_inner_q1.py, and
analyze_ba_ideal_causal_inner_prefix.py.

This closes the finite distance proof for the Toeplitz proof model. It does
not make the construction linear time. Implementing full-length Toeplitz
multiplication, or replacing it by a bounded-state inner with a sufficient
suffix law, remains open.

This result is a substantial finite-proof milestone. The earlier
independent random-block calculation at \(B=240\) retained only 54.043528
bits because its ensemble expectation paid for rare locally bad codes. The
fixed setup has an actual prefix profile that is nearly ideal and retains
190.521691555 outward bits. The certificate covers all parent messages and
does not rely on the expected BA spectrum.

The result is fully certified for the stated mathematical construction. It
does not certify a Toeplitz implementation, a linear-time bound, or the
RM2Sub-S19 variant. The source verifier is reproducible and outward-rounded;
the proof is not formalized in Lean.

### D.1. Return to RM2Sub-S19

The first bridge audit found an exact obstruction to a direct prefix-rank
port. RM2Sub-S19 observes 19 syndrome bits per 128-bit epoch. Across all
16,560 epochs, the syndrome map has rank at most 314,640. Therefore the
shortened outer contains a silent subcode of dimension at least 733,936.
RM2Sub acts as the identity on this subcode for every multiplier schedule.

This dimension bound is not a low-distance counterexample. It proves that
raw prefix ranks omit a necessary weight property. A complete RM2Sub proof
must certify the weight enumerator of the silent subcode and weight-sensitive
profiles of suffix-silent subcodes.

The live state remains promising. The selected RM2Sub \(A\) code gives a
uniform coset-moment bound from its exact spectrum. An optimistic live-only
diagnostic has more than 228,000 bits of single-word margin over 16,559
epochs. The next proof object is therefore a weighted syndrome-prefix
enumerator joined to the existing four-state RM2Sub transfer. See
`FINITE_K20_RM2SUB_FIXED_OUTER_BRIDGE.md`,
`audit_rm2sub_fixed_outer_bridge.py`, and
`rm2sub_fixed_outer_bridge_audit.json`.

The silent subcode's random-code benchmark has about 326,094 bits of margin.
Charging the existing \(B^2\) selected-BA spectrum factor in every row would
cost about 139,667 bits. This leaves a provisional 186,427-bit budget. The
comparison is encouraging but conditional: the structured syndrome map has
not been proved to generate the random-code benchmark.

### D.2. Random bounded-memory step baseline

`FINITE_K20_RANDOMSTEP_CONV_BASELINE.md` defines the bitwise random
bounded-memory linear state machine. At each position it samples one fresh
unrestricted linear map from the input bit and \(M\) state bits to the output
bit and \(M\) next-state bits. The maps are sampled once and shared by every
message.

With a full random linear outer, the exact two-state transfer covers every
nonzero message. At \(k=2^{20}\), \(N=2{,}119{,}680\), and the strict 11%
cutoff, the binary64 first-moment diagnostic gives 1,422.7138 bits at
\(M=10\) and 11,421.6932 bits at \(M=19\). Memory 9 fails. Before shortening,
the \(N/2\)-dimensional parent retains 157.6932 bits at \(M=19\).

The first repeated-BA calculation covers occupation one exactly after the
row-local and region permutations. Conditional on the BA constituent having
no nonzero word outside weights 23 through 217, it gives 21.5360 bits at
\(M=19\) and 41.6332 bits at \(M=23\). Occupations \(Q\geq2\) remain open, so
this is not yet a repeated-BA distance certificate.

The baseline establishes that bounded memory is not the source of the
RM2Sub obstruction. The structured syndrome interface is. RandomStepConv is
not yet a practical replacement: unrestricted maps require about 101.1 MiB
of setup at \(M=19\), and the current arithmetic is not outward rounded.

The repeated-random-constituent follow-up closes the complete finite
comparison model. It uses one admissible \([512,256]\) random constituent in
all 4,140 rows, shortens 44 complete information rows, and sets $M=30$.
For every accepted constituent, an outward all-occupation calculation proves

\[
 \Pr[d_{\min}<233{,}165]<2^{-108}.
\]

The exact outer event

\[
 A_w\le(13/2)\mathbb E[A_w]\quad(1\le w\le512)
\]

holds for one uniform generator with probability greater than
0.6365388449160215. Trying at most 28 generators and retaining the first that
passes the exact test gives a complete bound below $2^{-40}$ for setup abort
or distance failure. This uses one retained code repeatedly; it does not use
independent outer codes in different rows. The exact event test may enumerate
$2^{256}$ messages, and explicit RandomStepConv-M30 maps occupy about
242.83 MiB. The result is therefore a mathematical comparator, not the
implementation candidate. See
`FINITE_K20_REPEATED_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md`.

The intermediate RandomStepConv routes remain useful diagnostics.

- A selected repeated Golay--BA-3 spectrum closes occupations 2 through 64
  at (M=30), but its simultaneous pointwise selection factor is charged
  once per active row and destroys the dense endpoint. The search code and
  receipts are `evaluate_selected_ba_randomstepconv_g1_q2_64.py`,
  `selected_ba240_randomstepconv_g1_s30_q2_64_d11.json`,
  `evaluate_selected_ba_randomstepconv_g1_dense.py`, and the
  `selected_ba240_randomstepconv_g1_s30_*diagnostic.json` files.
- Minimum-distance deletion, information-set deletion, and fixed-basis
  deletion are rejected proof relaxations. They retain late isolated input
  bits and therefore allow the causal state to reset before most of the
  output. The supporting code and receipts are
  `evaluate_fixed_ba_randomstepconv_pivot_deletion.py` and the
  `fixed_ba240_randomstepconv_*` files. These failures are not low-distance
  counterexamples for the actual BA construction.
- Rank-only deletion for one repeated random ([240,120]) outer is also too
  weak. It is recorded by
  `evaluate_repeated_random_outer_randomstepconv_rank.py` and
  `repeated_random240_randomstepconv_g1_s19_rank_d11.json`.
- Replacing the fixed constituent by its exact random-code mean spectrum
  closes the dense (B=240) transfer decisively. The coarse and refined
  receipts are
  `repeated_randomspectrum240_randomstepconv_g1_s30_q65_8738_diagnostic.json`
  and
  `repeated_randomspectrum240_randomstepconv_g1_s30_q65_8738_refined.json`.
  At (B=240), however, the available simultaneous fixed-code selection
  factor costs too much at full occupation. Moving to (B=512) creates
  enough shortening slack to pay a rigorous fixed-code envelope and yields
  the completed certificate above.

### E. Independent-row Golay--BA-3

Sampling a BA row independently at each outer position makes the expected
spectrum factor and therefore simplifies the finite first moment. The same
\(B=240\) low-occupation ranges close: \(Q=1\) has 46.479 outward bits and
\(2\le Q\le64\) has 84.617 aggregate outward bits. Occupations
\(65\le Q\le8832\) remain unproved. Efficiently enforcing or authenticating
the conditioning event also remains open.

The transposed map measured 10.407942 ms, 2.699% faster than the same-binary
BCH256 baseline and 1.411% slower than repeated BA. This is useful performance
evidence, but the probability space violates the requirement that one outer
code be repeated. It is therefore not the final construction. See
`FINITE_K20_COMPARISON.md` and `FINITE_K20_INDEPENDENT_SETUP.md`.

### F. EBCH32--ParityFanout--BA

This power-of-two construction uses eight genuine \([32,16,8]\) constituents,
ParityFanout, and BA accumulators. At

\[
 (k,N,d)=(2^{20},2^{21},209715),
\]

the outward first-moment calculation certifies distance strictly above 10%
with failure below \(2^{-51}\). A bounded six-draw wrapper gives a displayed
combined bound below \(2^{-45}\), but its specified exact acceptance test can
require examining as many as \(6\cdot8192\cdot2^{128}\) messages. That is a
mathematical wrapper, not practical setup.

The 11% effort closed low occupations and much of the dense cover, but 149
dense cells remained pending in the last checkpoint. The construction also
uses independently sampled rows, contrary to the desired final interface.
Thus, this route worked as a complete 10% mathematical certificate and failed
as a practical, single-code 11% construction. See
`FINITE_K20_D10_CERTIFICATE.md`,
`FINITE_K20_EBCH32_PARITYFANOUT_SETUP.md`, and the historical checkpoints in
`FINITE_K20_PROOF_OBLIGATIONS.md`.

### G. One repeated exact EBCH128 constituent

For the exact \([128,64,22]\) extended BCH spectrum, occupation one has about
26.714 diagnostic bits and occupations 2 through 100 have about 28.450
aggregate diagnostic bits. Their combined diagnostic margin is about 26.335
bits. These numbers do not meet a 40-bit target and are not outward-certified.

The dense even-row reduction fails near the all-active endpoint. Omitting one
region loses too much margin, while the tested common-norm Hölder shortcut
again loses about 1.77 million bits. The remaining issue is an exact
high-occupation transfer that retains the all-one and mixed-tail structure.
The earlier description of this route as a complete 26-bit certificate was
corrected. It has no complete finite 11% certificate. See
`FINITE_K20_REPEATED_EBCH128_D11_STATUS.md`.

With RandomStepConv-M30 replacing RM2Sub, the fixed repeated constituent now
has a complete outward 11% certificate. At \(N=2119680\), it proves
\(d_{\min}\ge233165\) except with probability below
\(2^{-26.1921836158}\). A dispersed parity-pivot transfer closes every dense
occupation without sampling or testing an outer code. Occupation one is
dominant, and occupation two is next. The earlier aligned-pivot and coverage
relaxations remain recorded as rejected proof methods. See
`FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_CERTIFICATE.md` and
`FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_STATUS.md`.

Removing the 176 zero rows gives exact dimensions
\((k,N)=(2^{20},2^{21})\). At 11%, the ideal random-code first-moment margin
is only 188.379 bits, so the pointwise BCH comparison is unaffordable.
Relaxing the target to 10.9% gives 6,529.062 ideal-random bits and permits a
complete M20 certificate. The parity-pivot transfer covers occupations 2
through 16319. An exact-complement recurrence covers the final 65
occupations. The outward verifier proves failure below
\(2^{-7.58436084118}\). See
`FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_CERTIFICATE.md` and
`FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_MANIFEST.json`.

At M19, the corresponding complete diagnostic closes with 0.759 bits. That
result is not outward-certified. It shows that the sparse one-row event, not
the dense endpoint, limits further memory reduction at 10.9%.

The complete 10.9% argument is frozen as the reusable
bit-transpose/RandomStepConv template in
`BITTRANSPOSE_RANDOMSTEP_CONV_PROOF_TEMPLATE.md`. Its immutable reference is
commit `b2788fdb81e30cb78d68efc2ed06eb8af2b94407`. The template fixes the
construction interface, setup probability space, four occupation ranges,
outward certificate condition, and substitution obligations. Future outer,
routing, inner, distance, or memory variants must create new receipts and
manifests rather than overwrite the frozen artifacts.

### Single-random-constituent size proxy

The Random SPIN theorem resamples its random local injection independently
across blocks. The new proxy instead samples one uniform rate-half injection
and reuses it in every row. It retains the region-permuted bit transpose and
RandomStepConv-M30.

Occupation one admits an exact ensemble average without a resampling
assumption. At 10.9%, the nearest-binary64 margins are 18.7859 bits for
B=128, 55.3387 bits for B=256, and 127.9284 bits for B=512. Thus B=256 is the
first tested size above the 40-bit sparse target.

At B=256, memory 20 gives 35.7452 bits, memory 21 gives 42.1087 bits, and
memory 22 gives 46.5928 bits. Memory 22 is the safer first target because
memory 21 leaves only 2.1087 bits for every other occupation.

This is not a complete certificate. Reusing one injection requires the joint
law of the active local messages. The size scan initially left that law open.
See
`SINGLE_RANDOM_CONSTITUENT_SIZE_PROXY.md`,
`evaluate_single_random_constituent_q1.py`, and
`single_random_constituent_B128_256_512_q1_s30_d109.json`. The B=256 memory
scan is in the `single_random_constituent_B256_q1_s*_d109.json` receipts.

The exact-shell occupation-two diagnostic now closes with 64.1559 bits at
B=256 and M22. Equal nonzero local messages form the rank-one class and
dominate the result. Distinct nonzero local messages form the rank-two class
and give 99.7046 bits. Combining occupations one and two leaves 46.5928
bits, so occupation one remains the bottleneck. See
`SINGLE_RANDOM_CONSTITUENT_Q2.md`,
`evaluate_single_random_constituent_q2.py`, and
`single_random_constituent_B256_q2_s22_d109.json`. Occupations three and
above remain open.

The first all-occupation wrapper conditioned on the uniform spectrum bound
\(A_w\le8\mu_w\). One sample satisfies this event with probability at least
0.7655, but the distance bound fails at \(Q=L\) by 18,102 bits. The dense
failure is the repeated factor \(8^L\), not the ideal random-spectrum
transfer; removing that factor makes the endpoint pass by 6,474 bits.

`SINGLE_RANDOM_CONSTITUENT_ALLQ_RAMP.md` replaces the uniform event by a
three-band target. An elementary calculation gives this event at least
94.07% one-sample probability. Its factor-1.5 central class passes every
occupation, with 1,681.73 bits at \(Q=L\). Mixed fringe/central compositions
remain open and require a categorical simplex cover. The rejected wrapper
and receipt are `evaluate_single_random_constituent_allq.py` and
`single_random_constituent_B256_allq_factor8_s22_d109.json`. The three-band
event calculation is in
`evaluate_single_random_constituent_three_band_event.py` and
`single_random_constituent_B256_three_band_event.json`.

### H. Pure BCH constituents near length 256

The genuine parity-extended BCH parent has parameters
\([256,131,\ge38]\). Three-coordinate shortening gives
\([253,128,\ge38]\), but its rate is too high: even the ideal-random
benchmark misses the finite 11% target by about 12,113 bits. Six-coordinate
shortening gives \([250,125,\ge38]\); taking a fixed codimension-one subcode
gives the desired \([250,124,\ge38]\) rate. A deterministic
\([256,128,\ge38]\) Eq3 subcode also exists, but it is not a pure BCH code and
its full spectrum is unknown.

For the BCH250 option, a random-like modeled spectrum gives 59.587 bits at
\(Q=1\). The authenticated integer packing envelope gives only 13.777 bits at
\(Q=1\) and 19.044 aggregate bits for \(2\le Q\le100\). The pointwise dense
envelope fails by 446,724 bits at \(Q=L\); correcting only the total mass still
leaves about 32,500 bits of failure. These are proof-envelope failures. They
show that minimum distance and total dimension alone are too weak; they do not
show that the actual BCH code has bad distance after composition.

No authenticated full spectrum for the required near-256 constituent was
obtained. The route was therefore stopped as requested. See
`FINITE_BCH256_SHORTENING_OPTIONS.md` and
`shortened_bch256_parameter_audit.json`.

### I. One sampled BCH250 code with one repeated fanout wrapper

The repeated-wrapper experiment asked whether ParityFanout could make one
sampled BCH250 subcode sufficiently random-like while preserving one effective
outer code in every row. In the random-even source model, the central spectrum
is nearly random and the all-active RM2Sub endpoint has 177.405 diagnostic bits
before outer transfer losses.

The attempted broad-band transfer fails on thin mixed faces. For example, two
weight-36 rows and 8446 central rows still have a diagnostic margin of
-116.335 bits after 128 fanout layers. Other tested relaxations lose about
-228.861 bits or much more. The obstruction is the categorical-conditioning
relaxation, not roundoff and not a demonstrated bad codeword. A theorem for one
sampled wrapper would also require a higher-moment or concentration argument
for setup. No such theorem was completed. See
`ONE_SAMPLED_BCH250_PARITYFANOUT_RAMP.md`.

### J. BCH250-124 with independent row-local Fanout-56

This is the only complete finite 11% certificate produced in the workstream.
It uses one fixed \([250,124,\ge38]\) base subcode, but each of the 8576 rows
receives an independently sampled 56-layer ParityFanout wrapper and local
coordinate permutation. Its parameters are

\[
 B=250,\quad L=8576,\quad N=2{,}144{,}000,\quad D=235{,}840.
\]

The outward certificate proves construction failure below \(2^{-40}\); its
displayed combined margin is about 41.51 bits. A 52-layer candidate fails the
present tail-plus-central proof. That failure does not prove that 56 layers are
necessary.

The result is mathematically complete in its stated probability space, but it
does not use one effective repeated outer code: row-local wrappers change the
row code. Its 56-layer work is also too expensive. It remains valuable as a
proof benchmark showing what spectral regularization is sufficient. See
`FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_CERTIFICATE.md` and
`FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_MANIFEST.json`.

### K. Fanout-56 implementation evidence

The production-structure BCH256 proxy preserves the packed, tiled, fused
encoder organization. On one pinned Peach Ryzen 9 7950X core, its median times
are 10.650813 ms with zero fanout layers and 20.454284 ms with 56 layers. The
increment is 9.803471 ms and the ratio is 1.920443. The frozen production BCH
path is about 10.823872 ms, which supports the proxy's relevance.

At the measured per-layer cost, the 11 ms ceiling allows roughly one layer;
two layers are borderline. Meeting the ceiling with 56 layers would require an
approximately 28-fold reduction in the added fanout cost. Ordinary local
tuning is not a credible way to obtain that reduction. The exact proof-model
BCH250 encoder was not benchmarked end to end, so these measurements are a
cost rejection of the 56-layer mechanism, not an exact timing claim for that
encoder. See `FINITE_K20_BCH256_FANOUT56_PROXY_BENCHMARK.md` and
`bch256_fanout56_proxy_peach_7950x.json`.

A separate direct-scatter prototype was removed at the user's request because
it did not preserve the optimized production structure. No performance
conclusion from that prototype is retained.

### L. Frozen Structured SPIN instance

The frozen structured encoder is fast enough to remain the performance
reference: its measured transposed time is about 10.823872 ms. Its existing
distance receipt, however, substitutes a modeled even-floor outer spectrum for
the actual constituent spectrum. Its displayed 11% margin is 55.864642 bits,
but that value is conditional on the model and is not a proof for the original
structured code.

This distinction motivated the finite constituent search. The frozen instance
has useful implementation evidence and a conditional proof calculation, but
no authenticated finite 11% theorem for its actual outer code.

## Supporting material generated

All paths in this section are relative to
`workstreams/finite_asymptotic_theory/`. The lists distinguish proof artifacts
from diagnostics. A filename containing `*` denotes the complete local family
with that prefix or suffix. The certificate manifests give the authoritative
hash-bound lists for certified claims.

### Shared finite framework and work ledger

- Framework and interfaces: `FINITE_LENGTH_FRAMEWORK.md`,
  `ARBITRARY_LENGTH_WRAPPER.md`, `ASYMPTOTIC_SCALING.md`,
  `STRUCTURED_SPIN_THEOREM_TARGET.md`, and
  `STRUCTURED_TRANSLATION_FROM_RANDOM_BASELINE.md`.
- Finite planning and status: `FINITE_K20_COMPARISON.md`,
  `FINITE_K20_PROOF_OBLIGATIONS.md`, `FINITE_K20_CONDITIONING_WINDOW.md`,
  `FINITE_K20_DENSE_TRANSFER_TARGET.md`, `FINITE_K20_DENSE_CELL_COVER.md`,
  `FINITE_K20_TWO_TRACK_DENSE_PLAN.md`, and this attempt log.
- Partial-certificate index: `FINITE_K20_PARTIAL_CERTIFICATE_MANIFEST.json`
  and `build_finite_k20_partial_manifest.py`.

### Random-ensemble baselines

- Random SPIN audit and evaluators: `RANDOM_SPIN_PROOF_AUDIT.md`,
  `evaluate_random_spin_finite.py`, `certify_random_spin_linear_exponent.py`,
  and `optimize_random_spin_asymptotic_constants.py`.
- Random-outer theorem: `RANDOM_OUTER_RM2SUB_RAMP.md`,
  `RANDOM_OUTER_RM2SUB_CERTIFICATE.md`, and
  `RANDOM_OUTER_RM2SUB_CERTIFICATE_MANIFEST.json`.
- RM2Sub theorem components: `RM2SUB_ONE_ACTIVE_CONTINUUM.md`,
  `RM2SUB_TWO_ACTIVE_CONTINUUM.md`,
  `RM2SUB_FIXED_OCCUPATION_CONTINUUM.md`,
  `RM2SUB_UNIFORM_FIXED_OCCUPATION.md`, and
  `RM2SUB_DENSE_OCCUPATION.md`.
- RM2Sub generators and verifiers: `derive_rm2sub_one_active_continuum.py`,
  `derive_rm2sub_two_active_continuum.py`,
  `analyze_rm2sub_fixed_occupation_continuum.py`,
  `certify_rm2sub_uniform_fixed_occupation.py`,
  `analyze_rm2sub_dense_occupation.py`, `certify_rm2sub_dense_small.py`, and
  `certify_rm2sub_dense_compact_interval.py`.
- Theorem receipts: `rm2sub_one_active_continuum.json`,
  `rm2sub_two_active_continuum.json`,
  `rm2sub_fixed_occupation_continuum_d11.json`,
  `rm2sub_uniform_fixed_occupation_d11.json`,
  `rm2sub_dense_small_exact_d11.json`, and
  `rm2sub_dense_compact_interval_d11.json`.
- Finite one-active diagnostic: `evaluate_random_outer_rm2sub_one_active.py`
  and `random_outer_rm2sub_one_active_B256_d11.json`.

### One accumulator and structured permutations

- Theorem and verifier: `ACCUMULATOR_SPIN_ASYMPTOTIC.md` and
  `certify_accumulator_spin_asymptotic.py`.
- Pure-transpose obstruction: `BIT_TRANSPOSE_ACCUMULATOR_OBSTRUCTION.md` and
  `analyze_bit_transpose_accumulator_obstruction.py`.
- Region-shuffled analysis: `REGION_SHUFFLED_TRANSPOSE_ACCUMULATOR.md` and
  `analyze_region_shuffled_transpose_accumulator.py`.
- Preliminary BA scaling work: `BA_ASYMPTOTIC_EXPLORATION.md` and
  `analyze_ba_asymptotics.py`.

### One repeated Golay--BA-3 constituent

- Asymptotic theorem: `SINGLE_SAMPLED_BA_RM2SUB_D11.md` and
  `SINGLE_SAMPLED_BA_RM2SUB_CERTIFICATE_MANIFEST.json`.
- Asymptotic verifiers: `certify_golay_ba_concave_majorant.py`,
  `certify_golay_ba_rm2sub_joint_interval.py`,
  `certify_golay_ba_rm2sub_sparse.py`, and
  `certify_golay_ba_rm2sub_weight_coupled_fixed.py`.
- Asymptotic receipts: `golay_ba3_concave_majorant.json`,
  `golay_ba3_rm2sub_joint_interval_d11.json`,
  `golay_ba3_rm2sub_sparse_d11.json`, and
  `golay_ba3_rm2sub_weight_coupled_fixed_d11.json`.
- Supporting theorem: `WEIGHT_COUPLED_FIXED_OCCUPATION.md`.
- Finite low-occupation verifiers:
  `certify_golay_ba_rm2sub_finite_one_active.py` and
  `certify_golay_ba_rm2sub_finite_q2_64.py`.
- Finite outward receipts:
  `golay_ba3_rm2sub_finite_B240_q1_outward_w23_217_k20_d11.json` and
  `golay_ba3_rm2sub_finite_B240_q2_64_outward_w23_217_k20_d11.json`.
- Finite evaluators: `evaluate_golay_ba_rm2sub_finite.py`,
  `evaluate_golay_ba_rm2sub_finite_holder.py`,
  `diagnose_finite_k20_conditioning_window.py`,
  `diagnose_finite_k20_qL_column_holder.py`,
  `diagnose_finite_k20_band_pure.py`,
  `diagnose_finite_k20_band_compositions_qL.py`,
  `diagnose_finite_k20_band_cover_qL.py`,
  `diagnose_finite_k20_band_simplex_cover_qL.py`, and the
  `diagnose_finite_k20_three_band_renyi_*` family.
- Finite diagnostic receipts: the
  `golay_ba3_rm2sub_finite_B240_*` family, including the conditioning-window,
  column-Hölder, pure-band, composition, dense-probe, and three-band receipts.
- Block-size comparison data: `golay_ba3_B216_conditioned_spectrum_upper.json`,
  `golay_ba3_B240_conditioned_spectrum_upper.json`,
  `golay_ba3_B264_conditioned_spectrum_upper.json`, and the
  `golay_ba3_rm2sub_finite_B264_*` receipts.

The ideal causal-inner proof of concept is supported by
FINITE_K20_BA_IDEAL_CAUSAL_INNER_POC.md,
evaluate_ba_ideal_causal_inner_q1.py, and
ba3_B720_ideal_causal_g1_q1_d11.json.

The all-message prefix analysis is supported by
FINITE_K20_RANDOM_CONV_PREFIX_ANALYSIS.md,
analyze_ba_ideal_causal_inner_prefix.py,
sweep_ba_ideal_causal_inner_blocks.py,
ba3_B240_ideal_causal_g1_prefix_diagnostic.json,
ba3_B720_ideal_causal_g1_prefix_diagnostic.json, and
ba_ideal_causal_g1_block_sweep.json.

The frozen Toeplitz certificate and RM2Sub bridge are supported by
`TOEPLITZ_PREFIX_ARGUMENT.md`,
`FINITE_K20_BA240_REPEATED_TOEPLITZ_CERTIFICATE.md`,
`FINITE_K20_BA240_REPEATED_TOEPLITZ_MANIFEST.json`,
`certify_ba240_repeated_toeplitz_prefix_outward.py`,
`ba3_B240_repeated_toeplitz_prefix_outward.json`,
`FINITE_K20_RM2SUB_FIXED_OUTER_BRIDGE.md`,
`audit_rm2sub_fixed_outer_bridge.py`, and
`rm2sub_fixed_outer_bridge_audit.json`.

The random bounded-memory baseline is supported by
`FINITE_K20_RANDOMSTEP_CONV_BASELINE.md`,
`evaluate_random_outer_randomstepconv_g1.py`,
`random_outer_randomstepconv_g1_k20_d11.json`,
`random_outer_randomstepconv_g1_kparent_d11.json`,
`evaluate_ba_randomstepconv_g1_one_active.py`, and the
`ba240_randomstepconv_g1_s*_q1_d11.json` receipts for memories 19 through 23.
The complete repeated-random-constituent closure is supported by
`FINITE_K20_REPEATED_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md`,
`evaluate_repeated_random512_randomstepconv_g1.py`,
`repeated_random512_randomstepconv_g1_s30_allq_d11.json`,
`certify_repeated_random512_randomstepconv_g1_outward.py`, and
`repeated_random512_randomstepconv_g1_s30_allq_outward_d11.json`. The
hash-bound index is
`FINITE_K20_REPEATED_RANDOM512_RANDOMSTEP_CONV_MANIFEST.json`.

### Independent-row Golay--BA-3

- Construction interface: `FINITE_K20_INDEPENDENT_SETUP.md`.
- The shared finite comparison, proof-obligation documents, low-occupation
  verifiers, and receipts listed above support this route.
- `FINITE_K20_PARTIAL_CERTIFICATE_MANIFEST.json` records the exact proved and
  unproved components. It also records the hash of the external performance
  receipt used by `FINITE_K20_COMPARISON.md`.

### EBCH32--ParityFanout--BA

- Construction and certified result: `FINITE_K20_EBCH32_PARITYFANOUT_SETUP.md`,
  `FINITE_K20_D10_CERTIFICATE.md`, and
  `FINITE_K20_D10_CERTIFICATE_MANIFEST.json`.
- Constituent spectrum: `ebch32_16_delta8_spectrum.csv`.
- Setup and low-occupation verifiers:
  `certify_ebch32_parityfanout_ba_setup.py`,
  `certify_ebch32_parityfanout_ba_rm2sub_q1.py`, and
  `certify_ebch32_parityfanout_ba_rm2sub_q2_64.py`.
- Dense and combined verifiers:
  `certify_ebch32_parityfanout_ba_one_band_interval_cover_d10_outward.py`,
  `certify_ebch32_parityfanout_dense_cells_outward.py`,
  `certify_ebch32_parityfanout_ba_combined_d10.py`, and
  `certify_ebch32_parityfanout_ba_bounded_setup_d10.py`.
- Certified receipts:
  `ebch32_parityfanout31x33_ba3_B256_setup_outward.json`,
  `ebch32_parityfanout31x33_ba3_rm2sub_B256_q1_outward_d11.json`,
  `ebch32_parityfanout31x33_ba3_rm2sub_B256_q2_64_outward_d11.json`,
  `ebch32_parityfanout31x33_ba3_B256_one_band_interval_cover_d10_outward.json`,
  `ebch32_parityfanout31x33_ba3_B256_combined_outward_d10.json`, and
  `ebch32_parityfanout31x33_ba3_B256_bounded_setup_combined_outward_d10.json`.
- The incomplete 11% dense-cover record:
  `ebch32_parityfanout31x33_ba3_B256_three_band_cover_all_q_d11.json`,
  `ebch32_parityfanout31x33_ba3_B256_dense_cells_outward_partial.json`,
  `ebch32_parityfanout31x33_ba3_B256_dense_cover_partition_audit.json`, and
  the `ebch32_parityfanout31x33_ba3_B256_split_low_*` receipts.
- Diagnostic and audit programs: the
  `diagnose_ebch32_parityfanout_*` family,
  `audit_ebch32_parityfanout_dense_cover_partition.py`, and
  `check_ebch32_parityfanout_ba_row_min_distance.py`.

### One repeated EBCH128 constituent

- Status document: `FINITE_K20_REPEATED_EBCH128_D11_STATUS.md`.
- Diagnostic programs: `diagnose_ebch128_even_body_dense.py` and
  `diagnose_ebch128_qL_column_holder.py`.
- Diagnostic receipts: `ebch128_direct_rm2sub_s19_k20_d11_one_active_diagnostic.json`,
  `ebch128_direct_rm2sub_s19_k20_d11_one_active_three_state_diagnostic.json`,
  `ebch128_direct_rm2sub_s19_k20_d11_q2_100_diagnostic.json`,
  `ebch128_even_body_dense_k20_d11_diagnostic.json`,
  `ebch128_even_body_dense_high_k20_d11_diagnostic.json`, and
  `ebch128_rm2sub_s19_k20_d11_qL_column_holder_diagnostic.json`.
- RandomStepConv status:
  `FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_STATUS.md`.
- RandomStepConv diagnostic programs:
  `evaluate_ebch128_randomstepconv_q1_exact.py`,
  `evaluate_ebch128_randomstepconv_g1.py`, and
  `evaluate_ebch128_randomstepconv_coverage.py`.
- RandomStepConv receipts:
  `ebch128_randomstepconv_g1_s30_q1_exact_d11_diagnostic.json`,
  `ebch128_randomstepconv_g1_s30_allq_d11_diagnostic.json`,
  `ebch128_randomstepconv_g1_s45_allq_d11_diagnostic.json`, and
  `ebch128_randomstepconv_g1_s30_coverage_d11_diagnostic.json`.
- Closing RandomStepConv theorem and manifest:
  `FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_CERTIFICATE.md` and
  `FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_MANIFEST.json`.
- Closing diagnostic and outward receipt:
  `ebch128_randomstepconv_g1_s30_parity_pivot_d11_diagnostic.json` and
  `ebch128_randomstepconv_g1_s30_parity_pivot_outward_d11.json`.
- Power-of-two 10.9% theorem and manifest:
  `FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_CERTIFICATE.md` and
  `FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_MANIFEST.json`.
- Frozen reusable proof template:
  `BITTRANSPOSE_RANDOMSTEP_CONV_PROOF_TEMPLATE.md`.
- Power-of-two M20 witness and outward receipts:
  `ebch128_randomstepconv_g1_s20_pow2_parity_pivot_d109_diagnostic.json`,
  `ebch128_randomstepconv_s20_pow2_dense_complement_r64_d109_u074.json`, and
  `ebch128_randomstepconv_g1_s20_pow2_outward_d109.json`.
- Power-of-two M19 diagnostics:
  `ebch128_randomstepconv_g1_s19_pow2_q1_exact_d109_diagnostic.json`,
  `ebch128_randomstepconv_g1_s19_pow2_parity_pivot_d109_diagnostic.json`, and
  `ebch128_randomstepconv_s19_pow2_dense_complement_r64_d109_u074.json`.

### BCH constituents near length 256

- Option audit: `FINITE_BCH256_SHORTENING_OPTIONS.md` and
  `shortened_bch256_parameter_audit.json`.
- Spectrum-free envelope tools: `audit_shortened_bch250_125.py` and
  `build_shortened_bch250_subcode_envelope.py`.
- Parent diagnostics: the `shortened_bch250_125_{modeled,packing}_*` receipt
  families.
- Subcode sweep: the `shortened_bch250_*_subcode_constant_weight_envelope.json`
  and `shortened_bch250_*_subcode_rm2sub_s19_*_q1_diagnostic.json` families.
- Modeled power-of-two comparison:
  `ebch256_128_modeled_rm2sub_s19_k20_d11_one_active_diagnostic.json`.

### One sampled BCH250 code with one repeated fanout wrapper

- Analysis: `ONE_SAMPLED_BCH250_PARITYFANOUT_RAMP.md`.
- Core programs: `analyze_one_sampled_bch250_parityfanout.py`,
  `sweep_bch250_parityfanout_mixing.py`,
  `diagnose_bch250_parityfanout_three_band_qL.py`, and
  `diagnose_bch250_parityfanout_renyi_defect_edges.py`.
- Mixing sweep: `bch250_parityfanout_mixing_sweep.json`.
- Source-model and envelope families:
  `bch250_parityfanout31x33_l*_random_even_model.json`,
  `bch250_parityfanout31x33_l*_packing_expected_envelope.json`, and
  `bch250_parityfanout31x33_l*_spectrum_comparison.json`.
- Thin-face and endpoint diagnostics: the
  `bch250_parityfanout31x33_l*_cutoff*_edges_qL_d11_diagnostic.json`,
  `bch250_parityfanout31x33_l*_pure_qL_d11_diagnostic.json`,
  `bch250_parityfanout31x33_l*_renyi_*_diagnostic.json`, and
  `bch250_parityfanout31x33_l*_finite_*_d11_diagnostic.json` families.
- Inner reference: `rm2sub_dense_alpha1_d11_parityfanout_reference.json`.

### BCH250-124 with independent row-local fanout

- Initial certificate search: `FINITE_K20_BCH250_124_ROWLOCAL_FANOUT_STATUS.md`,
  `FINITE_K20_BCH250_124_ROWLOCAL_FANOUT_CERTIFICATE.md`, and
  `FINITE_K20_BCH250_124_ROWLOCAL_FANOUT_MANIFEST.json`.
- Final Fanout-56 certificate:
  `FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_CERTIFICATE.md` and
  `FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_MANIFEST.json`.
- Spectrum and finite verifiers:
  `certify_bch250_124_rowlocal_fanout_spectrum_arb.py`,
  `certify_bch250_124_shell_sensitive_q1_outward.py`,
  `certify_bch250_124_shell_sensitive_q2_31_outward.py`, and
  `certify_bch250_124_rowlocal_fanout_dense_intervals_outward.py`.
- Final certified receipts:
  `bch250_124_parityfanout31x33_l56_cutoff12_spectrum_outward.json`,
  `bch250_124_parityfanout31x33_l56_q1_outward_d11.json`,
  `bch250_124_parityfanout31x33_l56_q2_31_outward_d11.json`, and
  `bch250_124_parityfanout31x33_l56_cutoff12_dense_q32_8576_outward_d11.json`.
- Search and diagnostic programs:
  `evaluate_bch250_rowlocal_fanout_finite_cover.py`,
  `diagnose_bch250_rowlocal_fanout_finite_interval_cover.py`,
  `diagnose_bch250_rowlocal_fanout_poisson_edges.py`,
  `diagnose_bch250_rowlocal_fanout_pure_occupations.py`, and
  `diagnose_bch250_rowlocal_fanout_renyi_dense.py`.
- Diagnostic receipts: the `bch250_rowlocal_fanout_l*_*_diagnostic.json`
  family and the `bch250_124_parityfanout31x33_l*_*_diagnostic.json` family.
- Manifest builders: `build_bch250_124_rowlocal_fanout_manifest.py` and
  `build_bch250_124_rowlocal_fanout56_manifest.py`.

### Fanout-56 implementation evidence

- Standalone row-local benchmark:
  `FINITE_K20_BCH250_FANOUT56_IMPLEMENTATION_BENCHMARK.md`,
  `FINITE_K20_BCH250_FANOUT56_IMPLEMENTATION_MANIFEST.json`,
  `benchmark_bch250_rowlocal_fanout56.cpp`,
  `bch250_rowlocal_fanout56_smoke.json`, and
  `bch250_rowlocal_fanout56_benchmark_windows_i7_13700h.json`.
- Production-structure proxy:
  `FINITE_K20_BCH256_FANOUT56_PROXY_BENCHMARK.md`,
  `FINITE_K20_BCH256_FANOUT56_PROXY_MANIFEST.json`,
  `bch256_fanout56_proxy_peach_7950x.json`, and the
  `proxy256_fanout56/` source tree.
- Transpose-circuit prototype: `build_bch250_124_transpose_circuit.py`,
  `bch250_124_transpose_circuit.json`, `Bch250x124TransposeCircuit.h`, and
  `Bch250x124TransposeCircuit2.h`.

### Frozen inputs and companion evidence

The following inputs support the workstream but were not generated in this
directory.

- The frozen Structured SPIN construction is under
  `constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/`.
  Its relevant files are
  `constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/MAIN_CODE_FREEZE.md`,
  `constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/PROOF_STATUS.md`,
  `constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/frozen_source/SOURCE_MANIFEST.json`,
  and the construction's `receipts/` directory.
- The selected RM2Sub-S19 receipts are under
  `constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/receipts/min_state/`.
- The companion linear-time evidence was read from
  `C:/Users/peter/.codex/worktrees/3ef9/permute_conv/workstreams/linear_time_audit/`.
  Its principal files are `ASYMPTOTIC_DISTANCE_CERTIFICATE.md`,
  `LINEAR_OUTER_CERTIFICATE.md`, `OUTER_INTERVAL_RECEIPT.json`,
  `verify_outer_interval.py`, and `verify_outer_sparse_constants.py`. That
  companion directory is not replicated in this checkout.

## What the investigation established

1. Random-like outer spectra have ample finite margin at 11%. The central
   difficulty is proving a sufficiently tight spectrum property for one
   practical repeated constituent.
2. Minimum distance, dimension, and a packing envelope are not enough for the
   dense occupation proof. Shell shape and mixed compositions matter.
3. Independent row randomization can close the theorem, but Fanout-56 is too
   expensive and violates the desired repeated-code interface.
4. A complete 10% certificate can coexist with impractical setup. Setup must
   therefore be part of the construction theorem, not an afterthought.
5. The current BA implementation has only about 0.74 ms of online budget.
   Any robustifying mechanism must be much cheaper than 56 row-local fanout
   layers or must replace existing work rather than merely add to it.
6. None of the failed dense bounds is a proof that the corresponding code has
   distance below the target.
7. Random Toeplitz convolution gives 20.244610 all-message bits for
   independent BA rows at \(B=720\). Occupations \(Q\ge2\) contribute only
   \(2^{-41.491}\); rare occupation-one words cause essentially all loss.
8. Four fixed repeated-BA setups at each of \(B=240\) and \(B=720\) give
   190.521692 binary64 bits. The ensemble average is heavy-tailed. It is not
   representative of a typical selected BA code.
9. The exact fixed-code interface is the prefix-rank sequence. It can be
   computed by Gaussian elimination without enumerating the message space.
10. One frozen \(B=240\) repeated-BA setup has a complete outward 11%
    Toeplitz certificate with 180 claimed bits.
11. RM2Sub-S19 has a shortened silent subcode of dimension at least 733,936.
    A return to RM2Sub needs weight-sensitive syndrome-prefix profiles.
12. A fresh random linear bounded-memory step nearly matches the ideal random
    inner. With a full random outer, memory 10 already closes the finite 11%
    target; memory 19 has 11,421.6932 diagnostic bits. The corresponding
    repeated-BA occupation-one calculation has 21.5360 bits at memory 19 and
    41.6332 bits at memory 23 after the stated BA conditioning.

## Untried or materially changed routes

The following ideas are candidates, not completed attempts.

- Benchmark forward and transposed multiplication by the resulting random
  binary Toeplitz matrix. This inner is exact for the proof model but is not
  the current linear-time RM2Sub implementation.
- Build a weighted syndrome-prefix evaluator for the frozen repeated-BA
  outer and join it to the four-state RM2Sub transfer. The silent sector has
  enough provisional budget to absorb a \(B^2\) per-row selection factor.
- If the larger BA-3 constituent remains too weak, add a third accumulator
  inside the constituent. This is a genuinely changed BA-4 construction and
  must be costed against the approximately 0.74 ms budget.
- Test at most one to three very cheap post-BA mixing layers. The Fanout-56
  benchmark implies that this is the only additive fanout scale plausibly
  compatible with the current online budget. No distance certificate for such
  a small layer count is presently known.

These routes must be entered in this log only after their probability space,
proof result, and measured or estimated cost are explicit.

## Continuation rule

Before starting another finite route, record its exact constituent,
probability space, repeated-code semantics, length and shortening rule,
distance threshold, setup algorithm, proof target, and online cost model. On
completion, classify the outcome with the vocabulary above and link every
receipt. A numerical relaxation that fails must be logged as a failed argument,
not as evidence of a low-distance code.

## EBCH128x4 followed by repeated accumulator stages

`EBCH128X4_BA_REPEATED_CONCENTRATION_STATUS.md` records the return to a
single repeated structured constituent. The base is the direct sum of four
fixed extended BCH \([128,64,22]\) codes. Setup samples one independent
uniform 512-coordinate permutation before each accumulator stage. The sampled
constituent is reused in every one of the 4,096 outer rows.

The expected-spectrum recurrence shows that three accumulator stages are
already close to the random-linear spectrum in the bands relevant to the
finite transfer. An expected spectrum is not a high-probability statement.
The present sufficient concentration target is

\[
 \operatorname{Var}(A_w)\le512\,\mathbb E[A_w]
 \qquad(1\le w\le512).
\]

This target is not proved. Conditional on it, shellwise Markov--Cantelli caps
fail with probability at most \(2^{-42.255553}\). At RandomStepConv-M22 and
the 10.9% cutoff, exact shared-category binary64 diagnostics cover every
occupation \(Q<160\). The bottleneck is the all-defect \(Q=3\) composition,
with 46.114992 bits. Exact \(Q=1\) and \(Q=2\) margins are 71.741239 and
152.726349 bits. Sampled dense compositions have at least 559.784886 bits,
but a complete dense convex cover remains open.

The ordinary spectrum is insufficient to prove the variance target. The
second moment depends on joint pair types, equivalently on the triple of
weights of \(c,d,c+d\). The next proof object is therefore the exact pair-type
transition kernel for a uniform interleaver followed by an accumulator.

That kernel is now fixed by a four-state coefficient formula. The matrix
entry from accumulated pair state \(s\) to \(s'\) is
\(x_{s+s'}y_{s'}\). Extracting the input- and output-type coefficient from
the length-512 matrix power and dividing by the input multinomial gives the
transition probability. `verify_accumulator_pair_type_kernel.py` agrees with
exhaustive permutations on four representative length-8 types. The remaining
problem is to bound the three-stage length-512 operator, not to identify its
transition law.

An exact MacWilliams transform gives dual distance 22 for the local EBCH code
and hence for the four-block direct sum. Two independent uniform base-code
words therefore form a 21-wise independent \(\mathbb F_2^2\)-valued
coordinate process. All pair-type falling-factorial modes through degree 21
match the multinomial stationary law. Conditioning on distinct nonzero
messages introduces only the explicit mass
\(3\cdot2^{-256}-2\cdot2^{-512}\). The contraction proof should exploit this
mode annihilation; a generic second-singular-value estimate throws it away.

`analyze_accumulator_pair_chain_small.py` materializes the rank-two chain at
length 8. Its stationary residual is zero in binary64, its second weighted
singular value is 0.5174040542, and an RM(1,3) base has shell
variance-to-mean below one after two stages. This is interface validation, not
a length-512 extrapolation.

A length-512 one-word proxy shows that dual-distance mode annihilation is not
sufficient by itself. After projecting away polynomial degrees zero through
21, the generic L2 norm contracts by only 15.2928 bits after three stages.
The pair proof must exploit the target shell or the specific base pair
distribution. `analyze_accumulator_weight_chain_singular.py` and
`accumulator_weight_chain_B512_singular_probe.json` record this failed proof
interface.

The cap generator originally chose Cantelli for every shell above the failure
allocation. It now takes the better of Markov and Cantelli on every shell.
This reduces the weight-42 cap from 33 to 2 and improves the exact \(Q=1,2\)
margins, but it does not change the two-band \(Q=3\) bottleneck because the
tail majorant is attained at weight 79.

Two band refinements were tested and rejected. Splitting the tail by severity
gave at most 44.0334 bits for \(Q=3,\ldots,8\), worse than the unsplit bound.
Moving the tail boundary inward gave as much as 52.1069 sparse bits but
inflated the central majorant by roughly 32 bits per active row. Splitting the
new shoulder into its own category then failed on pure-shoulder sparse
compositions. None of these failures is evidence against the code.

A new triangle-convolution route avoids the unknown genus-two BCH enumerator. For
independent messages \(U,V\), the messages \(U,V,U+V\) are pairwise
independent and have the known ordinary spectrum. If the exact pair-kernel
shell probability is bounded by

\[
 f(\operatorname{wt}U)f(\operatorname{wt}V)
 g(\operatorname{wt}(U+V)),
\]

then Fourier expansion, Parseval's identity, and
\(\lVert\widehat g\rVert_\infty\le\mathbb E[g]\) bound its expectation by
\(\mathbb E[f(H)^2]\mathbb E[g(H)]\). The zero and diagonal pairs
are separated exactly and contribute the shell mean. A length-8 exact-kernel
test with binary64 convex factorization gives worst variance-to-mean bounds
3.0755, 1.0330, and 0.9416 after one, two, and three stages. A length-16
direct-sum test gives 19.2869, 2.7092, and 1.2642. In particular, two and
three stages beat the length-16 analogue target \(F=B=16\). This is the first
tested second-moment interface that uses only the ordinary spectrum and still
meets the desired form. The length-512 factorized kernel bound remains open.

The two-stage finite transfer was also tested under the same hypothetical
variance inflation \(F=512\). A proof event forcing all shell counts through
weights 23, 24, 25, or 26 to zero has respectively 41.5402, 41.2674,
40.9189, or 40.5002 bits of outer margin. This is a Markov event for an
unbiased sample, not rejection sampling. At RandomStepConv-M40, the matching
\(Q=1\) margins are only 32.4548, 34.2823, 36.1182, and 38.1036 bits. The
best outer-plus-\(Q=1\) union has 37.8528 bits. Thus two accumulator stages do
not meet the 40-bit target through this shell-cap interface. This is a failed
proof interface, not a lower bound on the construction's failure probability.

The inverse BA-2 sweep shows how sharp that variance lemma would need to be.
At RandomStepConv-M64 with all shells through weight 26 forced to zero, the
event-plus-\(Q=1\) margin is 40.2090 bits for \(F=1\). The largest tested
factor that still passes this necessary gate is \(F=2^{1.5}=2.8285\), with
40.0260 bits. This leaves almost no budget for \(Q\ge2\). The next tested
factor, \(F=2^{1.75}\), fails.

The triangle-convolution model now reaches length 24 after replacing dense matrix
powers with matrix--vector batches and replacing message-pair enumeration
with exact direct-sum convolution. Its one-, two-, and three-stage upper
bounds are 107.1594, 7.7268, and 1.5681. The corresponding exact ratios are
1.6389, 2.1517, and 1.2607. This strengthens the case for three stages. It
does not supply the length-512 factorization.

A proposed scalable shortcut was rejected. Projecting the four-symbol pair
kernel onto each of its three nonzero binary characters gives a valid
geometric-mean kernel inequality. Exact rational checks at lengths 8 and 16
verify the three underlying inequalities, including equality cases. Greedy
factor propagation through the combined inequality and the sharper
convolution average gives length-8 variance ratios 64.04, 474.28, and
8232.77. Joint optimization of every intermediate factor only lowers the
three-stage value to 8086.06. The projection majorant, rather than the greedy
optimizer, causes the loss; it cannot replace the genuine pair kernel.

The pair chain nevertheless admits an exact useful decomposition. Functions
of any one of the three character weights \(\operatorname{wt}(x)\),
\(\operatorname{wt}(y)\), and \(\operatorname{wt}(x+y)\) form a reducing
sector, because the common linear accumulator and its inverse commute with
all three characters. The complementary interaction sector is exactly the
information discarded by the rejected marginal reductions. At lengths 8 and
16, binary64 interaction-block norms are 0.279680 and 0.185857; both forward
and reverse leakage into the scalar sector are at roundoff scale. This gives a
clean representation target, but a generic L2 application still loses to the
initial density of a fixed rate-half code and does not prove BA concentration.

A genus-two MacWilliams LP gives a second route that does not assume the
unknown pair enumerator. Its variables are the pair-type counts of the fixed
base code. The exact ordinary primal and dual spectra fix all three character
marginals. The MacWilliams identity relates the two nonnegative pair tables.
The BA shell second moment is a linear objective in the primal table.

The length-16 direct-sum model is encouraging. Marginal constraints alone
permit variance inflation 22.9738, while adding dual-table nonnegativity
reduces every tested optimum to the exact value within binary64. The worst
exact ratio is 1.4622. This result is diagnostic because the LP is not
outward-rounded and the length-512 objective is not yet available.

Combining the LP with the existing character-projection envelope fails. That
envelope already gives variance inflation 115779.95 on the actual length-16
pair table after three stages. The LP cannot repair a pointwise bound that is
loose on the true code. The scalable route needs the exact second-order
input-output distribution of the accumulator, or a tighter shell-specific
majorant.

The primal/dual-distance route was strengthened with the optimal polynomial
moment bound. At twelve accumulator stages, a 43.7852-bit event gives primal
support (44,ldots,468) and dual distance 38. The exact Christoffel bound
still permits about (2^{145}) weight-44 words, and its M64 (Q=1) transfer
is (2^{65.34}). Thirteen stages raise the event margin to 45.1895 bits and
the dual distance to 44, but the M64 (Q=1) result remains (2^{53.25}).
This rejects every proof that uses only those two distances, including the
best degree-limited moment polynomial. It does not reject BA concentration.

A second scalable reduction uses only three one-word events. Two words of
weight (w) have difference weight at most (2\min(w,512-w)). An optimized
weighted geometric-mean bound on those events, followed by Triangle--Holder,
is exact but gives a worst variance factor (2^{256}) through every tested
depth from one to sixteen. Like the character-projection attempt, it throws
away the shared-interleaver overlap that the proof must preserve.

Supporting material consists of:

- `evaluate_ebch128x4_ba_expected_spectrum.py` and the
  `ebch128x4_ba0_{6,12,16}_B512_expected_spectra.json` receipts;
- `evaluate_ebch128x4_ba_variance_caps.py` and
  `ebch128x4_ba_variance_cap_requirements.json`;
- `evaluate_ebch128x4_ba_variance_proxy_q1q2.py`,
  `evaluate_ebch128x4_ba_variance_proxy_sparse.py`, and the seven
  `ebch128x4_ba3_variance9_*_m22.json` sparse receipts;
- `evaluate_ebch128x4_ba_variance_proxy_dense.py` and
  `ebch128x4_ba3_variance9_dense_probe_m22.json`;
- `probe_small_ba_spectrum_variance.py` and
  `small_ba32_spectrum_variance_probe.json`;
- `verify_accumulator_pair_type_kernel.py` and
  `accumulator_pair_type_kernel_small_exact.json`;
- `analyze_accumulator_pair_chain_small.py` and
  `accumulator_pair_chain_{B8,B16}_exact.json`; and
- `analyze_accumulator_pair_interaction_small.py` and
  `accumulator_pair_interaction_{B8,B16}_probe.json`; and
- `probe_genus2_macwilliams_lp_small.py` and
  `genus2_macwilliams_lp_{B8,B16}_probe.json`; and
- `probe_projection_macwilliams_lp_small.py` and
  `projection_macwilliams_lp_B16_probe.json`; and
- `analyze_accumulator_weight_chain_singular.py` and
  `accumulator_weight_chain_B512_singular{,_long}_probe.json`;
- `evaluate_ebch128x4_ba_variance_proxy_split_sparse.py`,
  `evaluate_ebch128x4_ba_variance_proxy_shoulder_sparse.py`, and their
  boundary/split/shoulder diagnostic receipts; and
- `probe_triangle_holder_pair_bound_small.py` and
  `triangle_young_pair_bound_{B8,B16,B24}_probe.json`;
- `ebch128x4_ba2_variance9_zero{23,24,25,26}_{caps,q1_m24,q1_m40}.json`
  for the two-stage forced-zero diagnostic; and
- `ebch128x4_ba2_variance9_q1_m{24,28,32,40,64}.json`,
  `ebch128x4_ba2_variance9_q1q2_m24.json`, and
  `ebch128x4_ba2_variance9_sparse_q3_*.json` for the ordinary-cap two-stage
  memory sweep; and
- `sweep_ebch128x4_ba2_variance_threshold.py` and
  `ebch128x4_ba2_variance_threshold_q1_m64.json`;
- `probe_pair_projection_holder_bound_small.py`,
  `pair_projection_young_B8_probe.json`, and
  `probe_pair_projection_global_small.py` with
  `pair_projection_global_young_B8_probe.json`, plus the exact checker
  `verify_pair_projection_kernel_inequality.py` and its B8 and B16 receipts;
  and
- `build_ebch128x4_ba_christoffel_caps.py`,
  `ebch128x4_ba{12,13}_christoffel_caps_probe.json`, and
  `ebch128x4_ba{12,13}_christoffel_q1_m64_probe.json`; and
- `evaluate_ebch128x4_ba_marginal_geometric_variance.py` and
  `ebch128x4_ba_marginal_geometric_variance_probe.json`; and
- the rejected limited-independence, minimum-distance, and column-Holder
  probes listed in the status document.

This is a conditional proof route, not a finite certificate. No performance
claim has been made for the three-stage implementation.

## One-shot high-probability random constituent

`SINGLE_RANDOM_CONSTITUENT_ONE_SHOT_STATUS.md` records the attempt to remove
the enumerated spectrum acceptance test from the repeated-random-constituent
comparator. Setup samples one uniform generator once and performs no test.
An explicit Markov--Cantelli spectrum event fails with probability at most
(2^{-42.26147}) for ([512,256]), so the sampler has more than the requested
40 bits of setup margin.

Conditional on the event, nearest-binary64 diagnostics close (Q=1) with
75.3684 bits, (Q=2) with 113.0370 bits, (Q=3\ldots8) with 54.6796
aggregate bits, and (Q=9\ldots16) with 175.1893 aggregate bits. The sampled
(Q=4096) three-band cover has 5035.5203 bits. A categorical conditioning
relaxation fails in the middle occupations because it pays the row-category
conditioning cost independently in all 512 regions. This is a failed proof
relaxation, not evidence that the code fails.

The ([1024,512]) fallback improves the sparse margins but exhibits the same
middle-occupation artifact under that rejected relaxation.

The (B=512) route subsequently closed. Monotonicity of the RandomStepConv
moment merges the symmetric low and high tails into one defect category. A
shared-category recurrence covers (Q=3,ldots,159) without repeating the
category-conditioning loss. A convex cover handles every integer composition
for (Q=160,ldots,4096).

The outward receipts give 51.6439589890 conditional bits. The exact spectrum
event gives 42.2614729204 bits. Their unconditional union has
42.2593129905 bits, above the requested 40-bit target. The construction uses

\[
 (k,N,D,M)=(2^{20},2^{21},228{,}590,22)
\]

and proves (d_{\min}\ge D\), hence relative distance at least 10.9%. Setup
samples one uniform (256)-by-(512) matrix and performs no spectrum test.
See `FINITE_K20_ONE_SHOT_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md` and its
manifest. This result supersedes the incomplete status in the preceding
paragraphs but retains the failed relaxations as proof-search history.

Supporting material for the closure consists of:

- `evaluate_single_random_constituent_shared_two_band.py` and the five
  `single_random_constituent_B512_shared_two_band_q*_s22.json` sparse
  receipts;
- `evaluate_single_random_constituent_shared_two_band_dense.py`, including
  the failed (Q=129) categorical probe and the successful sampled dense
  slices;
- `diagnose_single_random_constituent_shared_two_band_cover.py` and
  `single_random_constituent_B512_shared_two_band_dense_cover_s22.json`;
- `certify_single_random_constituent_{sparse,dense}_outward.py` and their
  outward receipts;
- `certify_single_random_constituent_combined.py` and
  `single_random_constituent_B512_combined_outward_s22.json`;
- `FINITE_K20_ONE_SHOT_RANDOM512_RANDOMSTEP_CONV_MANIFEST.json`, generated by
  `build_single_random_constituent_manifest.py`.

### Compact pairwise-systematic outer

The dense random generator can be replaced without changing the conditional
SPIN transfer. Identify a message row with
\(u\in\mathbb F_{2^{256}}\), sample one pair
\(a,b\gets\mathbb F_{2^{256}}\), and use

\[
 u\longmapsto(u,au+bu^2).
\]

For distinct nonzero \(u,v\), the determinant \(uv(u+v)\) is nonzero.
Consequently, their parity halves are exactly independent and uniform. The
shell indicators are pairwise independent, so every shell satisfies
\(\operatorname{Var}(A_w)\le\mathbb E[A_w]\). The systematic half removes
rank failure.

The exact shell mean is

\[
 2^{-256}\bigl(\binom{512}{w}-\binom{256}{w}\bigr)
\]

for \(w\le256\), and \(2^{-256}\binom{512}{w}\) for \(w>256\). The verifier
checks these identities and applies Markov--Cantelli to the original cap
vector. The cap-event margin remains 42.2614729204 bits. Reusing the outward
cap-conditional sparse and dense transfers gives 42.2593129905 bits overall.

This construction samples 512 outer bits and performs two
\(\mathbb F_{2^{256}}\) multiplications per row. No rejection test or shell
enumeration occurs. This is a proved operation count, not a benchmark.
Supporting material consists of:

- `FINITE_K20_PAIRWISE_SYSTEMATIC_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md`;
- `certify_pairwise_systematic_constituent_combined.py` and
  `pairwise_systematic_constituent_B512_combined_outward_s22.json`; and
- `FINITE_K20_PAIRWISE_SYSTEMATIC_RANDOM512_RANDOMSTEP_CONV_MANIFEST.json`.

### Fixed-base linearized-mixer fallback

The same pairwise-independence mechanism can preserve the fixed EBCH base.
Identify the 512 output bits of four EBCH \([128,64,22]\) blocks with
\(x\in\mathbb F_{2^{512}}\). Sample
\(a,b\gets\mathbb F_{2^{512}}\) and apply

\[
 x\longmapsto ax+bx^2.
\]

For distinct nonzero \(x,y\), the determinant \(xy(x+y)\) is nonzero. Thus
their images are independent and uniform. The restriction to the fixed
\([512,256]\) base code loses rank with probability

\[
 \frac{1+(2^{512}-1)(2^{256}-1)}{2^{1024}}<2^{-256}.
\]

The original shell caps and conditional transfer apply unchanged. The
outward certificate again gives 42.2593129905 bits overall at 10.9%
distance. This route uses a 1,024-bit mixer sampler and two
\(\mathbb F_{2^{512}}\) multiplications per row. It needs zero accumulator
stages, but it changes the requested BA-only construction.

Supporting material consists of:

- `FINITE_K20_EBCH128X4_LINEARIZED_MIXER_RANDOMSTEP_CONV_CERTIFICATE.md`;
- `certify_fixed_outer_linearized_mixer_combined.py` and
  `fixed_outer_linearized_mixer_B512_combined_outward_s22.json`; and
- `FINITE_K20_EBCH128X4_LINEARIZED_MIXER_RANDOMSTEP_CONV_MANIFEST.json`.

### Universal BA mixing fallback

A long-chain BA route now reduces the complete finite theorem to one precise
spectral lemma. Start with the direct sum of four fixed extended BCH
([128,64,22]) codes. Sample (t) independent uniform permutations of its
512 coordinates, applying the zero-state prefix accumulator after each one.
These permutations are sampled once. The resulting single constituent is
repeated in all 4096 outer rows.

The conditional lemma asserts that the second stationary-(L_2) singular
values of the one-word and ordered-distinct-pair operators for one such stage
are at most (63/1000). Under that lemma, chi-square contraction bounds the
first and second factorial moments of every shell count for every fixed
([512,256]) base code. The argument therefore needs neither the BCH
genus-two enumerator nor a random-code approximation.

The exact cap vector from the certified random-([512,256]) proof is held
fixed. BA-127 violates its union-bound target: the cap-event margin is only
36.7702 bits. BA-128 is the first passing stage count and has an outward
cap-event margin of 40.3474252482 bits. The existing outward sparse and dense
transfers depend only on this cap vector and already cover every
(1\le Q\le4096). Combining their 51.6439589890 conditional bits with the
BA-128 event gives 40.3468518020 bits overall at
(D=228{,}590), or relative distance at least 10.9%.

This is a conditional outward certificate, not a completed proof. Its only
mathematical gap is the singular-value lemma. Binary64 evidence gives the
one-word value 0.0626241431 at length 512. Exact pair kernels followed by
binary64 linear algebra at even lengths 4 through 20 place the complementary
interaction norm below (3/B); at length 20 it is 0.141808598. This suggests
proving a scalar-character reduction plus an interaction bound, but the
finite data are not extrapolated.

An outward Hilbert--Schmidt calculation now proves the weaker one-word bound
\(\sigma_2<0.11\) at length 512. It does not prove the \(0.063\) hypothesis
needed for BA-128 and does not control the pair interaction sector. Under a
matching pair bound of \(0.11\), the diagnostic cap calculation first passes
at BA-161 with 41.5690 event bits. This supplies a rigorously grounded
one-word fallback but leaves the decisive pair statement conditional.

Supporting material consists of:

- `evaluate_ba_pair_singular_mixing_sweep.py` and
  `ba_pair_singular_mixing_sweep_B512.json`;
- `evaluate_ba_pair_singular_random_cap_event.py` and
  `ba_pair_singular_random_cap_event_B512.json`;
- `certify_ba128_pair_singular_random_cap_event_outward.py` and
  `ba128_pair_singular_random_cap_event_outward_B512.json`;
- `certify_ebch128x4_ba128_randomstepconv_combined.py` and
  `ebch128x4_ba128_randomstepconv_combined_outward_s22.json`; and
- `analyze_accumulator_pair_interaction_small.py` and the even-length
  `accumulator_pair_interaction_B*_probe.json` receipts through \(B=20\);
- `certify_accumulator_oneword_singular_outward.py` and
  `accumulator_oneword_singular_B512_outward.json`; and
- `ba_pair_singular_random_cap_event_B512_lambda011.json`, the conditional
  BA-161 diagnostic under a common \(0.11\) one-word and pair bound.

The 128 stages are a proof fallback and are not performance competitive. The
short-stage optimization remains separate: BA-3 still needs a substantially
more constituent-specific second-order argument, and BA-2 appears fragile
under the present cap interface.

## Global two-sided regular Expand--Convolute alternative

The intended use of the expander work is now recorded separately in
`BLOCK_EXPAND_CONSTITUENT_ROUTE.md`. It replaces the repeated BA constituent
inside Structured SPIN. The remainder of this section records the distinct
end-to-end expander alternative.

The local constituent route now has an exact regional pair kernel. For input
nonzero pair counts \((p,q,r)=(n_{10},n_{01},n_{11})\), a region of length
\(\ell\), and output pair type \(m\), its labeled assignment count is

\[
 p!q!r!\binom{\ell}{m}
 [X^pY^qZ^r]\prod_{a,b}F_{ab}(X,Y,Z)^{m_{ab}},
 \qquad
 F_{ab}=\frac14\sum_{s,t=\pm1}s^at^b e^{sX+tY+stZ}.
\]

The fourteen independent region laws compose with the already verified
common-permutation accumulator pair kernel. This gives an exact finite
algorithm for every shell second factorial moment of Block Expand--\(t\).
An exact small-instance verifier agrees with direct
enumeration in 176 cases and 3,036 output-type entries. The remaining
obligation is a non-materialized length-512 outward implementation: the full
four-symbol state space has 22,632,705 types, so a dense transition matrix is
not a viable certificate representation.

The present finite gate separates the layer counts. Keeping the random-code
cap support 42 through 470, Block Expand--2 has \(2^{-32.4286}\) expected
positive mass outside the support and cannot reach 40 setup bits. The
corresponding values for \(t=3,4,5\) are \(2^{-51.6918}\), \(2^{-52.4060}\),
and \(2^{-52.6725}\). The separate kernel term is \(2^{-57.7008}\).

Conditional on \(\operatorname{Var}(A_w)\le2\mathbb E[A_w]\), recentered
Cantelli caps with a uniform \(2^{-51}\) shell budget inflate the
low/central/high two-band majorants by 0.3653/0.3378/0.3618 bits at \(t=3\),
0.1285/0.0964/0.1270 bits at \(t=4\), and 0.0627/0.0309/0.0625 bits at
\(t=5\). Exact analogues at \((K,B)=(4,8),(6,12),(8,16)\) give five-stage
maximum variance factors 1.3505, 1.2397, and 1.1659. This is evidence for the
factor-two target, not a length-512 proof.

The final cap calculation spends half of the available combined 40-bit budget
on the 429 positive shells. At \(t=5,F=2\), its per-shell allowance is
\(2^{-49.7455}\). The resulting low/central/high band majorants are 0.0082,
0.0018, and 0.0084 bits below the frozen random-code values. Therefore the
existing outward \(Q\ge3\) transfers apply unchanged. New exact-shell outward
calculations give 67.9441 bits at \(Q=1\) and 105.3578 bits at \(Q=2\). The
complete conditional distance failure is \(2^{-51.6439589890}\).

The mean and cap-event arithmetic is now outward. Exact integer regional
occupancy polynomials and five exact binomial accumulator transitions give
42.6283 outer-event bits under the factor-two variance premise. Combining
this with the conditional transfer gives 42.6254 bits. The shell means,
kernel term, zero-cap tail, caps, and all occupations are therefore closed.

The remaining proof is now isolated: prove

\[
 \operatorname{Var}(A_w)\le2\mathbb E[A_w]
 \quad(1\le w\le512).
\]

The regional pair identity and accumulator pair kernel define these moments
exactly. A non-materialized length-512 evaluator is still required.

The expander-code workstream supplies a direct finite alternative that avoids
the repeated-constituent spectrum problem. Its parent message length is

\[
 k_+=1{,}048{,}585,
\]

the smallest integer above \(2^{20}\) that is odd and divisible by five. A
binary two-sided regular expander has left degree 10 and right degree 5. It is
sampled by ten independent permutations, one for each region of
\(k_+/5=209{,}717\) intermediate coordinates. A memory-15 wrapped random
convolution follows the expander. All permutations and convolution
coefficients are sampled once and are shared by every message.

The existing outward first-moment verifier originally evaluated message
supports 1 through 16 exactly. At the present cutoff, the following
coefficient relaxation began at support 17 and limited the proof to 34.7518
bits. Extending the exact regional transfer through support 32 changes no
code parameter and gives

\[
 \Pr[d_{\min}(G_+)\le228{,}607]
 <2^{-50.2054082202}.
\]

The independent standard-library checker verifies that the frozen support
partition covers every nonzero parent support exactly once. The 192-bit Arb
verifier evaluates every probability bound outwards. The weight-one term is
the dominant contribution.

A deterministic wrapper fixes nine parent input coordinates to zero and
punctures eighteen fixed output coordinates. Thus

\[
 (k,N)=(2^{20},2^{21}),\qquad d_{\min}\ge228{,}590,
\]

so the relative distance is at least \(0.10900020599365234\). Shortening does
not decrease distance, puncturing loses at most eighteen, and the remaining
positive distance proves that puncturing preserves dimension. The wrapper
adds no failure probability.

This is a certified finite theorem and a linear-time operation count for both
ordinary and transposed encoding. It is
not a certificate for BCH--BA, RM2Sub, or structured SPIN: the outer map is
one global regular expander and the inner map is a wrapped random convolution.
It also supplies no decoder theorem. Its online constant has not been
measured; ten sparse edge updates per input coordinate plus the memory-15
recurrence could exceed the current 11 ms target.

Supporting material consists of:

- `EXPANDER_CODE_ALTERNATIVE_AUDIT.md`;
- `BLOCK_EXPAND_CONSTITUENT_ROUTE.md`,
  `evaluate_block_expand_accumulate_spectrum.py`, and
  `block_expand_accumulate_512_256_d14_t0_5_spectrum.json`;
- `verify_block_expand_pair_region_small.py` and
  `block_expand_pair_region_small_exact.json`;
- `evaluate_block_expand_cap_budget.py` and
  `block_expand_cap_budget_diagnostic.json`;
- `evaluate_block_expand_q1_cap_transfer.py` and
  `block_expand3_q1_cap_transfer_diagnostic.json`;
- `analyze_block_expand_pair_moments_small.py` and
  `block_expand_pair_moments_K4_B8_exact.json`,
  `block_expand_pair_moments_K6_B12_exact.json`, and
  `block_expand_pair_moments_K8_B16_exact.json`;
- `certify_block_expand5_F2_conditional_transfer.py` and
  `block_expand5_F2_conditional_transfer_outward.json`;
- `certify_block_expand5_mean_cap_outward.py` and
  `block_expand5_mean_cap_outward.json`;
- `build_block_expand5_conditional_manifest.py` and
  `FINITE_K20_BLOCK_EXPAND5_F2_CONDITIONAL_MANIFEST.json`;
- `FINITE_K20_EXPANDER_EC_D10_M15_CERTIFICATE.md`;
- `generate_expander_ec_k20_candidate.py` and
  `expander_ec_d10_m15_k20_parent_candidate.json`;
- `certify_expander_ec_k20_wrapper.py` and
  `expander_ec_d10_m15_k20_wrapper_outward.json`; and
- `build_expander_ec_k20_manifest.py` and
  `FINITE_K20_EXPANDER_EC_D10_M15_MANIFEST.json`.

The combined receipt and manifest record absolute paths and SHA-256 hashes
for the independent checker and outward verifier in the expander-code
worktree. Vendoring those verifier dependencies and benchmarking the encoder
are the two integration tasks. Neither is an open mathematical obligation.

## Pure sparse-mixer--accumulate reset

The current route is distinct from the earlier Block Expand-\(t\) experiment.
It does not place one expander before a chain of permutation--accumulator
stages. Instead, every proposed mixing stage is a fresh sparse linear map.
The encoder data flows as

\[
 u\longmapsto E_0u\longmapsto AE_0u
 \longmapsto E_1AE_0u\longmapsto AE_1AE_0u
 \longmapsto\cdots.
\]

where \(E_0:\mathbb F_2^{256}\to\mathbb F_2^{512}\) supplies the rate-half
embedding and later maps are square. Each map is sampled once, and the same
composite constituent is reused in every SPIN outer row.

The first proof ensemble fixes the degree of every output vertex. Each output
independently selects a uniform subset of the input coordinates and emits its
parity. This ensemble is a diagnostic baseline, not a frozen construction.
Left-regular and biregular replacements remain open alternatives.

The exact one-word and ordered-pair laws are now derived. For an input pair,
the distribution of one output pair is determined by three Krawtchouk ratios
at the weights of \(x\), \(y\), and \(x+y\). Independent output neighborhoods
make the pre-accumulator pair symbols independent. A four-state common
accumulator then gives the exact joint shell kernel. The base pair-type
multiplicities are multinomial because the message space is the full
\(\mathbb F_2^{256}\). Thus this route removes the unknown BCH genus-two
enumerator.

Four small instances match exhaustive ensemble enumeration with exact
rational arithmetic. Exact multistage models through \((k,n)=(5,10)\) have
worst positive-shell variance-to-mean ratios below \(1.70\). These checks
validate the finite kernel implementation but do not imply a length-512
variance bound.

The first full-size diagnostic tests the zero caps outside weights 42 through
470. It uses a bounded, efficient rank-tested setup with sixteen attempts.
Binary64 first moments give the following frontier:

| stages | output degrees | XORs per constituent | diagnostic outer bits |
|---:|---:|---:|---:|
| 1 | \(33\) | 16,895 | 44.6379 |
| 2 | \(17,3\) | 10,238 | 43.9504 |
| 3 | \(13,3,3\) | 9,725 | 42.4540 |
| 3 | \(15,2,3\) | 10,237 | 46.5549 |

These values control only kernel and nonzero extreme-shell events. They do
not prove simultaneous central-shell caps. The remaining central obligation
is a scalable length-512 bound on the pair kernel, with the prior
\(\operatorname{Var}(A_w)\le512\mathbb E[A_w]\) interface as the first
sufficient target.

Supporting material is contained in pure_expander_accumulate/:

- README.md;
- PAIR_KERNEL_AND_FIRST_GATE.md;
- verify_pair_kernel_small.py and four exact small receipts;
- analyze_multistage_pair_moments_small.py and three exact small-chain
  receipts; and
- the one-, two-, and three-stage extreme-shell sweep scripts and receipts.

### One-stage concentration reduction

The first concentration attack now uses the one-stage degree-33 candidate.
Fourier inversion converts its shell count into a weighted sum of indicators
that selected sparse-map rows XOR to zero. The resulting variance is one
quadratic form of an exact covariance operator. The covariance entries are
explicit Krawtchouk sums and depend only on the two row-subset weights and
their intersection. The operator is therefore invariant under row
permutations.

Odd right degree creates a deterministic parity mode: the all-ones message
maps through the accumulator to the alternating word of weight 256. This is
a real property of the ensemble and cannot be discarded. Moderate-size
calculations also show large signed cancellations between complementary
message-difference weights. Absolute termwise pair bounds are therefore not
a viable proof method.

A global covariance-operator norm is also insufficient. Two-row collision
events have probability \(\binom{256}{33}^{-1}\), whose base-two exponent is
only 138.1822, while a global norm proof of the factor-512 target would need
roughly 247 bits. The viable route is the standard symmetric-group block
decomposition of functions on row subsets, combined with the accumulator
Fourier energy in each block.

Exact integer diagnostics show that the constant-weight row walk has
nontrivial spectral mass below \(2^{-40}\) after 114 steps. For the central
output shell, the nonconstant accumulator Fourier energy in row-subset
weights 1 through 64 is about \(2^{-7.4676}\) of the total. Away from the
center it is much smaller. These facts make a low-mode/high-mixing split
plausible, but they do not yet bound the full quadratic form.

Additional supporting material in pure_expander_accumulate/ consists of:

- ONE_STAGE_VARIANCE_ROUTE.md;
- probe_one_stage_variance_scaling.py and its one-stage scaling receipts; and
- analyze_dual_walk_and_accumulator_energy.py with
  dual_walk_accumulator_energy_K256_B512_r33.json.

An explicit formula for every row-permutation block has now been derived.
Two small instances reproduce both exhaustive rational shell variance and
the full covariance eigenvalue multiset. The three-index character moments
also admit an exact rank-257 factorization through the radial row walk; 165
small rational identities pass.

A cheaper weight-level norm bound initially lost a factor at most 6.72 in
the smallest tests, but it fails at moderate sizes. At \((24,48,5)\) it gives
4611.74 against an exact variance factor 2.7906, and at \((32,64,7)\) it gives
15854.46 against 4.6473. This relaxation is rejected. The length-512 proof
must retain consistent irreducible sectors or the signed block quadratic
form. No realized-spectrum certificate is claimed yet.

The additional validation artifacts are
verify_subset_covariance_blocks_small.py,
subset_covariance_blocks_K3_B4_r1_exact.json,
subset_covariance_blocks_K4_B5_r3_exact.json,
probe_subset_block_norm_bound_small.py, and the
subset_block_norm_bound_*_probe.json receipts. The radial compression and
failed moderate-size relaxation are recorded by
verify_radial_triple_factorization_small.py,
radial_triple_factorization_K4_r3_m8_exact.json,
probe_radial_block_bound.py, and the radial_block_bound_* receipts.

The irreducible-sector projector energies now have an exact Johnson-scheme
formula. Exact rational tests reconstruct every row-subset level energy and
the full Parseval energy. The sector-norm variance bound improves the
weight-level bound in both tested instances. Its remaining input would be an
explicit accumulator correlation table indexed by row-subset weight and
intersection size. This exact sector-energy route is now the fallback to the
Walsh-conjugated gate below. The validation artifacts
are verify_sector_projection_energy_small.py and the
sector_projection_energy_* receipts.

### Walsh-conjugated low-shell gate

An exact basis change gives a cheaper route for the shells used by the
finite transfer. If \(W\) is the Walsh matrix and
\(f_w(z)=\mathbf 1\{\operatorname{wt}(Az)=w\}\), then the dual coefficient
vector is \(c_w=Wf_w\). Conjugating the covariance operator by
\(2^{-n/2}W\) preserves every symmetric-group sector. At rate one half, the
variance is the quadratic form of the conjugated operator against \(f_w\).

The level norm of \(f_w\) is an exact run count. For low shell weight \(w\),
it vanishes above level \(2w\). Applying the valid level-norm relaxation in
this basis gives 1.07079 at \((k,n,r,w)=(32,64,7,5)\), versus the exact
ratio 0.999994. At the matched next-size analogue
\((64,128,9,10)\), it gives 10.4285. Both values are nondirected binary64
diagnostics. The latter is below the factor-512 target but does not
extrapolate to length 512.

The two-band transfer can use the exact constituent size for the aggregate
central mass. The immediate variance target is therefore the defect band
42 through 79 and its complement. The first target shell uses primal levels
only through 84. The primary remaining implementation task is to evaluate
the needed conjugated block entries at length 512 without forming every
dense block, then add outward rounding. The full accumulator correlation
table remains a fallback.

For a fixed ordered message pair with one-row output matrix \(P\), the primal
pair kernel is \(P^{\otimes n}\). Its exact sector-\(j\) block is
\(\det(P)^j\operatorname{Sym}^{n-2j}(P)\). Aggregating these blocks over
message-pair types gives the Walsh-conjugated covariance directly and avoids
forming the full dual blocks. At \((4,8,3)\), direct aggregation over all 256
ordered pairs agrees with the conjugated exact covariance blocks to relative
binary64 residual below \(1.4\times10^{-15}\). The remaining compression is
to sum the symmetric-power entries over the 2.8 million target pair types
without constructing one dense matrix per type.

An exact diagonal insertion lemma removes the remaining high-sector
calculation. If \(s_{m,u}(P)\) is a diagonal symmetric-power entry, then

\[
 |\det P|s_{m,u}(P)\le s_{m+2,u+1}(P).
\]

Consequently, every sector diagonal at fixed primal level is bounded by the
maximum of sectors zero, one, and two. Positive semidefiniteness then bounds
all cross-level entries by the geometric means of those diagonals.

The target implementation streams all 2,862,209 message-pair types. It
centers sector zero locally and retains only positive local differences. It
retains only positive determinant contributions in sector one. Sector two is
already a positive sum. These operations give a cancellation-free upper
bound in exact arithmetic.

At \((k,n,r,w)=(256,512,33,42)\), the resulting nondirected diagnostic
variance factor is 31.0225. The required factor is 512, leaving a factor of
about 16.5. The exact reduction is proved; the numerical value is not yet an
outward certificate. Levels through 159 suffice for all defect shells 42
through 79 and 433 through 470.

Directed binary64 interval arithmetic now certifies the individual
sector-zero upper bounds 1.277282466 and 206.661452 at primal levels 80 and
84. Together with nondirected lower-level inputs, these endpoints give a
mixed shell-42 diagnostic factor 33.8177. A complete outward shell receipt
still requires all lower levels.

The pairwise-positive sector-zero relaxation does not scale to later shells.
Its outward bounds are approximately \(2.72\cdot10^6\) at level 90 and
\(6.97\cdot10^{12}\) at level 100. The locally centered signed diagnostics
remain 1.00000023 and 1.00000048. Complementing either message gives an exact
four-term orbit. At level 100, orbit grouping reduces the nondirected
positive and absolute mass to 1.0000032. Small exact models contain negative
orbit contributions, so orbit positivity is not assumed.

An exact Hoeffding expansion now replaces that failed relaxation. If
\(\delta_{x,y}=\det(P_{x,y})\), each fixed-order aggregate is a quadratic
form in the Schur power \((\delta_{x,y}^{\,j})\). The determinant matrix is a
covariance matrix, so every Schur power is positive semidefinite. The next
certificate task is to complete or upper-bound these radial quadratic forms
through level 159. Level 100 is now complete: the merged 2,048-bit Arb
receipt encloses its sector-zero diagonal as
1.000003554708789325353157978375139901... with radius below
\(2.14\cdot10^{-111}\). A generic Arb exponentiation edge case caused the
first order-512 calculation to return `nan`; direct interval multiplication
fixed it without changing the formula. The remaining work also includes the
outward sector-one and sector-two bands.

The transform inside each Hoeffding order now has a second exact form:

\[
 F_{j,\ell,p}(t)=2^{k-(n-j)}\sum_{s=0}^{n-j}
 \binom{n-j}{s}K_p^{(n)}(j+s)v_t(s+j-\ell).
\]

The small exact verifier checks all 1,980 instances of this identity. The
polynomial-convolution Arb implementation matches the scalar implementation
exactly on the complete small model. It is slower for a single target level:
81.4 seconds rather than 44.2 seconds through order 128 at
\((256,512,33,p=100)\). A dense batched transform was also slower, and a
512-bit run lost its enclosure at high orders; the default remains the
1,024-bit scalar transform.

Complement grouping remains a valid upper-bound route without an orbit-
positivity theorem. One may enclose each signed four-term orbit and retain
only its positive upper endpoint. A directed binary64 prototype bounds the
small \((4,8,3,p=2)\) diagonal by 0.8136835098268285, agreeing tightly with
the nondirected positive-orbit value 0.8136835098266602. The exact diagonal
is 0.61776065826416015625. The first target implementation was stopped after
33,153 of 2,862,209 pair types took 72 seconds, implying roughly 1.7 hours
per level. Canonicalizing under complements and message swap reduces the
calculation to 361,985 representatives. A second rounding fix constructs
each determinant from its exact integer numerator; subtracting rounded bias
products had produced the uselessly loose upper bound
\(7.5011\cdot10^{18}\). The corrected target receipt proves

\[
 D^{(0)}_{100,100}\le 2.45522037042628
\]

in 479.6 seconds. The bound is conservative relative to the independent Arb
value 1.0000035547087893..., but it is far below 512. The orbit inequality is
therefore viable. A compiled version should update symmetric-power
coefficients by recurrence and amortize the traversal across multiple levels.

A full-trace shortcut was rejected. The normalized Walsh transform preserves
trace, but the dual sector-zero block contains enormous low-level collision
energy. Nondirected long-double evaluation also becomes negative from
catastrophic cancellation. The trace does not control the needed primal band.

Supporting material is:

- `pure_expander_accumulate/probe_walsh_conjugated_level_bound.py`;
- `pure_expander_accumulate/verify_primal_schur_power_factorization_small.py`;
- `pure_expander_accumulate/primal_schur_power_K4_B8_r3_probe.json`;
- `pure_expander_accumulate/probe_primal_schur_diagonal_target.py`;
- `pure_expander_accumulate/evaluate_primal_schur_diagonal_shell_bound.py`;
- `pure_expander_accumulate/primal_schur_diagonal_K4_B8_r3_validation.json`;
- `pure_expander_accumulate/primal_schur_diagonal_K32_B64_r7_validation.json`;
- `pure_expander_accumulate/primal_schur_diagonal_K256_B512_r33_j0_2_p84_probe.json`;
- `pure_expander_accumulate/primal_schur_diagonal_K256_B512_r33_j0_p84_centered_positive_probe.json`;
- `pure_expander_accumulate/primal_schur_diagonal_K256_B512_r33_j1_condition_probe.json`;
- `pure_expander_accumulate/primal_schur_shell_bound_K256_B512_r33_w42_probe.json`;
- `pure_expander_accumulate/certify_primal_schur_diagonal_outward.py`;
- `pure_expander_accumulate/merge_primal_schur_diagonal_outward.py`;
- `pure_expander_accumulate/primal_schur_diagonal_K256_B512_r33_p80_84_outward.json`;
- `pure_expander_accumulate/primal_schur_diagonal_K256_B512_r33_j0_p90_p100_outward_probe.json`;
- `pure_expander_accumulate/probe_sector_zero_complement_orbits.py`;
- `pure_expander_accumulate/sector_zero_complement_orbit_K256_B512_r33_p100_probe.json`;
- `pure_expander_accumulate/verify_sector_zero_hoeffding_small.py`;
- `pure_expander_accumulate/sector_zero_hoeffding_K4_B8_r3_exact.json`;
- `pure_expander_accumulate/sector_zero_hoeffding_K4_B8_r3_polynomial_exact.json`;
- `pure_expander_accumulate/evaluate_sector_zero_hoeffding_radial_arb.py`;
- `pure_expander_accumulate/merge_sector_zero_hoeffding_order_receipts.py`;
- `pure_expander_accumulate/sector_zero_hoeffding_radial_K4_B8_r3_p2_arb.json`;
- `pure_expander_accumulate/sector_zero_hoeffding_radial_K256_B512_r33_p100_j8_arb.json`;
- `pure_expander_accumulate/sector_zero_hoeffding_radial_K256_B512_r33_p100_j64_arb.json`;
- `pure_expander_accumulate/sector_zero_hoeffding_radial_K256_B512_r33_p100_j160_arb.json`;
- `pure_expander_accumulate/sector_zero_hoeffding_radial_K256_B512_r33_p100_j256_arb.json`;
- `pure_expander_accumulate/sector_zero_hoeffding_radial_K256_B512_r33_p100_j512_arb.json`;
- `pure_expander_accumulate/sector_zero_hoeffding_radial_K256_B512_r33_p100_j512_only_arb.json`;
- `pure_expander_accumulate/sector_zero_hoeffding_radial_K256_B512_r33_p100_complete_arb.json`;
- `pure_expander_accumulate/sector_zero_hoeffding_radial_K256_B512_r33_p100_j128_polynomial_arb.json`;
- `pure_expander_accumulate/certify_sector_zero_complement_orbits_outward.py`;
- `pure_expander_accumulate/sector_zero_complement_orbits_K4_B8_r3_p2_outward.json`;
- `pure_expander_accumulate/sector_zero_complement_orbits_K256_B512_r33_p100_canonical_exactdet_outward.json`;
- `pure_expander_accumulate/probe_sector_zero_trace.py`;
- `pure_expander_accumulate/sector_zero_trace_K256_B512_r33_probe.json`;
- `pure_expander_accumulate/walsh_conjugated_level_bound_K4_B8_r3_w6_probe.json`;
- `pure_expander_accumulate/walsh_conjugated_level_bound_K24_B48_r5_w31_probe.json`;
- `pure_expander_accumulate/walsh_conjugated_level_bound_K32_B64_r5_w5_probe.json`;
- `pure_expander_accumulate/walsh_conjugated_level_bound_K32_B64_r7_w5_probe.json`;
- `pure_expander_accumulate/walsh_conjugated_level_bound_K32_B64_r5_w40_probe.json`;
- `pure_expander_accumulate/walsh_conjugated_level_bound_K32_B64_r7_w40_probe.json`; and
- `pure_expander_accumulate/walsh_conjugated_level_bound_K64_B128_r9_w10_probe.json`.

## One-stage sparse-EA outer: completed finite certificate

The one-stage candidate is now a proved finite SPIN variant. One setup
attempt samples 512 independent weight-33 row vectors in
\(\mathbb F_2^{256}\), applies the zero-initialized accumulator, and obtains a
\(256\)-to-\(512\) constituent. Setup makes at most 16 independent attempts
and accepts the first full-rank constituent. The accepted constituent is used
in all 4,096 outer rows.

The final variance proof does not use the failed diagonal relaxation. For
each ordered pair \((x,y)\), it retains the largest of the three character
biases at \(x\), \(y\), and \(x+y\). The other two biases give a pointwise
likelihood factor relative to a one-parameter reference law. Arb evaluates
the reference shell probabilities at 1,536-bit precision.

Direct outward subtraction closes the defect shells. The certified maximum
variance-to-mean factor is 245.145100 on weights 42 through 79 and 433 through
470. The high range is evaluated directly; no realized-spectrum complement
symmetry is assumed.

Central reference ratios are too close to one for binary64 subtraction. A
second verifier asks Arb for their absolute deviations from one. It then
bounds covariance excess using only nonnegative terms. This covers all
weights 80 through 432. The central variance factors can exceed 512, but the
normalized pair excess remains near \(2^{-61}\). Cantelli therefore adds only
about 3% relative cap slack near the center.

The integer cap for each positive shell is the smaller of its independent
\(2^{-51}\) cap and the frozen three-band envelope cap. This allocation keeps
the exact \(Q=1,2\) transfer sharp and preserves the frozen band hypotheses
for every \(Q\ge3\). After rank conditioning, the spectrum event has
41.3632966090 bits of margin. The conditional distance calculation has
51.6422972013 bits.

At

\[
 k=2^{20},\qquad N=2^{21},\qquad D=228{,}590,
\]

the complete outward union bound is

\[
 \Pr[d_{\min}<D\text{ or setup aborts}]
 <2^{-41.3621359294}.
\]

Thus the variant has certified relative distance at least
10.9000205994%. Its inner remains RandomStepConv-M22, and its routing retains
the uniform row-coordinate and region permutations. This result is not a
certificate for the frozen Structured SPIN outer, interleaver, or RM2Sub
inner.

The earlier Schur-sector program is also complete at the requested levels.
Its outward receipts cover sector zero through level 159, sector one through
level 159, and sector two from level 2 through 159. The assembled diagonal
relaxation closes shells 42 through 53. A lower-bound receipt proves that it
first exceeds 512 at shell 65. These facts explain the switch to the
dominant-character proof.

The central verifier now supports independent shell workers. An eight-worker
Peach run completed in 84.7095 seconds, compared with 771.2017 seconds for
the serial local run. Every shell row and the claim object matched exactly.

Supporting material is:

- `ONE_STAGE_SPARSE_EA_RANDOMSTEPCONV_CERTIFICATE.md`;
- `audit_one_stage_sparse_ea_certificate.py` and
  `ONE_STAGE_SPARSE_EA_CERTIFICATE_MANIFEST.json`;
- `certify_one_stage_sparse_ea_caps_outward.py` and
  `one_stage_sparse_ea_K256_B512_r33_spectrum_caps_outward.json`;
- `certify_one_stage_sparse_ea_transfer_outward.py` and
  `one_stage_sparse_ea_K256_B512_r33_randomstepconv_M22_d109_outward.json`;
- `pure_expander_accumulate/certify_dominant_character_likelihood_outward.py`;
- `pure_expander_accumulate/dominant_character_likelihood_K256_B512_r33_w42_79_outward.json`;
- `pure_expander_accumulate/dominant_character_likelihood_K256_B512_r33_w433_470_outward.json`;
- `pure_expander_accumulate/certify_dominant_character_deviation_outward.py`;
- `pure_expander_accumulate/dominant_character_deviation_K256_B512_r33_w80_432_outward.json`;
- `pure_expander_accumulate/verify_dominant_character_deviation_small.py` and
  `pure_expander_accumulate/dominant_character_deviation_K4_B8_r3_all_shells_exact.json`;
- `pure_expander_accumulate/dominant_character_deviation_peach_workers8_validation.json`;
- `pure_expander_accumulate/sector_zero_complement_orbits_K256_B512_r33_p1_159_rust20_outward.json`;
- `pure_expander_accumulate/primal_schur_sectors1_2_K256_B512_r33_p1_159_rust20_outward.json`;
- `pure_expander_accumulate/primal_schur_sector2_K256_B512_r33_p2_159_mass_rust20_outward.json`; and
- `pure_expander_accumulate/diagonal_route_obstruction_K256_B512_r33_w42_79_outward.json`.

The next smallest useful task is to decide whether to transfer this mechanism
to the lower-XOR multistage outer or replace RandomStepConv-M22 with a
practical inner first.

## Depth-two sparse-EA investigation

The first lower-XOR follow-up uses

\[
  C_2=A E_1 A E_0.
\]

Setup samples the two sparse maps independently once and reuses the pair at
all outer positions. The complete binary64 mean-spectrum calculation compares
three candidates with every integer cap in the closed degree-33 theorem.

| Degrees | XORs | Zero-cap tail margin | Minimum cap/mean headroom | Uniform variance factor allowed at 40 end-to-end bits |
|---:|---:|---:|---:|---:|
| \((17,3)\) | 10,238 | 43.9504 bits | 0.04213 bits | 53.81 |
| \((19,3)\) | 11,262 | 49.4771 bits | 0.04213 bits | 57.46 |
| \((17,5)\) | 11,262 | 45.9211 bits | 0.04213 bits | 57.04 |

All three mean spectra lie below every frozen cap. Therefore, a proof of the
same simultaneous cap event would reuse the existing conditional SPIN
transfer without a new occupation calculation.

The two-stage variance is not obtained by substituting a product of degrees
into the one-stage proof. The shared first map adds the variance of the
conditional mean. An exact small verifier checks the corresponding law of
total variance and a direct dual return-kernel formula against exhaustive
enumeration.

The equal-cost candidates \((19,3)\) and \((17,5)\) expose a useful tradeoff.
The former has 3.56 more bits of zero-cap-tail margin. The latter reduces the
dense squared-character average of the second map by 13.6583 bits. The next
proof step is a compressed target-size evaluation of the composed covariance.

Supporting material is:

- pure_expander_accumulate/TWO_STAGE_SPARSE_EA_ROUTE.md;
- pure_expander_accumulate/analyze_two_stage_cap_reuse.py;
- pure_expander_accumulate/two_stage_cap_reuse_B512_diagnostic.json;
- pure_expander_accumulate/verify_two_stage_dual_covariance_small.py; and
- pure_expander_accumulate/two_stage_dual_covariance_K2_B3_r1_1_exact.json.

The first square-mixer covariance diagnostics do not close. At shell 80,
the degree-5 square map has normalized positive pair excess
\(9.83582\cdot10^{-4}\) and variance/mean upper bound
\(1.11852\cdot10^{92}\). Degree 7 improves both values by a factor of 85.25,
to \(1.15379\cdot10^{-5}\) and \(1.31208\cdot10^{90}\), but remains far
from the factor-512 target. This rejects fixed-degree escalation under the
current positive-only relaxation; it does not reject depth two itself.

Invertible conditioning isolates the obstruction. For every fixed invertible
square \(E_1\), the dense uniform reference is invariant under \(AE_1A\).
However, conditioning the sparse degree-7 ensemble on invertibility couples
its rows, and the manageable first-stage comparison uses a biased
dominant-character reference. A new conditional contraction or signed-defect
lemma is still required. In a noncertified rank experiment, 39 of 256
degree-7 square samples were full rank.

Exact common-subexpression synthesis materially changes the cost comparison:

| Constituent | Raw XORs including accumulators | Optimized forward XORs including accumulators | Optimized transposed XORs including accumulators |
|---|---:|---:|---:|
| one-stage degree 33 | 16,895 | 9,477 | 9,733 |
| depth-two degrees \((17,7)\) | 12,286 | 9,300 | 9,556 |

Thus the optimized depth-two sample saves only 177 XORs, or 1.87%, relative
to the optimized sample from the already-certified one-stage ensemble. The
one-stage circuit has a generic-schedule live-intermediate upper bound of
2,025, so this is not yet a wall-clock performance result.

Additional supporting material is:

- pure_expander_accumulate/sample_optimize_sparse_matrix.py;
- pure_expander_accumulate/sampled_one_stage_r33_matrix_circuit.json;
- pure_expander_accumulate/sampled_first_stage_r17_matrix_circuit.json;
- pure_expander_accumulate/sampled_square_r7_matrix_circuit.json;
- pure_expander_accumulate/square_sparse_r5_w80_variance_probe.json; and
- pure_expander_accumulate/square_sparse_r7_w80_variance_probe.json.

The recommended next action is to generate a register-pressure-aware
implementation of the exact degree-33 circuit and benchmark it. Resume the
depth-two proof only if that implementation shows a material performance
problem.

## One-stage sparse-EA versus optimized BCH benchmark

The degree-33 circuit and the existing ExtendedBch256x128-Eq3 circuit were
implemented in one self-contained AVX2 benchmark. The sparse path computes
\(E^{\mathsf T}A^{\mathsf T}\) for a sampled \([512,256]\) map. The BCH
path computes two \([256,128]\) transposes. Both paths therefore map
\(2^{21}\) input blocks to \(2^{20}\) output blocks.

The generator verifies both linear circuits exactly. The executable also
checks each packed kernel against an independent dense reference before
timing. Four 31-sample benchmark processes ran sequentially on Peach CPU 0
in sparse, BCH, BCH, sparse order.

| Path | First median | Second median | Mean of medians |
|---|---:|---:|---:|
| sparse-EA | 8.658 ms | 8.747 ms | 8.703 ms |
| optimized BCH | 4.891 ms | 4.901 ms | 4.896 ms |

The sparse outer is 1.777 times slower and adds 3.807 ms to the isolated
outer transform. Its static transposed-XOR ratio is 1.659. The remaining
runtime gap is consistent with the sparse schedule's larger live set: 1,305
packed values, versus 448 for BCH. The compiled sparse function is 164,189
bytes with a 53,312-byte stack frame. The BCH function is 51,914 bytes with
a 19,328-byte stack frame.

This benchmark rejects the present one-stage implementation as a direct
performance replacement for BCH. It does not change the sparse-EA distance
certificate. It also does not provide an end-to-end encoder measurement.

Supporting material is:

- ONE_STAGE_SPARSE_EA_VS_BCH_PERFORMANCE.md;
- pure_expander_accumulate/build_one_stage_vs_bch_benchmark.py;
- pure_expander_accumulate/one_stage_sparse_ea_vs_bch_bench.cpp;
- pure_expander_accumulate/one_stage_sparse_ea_vs_bch_build.json; and
- pure_expander_accumulate/one_stage_sparse_ea_vs_bch_peach.json.

The next proof-compatible performance experiment would search jointly for
XOR count and transposed live pressure. The alternative is to return to the
depth-two proof, whose optimized circuit saved only 177 XORs in the first
sample.

## Sparse-EA freeze and switch to RM(4,9)

The sparse-EA lane was frozen without moving any artifacts, because its
certificate manifests use relative paths. `SPARSE_EA_FREEZE.md` records the
closed theorem, the unresolved depth-two contraction problem, and the
matched Peach performance result. `SPARSE_EA_FREEZE_MANIFEST.json` binds 23
canonical artifacts. `audit_sparse_ea_freeze.py` verifies those hashes and
then runs the complete semantic certificate audit; both checks pass.

The replacement outer candidate is one fixed RM(4,9) \([512,256,32]\)
constituent repeated in all 4096 outer rows. Its exact spectrum is
`scripts/rm512_256_spectrum.csv`, whose SHA-256 is
`995aab561da18f22074b5c6f5413882f492084510aa1cd19b848355f1fcd4ed7`.
The Gleason-form verifier passes, with total mass \(2^{256}\), minimum
nonzero weight 32, and complement symmetry. The outer is deterministic, so
the intended theorem has no outer setup-failure probability.

The first exact gate keeps the closed proof model's uniform row-coordinate
and transposed-region permutations and RandomStepConv-M22 inner. An outward
sum over every nonzero RM shell proves

\[
  Q_1 < 2^{-26.5608944032}.
\]

The dominant shell is weight 32, whose exact multiplicity is 52,955,952.
Thus this particular occupation-one bound does not meet a 40-bit target at
10.9%. This does not prove that the construction fails; it identifies the
first-moment transfer as insufficient. Binary64 witness searches at memories
23, 27, and 40 give 29.3433, 33.9054, and 34.5109 bits. Increasing memory
alone therefore appears unable to recover 40 bits under the same bound.

Supporting material is:

- `rm49_outer/README.md`;
- `rm49_outer/evaluate_rm49_q1_gate.py`;
- `rm49_outer/rm49_q1_randomstepconv_M22_d109_diagnostic.json`;
- `rm49_outer/certify_rm49_q1_outward.py`;
- `rm49_outer/rm49_q1_randomstepconv_M22_d109_outward.json`; and
- the M23, M27, and M40 diagnostic receipts in `rm49_outer/`.

The next decision is whether to retain 10.9% with about 26 bits from the
present Q1 bound, lower the distance until Q1 reaches 40 bits, or seek a
stronger treatment of the weight-32 shell. Only after that gate is selected
should the exact Q2 and high-occupation transfers be rebuilt.

## RM(5,11) deterministic outer

The next half-rate Reed--Muller constituent is
\(\operatorname{RM}(5,11)=[2048,1024,64]\), repeated 1024 times.  The outer
is fixed.  The probability space retains the independent uniform row and
region permutations and RandomStepConv-M22 maps.

Exact formulas give every nonzero constituent shell below weight 128.  An
outward replay over weights 64, 96, 112, 120, and 124 proves that their
complete occupation-one contribution is at most

\[
  2^{-83.3743314000}.
\]

The dominant shell is weight 64.  Thus the minimum-distance shell clears the
40-bit gate by more than 43 bits, unlike RM(4,9).

This is not a full Q1 certificate.  The unknown RM(5,11) enumerator remains
relevant from weight 128 through 380.  A total-mass bound is sufficient from
weight 384 onward.  Two generic replacements were tested and rejected:

- nonnegative Type-II Gleason enumerators can place nearly all relaxed mass
  at weight 128; and
- distance-64 Johnson packing exceeds the needed per-shell caps by hundreds
  of bits below weight 384.

Published cumulative RM bounds explain the correct asymptotic behavior, but
the versions checked here do not expose a usable finite constant at eleven
variables.  Sampling estimates are not accepted as upper bounds for this
fixed code.  Exact Kasami--Tokura--Azumi formulas below weight 160 would help
but would not cover weights 160 through 380.

Supporting material is:

- `rm511_outer/README.md`;
- `rm511_outer/RM511_Q1_ENVELOPE_STATUS.md`;
- `rm511_outer/verify_rm511_low_weight_formulas.py` and
  `rm511_outer/rm511_low_weight_exact.json`;
- `rm511_outer/evaluate_rm511_q1_low_gate.py` and its diagnostic receipt;
- `rm511_outer/certify_rm511_q1_low_outward.py` and
  `rm511_outer/rm511_q1_low_randomstepconv_M22_d109_outward.json`; and
- `rm511_outer/evaluate_rm511_packing_gap.py` and
  `rm511_outer/rm511_q1_packing_gap_diagnostic.json`.

The next useful task is a finite RM-specific cumulative envelope for weights
128 through 380.  Q2 and implementation work remain paused until complete Q1
closes.

## Direct combined-spectrum EA route

The EA lane was reopened with a different proof object.  Instead of requiring
simultaneous upper caps on every realized shell, the new route averages the
actual nonnegative SPIN transfer functional over one sampled constituent.
The constituent remains one fixed right-regular sparse map followed by an
accumulator, rank-tested once and repeated in every outer row.

At ((K,B,L)=(256,512,4096)), the full occupation-one diagnostic first clears
40 bits at degree 23, with 44.826 bits.  Degree 21 gives 39.174 bits.  This is
a substantial improvement over the degree-33 shell-cap requirement.

The equal-message sector of occupation two is the real obstruction.  Its
complete degree-23 Chernoff sum has only 19.174 bits; degree 29 has 33.088
bits; degree 33 has 42.707 bits.  Outer image weight 16 dominates.  Thus the
old low-shell caps were conservative for Q1 but were protecting against a
genuine repeated-row effect.

Increasing the rate-half constituent size improves this sector without
changing the leading raw XOR count (Nr).  Truncated positive sums identify
candidate degrees 31, 29, 27, and 25 at block sizes 1024, 2048, 4096, and
8192.  These candidates are not certificates: the reported larger-block
sums omit higher outer weights and all distinct-message pairs.

An exact rational small model now verifies both moments of

\[
F_f(C)=\sum_{x\ne0}f(\operatorname{wt}(Cx))
\]

against exhaustive enumeration of every sparse map at ((K,B,r)=(4,5,3)).
This validates the combined-functional pair kernel, not its large-parameter
bound.

Supporting material is:

- `ea_combined_spectrum/README.md`;
- `ea_combined_spectrum/EA_COMBINED_SPECTRUM_STATUS.md`;
- `ea_combined_spectrum/evaluate_ea_combined_q1.py` and both Q1 receipts;
- `ea_combined_spectrum/evaluate_ea_combined_q2_rank_one.py` and its four
  retained rank-one receipts; and
- `ea_combined_spectrum/verify_combined_functional_small.py` with
  `combined_functional_K4_B5_r3_exact.json`.

The recommended next proof target is the ([2048,1024]), degree-29 candidate:
outward-round and optimize the rank-one Q2 witness, then bound the
distinct-message rank-two functional.  Benchmarking should wait until those
two gates close.

This route is now paused by user decision.  The Q1 gain did not survive the
equal-message Q2 reuse test strongly enough, and the observed improvement of
about two degree points per block-size doubling was not compelling.

## Exact-spectrum replay at smaller message lengths

The authenticated extended-BCH and Reed--Muller spectra were replayed at
smaller \(k\) and smaller RandomStepConv memory. The outer remains one fixed
deterministic constituent repeated in every row. Only the routing and inner
maps are sampled.

Neither extended BCH \([32,16,8]\) nor extended BCH \([128,64,22]\) reaches
40 bits at occupation one in the tested range. RM(4,9) does. The selected
all-proof-regime candidate is \(k=2^{13}\), \(N=2^{14}\), 32 outer rows, and
memory 22. Its exact-spectrum binary64 diagnostics give 40.71735 bits for
occupation one and 84.27710 bits for occupation two. The random rate-half
first-moment reference gives 59.9674 bits at this length; it falls below 40
bits at every tested \(k<2^{13}\) except that no claim of universal
impossibility is inferred from this random benchmark.

The result is not a distance certificate. Attempts to cover occupations 3
through 32 by minimum-distance deletion, banded Bernoulli domination, and an
earliest-pivot reduction were vacuous. The first loses the exact spectrum;
the third reduces the input to one pulse and exposes a \(2^{-M}\)-scale
inner-annihilation floor. The first banded calculation also charged an
all-zero input atom with the Chernoff distance factor. Splitting that atom
correctly did not close the bound: live low-band compositions remain vacuous,
with a positive 734.161-bit bound already at occupation three.

The random comparison was corrected from a global random code to one uniform
full-rank constituent sampled once and reused in every row. At block size 512,
\(k=2^{13}\), and memory 12, Q1 has 91.019 bits and exact reuse-aware Q2 has
180.187 bits. For Q2, equal local messages share the same random image and are
handled separately from distinct messages. Occupations at least three remain
open because the rank and binary relation type of each active local-message
tuple must be retained. This 512-bit result is only a matched control for
RM(4,9), not evidence for the smaller BCH blocks. At the same \(k\) and memory,
random \([32,16]\) and \([128,64]\) have Q1 margins -10.529 and 6.180 bits;
the corresponding BCH margins are -10.134 and 7.181 bits. Thus neither small
random block closes, and the exact BCH constituents are slightly better.

Supporting material is in `small_k_replay/`, led by
`SMALL_K_EXACT_SPECTRUM_STATUS.md`. The next decision is whether block size
512 is acceptable. If not, the useful search is for a stronger authenticated
spectrum near block size 128; continuing the random-512 all-occupation proof
would not address the desired size regime.

### Controlled 10% margin curves

A matched curve experiment now fixes
\(M(k)=\lceil\log_2 k\rceil+2\) for every constituent and uses
\(D=\lceil0.10N\rceil\). The curve quantity is the diagnostic Q1 margin, not a
full certificate margin. Each random control has the same block length and
dimension as its structured counterpart and samples one full-rank map reused
across rows.

The best Q1 margins over \(2^8\le k\le2^{20}\) are 4.107 bits for BCH
\([32,16,8]\), 30.014 bits for BCH \([128,64,22]\), and 38.925 bits for
RM(4,9). Their matched random controls reach 0.866, 24.885, and 129.461 bits.
Thus no authenticated structured constituent reaches 40 Q1 bits under this
neutral schedule. RM is closest; the large gap to random at block size 512 is
a same-size spectrum comparison, while both smaller BCH constituents beat
their own random controls.

Supporting material is `small_k_replay/MATCHED_CONSTITUENT_CURVES.md`,
`evaluate_matched_constituent_curves.py`, and the JSON/CSV receipts named
`matched_constituent_k_margin_d100`.

### Complete rate-half family Q1 sweep through block length 1024

The neutral (M(k)=\lceil\log_2 k\rceil+2), 10% Q1 sweep now covers every
accepted nontrivial rate-half BCH-derived and RM constituent below the
requested block-length cap. The BCH-derived rungs are lengths 8, 32, 64, and
128. The length-64 rung is the published Philips shortened-XBCH code; its
complete (2^{32})-word enumerator was reconstructed locally and matches the
frozen spectrum hash. The RM rungs are RM(1,3), RM(2,5), RM(3,7), and RM(4,9)
at lengths 8, 32, 128, and 512. Random full-rank controls cover every doubling
from length 8 through 1024, with one sampled constituent reused in all rows.

No structured curve reaches 40 diagnostic Q1 bits. The maxima are -1.282,
4.107, 11.916, and 30.014 bits on the BCH-derived ladder, and -1.282, 4.107,
13.587, and 38.925 bits on the RM ladder. The random ladder first clears 40
bits at block length 256, where its maximum is 58.745 bits. Lengths 512 and
1024 reach 129.461 and 274.879 bits.

The BCH-derived sweep stops at length 128 because the required complete
([256,128]) or larger rate-half spectrum is unavailable. The RM sweep stops
at length 512 because rate-half RM codes occur at the odd-(m) middle orders;
there is no rate-half member at length 1024, and the next is length 2048.

These are binary64 occupation-one diagnostics. They do not prove a 40-bit
distance certificate. For random constituents, the exact expected spectrum
supports the joint Q1 first moment but does not certify one sampled spectrum
simultaneously and does not solve reused-map relation types at higher
occupation.

Supporting material is:

- `small_k_replay/RATE_HALF_FAMILY_CURVES.md`;
- `small_k_replay/SPECTRUM_SOURCES.md`;
- `small_k_replay/evaluate_rate_half_family_curves.py`;
- `small_k_replay/rate_half_family_k_margin_d100.json` and `.csv`; and
- the authenticated spectra in `small_k_replay/spectra/`.

The recommended next step is to select the acceptable block-size regime. If
256 bits is acceptable, it is the smallest random control that clears the Q1
screen; the proof task is then a structured spectrum or combined-functional
event near ([256,128]), followed immediately by a reuse-aware Q2 bound.

### RM2Sub calibration and persistence-matched replay

The random-inner family scan was replayed with RM2Sub after separating state
dimension from epoch length. At fixed persistence exponent
\(s+\log_2t=26\), epochs of 64 bits outperform the tested 128- and 256-bit
epochs on exact BCH \([32,16,8]\) and \([128,64,22]\) spectra. The retained
epoch size is \(t=64\).

For \(k=2^e\), the neutral replay uses \(s(e)=\max(7,e-4)\). For \(e\ge11\),
this matches the earlier RandomStepConv schedule \(M(e)=e+2\) at the level of
the first-order live-state reset exponent. A direct \(k=2^{13}\) comparison
also checks occupation two for both BCH constituents. RM2Sub stays within a
few bits of its persistence-matched random-step control; the comparison is
not an equality of transition laws.

At \(k=2^{20}\), \(t=64,s=16\), and 10% distance, the exact BCH
\([128,64,22]\) Q1 diagnostic is 32.338 bits. RM(4,9) gives 38.941 bits. Thus
neither produces a 40-bit \(k=2^{20}\) certificate. Across the full scan,
RM(4,9) exceeds the 40-bit Q1 screen for \(2^{14}\le k\le2^{18}\), peaking at
41.462 bits for \(k=2^{16}\). This is a useful inner-model result, not a full
distance certificate: arithmetic is binary64 and occupations at least three
remain open.

The complete record is
`small_k_replay/RM2SUB_CALIBRATION_AND_FAMILY_REPLAY.md`. The principal
machine receipts are `rate_half_family_k_margin_d100_rm2sub_t64_matched.json`
and `.csv`; `rm2sub_calibration_replay_audit.json` records a passing exact
audit of the selected A/B maps and receipt structure.

### 10% primary parameter tranche with larger RM2Sub epochs

The optimization target was relaxed to 10% distance. Distances near 11% are
now bonuses. The study must ultimately return the fastest fully certified
candidate and the fastest conditional candidate under a narrow, stated outer
spectrum assumption.

The first tranche generated and exactly enumerated the missing RM2Sub maps for
`t=128` and `t=256`. It evaluated exact BCH `[128,64,22]`, exact RM(4,9), and
matched random Q1 controls at native lengths. Only RM(4,9) at `k=2^16` crosses
the 40-bit Q1 screen. BCH128 remains below 38 bits.

The corrected RM(4,9) Q2 screen retains `t128_s13` with 85.947 bits and all
tested `t256` states. The earlier -11.587-bit result for `t128_s13` was an
invalid elimination: two aggregate margins from disjoint tilt grids were
compared after aggregation. The valid calculation minimizes pointwise over
one complete grid before summing the spectrum.

A nested `t=128` family removes the resampled-map confound. Its `s=13` and
`s=14` members give 85.819 and 87.181 Q2 bits, so the apparent sharp state
threshold disappears. At matched persistence exponent 20, `t64_s14`,
`t128_s13`, and `t256_s12` give Q1/Q2 pairs 42.583/87.424,
42.147/85.947, and 41.406/83.721. Larger epochs reduce the screened margins
in this controlled comparison. No all-occupation or outward claim has been
made.

Supporting material is:

- `small_k_replay/TEN_PERCENT_PARAMETER_STUDY.md`;
- `small_k_replay/evaluate_rm2sub_primary_tranche.py`;
- `small_k_replay/summarize_rm2sub_primary_tranche.py`;
- `small_k_replay/rm2sub_primary_tranche_d100.json` and `.csv`;
- `small_k_replay/rm2sub_primary_frontier_summary_d100.json` and `.csv`;
- `small_k_replay/rm2sub_primary_tranche_audit.json`;
- the regenerable `small_k_replay/nested_rm2sub/t128_chain0/` calibration
  bundle, preserved in commit `97c7a7e` and omitted from the Overleaf tip;
- `small_k_replay/RM2SUB_PERSISTENCE_UNCONFOUNDED.md`; and
- the added exact maps in
  `small_k_replay/rm2sub_calibration_constituents/`.

The next proof step is a representative all-occupation screen for `t64_s14`,
`t128_s13`, and `t256_s12`, followed by a complete binary64 replay only for
survivors.

### Fixed-RM RM2Sub occupation ladder at $k=2^{16}$

The representative occupation screen is now complete through $Q=8$ for the
three persistence-20 configurations. One fixed RM$(4,9)$ constituent is
repeated in all 256 rows. The proof averages only over the row-coordinate
permutations, region permutations, and nonzero epoch scalars.

A new positive two-colour recurrence tracks regular Bernoulli-reference rows
and forced all-one rows without subtraction. A direct six-position enumerator
checks all 15 type pairs with maximum absolute discrepancy
$1.6654\times10^{-15}$. An exact-spectrum Holder change of measure preserves
the fixed-code model; it does not replace the repeated RM constituent by
independent random codes.

The $Q=3,\ldots,8$ margins for `t64_s14` are 82.414, 110.168, 135.423,
166.847, 174.884, and 217.458 bits. The corresponding `t128_s13` and
`t256_s12` margins are smaller at every occupation. Extending `t64_s14`
gives a positive 30.683-bit margin at $Q=29$ and a vacuous bound at $Q=30$.
A fine local tilt grid confirms the last transition. The
dominant face has no all-one row.

This is an exploratory baseline, not a certificate. Arithmetic is nearest
binary64, and occupations 30 through 256 are not yet covered by nonvacuous
bounds. The exploration does not impose a 40-bit acceptance threshold.
The controlled account is `small_k_replay/RM2SUB_Q_LADDER.md`; its scripts and
JSON receipts are in the same directory.

The first bridge attempt applied Holder once to the full RM spectrum and the
finite three-state dense transfer. It fails badly at $Q=29$ because the
change of measure loses too much of the uniform-reference event exponent.
This is a failed proof bound, not a distance counterexample.

A jointly optimized nine-band probe localizes the obstruction. At $Q=29$,
every pure face of weight at least 64 is positive; the 32--47 and 48--63 faces
have -765.880 and -141.613 bits. The minimum-weight dense-face margins are
-29.391, 4.780, 39.057, and 73.434 bits at $Q=52,53,54,55$. The
exact sparse transfer restricted to the minimum-weight shell has more than
1,200 bits at $Q=29$ and $Q=30$ and more than 2,100 bits at $Q=54$.

The first mixed-composition step is now closed diagnostically. Conditional on
the total number of surviving reference bits, permutation symmetry collapses
two Bernoulli row types to the established forced-support recurrence. A
direct enumerator checks the identity to $4.8850\times10^{-15}$.

Using reference probabilities 0.25 for weights 1--47, 0.5 for weights
48--416, and 0.75 for weights 417--512, every minimum-band-plus-partner family
has more than 1,170 bits over $30\le Q\le52$. Their union has 1172.882 bits.
This is strong evidence that the earlier gap was caused by a common-reference
proof relaxation, not the RM2Sub parameters.

The unresolved sparse cases use three or more spectrum bands. The dense union
for $53\le Q\le256$ also remains open. The exact account and supporting
receipts are in `small_k_replay/RM2SUB_TWO_BAND_BRIDGE.md`.

### Refined-band closure of the RM2Sub middle and dense occupations

The preceding open cases are now closed as binary64 diagnostics. The key
change is a better partition of the same exact RM spectrum; the construction
does not change. The final groups and Bernoulli references are

- weights 1--95 with probability 0.25;
- weights 96--416 with probability 0.5; and
- weights 417--512 with probability 0.8677722630069483.

Moving weights 48--95 into the low group does not increase that group's
pointwise density envelope. The resulting three-binomial reduction covers
every group composition for every occupation from 30 through 256. The audit
checks 2,857,249 compositions and 227 occupations. Their union has 1228.638
diagnostic bits, with occupation 30 weakest.

The former dense failure was proof slack. Under the old coarse groups, the
Q=224 margin was -3665.254 bits. The refined groups give 6998.852 bits at
Q=224 and 1611.249 bits at Q=256. All nine pure spectrum bands also close at
Q=256; the weakest pure face has 2015.626 bits.

Supporting material:

- `small_k_replay/RM2SUB_REFINED_BAND_BRIDGE.md`;
- `small_k_replay/probe_rm2sub_three_group_bridge.py`;
- `small_k_replay/consolidate_rm2sub_refined_band_bridge.py`;
- `small_k_replay/rm2sub_refined_band_bridge_q30_q256_diagnostic.json`;
- `small_k_replay/rm2sub_refined_band_bridge_audit.json`;
- `small_k_replay/rm2sub_fixed_rm_dense_pure_bands_q256_probe_d100.json`; and
- `small_k_replay/rm2sub_central_split_95_q256_probe_d100.json`.

The remaining finite-proof work is outward rounding and combination with the
existing Q=1--29 receipts. No additional band argument is currently needed.

The diagnostic combination is now complete. It selects the exact-spectrum
Q1 receipt, the fixed-spectrum Q2 receipt, the two-colour ladder for Q=3--4,
the three-group bridge for Q=5--29, and the refined bridge for Q=30--256.
The full Q=1--256 union has 42.577982 bits. It clears the 40-bit line by
2.577982 bits, with Q=1 dominant. The remaining task is directed outward
rounding; the present result is not yet a formal certificate.

The combined receipt and audit are
`small_k_replay/rm2sub_full_occupation_q1_q256_diagnostic.json` and
`small_k_replay/rm2sub_full_occupation_q1_q256_audit.json`.

### Outward closure of the fixed-RM small-k instance

The finite proof-model instance is now closed with outward arithmetic. The
construction is unchanged from the final diagnostic: one fixed RM(4,9)
([512,256,32]) constituent is reused in all 256 rows; routing uses
independent uniform row-coordinate and region permutations; and the fixed
audited RM2Sub inner has (t=64,s=14) with an independent nonzero scalar per
epoch.

At (k=2^{16}), (N=2^{17}), and bad output weight at most 13,107, the
final receipt proves

\[
  d_{\min}\ge13{,}108
\]

except with probability below (2^{-42.5779817562755}). The exact outward
margin is 42.577981756276 bits. Occupation one remains dominant; adding all
occupations two through 256 costs only about (3.32\times10^{-10}) bits.

The occupation cover is:

- Q=1: 42.577981756608 bits;
- Q=2: 74.611094272379 bits;
- Q=3--4 union: 83.102187705873 bits;
- Q=5--29 union: 81.502991282228 bits;
- Q=30: 1228.637966758792 bits; and
- Q=31--256 union: 1293.444358402287 bits.

The dense witness representation was reduced from roughly 800 MB of raw JSON
to 2,856,753 signed bytes, one exact-tenth Chernoff witness per canonical
composition. The dense outward checker uses Arb-generated coefficient
endpoints, explicit positive-operation rounding inflation, upward-rounded
matrix products, and exact power-of-two scaling. The final manifest verifies
all six component checkers and receipts and checks exact Q=1--256 coverage.

Three attempted dense compressions were rejected before the packed format was
selected. A single Bernoulli-half envelope for all regular non-all-one RM
words loses tens of thousands of bits at high occupation. A scalar row-sum
envelope destroys the state-activation information and is also vacuous. One
Chernoff witness shared by every composition works through Q=128 but becomes
vacuous in the denser range. None is a construction failure. The packed table
retains one witness per composition without retaining the bulky diagnostic
rows and reproduces or improves every archived binary64 occupation margin.

Supporting material:

- `small_k_replay/RM2SUB_RM49_FINITE_CERTIFICATE.md`;
- `small_k_replay/RM2SUB_RM49_FULL_CERTIFICATE_MANIFEST.json`;
- `small_k_replay/certify_rm2sub_rm49_full_distance.py`;
- `small_k_replay/rm2sub_rm49_t64_s14_full_distance_outward.json`;
- the six occupation-range checkers and outward receipts; and
- `small_k_replay/RM2SUB_RM49_Q31_Q256_INPUT_MANIFEST.json` plus the compact
  dense witness receipts.

This closes one finite proof-model instance, not the frozen Structured SPIN
(B=256, t=128, s=19) implementation and not an asymptotic family. The next
work is the planned parameter study: use this certified engine to compare
exact-spectrum BCH/RM constituents and RM2Sub schedules across smaller and
larger message lengths, while keeping performance metadata beside margin.
