# Three-quarter-rate GF256 SPIN goal

The target is nominal rate 3/4, 10% symbol distance, and at least 40 bits
of whole-code setup-failure margin. Start near K=2^20. Encode sixteen
GF256 lanes in each 128-bit application record, with one global regional
routing stage and no general large-field multiplications.

The target now has a [complete K20 certificate](GF8_PUNCTURED_WHOLE_CERTIFICATE.md)
and a [matching research encoder](GF8_PUNCTURED_NATIVE.md): nominal rate3/4,
10% symbol distance, and116-bit whole-code margin. The winning construction
uses outer[128,84] and retains seven of eight full-Cauchy inner outputs.
Its requested K20 rate is256/343; full transpose takes5.57--5.82ms in three
serial101-call launches. Release, AddressSanitizer, and UndefinedBehaviorSanitizer
checks pass. This closes the goal's proof and validated-implementation gates;
it does not establish an optimal rate or replace the binary library.

The original, unclosed attempt below used an unpunctured [128,96] outer.
For that outer and the current natural row multiple 64, a request
for K=1048576 uses L=10944 rows, padded dimension 1050624, and N=1400832.
The requested rate is 128/171, just below 3/4 because of shortening. The
distance target excludes output weights through 140083. The full padded
code has nominal rate 3/4; a certificate for it also covers the shortened
request.

The unpunctured attempts below are historical evidence, not the new certificate.
The prior [65-bit cubic certificate](GF8_PACKET_CERTIFICATE.md), the older
higher-rate certificates, and all production profiles remain unchanged.

## Work and acceptance gates

1. Derive valid packet bounds before using a new map in the existing
   regional/source calculation. Preserve the conditioning order and
   distinguish actual support from the rank of source forms.
2. Screen sparse and intermediate classes, not only a single active row.
   Keep failed bounds as evidence about the analysis, not counterexamples.
3. For a surviving candidate, cover every class with outward arithmetic,
   check the joins, and obtain an independent mathematical audit.
4. Implement a faithful research encoder, check scalar references and
   sanitizers, then measure whole calls serially. Count routing, arithmetic,
   conversions, and stores. Do not substitute instruction counts for time.

The first parallel branches are stronger feedthrough, a larger state per
packet, and inexpensive circuits for any surviving map. Rate 3/4 is the
next target, not an asserted upper limit.

## Starting evidence

The previous [higher-rate investigation](GF8_HIGH_RATE_INNER.md) found
that enlarging the outer, modest rank-history refinement, and pairwise
premixing did not close 3/4. A six-symbol state reduced the proof gap but
did not eliminate it. Removing return events as an explicitly unrealistic
diagnostic greatly improves sparse scores. This motivates strengthening
the joint activation/return behavior.

The [paired cubic implementation](GF8_CUBIC_NATIVE.md) measures about
.237 / .94 / 5.71--5.86ms at K16/K18/K20. This is an implementation control
at a different rate, not a timing claim for a three-quarter-rate code.

## Dense-class control

The existing directed dense theorem, with outer geometry [128,96] and
K20 padding, covers inverse-input weights 392233 through 1400832. These
are densities at least 7/25. The replay gives 17894 bits, 144 cells, and
192-bit Arb precision. Its witness digest is
`2cc3870f9307ad8ec18a0d30107843abcc286f6f4f14baa625a910af0297d8cb`.

Reproduce from the repository root:

```text
python research/workstreams/char2_field_spin/theory/gf8_outer_rate_dense.py --dimension 96 --k-log 20 --first 7/25 --cell-bits 128
```

This is a complete dense component, not a whole-code certificate. The
existing helper labels its original width16/state8 construction. Do not
copy that label onto a new code. A new certificate must explicitly verify
the dense theorem's inner hypothesis and bind its own geometry and maps.
The remaining sparse/middle cover must join no later than weight 392233.

## Progress log

- Goal initialized. Three construction/proof branches started; the dense
  component was independently replayed. No new encoder, benchmark, commit,
  or production change has been made in this goal round.

- Stronger feedthrough and non-self-orthogonal state maps were derived and
  checked. The [full-state branch](GF8_FULLSTATE_PROOF.md) tests state sizes
  8 and 16, exact output-zero branch handling, retained source-form rank,
  and an outer [256,192] control. None closes the intermediate sparse gap.
- Two local proof refinements passed independent audits and exhaustive
  small-field tests: [joint live-kernel rank](LIVE_KERNEL_RANK.md) and
  [partial-live return moments](PARTIAL_LIVE_RANK.md). They preserve the
  conditioning order and improve valid bounds, but their observed gains
  do not close the three-quarter-rate target.
- An extra positive mark enforcing total input weight at least 33A did
  not help the sampled sparse classes: its optimum was the neutral mark.
  Direct fixed-word first moments were substantially worse. These controls
  rule out two easy replacements for the current rank-potential method.
- The [full-Cauchy packet](GF8_MDS_PACKET.md) has MDS birth and return maps.
  Explicit 8/4, 8/8, and 8/16 maps pass minor, scalar, adjoint, and batched
  chronology checks. An exact factorization removes a separate feedthrough
  map in the 8/8 and 8/16 schedules. Proposed arithmetic is 31 and 58
  GFNI512 instructions per sixteen records, respectively. These are
  checked circuit counts, not measured whole-encoder times.

### Geometry gates from the first goal round

All entries below are floating individual-class scores at outer [128,96],
K20, and 10% distance. A negative score means that this upper bound fails
to prove the class; it is not a counterexample to the code. Passing one
class does not establish the requested whole-code 40-bit margin.

| Full-Cauchy packet | 64 active rows | 128 active rows | 256 active rows |
|---|---:|---:|---:|
| width 8, state 8, refined returns | -875.98 | -2312.06 | not run |
| width 8, state 16, refined returns | +355.66 | -1041.33 | not run |
| width 16, state 8, conservative returns | -873.60 | -2308.48 | not run |
| width 16, state 16, conservative returns | +349.40 | -1065.00 | not run |
| width 8, state 32, conservative returns | not run | +1315.19 | -1131.11 |

All listed sparse optimizations converged. The generic width-16 or
state-32 screens use the full-Cauchy algebraic family; they do not claim
a packed implementation for those geometries. The independent geometry
helper agrees exactly with the existing width-eight conservative transfer.

At state 32, middle densities 0.03, 0.05, and 0.07 still give scores of
-14880.97, -17188.51, and -14180.13 bits. At the 0.05 witness, the marked
return entry is only about 5.2e-25 and the continuation entry is 0.999958.
The source contribution is about -101295 bits, while the output-weight
mark costs +118462 bits. More return suppression alone barely changes
that witness. See the [feedthrough analysis](GF8_FEEDTHROUGH_PROOF.md).

For contrast, an explicitly unrealistic diagnostic deleting every return
from the width-eight/state-eight sparse transfer gives +18061, +36383,
and +73308 bits at A=128,256,512. This is not a valid code bound. Together,
the controls show two distinct gaps: returns matter in the sparse cover,
while the current middle cover loses even when they are negligible.
The diagnostic sets the live-to-zero entry to zero for every nonempty
input weight, leaving all other entries unchanged. Its optimized natural
log marks `(logtilt,logtau)` are `(6.5691544191,5.5412635452)`,
`(6.5692355433,5.5412635452)`, and `(6.5694969332,5.5412635452)`.

Reproduce the independent geometry gates:

```text
python research/workstreams/char2_field_spin/theory/gf8_mds_geometry_screen.py --width 16 --state 8
python research/workstreams/char2_field_spin/theory/gf8_mds_geometry_screen.py --width 16 --state 16
python research/workstreams/char2_field_spin/theory/gf8_mds_geometry_screen.py --width 8 --state 32 --sparse 128 256 --density .03 .05 .07 .10
```

### Decision after the first round

Keep the factored full-Cauchy 8/8 packet as the most attractive local
implementation candidate, not as a proved three-quarter-rate code.
Prioritize a joint source-support/active-excursion bound that can bridge
the sparse and middle regions. The current separation forgets useful
relationships between which outer rows are active and how state activity
continues. Test any proposed refinement against the failed A=128/256 and
density-0.05 controls before running a broad search or writing a kernel.
If that coupling cannot be recovered, change the construction to make it
available. Simply enlarging the state or doubling the packet width has
not justified its additional work in these screens.

At the end of that first round, the goal remained active. No whole-code three-quarter-rate certificate,
new native benchmark, commit, or production change resulted from this
round. All prior complete certificates and native controls are preserved.
The final serial validation passed all 44 tests across the twelve new
lemma, geometry, map, feedthrough, full-state, and diagnostic test suites.

## Second round: local puncturing closes the target

The [excursion diagnostic](GF8_EXCURSION_DIAGNOSTIC.md) showed many live
segments in the failed unpunctured witnesses. Joint active-row/input-weight
conditioning recovered only small amounts of slack. We preserved those
results and changed the rate allocation instead: outer[128,84] supplies
minimum row weight45, and a fixed seven-of-eight local puncture raises the
complete construction's nominal rate to3/4.

The directed components cover A=1..713, W=32113..722534, and
W=722535..1605632, with116/2782/2539-bit component margins. The exact union
retains116bits. Independent audits verify the selected Cauchy minors,
coverage joins, matrix digest, and the full256^200704 omitted-output
completion factor. Both the independent agent and root replayed the whole
certificate. See the [proof statement](GF8_PUNCTURED_WHOLE_CERTIFICATE.md).

The native implementation reuses the factored full-Cauchy8/8 maps and
processes four packets at a time. Dense scalar references, all-byte local
tests, shortening, schedules, aliasing, and large-size checks pass.
Inlining improves the kernel; scratch-stride padding loses at K20 and is
retained only as a control. Three final full-call K20 launches give
5.816/5.571/5.815ms against5.880/5.907/5.824ms for the lower-rate cubic
control. The new retained output is13.38% shorter. This is a rate gain at
comparable field-encoder cost, not a speed win over binary SPIN.

No commit, push, production change, or replacement of earlier proofs is
part of this checkpoint. Recommended next turn: preserve this winner,
then test outer dimensions above84 with the same seven-output packet,
using sparse and middle proof gates before any kernel change.

Final validation passes28 proof/diagnostic tests, the independent matrix
and packing audit, the root's complete directed replay, and release plus
AddressSanitizer/UndefinedBehaviorSanitizer checks. The goal is complete.
The first concrete follow-up candidate is outer[128,85], nominal rate85/112;
it has no certificate yet and must remain separate from this winner.
