# Repaired cubic inner: directed middle and dense components

The K20 GRS[128,83] candidate now has complete middle and dense component
certificates for 10% symbol distance. The middle component has a
3303-bit failure margin; the dense component has 212483 bits. These
components alone are not a whole-code claim. They are now joined with
matching sparse coverage in the [65-bit whole-code certificate](GF8_PACKET_CERTIFICATE.md).

## Fixed candidate

The scalar field is AES GF256. Each padded source row is a GRS[128,83]
codeword. Its coordinates enter 128 separate regions, each independently
uniformly permuted. Every input coordinate and every state coordinate
uses an independently sampled uniform nonzero field multiplier. Setup
is sampled once and fixed for subsequent encodings.

The inner has packet width eight and state dimension four. Its maps are
the eight phases of `cubic_matrices` in
[the exact map reference](theory/gf8_packet_maps.py). Each phase uses four
residue nodes `(16,17,172,173)` and eight ordered input labels from the
embedded-GF4 geometry. For input label u and residue node v, set

\[
A_{u,v}=(u+v)^3+1,\qquad C_{v,u}=(u+v)^3+3.
\]

Field addition of byte labels is XOR. The scalar eta is byte label two.
The map audit checks superregularity of both A and C and the identity
`CA=0`. The inner computes `Y=X+2AR` and `Q'=R+CX`, where R is the
diagonally refreshed incoming state and X is the diagonally scaled input.
It starts with zero state and discards the final state.

The audit exhausts all support masks of the returned map `I+2AC`. Its
exact minimum output weights for exact input weights one through eight
are `[8,6,5,4,3,2,2,1]`, uniformly across all phases. The component replay
uses these minima and the generic joint return/rank cover justified in
[the width investigation](GF8_PACKET_WIDTH_SCREEN.md).

Natural padding gives

\[
K_{\rm req}=1048576,\quad L=12672,\quad K_{\rm pad}=1051776,
\quad N=1622016.
\]

The requested rate is `64/99`, approximately 64.646%, and the nominal
rate is `83/128`. The bad output threshold is `162201`; excluding that
threshold proves minimum distance at least `162202`, hence at least 10%.

## Complete integer coverage

| Component | Input-weight interval | Checked cells | Certified margin |
|---|---:|---:|---:|
| Middle | 48661–389283 | 210 | 3303 bits |
| Dense | 389284–1622016 | 152 | 212483 bits |

The middle interval is `[.03N,.24N)`. The replay aggregates the local
two-state packet matrices as exact fractions, then converts their entries
to outward-rounded decimal intervals with 80-digit precision. The
source envelope deliberately remains the older relaxed one, not the
stronger floating-screen truncation. An exact modal-conditioning factor
is used for the regions. Frozen rational marks are reused without an
optimizer; deterministic bisection verifies every integer weight.

Its upward-rounded failure bound is

`3.5293148762285843312569814012058950504033303692153197007186574154557181842827022E-995`.

The middle partition hash is
`a6bb9147e34372abfd664a6c75d5ea2bcbdec8c6fb3dbcf7bd99d37828216c9f`.

The dense theorem does not depend on packet width or state size. For a
fixed output word, causal identity feedthrough determines a unique
inverse-inner word before the independent input diagonals are exposed.
The existing source-membership and output-ball arguments therefore apply
to this changed inner. The dense wrapper reuses only that directed
arithmetic, checks the complete finite outer geometry, and constructs
fresh width-eight/state-four metadata. It does not relabel an old
inner-specific certificate.

The dense bound is at most `2^-212483`. Its witness hash is
`efcdb9e1ded74c57ee496378e5665d96820755cf5b3aa2079aab01e5a2628b5b`.

Every nonzero outer row has at least 46 nonzero coordinates. Sparse
coverage through 1057 active rows therefore bridges into the middle
interval, since `46*1058=48668>48661`. This statement identifies the
coverage requirement. The separate whole-code replay verifies that sparse
class and the probability sum.

## Reproduction and checks

```text
python -B research/workstreams/char2_field_spin/theory/gf8_packet_middle_certificate.py
python -B research/workstreams/char2_field_spin/theory/gf8_packet_middle_certificate.py --dense
python -B research/workstreams/char2_field_spin/theory/gf8_packet_middle_certificate_tests.py
```

All six tests pass. They verify fresh geometry and map metadata, exact
packet aggregation, agreement with the floating bound after removing its
source refinement, interval domination of individual weights, complete
middle/dense replays, and rejection of invalid or inexact witnesses.

The sparse replay and exact aggregate are recorded in the whole-code
certificate. Native performance measurement remains the next gate;
these certificates do not imply a measured speedup.
