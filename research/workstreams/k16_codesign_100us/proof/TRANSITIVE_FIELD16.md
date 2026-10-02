# Exact transfer to field16 randomizers

The retained small RS16 code has a fresh complete bound at K = 65,536 and
N = 131,072. The replay in `fresh-rs16-small-p256.json` gives

\[
\Pr[d_{\min}\le 13107]
\le 2^{-49.7217337682000108916446781376602835635}.
\]

The probability is over independently sampled setup maps and shuffles.
Outside this event, the minimum distance is at least 13,108, above 10% of N.
The replay covers every nonzero-message occupancy from 1 through 512.
It uses 256-bit outward Arb arithmetic and freshly regenerated fixed-map data.

The construction consists of 512 independent 128-to-256-bit outer groups.
Each group has four parallel GF(16) RS[16,8] rows. Each aligned 16-bit symbol
receives an independent uniform GL(16,2) map. Each group then independently
shuffles its 64 four-bit packets. Every region independently shuffles the
512 group slots. The retained physical t64/s16 inner starts in state zero.
For input x and entering state a, step i emits x + A a and updates the state
to M_i a + A^T x. Here A is the selected fixed expansion, and each M_i is an
independent uniform GL(16,2) map. State continues across regions without a
final flush.

The following transfer permits cheaper randomizers in the outer, the inner,
or both, without changing this complete bound.

Fix a binary basis of F = GF(2^16). For c in F, let L_c be the binary matrix
of multiplication by c. Independently sample c from F minus zero for every
replaced outer symbol and every replaced physical inner step. At each such
location, use either U L_c V or U L_c^T V, where U and V are arbitrary fixed
invertible binary matrices for that location. All these choices define linear maps on
the 16 code coordinates; they do not multiply the carried 128-bit payloads.

For any fixed nonzero x, the map c ↦ L_c x is a bijection of F minus zero.
The same statement holds for c ↦ L_c^T x. Indeed, equality at distinct c and
d would imply L_(c+d)^T x = 0. Multiplication by c+d is invertible, so this
contradicts x ≠ 0. Both families fix zero. For fixed nonzero x, V x is
nonzero, and multiplication by U permutes the nonzero outputs. Thus every permitted family
has exactly the fixed-input distribution of a uniform GL(16,2) matrix.

Fix a complete message. Independent outer randomizers therefore give the
same joint distribution of its randomized symbols. The subsequent shuffles
also have the same distribution. To compare inner outputs, condition on
this outer output and all preceding inner randomness. The next input x and
entering state a are then fixed. The emitted value x + A a is independent
of the current randomizer. The conditional law of M_i a + A^T x is unchanged
by the replacement. Induction over physical steps proves equality of the
entire output distribution for the fixed message.

Let Z be the number of nonzero messages whose output weight is at most
13,107. Summing the preceding equality over messages preserves E[Z].
The retained proof bounds Pr[Z ≥ 1] by E[Z], so its complete first-moment
bound transfers unchanged. No independence between different messages is
required. The random matrix families need not have the same distribution.

For transposed encoding, sampling the forward family L_c^T is convenient:
the encoder applies L_c. The existing tower identifies F with
GF(256)[u]/(u²+u+0x20), where GF(256) has modulus 0x11b. The absolute trace
of 0x20 is one, so the quadratic is irreducible. Multiplication has the
three-product circuit

\[
p_0=c_0x_0,\qquad p_1=(0x20\,c_1)x_1,\qquad
p_2=(c_0+c_1)(x_0+x_1),
\]

with output (p_0+p_1, p_0+p_2). The products are byte-field products.
`k16_design/field_symbol_randomizer.py` supplies the retained `mul16` and
`coefficients16` reference routines. This arithmetic count predicts no
particular end-to-end timing.

The native outer layout uses a fixed coordinate permutation Q. It treats
the sixteen coordinates as a four-by-four array and transposes that array;
its index permutation is q(c)=4(c mod 4)+floor(c/4), for 0≤c<16.
Thus Q = Q^T = Q^-1. If B_c = Q^T L_c^T is the forward symbol family, then the
transposed encoder applies B_c^T = L_c Q. This is the permitted case U=Q^T
and V=I. The older input representation supplies Q times the physical
sixteen-coordinate input. A native-layout kernel may absorb this fixed
permutation into its symbol map or its loads. Correctness testing must
compare these coordinate conventions explicitly. The fixed-input law is
uniform nonzero in either convention, so the first-moment transfer applies.

The transfer requires independent uniform nonzero scalars at the replaced
locations. Reusing a scalar across locations or removing refreshes changes
the fixed-message process and is not covered. A seeded benchmark draw also
does not turn this ideal-ensemble bound into a guarantee for that seed.

The new t128 screens are separate constructions. Their physical expansion,
feedback, and step count differ, so the t64 numerical certificate does not
transfer to them. `screen_t128.py` records their explicit maps and recomputes
their bounds. It keeps its whole-code certificate flag false.
