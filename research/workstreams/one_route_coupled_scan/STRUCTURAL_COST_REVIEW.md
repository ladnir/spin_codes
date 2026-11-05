# Structural cost review of the persistent inner

The current Cauchy-coordinate inner is a useful proved reference, but its
measured cost remains far above the historical 3.798 ms GF256 encoder.
Completing the one-diagonal proof would certify an intermediate checkpoint;
it would not establish a competitive implementation. This review records
operation counts and two construction alternatives. It adds no implementation,
benchmark, numerical proof run, or whole-code certificate.

## What is measured

The [native record](NATIVE_PERSISTENT.md) supplies the following K20 results.
Each record contains four GF32 lanes. Setup and allocation are excluded;
the runs used Peach CPU2 and were serial. Component timings are not additive.

| Construction | Contiguous inner | Routed inner | Full encoder |
|---|---:|---:|---:|
| Two diagonals, scalar B64/state64 | 10.401 ms | 19.353 ms | 22.973 / 24.246 ms |
| One diagonal, scalar B128/state32 | 4.356 ms | 12.086 ms | 18.529 / 18.162 ms |
| One diagonal, packet-four B64/state128 | 9.311 ms | 12.075 ms | 16.132 / 16.015 ms |

The later full-Cantor two-diagonal controls reach approximately 22 ms.
A matched packet control/single/single/control comparison gives
22.250/17.034/16.934/22.536 ms. Thus 16--17 ms is the appropriate current
one-diagonal range, not a claim based only on its fastest sample.
The [one-diagonal construction](SINGLE_DIAGONAL.md) remains uncertified.

The [historical 3.798 ms reference](../char2_field_spin/HD_SHIFTED_CLOSED_CERTIFICATE.md)
has a different outer, route, inner, and alphabet distribution. It was not
rerun in the same native session. Its complete certificate gives 35% distance
with 55.092738508223 bits of setup-failure margin over GF256, applied in
parallel to sixteen byte lanes.

Compact random keys and outer-store tiling have already received bounded
native comparisons. Neither gave a reliable full-encoder improvement.
Their isolated savings are not grounds for another tuning loop.

## Work intrinsic to the current representation

In [additive coordinates](ADDITIVE_MDS.md), the one-diagonal candidate computes

    V = D_delta Z,
    Y = X + r .* V,
    Z_next = V + Y.

Here D_delta is interpolation on a b-point additive subspace followed by
evaluation on a disjoint coset. The two-diagonal control applies this map
twice and uses two random field products per output.

For the normalized radix-two butterfly representation, write b=2^h.
The inverse transform at coset zero has h*b/2-(b-1) nonzero twiddle products.
At each level, precisely the first node has zero twiddle. The forward
transform on a disjoint coset has h*b/2 nonzero twiddle products. Therefore
one D_delta has the following algebraic count:

    products(D_delta) = h*b - b + 1.

| b | Products per map | Products per output |
|---:|---:|---:|
| 16 | 49 | 3.0625 |
| 32 | 129 | 4.03125 |
| 64 | 321 | 5.015625 |
| 128 | 769 | 6.0078125 |

This counts scalar butterfly multiplications, not machine instructions.
Zero twiddles are excluded; multiplication by one, SIMD duplication,
register shuffles, loads, and stores are not separately priced. GF32
twiddles cost more than GF256 twiddles in the existing implementation.
The full Cantor chain reduces the frequency of those expensive twiddles,
but does not remove the two transforms in each D_delta.

The prepared [GF32 multiplier](../../../spin/experiments/char2_field_spin/WideDiagonal32.h)
uses four GFNI products, three word rotations, four coefficient shuffles,
and two Boolean operations per four records. It reads 64 coefficient bytes
for those four records. The one-diagonal candidate consequently retains
32 MiB of random coefficient material at N=2^21, besides its transform work.

The global MDS reset proof also imposes a structural alphabet requirement:
its scalar output columns use N distinct field points. Replacing GF32 by
GF65536 at K20 does not preserve that interface. Changing state coordinates
alone moves conversion work; it does not automatically remove the cost of
injecting ordinary output coordinates into the persistent state.

These counts and measurements do not prove a machine-level lower bound.
They do show why closing the present proof, shrinking one table, or tuning
one transform stage is insufficient evidence for a 3.8 ms path.

## Primary proposal: amortize full-field randomness on a smaller state

Use the exact outer, paired route, and fixed inner maps of the historical
GF256 reference. Change the scalar field and the refresh coefficients only.
Let E=GF256 and F=GF(2^32), with the existing embedding of E into F.
One 128-bit record carries four parallel F symbols under the same setup.
The proposed scalar code is F-linear and genuinely mixes each group of
four byte lanes.

The retained outer has 8192 groups. Each group contains four RS[64,32]
constituents and independent invertible 4-by-4 Toeplitz mixers over E,
one per output column. Its 256 regions retain the alternating adjacent-pair
matchings, independent uniform pair permutations, and independent fair
pair swaps. This is the old paired route, not the new unshuffled
four-record packet route.

Retain 64-symbol input batches, a 16-symbol state, the four RS[16,4]
emission maps A, C=A^T with CA=0, and the fixed E-valued map H16. Compute

    Y_t = X_t + A Q_t,
    Q_(t+1) = D2_t H16 D1_t (Q_t + C X_t),
    Q_0 = 0.

Sample every diagonal entry of D1_t and D2_t independently from F minus
zero. Setup remains fixed for every message; the final state is unrestricted.
Each refresh is invertible. Sampling from all of F, including zero, would
change the construction and invalidate that persistence property.

This DHD refresh is not the two-diagonal additive feedback map. It is not
uniform on all states for each fixed nonzero input. Its exact support law
comes from the MDS matrix H16: D1 makes the active input amplitudes uniform
and nonzero, H16 determines their output support, and D2 makes amplitudes
uniform and nonzero conditional on that support. The support kernel in
[the grouped-subfield analysis](../char2_field_spin/HD_GROUPED_SUBFIELD.md)
therefore has a well-defined q=2^32 version. This local interface is not
itself a whole-code theorem.

The structural saving is explicit:

- There are 2*16/64=0.5 general-field products per output, instead of one.
- One H16 application serves 64 outputs. Its normalized butterfly count
  is 49/64 products per output, instead of 769/128 for D128.
- Prepared state-diagonal material is 16 MiB at K20, instead of 32 MiB.
- A, C, the outer mixers, and the old paired-route work remain. They must
  not be omitted when estimating full cost.

This is a proposal, not a timing prediction. The old A/C/H geometry was
already explored by the [cheap-lane project](../char2_field_spin/CHEAP_LANE_MIXING_20261016.md).
That project used restricted refresh coefficients, with only 65280 possible
values in one combined coefficient, and still has nonquiet proof gaps.
The proposed distinguishing change is full F* randomness in both state
diagonals, while retaining the inexpensive E-valued outer. It does not add
full-field random multipliers to every outer coordinate.

### Exactly what the quiet-word proof preserves

Fix a realization of the old E-valued outer O and old paired route P.
Let M impose C X_t=0 on X=P O m for every batch except the last. Since the
refresh is invertible, these equations characterize histories whose incoming
state is always zero. Their emitted word is P O m, independently of all
refresh coefficients. M is an E-matrix, so ker_F(M) is the scalar extension
of ker_E(M). Their symbol-Hamming distances are equal: expansion in an
E-basis gives a nonzero component with contained support, and scalar
embedding gives the reverse inequality.

Consequently the [retained quiet proof](../char2_field_spin/CHEAP_LANE_MIXING_GOAL.md#second-round-the-quiet-class-closes-at-40)
also bounds this proposal's quiet bad event at 40% by

    8394752 / 2^100 < 2^-76.

The probability is over the same ideal outer mixers, pair permutations,
and swaps. No new numerical replay was performed for this review.
This transfer does not apply if those E-valued mixers are removed, if the
route is replaced by the current packet-four route, or if A/C or the batch
ordering changes. The present fixed-RS persistent candidate has no such
outer mixers; its distribution is therefore not covered by this quiet receipt.

The [quiet-subcode obstruction](../char2_field_spin/CHEAP_LANE_MIXING_GOAL.md#an-actual-distance-ceiling)
also persists: the unchanged subcode gives distance at most 255/512,
or 49.8046875%, for every setup. That leaves room for 35--45%, but proves
none of the nonquiet classes at those distances. In particular, the quiet
receipt stops at 40%; no 45% guarantee transfers. The separate certified
GF128 45% persistent profile uses different parameters and does not certify
this GF32 proposal.

The next discriminating proof task is a conditional-support DHD bound for
nonquiet trajectories, combined with the actual subfield outer and paired
route. The new shortening and guarded-moment ideas may help, but their
route hypotheses must be re-established for this retained distribution.
Do not relabel a GF256 outer point-density bound as uniform-F density.

## Alternative: split the rate between outer and systematic inner

A different construction could expose its state directly as parity.
For an input batch of b symbols and state dimension s, choose C so that
[I;C] is an MDS systematic map. Let T_t have the existing exact-uniform
feedback property: T_t q is uniform in F^s for every fixed q!=0.
Set

    Q_(t+1) = C X_t + T_t Q_t,
    emitted batch = (X_t, Q_(t+1)),
    Q_0 = 0.

Unlike DHD, this uniform feedback may produce zero. A fixed input/parity
trajectory receives probability q^-s at each nonzero preceding state;
a zero preceding state instead requires Q_(t+1)=C X_t. Thus the state is
observable in the emitted word, and quiet transitions have an MDS counting
interface. These facts do not supply a finite whole-code bound.

To retain total rate one half, the outer rate must be (b+s)/(2b). There is
one route, on the outer output only. However, a message confined to one
outer block has at most B nonzero systematic outputs. Even if every parity
coordinate is nonzero, it yields the deterministic ceiling

    d_min/N <= s/(b+s) + B/N.

Therefore b32/s8 and b32/s16 cannot reach 35% at large N: their parity
fractions are only 20% and one third. These are not useful candidates to
retry for this target. An admissible first shape is b48/s32 with outer rate
5/6, for example RS[96,80]. It routes 60% as many records as a rate-one
inner and uses 0.8 feedback field products per final output. Its parity
fraction is 40%, which itself prevents a 45% target asymptotically.

The non-power-of-two outer and parity maps add real arithmetic, and a
late-activation/quiet-subcode proof remains necessary. This alternative has
a clear route-volume saving but no established path to 3.8 ms.

## Priority

Preserve the one-diagonal code, partial proofs, and 16--17 ms measurements
as an intermediate research checkpoint. For the performance goal, first
test whether the small-state full-GF32 DHD interface closes a meaningful
nonquiet class without restoring expensive outer randomness. Only then
consider a bounded native implementation comparison. Do not equate a
completed one-diagonal certificate with a completed performance objective.
