# Gamma30 at 40% symbol distance

The [v31 index](theory/hd_field128_gamma30_d40_partial_cover_v31_20261010.json)
covers every active-group count from 1 through 8192. Its 274 event bounds
have a directed union ceiling of `2^(-49.171469646557)`. The separate
[completion audit](theory/hd_field128_gamma30_d40_completion_v31_20261010.json)
passes: the 40% rung is closed with 49.171469646557 bits of ideal-setup
failure margin. The frozen index intentionally retains
`whole_code_certificate=false`; the separate completion manifest records
the whole-code conclusion.

## Construction and probability space

Let `E=GF(256)` and `F=GF(2^128)`, with a fixed field embedding `E -> F`.
The code has source dimension `K=2^20` and output length `N=2^21`, both
measured in `F` symbols. Its 8192 outer groups each encode 128 source
symbols using four fixed RS[64,32] constituents. Independently for each
group and each of its 64 columns, sample a uniform invertible 4-by-4
Toeplitz matrix over `F`.

Use the alternating pair route from the
[lift proof](HD_FIELD128_LIFT_PROOF.md): each of 256 regions has an
independent uniform permutation of its 4096 pairs and independent fair
pair swaps. Let `A` be the 64-by-30 Lagrange evaluation matrix from
nodes 128 through 157 to nodes 0 through 63 in `E`, and put `C=A^T`.
Fix a scalar `gamma` of degree 16 over `E`. Let `H_30` interpolate from
nodes 0 through 29 to nodes 30 through 59 in `E`.
For each chronological packet, use

    Y_t = X_t + A Q_t,
    Q_(t+1) = D2_t H_30 D1_t (Q_t + gamma C X_t).

The diagonal entries of every `D1_t,D2_t` are independently sampled
from `F` minus zero. Outer maps, route choices, swaps, and diagonals are
mutually independent. They form one setup `Theta`, fixed across messages.
The initial state is zero; state persists across regions; the final
state is unrestricted. The
[gamma return proof](HD_FIELD128_GAMMA_RETURN.md) establishes the
required local bounds for this fixed scalar, without a bad-gamma event.

Write `C_Theta:F^K -> F^N` for the resulting linear encoder. The outer
maps are injective. For fixed setup, recover each packet by
`X_t=Y_t-AQ_t` and then compute the next state. Thus the routed inner
map is invertible, and the code has dimension `K` and rate exactly 1/2.

## Completion statement

Set `L=floor((2/5)N)=838860`. The linked proofs and audited event bounds give

    Pr_Theta[there exists m != 0 with wt(C_Theta(m)) <= L]
        <= 2^(-49.171469646557) < 2^-40.

Consequently, with probability at least `1-2^-40` over ideal setup, the
code has minimum symbol Hamming distance at least `838861`, which exceeds
`(2/5)N=838860.8`.

To see the coverage implication, classify a nonzero message by the number
`r` of its nonzero outer groups. Injectivity makes `r` an integer in
`1..8192`; `r=0` contains only the zero message. The indexed weight
ranges partition these integers. Each event bounds the existence of a
forbidden codeword among all supports in its range, including every
feasible ordered even/odd matching profile. The terminal event covers
`6678..8192` simultaneously and contributes once, including `r=8192`.
Taking a union over the 274 events proves the displayed bound. Their
weight ranges are disjoint; the events in setup space need not be.

The [all-pair proof](HD_FIELD128_ALL_PAIR_COMPOSITION.md) keeps support
multiplicities outside the fractional power and includes the projective
divisor `q-1` once in each complete bound. The chronological transfer
starts from zero state and sums both terminal state classes. The
[integer-band argument](HD_FIELD128_MIDPOINT_WEIGHT_BANDS.md) evaluates
every original integer weight, without interpolation. Complementation
changes only support counts, not the encoder, regional matrices, or
outer prefactor. Where witnesses use different marks or fractional
exponents, the [profile-minimum argument](HD_FIELD128_PROFILE_MIN.md)
takes the minimum of complete bounds for the same class. It neither
mixes their matrix entries nor introduces another setup distribution.

## Pinned record and scope

The index SHA256 is
`b98f340ba669ae9712be7ff9ac401af5bcccb7e83d93bbccf5025f37173664f0`.
Its [explicit selection](theory/hd_field128_gamma30_d40_v31_selection_20261010.json)
has SHA256
`065dc8d8432cf70d1e4e02eb7cbadf3efc326155db78207f5192b598808bc6c9`.
The index preserves all 265 prior entries and adds two composite bands
and seven selected weights. Its 384-bit positive aggregation records
the rational logarithmic ceiling `-49171469646557/1000000000000`.
The separate completion manifest has SHA256
`99d7fe81ccd96e6b0bdb44cca4368ad2db8325a64e3f68c52294a10356a0de0c`.

The completion audit reconstructs the index, preserves its 265-entry parent
prefix, checks complete coverage and the terminal event, and independently
recomputes the positive union at 384-bit precision. It checks byte bindings;
it does not generate coefficients again. Its conclusion inherits the
reviewed numerical producers, saved-data audits, and recorded remote
archive verification. Historical optimization ancestry is not a proof
premise for a fresh bank built from exact positive witness literals.

Relative to the literal 35% construction, this variant replaces four
RS[16,4] emission blocks and state dimension 16 by one full RS[64,30]
emission map, state dimension 30, and feedback `gamma A^T`. Sharper
counting, paired composition, and retuned proof witnesses also improve
the bounds without further changing this gamma30 encoder. This is not
a claim that a field lift or parameter retuning alone establishes 40%.
It supplies no timing improvement, concrete pseudorandom-seed guarantee,
bit-Hamming distance claim, or certificate for shared setup randomness.
