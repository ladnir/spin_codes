# Randomized outer outputs: the new high-distance baseline

Independent outer-output scalars restore the proof mechanism used by the
successful binary packet code. Their measured cost is about 0.48 ms at K20.
The new forward prototype takes 7.02 ms, compared with 6.53--6.55 ms for
the same executable's identity-diagonal control. The sub-5-ms target and
the complete 40%/40-bit distance certificate remain open.

This is a construction change, not a new implementation of the fixed-outer
ensemble. The fixed-outer sources, proofs, and measured FFT/cube controls
remain unchanged. Production profiles and the certified high-rate winner
are not modified.

## Construction

Work over F=GF65536, with the AES GF256 subfield. One 128-bit application
record contains eight parallel F-symbol lanes. Every lane uses the same
scalar code. At K=2^20, there are 8192 outer rows and N=2^21 output records.

1. Encode each 128-symbol row with the fixed RS[256,128,129] evaluation
   circuit. Its coefficients are in GF256; its messages are in F.
2. Independently multiply each of the 256 coordinates of every row by a
   uniform nonzero F scalar. Sample these scalars once during setup.
3. Collect outer coordinate j into region j. Independently permute the
   8192 records of each region, using the existing regional route.
4. Process packets of 64 records with the shared beta64 inner. Its state
   has four F symbols and persists across region boundaries. Initial
   state is zero; final state is free.

The inner is unchanged:

    Y_i = F_64 X_i + A Q_i,
    Q_(i+1) = M_i Q_i + C X_i.

The eight width-eight feedthrough blocks, A, and C are fixed GF256 maps.
The mixer is the existing independent invertible rank-one update
M_i=I+u_i v_i^T over F, not a uniform GL(4,F) matrix. Sample nonzero u_i
uniformly and v_i uniformly subject to v_i^T u_i != 1. Outer scalars,
regional permutations, and packet mixers are independent.

There is one routing stage, no pre-accumulator, no extra local shuffle,
no puncturing, and no terminal flush. The local feedthrough is still needed:
randomizing the outer does not remove the deterministic quiet-subcode
obstruction for an identity feedthrough with width 64 and state dimension 4.

The mathematical claims use the specified independent setup ensemble.
The native cost prototype uses separate reproducible MT19937 streams;
its timing and arithmetic checks are not a certificate for a short-seed
setup distribution.

## What this changes in the proof

For a fixed outer word, independently scaled nonzero coordinates become
independent uniform nonzero F values. Correlations between different
messages remain; the proof does not assume their independence.

For every fixed nonzero row output x, MDS shortening bounds the expected
number of messages mapping to x by (65535)^(-128). Thus a randomized row's
expected counting measure is dominated by a scalar multiple of uniform
F^256. This is the whole-row comparison used in the binary proof, now at
the actual routed F-symbol granularity. It holds for every fixed route.

Moreover, the exact expected row measure depends only on its support size.
An additional within-row permutation would not change that measure. The
binary predecessor needed its local shuffle after splitting a larger
randomized symbol into smaller routed packets; this construction does not
perform that split.

The proof still must cap conditional failure before averaging rare routes.
The old large-kernel route example survives every nonzero diagonal scaling.
Randomizing values does not make an unconditional word first moment safe.

See the [class proofs](HD_RANDOMIZED_OUTER_PROOF.md) for exact counting
measures, directed bounds, and the routing-conditioned composition. The
[randomizer audit](HD_RANDOMIZED_OUTER_OPTIONS.md) explains why a shared
row scalar does nothing, and why GF256-only coefficients cannot inherit
the full-field dense bound. One coefficient per row may be fixed to one
without changing the random code-image distribution. The prototype keeps
all coefficients independent; the potential arithmetic saving is only 1/256
of the new layer.

## Verified proof ranges

Distance counts nonzero F symbols. Set D=838860, just below 40% of N.
A quiet word has zero state at every emitted packet; the final successor
state is unrestricted. The following are setup-failure bounds, unioned over
all messages and all active-row sets in each stated class.

| Class | Active outer rows | Failure upper bound |
|---|---:|---:|
| All words of weight at most D | 1 | 2^-391.9739 |
| All words of weight at most D | 7548--8192 | 2^-326.6848 |
| Quiet words of weight at most D | 5940--8192 | 2^-3832.4504 |
| Quiet words at any weight | 1--953 | 2^-142.1972 |

The [joint quiet/output refinement](HD_RANDOMIZED_OUTER_JOINT_QUIET.md)
extends the original dense quiet tail from 6566 down to 5940. It expands the
output-weight tilt into zero-coordinate constraints, then counts the
remaining feedback rank in each tile. Convexity of the resulting occupancy
bound certifies the interval from two directed endpoints, not interpolation
between sampled values.

The [weighted-state refinement](HD_RANDOMIZED_OUTER_STATE_PROFILE.md)
also closes the complete one-active-row class, including nonquiet histories.
It tracks a flat-density cap separately from arbitrary concentrated mass,
then uses the exact regional coefficient rather than a Bernoulli marker.
The reported 391.97-bit bound unions all 8192 choices of that active row.

These directed bounds do not cover quiet rows 954--5939 or nonquiet rows
2--7547. They must not be reported as a whole-code margin. The shared
eight-lane implementation needs no eight-fold failure union: a nonzero
record contains some nonzero scalar lane, whose support is contained in
the record support.

## Measured implementation cost

The randomizer is fused into the existing sequential region-staging pass,
before random route reads. It adds no full-code memory pass. The prepared
coefficients occupy 8 MiB at K20, plus 4 MiB of raw coefficients retained
by this research plan for its independent oracle. The identity control
allocates the same plan but compiles out all coefficient loads and products.

Peach Ryzen 7950X, CPU2, 4.5-GHz userspace request, boost disabled; GCC 15.2,
Release AVX512/GFNI, -O3 -mtune=znver4 -fno-lto. Each full measurement uses
three warmups and 101 calls. Calls regenerate the source outside timing.
Setup, allocation, checksums, and validation are excluded. Required outer
work, staging, routing, arithmetic, and output stores are included.
All native measurements ran serially under the shared benchmark lock.

| Full forward encoder | First median (p10--p90), ms | Reverse-order repeat, ms |
|---|---:|---:|
| Identity-diagonal control | 6.527616 (6.504613--6.568251) | 6.549697 |
| Independent GF65536 scalars | 7.015246 (6.990029--7.046454) | 7.023822 |

The paired overhead is 0.488 ms (7.47%) in the first launch and 0.474 ms
(7.24%) in the reverse-order repeat. A separate 51-call phase measurement
gives the following diagnostic split; phase medians need not add exactly.

| Phase | Identity, ms | Randomized, ms |
|---|---:|---:|
| Outer | 2.306536 | 2.287190 |
| Staging, routing, and inner | 4.224035 | 4.732995 |
| Full instrumented call | 6.535671 | 7.021056 |

The outer circuit is identical in both variants. The randomization cost
appears in the second phase, where it is implemented. The new cost is about
2.37 times faster than the retained 16.6-ms RS baseline, but has no complete
distance guarantee yet; this is not a certified replacement comparison.

The release and ASAN/UBSAN correctness runs pass. They check all 65535
nonzero scalars on eight lanes and both extension components; every packet
basis direction; and independent complete encodings at K8192, K16384,
and K24576. Both FFT and cube schedules, 16/32-bit routes, exact aliasing,
offset pointers, and guards are covered. The oracle uses schoolbook
extension-field multiplication, not the native multiplier being tested.

## Reproduction and source identity

New implementation files:

- [HdBlockwiseRandomizedOuter.h](../../../spin/experiments/char2_field_spin/HdBlockwiseRandomizedOuter.h)
- [hd_blockwise_randomized_forward.cpp](../../../spin/experiments/char2_field_spin/hd_blockwise_randomized_forward.cpp)

Build the isolated research target, then validate before timing:

```text
cmake -S spin/experiments/char2_field_spin -B build/hd-blockwise -DCMAKE_BUILD_TYPE=Release
cmake --build build/hd-blockwise --target spin_char2_hd_blockwise_randomized_forward
build/hd-blockwise/spin_char2_hd_blockwise_randomized_forward --check
build/hd-blockwise/spin_char2_hd_blockwise_randomized_forward --bench 101 full 20 fft both
build/hd-blockwise/spin_char2_hd_blockwise_randomized_forward --bench 51 phases 20 fft both
```

Pin a core and acquire the shared benchmark lock for comparative timings.
The prototype requires AVX512 and GFNI; unsupported hosts return 77.
`both` runs the two variants serially. No production API is added.

SHA256 identities of the measured files:

| File | SHA256 |
|---|---|
| Randomized header | 2df96821138ef335e127feab6fe85a364baee2f9bff8ba5390bc7258afa1d306 |
| Driver | 761021724b3a37cfca7cadf40536f4279eede0a4a8a33763f20e788ac0df018b |
| Peach release executable | 513a1f5de977136b3969feed13c8c6b5d10db9a31abc2c628b8fd161821e9482 |

The identity checksum is 5184726652255349648, matching the prior fixed-outer
FFT control. The randomized checksum is 15459768057287540600.

## Coarse state screen and next step

The new [floating diagnostic](theory/hd_randomized_outer_state_screen.py)
keeps a two-state zero/nonzero envelope. It bounds touched-tile cosets by
their zero-coset weight moments and charges at least 61 minus the touched
positions for a nonzero old state. It uses the actual rank-one return law.
It averages input values before fractional powers and retains the regional
occupancy normalization. Both one-step and chronological two-step versions
are included; the latter composes before taking entrywise powers.

At the nine sampled occupancies 1,16,128,512,1024,2048,4096,6144,7547,
the coarse bounds are uninformative. Replacing the mixer by an idealized
uniform nonzero-state refresh in the screen does not close those samples
either. These failures are failures of the chosen upper bound, not low-weight
codewords or evidence that increasing state size is required.

The subsequent three-component operator carries zero mass, a uniform
nonzero-state point-density cap, and arbitrary concentrated nonzero mass.
New births enter the concentrated component. Empty-packet refreshes contract
that component by a factor at most 1/q per packet and feed the flat-density component.
This retains mixing through empty gaps without asserting uniform births or
uniform states after output conditioning.

For one active row, each region has exactly one potential input position.
The exact coefficient counts its 128 possible packets; collisions are
impossible. At the same weight mark, the improved density bound with a
Bernoulli marker is still uninformative, whereas exact regional extraction
closes the class at 391.97 bits. This is direct evidence of coefficient-bound
slack, not a reason to enlarge the inner state.

An ablation separates the two improvements. With exact regional extraction,
even the earlier two-state envelope gives a floating log-bound of -204.44
for the actual rank-one mixer. The three-component envelope improves it
to the directed -391.97 result. For one row, exact routing is the primary
repair; retaining state density adds about 187.5 bits. The two-state ablation
is a floating diagnostic, not a separately issued certificate.

The subsequent [sparse ladder](HD_RANDOMIZED_OUTER_SPARSE_LADDER.md) now
closes one through five active rows with a 56.51-bit union margin. It combines
exact regional coefficients, chronological four-packet blocks, and a bound
on a fixed state's membership in a random span of feedback columns.
Six active rows remain unclosed even with sixteen-packet blocks.

The [quiet coefficient audit](HD_RANDOMIZED_OUTER_QUIET_EXACT.md) extends
quiet coverage to 1 through 962 and 5865 through 8192. It also separates
marker-grid loss, coefficient loss, and local-rank loss. At four specified
interior witnesses, even perfect residual ranks cannot close the existing
packet-averaged envelope. This is a limitation of those witnesses, not of
the construction.

The [occupied-tile moment](HD_RANDOMIZED_OUTER_FLAT_MOMENT.md) retains state
constraints lost by the earlier untouched-tile analysis. Its
[composed diagnostic](HD_RANDOMIZED_OUTER_LOCAL_REFINEMENTS.md) still leaves
the selected middle classes open. The next proof work retains finer packet
layouts before fractional powers and the outer's nonzero-support restriction.
The implementation and its measured source identities remain unchanged.

The [tile-conditioned calculation](HD_RANDOMIZED_OUTER_QUIET_TILE_CONDITIONING.md)
recovers about 26,584 bits at one fixed middle witness by delaying the
profile average. Re-optimizing that occupancy still leaves a positive
log bound of about 5439. The [support-mark audit](HD_RANDOMIZED_OUTER_QUIET_SUPPORT.md)
and [FFT16 probe](HD_RANDOMIZED_OUTER_FFT16_PROBE.md) likewise supply no
new whole-code certificate. The next experiment should retain actual
input and output weights together; it should not repeat a wider scan of
the same coarse envelope. There is no basis yet to optimize or promote
a replacement high-distance implementation.
