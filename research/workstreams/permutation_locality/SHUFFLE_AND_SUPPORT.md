# Shuffle screening and a spectrum-only support bound

2026-09-23. The best confirmed candidate remains the uniform 16-row streaming
kernel, about 7.86 ms versus a 9.57 ms production control. This iteration
rejects three implementation directions and supplies a proof component for
counting grouped BCH words. Neither the 2x target nor a new full-distance
certificate is achieved. Production code and the paper remain unchanged.

## Performance screens

Measurements use K=2^20, 128-bit elements, the existing BCH and IMT maps,
CPU 15 on Peach, three warmups, and serial execution under the shared locks.
Setup, allocation, initialization, and correctness checks are excluded.
These are selection screens, not repeat confirmations of close rankings.

The 101-call shuffle screen gives:

| Encoder | Full time (ms) |
|---|---:|
| Production distribution and kernel | 9.724 |
| Uniform 16-row groups, streaming | 7.847 |
| XOR lane shifts, existing tiled kernel | 8.748 |
| XOR lane shifts, SIMD streaming kernel | 7.976 |

The XOR candidate replaces each independent uniform lane permutation by
`lane -> lane XOR a`, where a is sampled independently from {0,...,15} for
each group in each region. Group positions and shared coordinate permutations
retain their previous distributions. A fixed row still has the original
one-row marginal, but joint row positions change.

The SIMD path fills its lane buffer sequentially, shuffles four 128-bit lanes
per AVX-512 register, and reorders the register stores. It does not improve
the complete encoder. A separate phase diagnostic gives 3.202 ms for uniform
inner/routing and 4.675 ms for repack/BCH. The XOR counterpart gives 3.326
and 4.712 ms. These means include warmups and instrumented code; they are
not replacements for the uninstrumented timings.

The XOR feedback census enumerates all 4,336 translation orbits of 16-bit
supports. For each orbit it enumerates all eight aligned windows and all
16 shifts. Exactly 32 nonempty supports can have zero feedback. Each has
weight eight and zero feedback with probability 1/128. Every other nonempty
support has zero probability. This worsens the worst fixed-support probability
from 1/51,480 under uniform lane permutations. Averaging over all weight-eight
supports recovers the uniform value; such an average cannot replace the bound
for a fixed message. A full-rank zero-state check passes at seed 17, but that
does not certify distance.

Direct BCH loads were mechanically specialized for column strides 16 and 20,
preserving the original circuit. Stride 20 adds one cache-line gap per column.
They avoid repacking but take 10.478 and 10.679 ms, versus 7.841 ms for the
repacking kernel in the same screen. The production control takes 9.588 ms.
Keep the repack for this implementation.

## Split groups

A larger transfer need not imply that all group inputs occupy one epoch.
The split family uses shared groups of G=32,64,128 rows. Divide each region
into G/16 bands. Each band independently permutes the group labels. A fresh
uniform permutation of G lanes in each group and region assigns 16 lanes to
each band. Thus every row still contributes once per region, with its original
one-row marginal. Tested band boundaries align with epochs.

The reverse kernel buffers one region's values, then streams a complete G-row
transfer after processing the final chunk of a group. This admits several
feedback constraints per group per region: the previous one-epoch constraint
count no longer applies. This is not a positive distance proof. The complete
zero-state matrix has full rank for G=32 at seed 17.

The 51-call screen gives:

| G | Same-route tiled kernel (ms) | Split streaming kernel (ms) |
|---:|---:|---:|
| 32 | 9.034 | 8.733 |
| 64 | 9.797 | 8.742 |
| 128 | 9.569 | 9.280 |

The matched production control is 9.629 ms and the prior 16-row streaming
candidate is 7.899 ms. The tested buffering schedule loses to that candidate.

## Bounding joint supports from an ordinary spectrum

Shared coordinate permutations require counts of tuples of BCH words with
small union support. The following bound uses only the ordinary spectrum.
It avoids requiring an exact higher-order weight enumerator, although its
tightness for a full SPIN proof is not yet established.

Fix a binary linear code C of length n and minimum distance at least d>0.
Let A_w count its nonzero words of weight w. For integers g>=h>=1, let
T_{g,h}(u) count ordered g-tuples of C-words whose span has dimension h and
whose union support has size at most u. Zero words are permitted in the
g-tuple. Define the nonzero enumerator and a cumulative coefficient sum:

    P(z) = sum_{w>=1} A_w z^w,
    S_h(W) = sum_{j<=W} [z^j] P(z)^h.

Here [z^j] selects the coefficient of z^j. Thus S_h(W) counts ordered
h-tuples of nonzero words with total weight at most W, including dependent
tuples. Including dependent tuples makes it an upper bound on the number
of such ordered bases.

For an h-dimensional subcode V with union support U, choose an ordered basis
uniformly among all its ordered bases. This is auxiliary randomness used for
counting; the code C and subcode V are fixed. Each basis vector is marginally
uniform over V's nonzero words. Each used coordinate is one in exactly
2^(h-1) words. Consequently the expected total basis weight is

    h * 2^(h-1) * |U| / (2^h-1) <= mu,
    mu = h * 2^(h-1) * u / (2^h-1).

Every basis has total weight at least hd. For any integer W>=hd with W+1>mu,
the fraction of bases with total weight at most W is at least

    p = (W+1-mu) / (W+1-hd).

Indeed, a basis above W has integer weight at least W+1. Lower-bounding the
weight by hd on the remaining bases proves this inequality by averaging.
If mu<hd, no such subcode exists.

There are H_h=product_{j=0}^{h-1}(2^h-2^j) ordered bases of V. At least

    L_h(W) = max(h!, ceil(H_h*p))

of them have weight at most W. The h! term follows because one such basis
exists and all its reorderings have the same weight. An ordered basis belongs
to exactly one subcode. Finally, each h-dimensional subcode has exactly
N_{g,h}=product_{j=0}^{h-1}(2^g-2^j) ordered g-tuples that span it. Therefore

    T_{g,h}(u) <= floor(N_{g,h} * S_h(W) / L_h(W)).

Minimize over admissible W. Coefficientwise upper bounds on A_w may replace
the exact coefficients: all polynomial operations and cumulative sums are
nonnegative. A valid lower bound on d remains necessary. The usual binary
Griesmer sum also implies T_{g,h}(u)=0 when
u < sum_{j=0}^{h-1} ceil(d/2^j).

For h=1, use the exact identity

    T_{g,1}(u) = (2^g-1) * sum_{w<=u} A_w.

`joint_support.py` implements these bounds using integers and rational
arithmetic. Its tests exhaustively enumerate 16,640 tuples in five small
linear-code cases and check 103 rank/support inequalities. They also check
coefficientwise inflation and the rank-one identity. These tests corroborate
the derivation; they do not replace it or establish a SPIN certificate.

The retained BCH shell caps are available under
`research/workstreams/bch_rm2sub_bridge/generated/`. The existing
`verify_progress.outer_caps()` path authenticates their dependencies. This
iteration has not fed those caps into the new bound or reconstructed that
authentication chain. Their later use must retain that verification.

## Validation and next step

All timed candidates first pass full-output comparisons, including suffix
preservation, against the library with the same supplied route. Dense and
adjoint checks pass for the XOR family at K=2^14 and the split G=32 family
at K=2^15. ASan/UBSan passes for the new XOR and strided kernels at K=2^14,
and split kernels at K=2^15,2^16,2^17 for G=32,64,128 respectively. All four
existing library tests pass under sanitizers. The remote log is
`/tmp/spin-locality-JjY7yR/sanitize-next.log`.

Replay scripts are `run_xor.sh`, `run_strided.sh`, `run_split.sh`, and the
extended `run_sanitize.sh`. `xor-census` runs the complete local census.
Run `python joint_support.py` for the exact small-code checks. Raw timing
receipts remain outside version control. The final release executable hash is
`e88da74236b3c32732447df0c7a7e891594233acee77cfe3d8430c5e5dc2f872`;
the production library hash is unchanged.

Next, apply the support bound to authenticated BCH caps and test whether it
can control one active 16-row group. Then combine it with the grouped inner
process, rather than transferring an independent-row certificate. On performance,
the repack/BCH phase is now the larger measured cost. Further progress toward
2x needs a materially different arithmetic or data-layout schedule, not the
tested cheaper lane shuffle or wider buffered transfers.
