# Near-half distance with the retained Toeplitz-window design

Started 2026-10-21. The user selected this design for the high-distance task
and [parked the small-state investigation](PARKED_SMALL_STATE.md). Preserve
the existing [45% certificate](STRATIFIED_CLOSED_CERTIFICATE.md) unchanged.

**Current checkpoint: >49.1% symbol distance, with setup failure below
2^-40, at K=2^20 and rate 1/2 over GF(2^128).** Use the same RS[256,128]
outer and stratified route, with Toeplitz memory 1024 and batch 256.
The [native follow-up](NEAR_HALF_NATIVE.md) records the exact implementation
and an explicitly proved coordinate variant, including their measured costs.

## Construction and scope

Fix F=GF(2^128), q=|F|, K=2^20, N=2K, and B=256. Every 128-symbol message block
uses the same RS[256,128] outer with GF256 evaluation points. Region j
collects coordinate j from all 8192 outer blocks, in an independent uniform
order. Concatenate the 256 regions. This remains one global permutation.

For memory L and batch width b dividing L and N, let S_t be the previous
L output symbols, initially zero. Sample independent b-by-L Toeplitz maps
T_t over F, with b+L-1 independent uniform diagonals per map, and compute

    Y_t = X_t + T_t S_t.

Append Y_t to the state and discard its oldest b symbols. Emit N symbols;
the final state is free. All setup is sampled once and fixed for all messages.
The inner is invertible for every setup. The rate is exactly 1/2.

There is no requirement that b<=B. The triangular Toeplitz argument makes
T_t S jointly uniform on F^b for every fixed nonzero S. For a fixed output,
states are determined by that output; independence across maps therefore
gives the retained inverse-input distribution. The finite proof applies
directly, without an additional asymptotic batch-ratio assumption.

The working target was 49% distance and at least 40 bits of setup-failure
margin, not another 35% checkpoint. The completed result below exceeds it.
The full-field volume threshold at rate one half is approximately 49.2189%;
this is a limit of the current
random-volume proof, not a universal upper bound on code distance.

## Why the original memory is insufficient for this bound

Write n=N/b and m=L/b. Classify output activity by the number c of intervals
and their total length v in batches. The random inverse coordinates have
size a=b(v-c). The free output positions number f=bv-L(c-1).

At (L,b)=(256,64), take c=128 and v=16512. Then a=N/2 and f=1024256,
below the 49% bad-weight threshold e=1027604. The output-volume bound is
therefore q^f. The retained route moment includes

    H(theta) = E exp(theta (Z-128)_+),  Z~Bin(256,1/2).

Its mean satisfies

    8192 E(Z-128)_+ = 2^19 binom(256,128)/2^256 > 24320 = a-f.

Jensen's inequality makes theta(f-a)+8192 log H(theta) positive for every
theta>0, before charging interval placements. Thus no marker refinement
can certify this class using the current envelope and these parameters.
This is not a low-distance codeword or a proof that the construction fails.

Increasing L strengthens the forced-zero tails. Increasing b as well can
amortize the larger window while preserving the same construction.

## Numerical approach

The new bounded driver is `theory/climb_stratified_batch.py`. It imports the
three retained certificate helpers without changing them. It enables the
previously proved entropy output-volume bound, and optionally tests rational
markers 1+i/32 for i=1,...,32, plus 4 and 8. Markers affect only the proof.

Every accepted rectangle is replayed with exact rational route moments and
outward-rounded Decimal arithmetic. An independent integer checker verifies
adjacent, nonoverlapping v-intervals for every feasible c. It checks the exact
number of feasible pairs, not just a sum of rectangle areas. Removed,
duplicated, and boundary-modified boxes must fail this partition check.
Floating scores and individual passing classes are never certificates.

Complete-cover/replay jobs are serial, capped at 50 seconds internally,
55 seconds externally, and 1 GiB virtual memory. Selected-point diagnostics
use the external caps. Runs use Python -B and write no data.
The 80- and 120-digit replays are separate jobs. No native kernel or paper
is changed, and no encoder timing is inferred from certificate runtime.

## Initial bounded diagnostics

These are failed per-term cover attempts or selected-point diagnostics,
not complete upper bounds and not distance counterexamples.

| L,b | Distance | Check | Floating log2 bound |
|---|---:|---|---:|
| 256,64 | 48% | First rejected class v17641,c479; original markers, target -68 | -65.57 |
| 512,128 | 49% | Near-half class v8248,c56; refined markers | +848.25 |
| 768,256 | 49% | First rejected class v4487,c158; refined markers, target -64 | -60.66 |
| 1024,256 | 49% | Near-half class v4125,c29; refined markers | -635.47 |
| 1536,512 | 49% | Near-half class v2070,c22; refined markers | -1535.25 |

## Closed 49% checkpoint

At (B,L,b)=(256,1024,256), the finite whole-code setup-failure probability
for output weights at most 1027604 is below 2^-40. Thus the minimum symbol
distance exceeds 49% except with that probability. This covers all nonzero
messages chosen after the complete setup, not only a fixed input.

Both 80- and 120-digit directed replays passed on the final source. The
851 boxes cover exactly 6714983 feasible (v,c) pairs. Every per-pair bound
is at most 2^-64, and the deliberately conservative rectangular class count
is below 2^24. The worst directed per-pair log2 upper bound is
-73.5691515254769370487277493838943732495484601135630732647195858891857609567.
This supplies spare margin, but the claimed bound remains 40 bits.

Partition SHA-256:
`c971f0f6a77e0e084f09e31b143370127fcf1a3504cbb2fd07859b6b3b787bc7`.
The removed-box, duplicated-box, and modified-boundary controls all rejected.
Final-source total replay times were 17.28 and 17.39 seconds respectively;
these are verifier times, not encoder timings.

From the repository root, run the following separately (not concurrently),
with an external 55-second/1-GiB limit:

```text
python -B research/workstreams/one_route_coupled_scan/theory/climb_stratified_batch.py --distance 49/100 --L 1024 --batch 256 --fine-markers --check-controls --precision 80 --seconds 50
python -B research/workstreams/one_route_coupled_scan/theory/climb_stratified_batch.py --distance 49/100 --L 1024 --batch 256 --fine-markers --check-controls --precision 120 --seconds 50
```

Retained source SHA-256 manifest:

| File in theory/ | SHA-256 |
|---|---|
| climb_stratified_batch.py | 825f8d813e1abf472583089b4cf3520f54dd2a84fa4bced11e4f3f8f616e4bf2 |
| stratified_batch_certificate.py | c7e21a4a3e8a70721f0cfec915feb488107c0144e79aa05fdb05f661ce32e3c8 |
| stratified_certificate.py | a7a28734f9ad271c409eaa2623e0602451418c2bb0852ff626663df12c71ddfc |
| stratified_fractional.py | 7f129e764ec02f5ba3848f68ce4ac57d7746c82c947fd642e6fb4d616c314590 |

The last three helpers are unchanged from the 45% checkpoint. Independent
analytic review checked the larger-batch interface and finite count;
independent static review checked the new wrapper, partition checker, and
failure controls. Agreement of two precisions is corroboration, not a
substitute for the outward-rounding argument.

## Closed 49.1% checkpoint, with the same encoder

Without changing B, L, b, or any setup distribution, raise the threshold to
e=floor((491/1000)N)=1029701. The whole-code setup-failure probability remains
below 2^-40. Thus the minimum symbol distance exceeds 49.1% except with that
probability. Both 80- and 120-digit directed replays passed on the source
manifest above.

This cover has 926 boxes and the same 6714983 feasible pairs. Its partition
SHA-256 is
`b556edf695389e6b7f028b39c7704fd8a7de9175099cefa417af51a2a508a953`.
The worst directed per-pair log2 upper bound is less than -72.23754810968,
again stronger than the uniform -64 acceptance threshold. The conservative
union argument still claims 40 bits. Total verification times were 20.05
seconds at 80 digits and 19.79 seconds at 120 digits. Both runs passed the
corrupted-cover controls.

Reproduce by replacing `--distance 49/100` with `--distance 491/1000` in the
two commands above. Root ran these jobs serially on Python 3.13.7 with
NumPy 2.5.2. An immediate-timeout control rejected with no certificate.
An exact integer check also confirmed the original-parameter Jensen
obstruction. Local and remote copies match the complete source manifest.

The next task was performance engineering for this fixed high-distance
construction, not renewed work on the parked 35% proposal. That work is now
recorded in the [native follow-up](NEAR_HALF_NATIVE.md). Neither the
existing 3.798 ms GF256 timing nor any smaller-field window timing applies
to this new profile. Keep the proof certificate intact while measuring and
optimizing the actual Toeplitz products.

## Cost boundary

The baseline stores 319/64=4.984375 sampled coefficients per output.
The 1024/256 proposal stores 1279/256=4.99609375, with a fourfold larger state
and batch. The 1536/512 alternative stores 2047/512=3.998046875.
These are setup counts, not encoding times. Direct multiplication still
costs L full-field products per output. Fast middle products or transforms
need separate implementation work; larger transforms do not automatically
inherit the existing 256-point GF256 kernel.
