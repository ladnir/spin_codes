# Optimizing the 35% distance point

## Target and preserved baseline

The target is full forward encoding at K=2^20, rate 1/2, and at least 35%
GF65536-symbol distance with 40 bits of whole-code setup-failure margin.
Each 128-bit record carries eight SIMD field lanes under the same scalar code. Setup is fixed
before encoding and excluded from timing. The terminal state is unrestricted.

Commit `412e1102` preserves the [40% checkpoint](HD_OPTIMIZED_CHECKPOINT.md)
and the preceding characteristic-two research. No production or paper profile
changes in this workstream. Lowering the claimed distance alone does not
reduce work; the experiments below change the outer or inner maps.

At N=2^21, the forbidden output weights are 0 through 734003. Excluding these
weights gives minimum distance at least 734004, which exceeds 0.35N.
The screen and certificate APIs now accept this integer boundary explicitly;
their defaults remain the old 40% boundary 838860.

## Candidate family

Use independently randomized RS[L,L/2] outer rows, for L in {128,256}.
Each output coordinate receives an independent uniform nonzero GF65536
scalar. Each of the L regions independently permutes its N/L row entries.
The route distribution remains the existing uniform regional distribution.

The inner processes 64-symbol packets with a 16-symbol state Q, initially zero.
Choose h in {1,2,4}. The emission map A is the direct sum of 16/h copies of
RS[4h,h], in the normalized additive-polynomial basis. Set C=A^T. Each packet
has fresh independent nonzero diagonal entries in D1 and D2, sampled at setup.
The fixed matrix H interpolates the 16 input nodes 0,...,15 at nodes 16,...,31.
Compute

    Y = X + A Q,
    Q_next = D2 H D1 (Q + C X).

This is the [common-refresh construction](HD_COMMON_REFRESH.md). It omits
the random dot product of the 40% winner. For h=4, this is the existing local
proof, with a smaller forbidden-weight threshold and optionally smaller outer.
For h=2, each local map is RS[8,2]; one emission FFT level and one feedback
FFT level disappear. For h=1, each local map is repetition over four positions.
The latter saves little time and has no complete certificate.

All these local maps satisfy CA=0: products of two basis polynomials have
degree at most 2h-2, below 4h-1. Thus I+AC is invertible. The common-refresh
return analysis and diagonal density bound apply.

For h=2, the return bound admits a simple refinement. Within a local block,
the columns of C are (1,x), for x in {0,...,7}, and A=C^T. For a coordinate
set D, put c=|D| mod 2 and tau=sum_{x in D} x. The determinant lemma reduces
the nullity of (I+AC)[D,D] to that of

    I + C_D A_D = [[1+c, tau], [tau, 1+tau^2]].

Its determinant is 1 when |D| is even, and tau^2 when |D| is odd. In the
singular case its rank is one. Thus the maximum principal nullity is zero
at even sizes and at most one at odd sizes, rather than the generic cap two.
The exhaustive scan of all 256 subsets agrees. Independent direct ranks also
check every input-support/output-zero-set pair for the h=1 and h=2 maps.

This refinement improves the 2048-row single-packet regional screen from
about -37 to -2018 log2 units. These are floating-point proposal bounds,
not setup-failure margins. The gain comes from the proof; the native map
and its timing are unchanged.

The exact operators are in `theory/hd_common_refresh_packets.py`.
Floating-point screens propose witnesses; only directed replay and complete
coverage can certify distance. Native setup uses reproducible MT seeds for
benchmarking, not as a substitute for the ideal setup distribution in the proof.

## Performance results

All timings below are serial 501-call medians on Peach, pinned to CPU 2,
with precomputed setup and streaming output. No benchmark runs concurrently.

| Construction | Full K20 time | Proof status at35% |
|---|---:|---|
| Preserved40% winner, mode45 |5.459 /5.464 ms|Complete40% certificate,56.659 bits|
| Common refresh, RS256 outer, h=4 |5.156 ms|Local proof; full coverage pending|
| Common refresh, RS128 outer, h=4 |5.168 /5.259 ms|Complete 256/384-bit certificate,64.7737 bits|
| Common refresh, RS256 outer, h=2 |4.750 /4.739 ms|Complete 256/384-bit certificate, 76.7730 bits|
| Common refresh, RS128 outer, h=2 |4.769 ms|No complete certificate yet|
| Common refresh, RS256 outer, h=1 |4.709 ms|Failed bound at2048 active rows|

The certified two-check candidate removes about 13% of encoding time, or
improves throughput by about 15%, at the lower distance target. The smaller
outer has not established a speed advantage. Separate phase measurements for
h=2 give 2.455 ms outer and 2.305 ms routed inner for RS256, versus 2.520 and 2.286 ms
for RS128. These are separately timed phases, not an exact decomposition of
the full run. They do not support attributing the smaller outer's cost to routing.

Reducing next-region prefetching is counterproductive. Serial 501-call medians
with 16, 8, 0, 8, and 16 prefetched cache lines per packet were respectively
4.801, 8.336, 9.302, 8.343, and 4.734 ms. All checksums match. The retained
kernel prefetches all 16 lines; the temporary compile-time tuning switch was
removed after this test.

The initial RS128/h4 cached-output control took 5.578 ms, versus 5.192 ms with
streaming stores. That exploratory pair used a different payload seed; subsequent
runs use exactly the preserved control's payload generator.

## Validation and next decisions

The new native driver is `spin/experiments/char2_field_spin/hd35_forward.cpp`.
Its small-size checks cover both outer sizes and all three local ranks, with
independent dense polynomial arithmetic, alignment guards, and exact aliasing.
Release and ASAN/UBSAN pass. K20 independent dense checks pass for h=2 and 4
with both outer sizes. The observed K20 checksums are:

| Outer / local rank | Checksum |
|---|---:|
|128 /4|14635816923942864182|
|256 /4|15132924349672518052|
|128 /2|2199020749081728522|
|256 /2|986088359013200724|

The full Python regression run passes 456 tests. These include ten
35% tests covering exact local operators, self-orthogonality, counting
normalization, directed distance shifts, and the single-profile global bound.
All 45 pre-existing native variants still pass Release and ASAN/UBSAN after
adding the separately instantiated smaller local maps.
The rebuilt 40% winner takes 5.417 ms, versus 5.415 ms for its frozen binary
in the same serial run. The additional templates do not introduce a measured
regression in that control.

The RS128/h4 checkpoint has 44 pinned witnesses in
`theory/hd35_common128_witnesses.py`. A complete 256-bit replay covers all 16384
active-row counts. It gives ideal union log2 bound -76.31543256667072 and
diagonal density loss 11.54164838312470, for final log2 bound -64.77378418354602.
The probability is over the ideal setup distribution for the full code,
not over a sampled message. Independent 384-bit replay confirms the same
coverage and margin; the final enclosing ball has radius below 7e-108.

The selected RS256/h2 checkpoint has 39 pinned rational witnesses in
`theory/hd35_common256_witnesses.py`. They use exact regional coefficients
of single-packet profile moments. The complete 256-bit replay covers all
8192 active-row counts, with no rejected witnesses or missing classes:

    ideal union log2 upper bound = -88.31470859712889
    diagonal density loss       =  11.54164838312470
    final union log2 upper bound = -76.77306021400419

Therefore the ideal sampled code has minimum distance at least 734004,
except with probability at most 2^-76.77306021400418. This is a statement
about all nonzero messages simultaneously. A separate 384-bit replay confirms
all 8192 classes and the same margin, with final ball radius below 7e-109.
The proof covers the scalar GF65536 code; using the same map
on eight SIMD lanes preserves the same distance guarantee without a union
bound over lanes.

The native executable SHA256 is
`81bc8f93fe7f1dc1380a4c91ab9d5679975e32b9b6a8ba668d3a60a85b722785`.
The frozen 40% control remains
`69329fefb8f3e4020488411adfcaa8c0aa694eaf43fdc8b29aaa9e88c57e46c3`.

Preserve the RS128/h4 certificate as a conservative checkpoint alongside
RS256/h2, the measured 4.74-ms winner at 35%. The one-check version saves only
about 0.04 ms. Its earlier conservative screen failed by more than 22000 log2
units; this is not a counterexample or a limit for tighter proofs.

The first h2 audit transferred rational search hints from the h4 certificate.
It recomputed all operators but left middle-class coverage gaps. That audit
and a paired-packet repair were stopped after the exact principal-nullity
bound made a cheaper proof viable. They are not complete certificates.

The successful replacement search uses exact regional coefficients of single-packet
profile moments. Apply fractional powers before averaging profiles. Retaining
the regional coefficient avoids the loose global marker bound in the middle
classes. This search checks every returned class with directed arithmetic;
no floating-point screen contributes to certified coverage.

## Reproduction and scope

From the theory directory, replay the selected certificate with
`python hd35_common256_witnesses.py --precision 384 --workers 4`.
The complete local regression command is
`python -m unittest discover -p 'hd*_tests.py'`.

Build the existing research CMake target `spin_char2_hd35_forward` in Release.
Run `--check` for small-size correctness and `--large-check` for the K20 dense
oracle. On the benchmark host, pin one process to CPU 2 and run
`--bench 501 full 20 nt 256 2`. Do not run benchmarks concurrently. Setup
and input generation are outside the timed interval; outer encoding,
routing, inner encoding, and output stores are inside it.

This is a certified research checkpoint, not a production promotion or a
new timing for the 40% code. Further meaningful gains need a new change to
the maps or data layout; the shorter outer and reduced prefetch did not help.
The result does not establish a global performance limit.
