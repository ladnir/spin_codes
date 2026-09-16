# Follow-up optimization and inner-design direction

The additional experiments did not establish a repeatable improvement over the
[17.996 ms baseline](PERFORMANCE.md). The supported implementation, generator,
and manifest were restored byte-for-byte to that tested version. The certificate
inputs are unchanged. No new inner design was adopted.

## What was tested

All measurements used the same Peach host and CPU-15 affinity, sequentially.
The benchmark lock rejected one attempted launch; no timing was collected from that attempt.
Most screening runs used 51 trials. Instrumented runs additionally measured phases;
these are diagnostics, not replacements for the uninstrumented performance report.

| Candidate | Screening time at K=2^20 | Outcome |
|---|---:|---|
| Interleaved two-row outer, instrumented | 16.988 ms | Initial apparent win; not confirmed in the final build |
| Interleaved four-row AVX-512 outer, instrumented | 17.693 ms | Slower than the two-row candidate |
| Specialized inner expansion, instrumented | 16.922 ms | No substantial gain over the interleaved control |
| Write-prefetch, 32–128 entries ahead, instrumented | 18.457–18.815 ms | Slower |
| One-byte bucket IDs with per-bucket cursors, instrumented | 18.244–18.350 ms | Saved 8 MiB of setup but slowed the inner/routing stage |
| Outer gathers directly from buckets, without tile scatter | 24.409–26.623 ms | Slower, with or without lookahead prefetch |

The specialized expansion reconstructed the selected quadratic polynomial and
shared its XOR expressions. Symbolic evaluation matched all 128 selected columns.
It used 296 XORs for the state expansion, but fewer XORs did not translate into
a meaningful end-to-end improvement. Register pressure and scheduling remain
possible explanations; this experiment does not isolate their effects.

The final uninstrumented interleaved candidate and original baseline were tested
in alternating runs, each with 101 trials:

| Version | Three run medians | Median of medians |
|---|---|---:|
| Original baseline | 17.615, 17.651, 17.526 ms | 17.615 ms |
| Interleaved candidate | 17.979, 17.856, 18.093 ms | 17.979 ms |

Consequently, the early approximately 17 ms result is not a claimed optimization.
Differences between instrumented and final builds also caution against promoting
phase-screening results without checking the actual release binary.
Both candidates passed the full dense-reference and regression tests.
The original build's retained sanitizer results still apply to the restored sources.

[optimization_records](optimization_records/) retains compact diagnostic measurements,
the candidate's test log, and final binary hashes. Failed experimental kernels
are not supported configurations and were removed from the working source.

## Where to investigate the inner

The interleaved candidate's phase samples were approximately 7.8 ms for inner
computation plus initial bucket writes, 6.3 ms for tile routing, and 2.9 ms for
the outer. The first figure is not pure inner arithmetic. Removing all inner
arithmetic would not remove its memory traffic.

A change of state coordinates is a useful first design experiment because it
can preserve the complete encoder, rather than just its image spectrum.
Write the forward inner as

    y_i = x_i + A q_i,
    q_(i+1) = M_i q_i + B x_i,      q_0 = 0.

Here A expands the 19-bit state, B maps an input block into that state, and M_i
is multiplication by the selected nonzero field coefficient. In the current
coordinates B=A^T. For any invertible binary matrix U, substitute q_i=U z_i:

    A' = A U,
    B' = U^(-1) B,
    M_i' = U^(-1) M_i U.

These three simultaneous changes preserve every output for every input and
every existing setup. They also preserve the image of A, the kernel of B,
and B'A'=0. Thus an exactly equivalent implementation can use the existing
certificates without a new distance estimate. In general B' is not A'^T;
replacing it by A'^T would not be the same change of coordinates.

The proposed search should score expansion, feedback, and field-transition
costs jointly. Candidate matrices can be precomputed during setup. Compare
each candidate against the original dense oracle and measure the complete
transposed encoder, not just its gate count.

If equivalent representations give little improvement, the next tier is a
joint search over RM2Sub maps and (t,s), scored by both runtime and certified
margin. That changes the construction and needs fresh spectrum/margin checks.
This note proposes that search; it does not claim a faster inner exists.
