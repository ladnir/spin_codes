# Riffle transpose co-design performance frontier

This note records implementation probes for the `K=2^20`, `N=2^21`
GF(2^128)-symbol transpose.  These are construction candidates, not additions
to the proved theorem.  Timings use one pinned CPU core and the synthesized
913-XOR extended-BCH transpose.

## Fixed-half adjacent recurrence

Split every 128-coordinate BCH word into fixed 64-coordinate halves
`E=(E_L,E_R)`.  Both halves of the committed cyclic encoder have binary rank
64.  For a permutation `pi` of the `B=N/64` physical symbol groups, define the
forward recursive chain in the order

```
pi(0), pi(1), ..., pi(B-1).
```

The transpose is evaluated in place from right to left:

```
state = 0
for i = B-1 ... 0:
    result = E^T(input[pi(i)], state)
    input[pi(i)] = result
    state = result
```

The 64-symbol state stays in a hot 1 KiB temporary.  Every completed result is
written immediately to its final slot in the caller's input buffer.  There is
no N-symbol permutation workspace and no inverse permutation.  The physical
buffer remains in outer-coordinate order, so outer compression is sequential.

The fastest random-chain implementation precomputes `pi(i)=Feistel(i)` for
the 32768 group indices (128 KiB) and prefetches all sixteen cache lines of
group `pi(i-2)` while evaluating the current BCH call.

## Measured frontier

Representative medians, in milliseconds:

| Candidate | Time |
|---|---:|
| Identity-order adjacent chain | 10.83 |
| Random 1 KiB group chain, no software prefetch | 13.78--15.4 |
| Random 1 KiB group chain, full-line lookahead 2 | 11.17 |
| Fixed split plus coordinate-level precomputed-Feistel pull | 18.85--20.4 |
| Existing random-split, coordinate-level construction | 25--31 |

The full-line, two-call lookahead leaves the strongest tested group-randomized
candidate only about 0.34 ms (3 percent) above the identity-order performance
ceiling in the clean chunk-sweep run.

Keeping adjacent chain nodes in randomly permuted physical chunks gives the
following performance/randomization curve:

| Nodes per permuted chunk | Bytes per chunk | Time (ms) |
|---:|---:|---:|
| 1 | 1 KiB | 13.78 |
| 2 | 2 KiB | 13.41 |
| 4 | 4 KiB | 12.57 |
| 8 | 8 KiB | 11.99 |
| 16 | 16 KiB | 11.36 |
| 32 | 32 KiB | 10.88 |
| 64 | 64 KiB | 11.48 |
| 128 | 128 KiB | 11.09 |
| 256 | 256 KiB | 11.12 |
| Identity | 32 MiB | 10.83 |

The staged random-chain control copied two 1 KiB groups into a 2 KiB local
word before every BCH call and was substantially slower.  The winning kernel
feeds the physical current group and hot state directly into a generated
two-pointer transpose circuit.  This confirms that the permutation should be
fused into the recurrence's address schedule rather than implemented as a
value-moving push or pull pass.

## Small-inner co-design probe

The 64-symbol group is an implementation choice, not a requirement of the
outer EBCH code.  Two ways of reducing the permutation atom were tested on the
same `N=2^21` block buffer:

1. retain the 64-wide systematic BCH state and assemble each step from
   independently permuted packets;
2. use a genuinely smaller systematic inner code and state while retaining
   the `[128,64,22]` outer.

The packetized 64-wide circuit is correct and keeps the data buffer in place,
but scattered pointer pressure is expensive.  A sequential isolated sweep of
the exact 930-XOR systematic step measured:

| Packet width | Packets per 64-wide step | Median (ms) |
|---:|---:|---:|
| 64 | 1 | 7.83 |
| 32 | 2 | 10.08 |
| 16 | 4 | 11.90 |
| 8 | 8 | 14.20 |

These rows record the original compiler snapshot. A 2026-08-13 rebuild of the same source with the matched `/O2` configuration measured 7.672 ms for packet width eight. The rebuilt executable passed the same direct reference check. Use the rebuilt result for current comparisons.

The historical small-inner probe uses a systematic `[16,8,4]` inner code.  Its right
projection has rank eight, its distance is exhaustively checked over all 255
nonzero messages, and a generated/common-subexpression step needs 13 XORs for
`P^T`, eight current merges, and seven accumulator XORs: 28 field XORs per
eight-symbol group.  The Feistel-derived group order is precomputed; each hot
access is one contiguous 128-byte group.  The same isolated run measured
`4.58 ms` (minimum `3.89 ms`), versus `7.83 ms` for the 64-wide systematic
chain. This was once a prospective co-design direction. It is now rejected.
The deployment construction keeps both
BCH codes at `[128,64,22]`; only the permutation packet width may change.

The implementation probes are in
`libOTe/Tools/RiffleCode/SystematicCode16x8.h`,
`libOTe/Tools/RiffleCode/RiffleSmallInnerChain.h`, and
`libOTe_Tests/RifflePacketInnerChain_Bench.cpp` in the
`codex/permute-conv-code` libOTe worktree.

## Optimized packet-width-eight checkpoint

The consolidated implementation and construction lessons are recorded in `explorations/riffle_packet8_implementation_lessons.md`.

The packet-width-eight implementation now specializes only the generated hot path. Both local codes remain `[128,64,22]`.

The inner circuit receives eight explicit packet pointers and prefetches one node ahead. The outer circuit reads unpunctured words in place and writes parity-adjusted outputs directly.

A 101-trial standalone run measured 7.765 ms for the inner, 4.017 ms for the outer, and 11.970 ms end to end. The matched Tungsten time was 9.608 ms.

An interleaved comparison alternated the old and optimized Riffle paths. It checked exact equality after each map. The optimized path used 80.5% of the old path's time.

The optimized packet-width-eight path is only 0.885 ms slower than the rebuilt 11.085 ms packet-width-64 historical encoder. The remaining packet-width cost is therefore modest.

The next locality pass rejected immediate outer-block fusion. Outer blocks become complete in pseudorandom order, so the fused path measured 15.566 ms against 11.793 ms for separate phases.

Three exact implementation changes were retained:

1. Write-intent prefetch requests ownership of packets before their read-and-overwrite step.
2. The graph correction is folded into each BCH addend, which removes a second pass over the 16 MiB output.
3. A 1,024-seed exact search reduced the packet circuit from 803 to 799 transpose XORs.

The final 101-trial standalone run measured 7.768 ms for the inner, 3.636 ms for the outer, and 11.670 ms end to end. Tungsten measured 9.767 ms. An interleaved comparison against the preceding optimized Riffle path measured 11.543 ms versus 12.063 ms, a 4.3% reduction.

The final packet-width-eight path is 0.585 ms above the rebuilt 11.085 ms packet-width-64 checkpoint.

## Proof obligations exposed by the fast design

1. Replace the fresh per-block random 64-of-128 split by the fixed-half joint
   kernel.  Both projections are invertible, so both marginals are uniform
   nonzero and the conditional turnoff-support law remains uniform.  The joint
   `(output weight,state weight)` distribution must be bounded or certified.
2. Replace the coordinate interleaver by a permutation of 64-symbol groups (or
   select a larger chunk point from the table).  The outer placement argument
   must work with group-occupancy profiles rather than individual-coordinate
   hypergeometric placement.
3. Decide whether the proof retains the per-block transitive mixer `A_i`.
   The timings above omit that mixer.  Retaining it gives the cleaner immediate
   local proof; removing it requires the no-Singer/identity-kernel proof lane.
4. If a mixer remains necessary, benchmark it fused immediately after the BCH
   transpose and before the result is written to its physical input slot.
