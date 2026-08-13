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

The more promising probe uses a systematic `[16,8,4]` inner code.  Its right
projection has rank eight, its distance is exhaustively checked over all 255
nonzero messages, and a generated/common-subexpression step needs 13 XORs for
`P^T`, eight current merges, and seven accumulator XORs: 28 field XORs per
eight-symbol group.  The Feistel-derived group order is precomputed; each hot
access is one contiguous 128-byte group.  The same isolated run measured
`4.58 ms` (minimum `3.89 ms`), versus `7.83 ms` for the 64-wide systematic
chain.  This is not yet an end-to-end Riffle timing and the small inner has not
yet been certified by the distance ledger.  It is the current co-design
candidate to prove first and then grow to widths 16 or 32 if the proof margin
allows.

The implementation probes are in
`libOTe/Tools/RiffleCode/SystematicCode16x8.h`,
`libOTe/Tools/RiffleCode/RiffleSmallInnerChain.h`, and
`libOTe_Tests/RifflePacketInnerChain_Bench.cpp` in the
`codex/permute-conv-code` libOTe worktree.

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
