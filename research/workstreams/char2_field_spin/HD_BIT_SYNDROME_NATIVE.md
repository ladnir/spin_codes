# Native cost of the four-XOR outer

The new [bit-syndrome outer](HD_BIT_SYNDROME_OUTER.md) runs at about 3.4 ms
in the first matched K20 comparison. It has no whole-code distance
certificate. The [reviewed column majorant](HD_BIT_SYNDROME_KERNEL.md) is
the next input to a full-state proof calculation, not an existing margin.

## Construction and implementation

Keep sixteen RS[16,8] constituents per 128-input group, their grouped
transpose, and independent nonzero GF256 diagonal randomizers on each
sixteen-symbol column. Replace the two-sum matrix H by H*, which adds the
four binary-index syndrome terms defined in the proposal. The existing
common-scalar normalization preserves each sampled diagonal-sandwiched map.

The alternating paired route, independent pair swaps, and h4 inner are
unchanged. In particular, the inner has its original single 16-symbol
state, not either two-state coupling proposal. The measured implementation
uses the [predecoded swap masks](HD_DECODED_ROUTE.md), unscaled 16-bit route
bases, and original outer load schedule.

The isolated sources are in `bit_syndrome_native/`. Flag
`HDGROUP_BIT_SYNDROMES=1` requires the experimental G16 mixer. The default
is zero. Four shared XORs modify the output corrections; the eight fixed
field products and sixteen explicit output stores remain unchanged.
The scalar oracle independently constructs the dense matrix using the
parity of the binary dot product of the row and column indices.

## Validation and timing

Both variants pass 8,192 SIMD mixer basis directions, sixteen random
columns, the involution check, all 30,720 input/state packet directions,
pair-wrap checks, small dense cases, and the K20 dense oracle. The new
variant also passes ASAN/UBSAN small-size and basis checks. Its K20 checksum
is 6721820652671321300; the old G16 control remains 5924045154917714530.

The following full-encoding medians use Peach CPU2, GCC15, 501 calls,
precomputed setup, 128-bit records, and streaming output. Both benchmark
locks are held and all runs are serial. Setup, allocation, and input
generation are excluded. The shared host is noisy.

| Round | G16 with predecoded route, ms | Bit-syndrome with predecoded route, ms |
|---:|---:|---:|
| 1 | 3.350616 | 3.335718 |
| 2, reversed order | 3.284903 | 3.410037 |
| 3 | 3.478614 | 3.621712 |

The medians of run medians are 3.350616 and 3.410037 ms. The changing
paired differences and final slower round do not support a precise
small-overhead claim. These results establish a useful cost candidate,
not a stable upper runtime bound or completion of the 3.4-ms proof goal.

Assembly retains 39 GFNI products per column: fifteen input-diagonal,
eight fixed, and sixteen output-diagonal products. It has no outlined
mixer/store helper and unchanged RS8 leaves. The compiler combines the
four extra XORs into one additional XOR and two ternary-XOR instructions.
The new register pressure spills only the uniform 0x0c multiplier constant,
not payload: one 64-byte store per coset and one reload per column. This
adds 2.25 MiB of logical L1 traffic at K20. No further native tuning is
needed before testing the proof.

## Binding and next step

Remote directory: `/tmp/hd-bit-syndrome-20261009-SA2rEj`.
The local build script, checks, timing log, and hash records bind this run.

| Binary | SHA256 |
|---|---|
| G16 control | `8065c7f4cc72d80029327ffc622f782bb411c1b57475bab8f9117a9a4fc1da74` |
| Bit-syndrome | `e37da61ac3d01e47bf3f9c89d1e47c5cc3b68def4a16b5908499dae868eb6fbc` |

The next gate is the full r=18, arity-eight, Fmax=1 class using the
original one-state regional bounds. It must retain every labelled source
set across all sixteen columns and include both quiet and nonquiet words.
A successful FLOAT gate would justify further coverage and directed
replay, but would not itself establish the requested whole-code margin.
The certified Toeplitz checkpoints remain unchanged. Nothing was promoted,
committed, pushed, or inserted into the paper.
