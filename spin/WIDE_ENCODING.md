IMT wide forward encoding composes the outer tile permutation with the final
gather during workspace construction. Encoding writes BCH outputs directly
into their natural positions, then reads through the composed index table.
The workspace retains one 32-bit index per encoded coordinate and no longer
needs a temporary permutation tile. Both packed-24 and full-32 setup indices
are supported, including a partial final tile.

On GNU-compatible builds with AVX-512 support enabled, 256-bit encoding can use
explicit register schedules for BCH and the T128S19 feedback map. The helpers
declare their AVX-512F/VL target and use all 32 vector registers. Workspace
construction checks CPU support once. An explicit AVX2 backend, an unsupported
CPU, or an AVX-512-disabled build uses the original C++ circuits. The 512-bit
encoder and the T64S12R2 feedback circuit also retain their existing arithmetic.

The generated assembly lives in ordinary functions, with private aligned
scratch. Input and output must be disjoint: a helper may read an output it
wrote earlier to recover a dependency. The internal caller satisfies that
requirement. Encoding allocates no memory; routing preparation happens when
the workspace is created or prepared for a new code.

These changes preserve the encoded map, descriptors, seeds, and wire format.
They do not change the code construction or its distance claims.

The `Code` byte-view and typed forward interfaces take an optional per-call
`OutputStores` argument. The default, `Streaming`, requests non-temporal output
stores where supported by the backend and alignment. `Cached` forces ordinary
output stores. Both choices work at the public 16-byte alignment minimum, with
a cached fallback wherever streaming stores cannot be used. The selection is
made outside the hot loop. Reuse the same workspace while changing the policy;
it does not change scratch allocation, the record layout, or the code descriptor.
The result is ready to read on return, including after streaming stores.

Choose the policy according to when the consumer will read the output. A large
matrix written across many encodings before hashing has different cache behavior
from one codeword repeatedly encoded and read immediately. Alignment alone does
not express this distinction. Transpose and bit-packed forward keep their
existing behavior. The old packet profile's large direct route can use output
as intermediate storage; the explicit preference applies to that output too.

Packet SPIN also supports 256- and 512-bit forward records through `Code`.
This includes the fixed-K16 `PacketRsT64S15K16` profile described in
[PAIRED15.md](PAIRED15.md). Its outer loop shares physical input loads across
lanes and uses its own exact forward circuit and generators. The following
prefetch and generator details describe the preceding s20 profile.
The packet kernels retain their fixed 128-bit GFNI schedules. Outer loads read
each lane from interleaved records directly. Inner evaluation advances every lane
over the same 64-coordinate block, then combines their outputs into complete
cache lines. This uses bounded local scratch instead of gathering the entire
message and interleaving the entire codeword in separate passes. There are no
FLOCK-specific dimensions or layouts in these kernels.

The 256-bit packet outer kernel prefetches the next 256-coordinate input group
into L2. This overlaps cold input reads with GFNI work without displacing the
current group's L1 working set. Prefetch addresses stay within the input;
the 512-bit path keeps its existing schedule.

Packet wide workspaces own one outer scratch plane per 128-bit lane. Compatible
workspaces can be rebound to different seeds without allocation. The portable
path implements the same layout, and output buffers retain the 16-byte alignment
contract. Wide transpose remains unsupported. Regenerate the packet forward
sources with `python3 tools/generate_packet_forward.py`; pass `--check` to check
them without editing files.

Regenerate the two schedules from the checked-in wide circuit:

```sh
python3 tools/generate_bch_register.py
python3 tools/generate_feedback_register.py
```

Both generators track symbolic input masks through every emitted instruction
and check the resulting linear map. The API and known-answer tests compare
wide lanes against independent paths. `spin_wide_routing` additionally checks
packed and 32-bit routes, partial tiles, explicit AVX2 and automatic dispatch,
and dense/zero/sparse inputs with reused scratch. Run `ctest` serially in
builds with `SPIN_ENABLE_AVX512` both enabled and disabled.
The API suites also compare both output policies and the default at 64-byte and
16-byte-only alignment, alternate policies on dirty scratch, and reject invalid
policy values before writing any output.
