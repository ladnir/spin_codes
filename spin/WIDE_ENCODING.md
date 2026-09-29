Wide forward encoding now composes the outer tile permutation with the final
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
