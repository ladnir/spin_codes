# Frozen paired-s15 K16 profile

This private module implements the selected mode-52 construction from
`spin/experiments/k16_codesign_100us`. It accepts only K = 65,536 and has
rate 1/2. The shipped library does not include or link research sources.

The setup is seed-compatible with the frozen harness. Routing uses
`Words(routeSeed)`, the native GF(2^16) symbol maps reset the stream to
`Words(routeSeed ^ 0x75a1dc09)`, and the GL15 updates reset it to
`Words(innerSeed ^ 0x3f625a92)`. The harness uses the same seed for both
arguments; the library also supports distinct route and inner seeds.
The rejected zero scalars and singular GL15 draws consume the same words.

`Paired15Fast.cpp` preserves the selected transposed inner loop and uses
`Paired15TransposeOuter.h` for the selected outer. Forward encoding applies
the adjoint RS circuit, the adjoint of each sampled symbol map, and the
transpose of each stored inner update. The fixed SIMD state permutation is
`(0,1,2,7,3,4,5,6,8,15,10,12,11,13,14,9)`; slot 9 is always zero.

Forward accepts interleaved 128-, 256-, and 512-bit records. Wide input loads
read every record once and share their register rearrangement across lanes.
Scratch has one contiguous routing plane per 128-bit lane. The inner writes
wide results after each 64-coordinate step, using bounded stack storage.
Encoding allocates nothing. Scratch never aliases either public buffer, and
the whole message is consumed before any output is written.
The trailing `stream` argument to `forwardFast` defaults to true and enables
non-temporal final output stores when 64-byte aligned, followed by one store
fence. Passing false forces cached output at every supported alignment. The
dispatch occurs once per call; there is no per-store branch. Scratch and
transposed encoding retain cached stores. `forwardFastCached` remains an
internal comparison entry point with the same always-cached final output.

`import_frozen.py` refreshes the retained circuits from the frozen research
files. `generate_forward_outer.py` transposes the explicit RS circuit and
inverts its byte permutations. These are development tools; neither runs
when building the library. The scalar path independently constructs the RS
matrix by interpolation and uses the literal compact s15 coordinates.
Both generators accept `--check` to verify checked-in output without writing.
Generated sources record SHA256 hashes of their immediate source circuits,
with newlines normalized to LF for reproducibility across Git checkouts.
