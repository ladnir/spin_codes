# BCH256 Fanout-56 performance proxy

This directory contains the only active Fanout-56 performance proxy. It is a
working copy of the frozen 256-by-8192 packed, tiled, fused transpose encoder.
The copy changes only schedule multiplicity and the number of row-local
fanout layers applied by the benchmark mode.

The remote build workspace supplies unchanged frozen dependencies, including
`RiffleRm2SubS19.h` and the BCH256 transpose headers.

Run `RiffleFieldCheckpoint_Bench` with the mode `fanout56-proxy`. The mode
checks staged zero-layer and 56-layer oracles before timing either path.

This proxy is not the exact BCH250-124 encoder from the finite 11% proof. It
must not be cited as an exact end-to-end implementation result.
