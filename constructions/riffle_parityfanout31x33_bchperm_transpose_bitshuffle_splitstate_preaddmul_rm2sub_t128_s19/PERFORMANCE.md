# Performance

The transposed ParityFanout-31x33 implementation evaluates the two outer
blocks consumed by one BCH pair in lockstep. This preserves the exact matrix
while exposing independent loads and XOR chains to the processor.

On `peach48.devcore4.com`, pinned to CPU 15, the 21-trial end-to-end median
is 10.823872 ms. The prior integrated state-size-16 baseline is 10.821225 ms.
The observed difference is 0.002647 ms.

The correctness harness checks the optimized five-epoch RM2Sub transpose
against an independent dense implementation. It then checks the complete
fused encoder against a materialized inner, independent inner-to-outer route,
independent ParityFanout transpose, and BCH transpose.

The full samples and source hashes are recorded in
`benchmark/parityfanout31x33_s19_integrated_peach_7950x.json`.
