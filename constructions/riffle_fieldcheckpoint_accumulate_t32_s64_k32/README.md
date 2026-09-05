# Riffle FieldCheckpointAccumulate t=32 s=64 K=32

This candidate replaces a dense randomized inner map at every step with 64
parallel accumulators and one randomized field checkpoint every 1024 output
bits. The checkpoint uses one nonzero element of \(\mathrm{GF}(2^{64})\).
The deployed transposed evaluator therefore performs ordinary field
multiplication only once per checkpoint, not once per output step.

The inner admits an exact two-state distance analysis. At every checkpoint,
the state is either zero or uniform over the nonzero 64-bit states. This
closure holds while retaining a generating variable for output weight.

The one-active-outer-block calculation has 81.2667 bits of margin at 9%
relative distance. Its 40-bit frontier is approximately 14.049%. A later
full-occupation calculation found that this sparse result is misleading:
near 1024 active outer blocks, the 9% first-moment bound misses by about
16039 bits. The loss is too large to come from the factor-two spectrum
envelope.

See `CONSTRUCTION.md` for the map, `PROOF_STATUS.md` for the exact transfer
formulas, and `CANCELLATION_ANALYSIS.md` for the weight-two audit.

At 6% relative distance, the regular-support envelope now covers every
active-block count. The low range has at least 55.9008 bits of margin, the
exact middle range has 834.4399 bits, and an analytic dense-subspace bound
takes over at 5467 regular active blocks. The layer containing exactly one
all-one outer word also closes through the middle range with 1161.5591 bits.

These are ordinary floating-point diagnostics for a modeled outer spectrum.
The remaining full-certificate work is to cover two or more all-one outer
words below the dense threshold and then repeat the retained calculation
with outward-rounded arithmetic. See `proof/FULL_CERTIFICATE_PLAN.md`.
