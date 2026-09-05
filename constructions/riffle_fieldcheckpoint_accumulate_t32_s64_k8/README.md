# Riffle FieldCheckpointAccumulate t=32 s=64 K=8

This variant refreshes the 64-bit accumulator state every 256 inner-output
bits. It retains the outer code, coordinate permutations, transpose, region
permutations, and 32-bit implementation step of the K=32 parent candidate.

The shorter epoch removes the regular middle-density obstruction at 9%
relative distance. Under the modeled outer spectrum, the floating-point
first-moment calculation covers every regular active-block count. The low
range has 55.9507 bits of aggregate margin. The exact middle and high ranges
have more than 19000 bits of margin. An analytic dense-subspace bound takes
over at 7242 regular active blocks.

This result is not a full certificate. The modeled spectrum is conditional,
configurations containing all-one outer words remain partially open, and the
calculation is not outward rounded.

The transposed evaluator performs approximately four times as many field
multiplications as the K=32 parent. The accumulator XOR work is unchanged.

See `CONSTRUCTION.md` and `PROOF_STATUS.md` for the precise scope.
