# Riffle FieldPair512-P2-ScalarStepConv g=4 sigma=18

Status: ACTIVE_ENTROPY_OPTIMIZATION

This candidate minimizes setup entropy while preserving the probability laws
used by the current first-moment proof. It changes both random linear-map
families and reduces the inner state from 20 to 18 bits.

Each of the 4096 data groups receives 512 independent setup bits, represented
as two elements of \(\mathbb F_{2^{256}}\). Each inner step receives 22 fresh
setup bits, represented as one element of \(\mathbb F_{2^{22}}\). Setup also
samples the inherited global packet permutation and two local parity-block
permutations. All sampled objects are frozen before encoding.

The resulting setup entropy is approximately 22,840,639 bits, or 2.723 MiB.
The previous dense-matrix description used approximately 848,105,407 bits.

The sigma-18 exact-saddle diagnostic sums to `-43.5986` over supports 50
through 3200. The omitted tails are far below this range in the sampled
Chernoff scan. A complete tail cover and outward-rounded evaluation remain
necessary before claiming the required `2^-40` first moment.

See `CONSTRUCTION.md`, `GOAL_01_ENTROPY_LEDGER.md`, and `REPORT.md`.
