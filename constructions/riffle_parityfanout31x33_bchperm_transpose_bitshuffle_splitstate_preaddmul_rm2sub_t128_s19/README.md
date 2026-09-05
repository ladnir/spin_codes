# Structured SPIN (B=256, t=128, s=19)

**Status:** frozen main implementation candidate. See `MAIN_CODE_FREEZE.md`.

Legacy exploration name: `Riffle ParityFanout-31x33 s=19`.

This candidate targets relative distance 0.11 at message length `2^20`
and output length `2^21`. It retains the BCH-per-block, coordinate
permutation, bit transpose, region permutation, and RM2Sub inner structure.

The construction change is one independently sampled sparse rank-one map per
outer block. Setup samples disjoint coordinate sets

\[
|S|=31,\qquad |T|=33.
\]

The forward outer map is

\[
y_t=x_t+\sum_{s\in S}x_s\quad(t\in T),
\qquad y_i=x_i\quad(i\notin T).
\]

It is invertible because `S` and `T` are disjoint. Both set sizes are odd,
so the map breaks the even-weight constraint and moves the all-one word into
the ordinary randomized-shell spectrum. There is no exceptional all-one
proof class.

The transposed encoder first XORs the 33 target PCG blocks into one temporary
and then XORs that temporary into the 31 source blocks. This costs exactly 63
block XORs per outer block, or 516,096 block XORs per full call. The paired
implementation evaluates two independent fanouts in lockstep before the
existing paired BCH transpose.

Current results:

- nearest-binary64 first-moment margin at distance 0.11: 55.864642 bits;
- 21-trial Peach end-to-end median: 10.823872 ms;
- dense-inner and independent staged end-to-end correctness: PASS.

The proof result still inherits the modeled outer spectrum and nearest-
binary64 arithmetic. The same fine tilt grid gives 66.270087 bits at distance
0.09 and 61.069584 bits at distance 0.10. Outward-rounded numerical hardening remains the final
certificate gate.
