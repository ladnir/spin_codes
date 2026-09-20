# Riffle TransposeBitShuffle-RandomStepConv g

This construction uses the same random rate-half outer blocks and one-lap
RandomStepConv inner encoder as Riffle
TransposePacketShuffle-RandomStepConv g. Its permutation law is different.

Place the \(L\) length-\(B\) outer outputs in an \(L\)-by-\(B\) bit matrix and
transpose it. Setup independently permutes all \(L\) individual bits in each
of the \(B\) transposed rows. The encoder then partitions every permuted row
into consecutive \(g\)-bit convolution inputs. There is no additional packet
permutation.

The bit-shuffle mechanism randomizes collisions between active outer blocks
independently in every row. It performs \(B\) permutations of \(L\) items,
whereas the packet construction performs \(B\) permutations of \(L/g\)
items. For one active outer block the two constructions induce the same input
support distribution. They differ when two or more outer blocks are active:
the packet construction reuses fixed block groups, while this construction
forms new random collisions in every row.

For \(a\) active outer blocks, their candidate bits in one row form a uniform
\(a\)-subset of the \(L\) bit positions. Let \(M_r(z)\) be the tilted inner
transition of a packet containing \(r\) candidate bits. The exact averaged
row transition is

\[
  \frac{[x^a]\left(\sum_{r=0}^g {g\choose r}x^rM_r(z)\right)^{L/g}}
       {{L\choose a}}.
\]

The full checker computes these matrix-polynomial coefficients with
outward-rounded mantissa-plus-exponent arithmetic. High occupations use the
reversed polynomial, so their cost depends on \(L-a\), without changing the
ordered matrix products.

This random ensemble now has complete distance certificates for \(g=8\) at
the three cells in `TRADEOFF_G4_G8.md`. The \(g=4\) cells are certified by the
stronger packed-support upper bound and independently reproduced by the exact
bit-shuffle diagnostic.
