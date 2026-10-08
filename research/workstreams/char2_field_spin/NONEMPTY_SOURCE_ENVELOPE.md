# A valid but unsuccessful nonempty-source envelope

The exact source program requires every selected outer row to have a nonzero systematic part.
Dropping that condition makes some middle-class bounds much weaker.
This note records a cheaper envelope that retains the condition, but loses too much information elsewhere.
It is not the next certificate method.

Fix the fresh GF16 packet construction in [RANDOM_PACKET_GATE.md](RANDOM_PACKET_GATE.md).
Let a be the number of selected source rows, and let R_h(z) be the normalized regional moment matrix for h nonzero input positions.
The state indices are zero and nonzero; all entries are nonnegative.
For positive v and u, define entrywise

    M(v)_(i,j) = max_(0<=h<=a) R_h(z)_(i,j) / v^h.

Every systematic region matrix is bounded by v^h M(v), and every parity region matrix by u^h M(u).
Matrix multiplication preserves these entrywise inequalities.
Summing over all source values then bounds the exact source moment by

    ((1+15v)^k - 1)^a
      * ((1+15u)/16)^((b-k)a)
      * e_0 M(v)^k M(u)^(b-k) 1.

The subtraction of one occurs separately for each row.
It excludes the zero systematic vector, not merely a zero total systematic weight.
The parity factor uses the independent uniform parity image of each fixed nonzero source row under its random outer Toeplitz map.
The complete class bound additionally includes binom(L,a), division by 15 for scalar multiples, and the output-weight Chernoff factor z^-e.

This inequality holds for arbitrary nonnegative regional matrices.
The implementation evaluates it in floating logarithms; it is not a directed numerical certificate.
Four tests check the row enumerator, entrywise domination, and comparison with the exact positive source program, including arbitrary positive matrices.

## Screening result

For the rate-5/8 [64,40], width16/state5 backup at K20, optimize a common v=u at the exact-source witnesses:

| Active rows | z | Exact source margin | Envelope margin |
|---|---:|---:|---:|
| 128 | 0.9703716207340908 | 208.435 bits | -3913.371 bits |
| 256 | 0.9413912472726532 | 478.169 bits | -7596.050 bits |

Negative margins describe the evaluated upper bound, not construction failure.
Allowing the worst regional weight independently for each matrix entry discards crucial weight dependence.
Retaining nonempty source rows alone does not compensate for that loss.
Do not replace the exact source program with this envelope or treat it as a promising route to full coverage.

The backup also has a substantial arithmetic disadvantage under the current circuits.
Splitting a 24-by-40 Toeplitz parity map into 16/8-sized blocks gives
2*81+3*54+27=351 GFNI512 products per four source rows.
Splitting each width16/state5 inner map into four-state and one-state pieces gives 52 scalar products per map.
Both maps therefore use 26 GFNI512 products per packet.
At rate 5/8, the combined count is 351/160 + 26/10 = 4.79375 per source symbol.
Two XOR-only reorganizations reduce this count. At eight-symbol block
granularity, reverse the five input blocks to make the 3-by-5 outer block
matrix Hankel. A symmetric 3-by-3 block product uses six block products;
the remaining 3-by-2 piece uses five. This gives 11*27=297 GFNI512 products
per four source rows. For each inner map, split sixteen coordinates into
5+5+5+1. Three symmetric 5-by-5 Hankel products and a five-term tail use
3*15+5=50 scalar products. Both inner maps then cost 25 GFNI512-equivalents
per packet, assuming suitable packing of the tails.

The improved total is 297/160 + 25/10 = 4.35625 per source symbol, about
2.17 times the regional 177/88 estimate. These are explicit circuit
estimates, not lower bounds or runtime measurements. Additional lane and
coefficient handling is not free. Even the improved count does not justify
another native implementation or a lengthy backup certificate campaign.

## Replay

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_source_envelope.py
python -B research/workstreams/char2_field_spin/theory/source_envelope.py --active 128 --z 0.9703716207340908
python -B research/workstreams/char2_field_spin/theory/source_envelope.py --active 256 --z 0.9413912472726532
```
