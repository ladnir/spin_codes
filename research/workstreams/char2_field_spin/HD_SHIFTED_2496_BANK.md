# One fixed bank centered at weight 2496

This experiment addresses the open weight interval `2154..2827` for the
[swapped-pair Toeplitz construction](HD_SHIFTED_PAIR_SWAPS.md). Weight counts
active scalar-code input groups. The experiment retains the
`invertible-toeplitz` outer bound; it does not change the encoder or setup
distribution.

The first proposed mark, `2^(-2496/1000)` evaluated in binary64, sufficed.
The fractional exponent `alpha` is the exact rational value of binary64 `0.1`. Fifty marker
pairs were fitted once at weight 2496: full-pair counts spaced by 32,
plus distances 1, 2, 4, 8, and 16 from each endpoint. A separate
19-evaluation scalar optimization fixed one common positive state
potential. No neighboring weight was fitted, and no second mark was tried.

The frozen exact parameters are

```text
z     = 798341867030649/4503599627370496
alpha = 3602879701896397/36028797018963968
v     = (1, 2576174089033825/4611686018427387904)
```

The grid is 128. All marker literals are retained in
[`hd_shifted_potential_band_2496_choice0_2154_2827_selected.json`](theory/hd_shifted_potential_band_2496_choice0_2154_2827_selected.json),
whose SHA-256 is
`e8d485bf47138fb9caaea288bc26fccaab925a75b2a0b66d93ca21a3208d92d8`.

## Screen and replay

The FLOAT screen evaluated every integer weight in `2154..2827`, using
log-factorial cyclic-run counts. Its contiguous passing component around
2496 is `2275..2683`: 409 weights, each below the chosen -60-bit threshold.
The FLOAT union bound is -70.8309761316041 bits. These values alone are
not certificates. The complete screen took about 42.5 seconds on the
local proof runner, including central fitting.

Fresh 256- and 384-bit replays give the same directed log2 union ceiling
`-41223743364519/500000000000` (-82.447486729038 bits). It covers all
630,556,572 feasible matching-profile triples in the 409-weight interval.
Every per-weight rational ceiling also agrees between the two precisions.
Both replays normalize every required histogram before taking atlas minima
and sum a bound for every selected integer weight. The interval remains a
partial family, not a whole-code certificate.

The saved receipts and SHA-256 hashes are:

- [256-bit receipt](theory/hd_shifted_potential_band_2496_2275_2683_directed256.json):
  `903c279770cb86ee223127f0f0af51fdced1fdf6a59390c5951602bfd02b647d`.
- [384-bit receipt](theory/hd_shifted_potential_band_2496_2275_2683_directed384.json):
  `e29de1ffd22bd891ab938475525f30cee4adc095b42c786246cacdecfc4d19d6`.

At the lower endpoint 2275, directed replay gives -201.796888489467 bits,
whereas the FLOAT screen gives -70.831435675500 bits. The final union is
instead dominated by the upper endpoint 2683. This comparison does not
certify any weight outside the selected interval.

From the theory directory, reproduce the frozen-witness replays with

```text
python -B hd_shifted_potential_band.py --proposal hd_shifted_potential_band_2496_choice0_2154_2827_selected.json --max-weights 409 --precision 256
python -B hd_shifted_potential_band.py --proposal hd_shifted_potential_band_2496_choice0_2154_2827_selected.json --max-weights 409 --precision 384
```

The bounded proposal driver is `hd_shifted_2496_bank.py`. Its choice list
contains at most three marks; this experiment used only choice 0. Its
tests check the selection boundary, frozen scope, and optional matching
cross-precision receipts without launching proof computations. All five
local tests pass; the combined related suite has 30 passing tests.

The next step is disjoint ledger integration. This bank alone does not
cover `2154..2274` or `2684..2827`; those intervals need separate evidence.
