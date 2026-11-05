# Two banks adjacent to the verified 2496 interval

These experiments address `2154..2274` and `2684..2827` for the
[swapped-pair Toeplitz construction](HD_SHIFTED_PAIR_SWAPS.md). Weight counts
active scalar-code input groups. Both banks retain the explicit
`invertible-toeplitz` outer bound and the same encoder and setup law.

The source is the frozen [2496 proposal](theory/hd_shifted_potential_band_2496_choice0_2154_2827_selected.json).
For each new center, the proposer scales `log(z)` by the ratio of center
weights and preserves the exact source `alpha`. It fits coefficient
markers once, with full-pair spacing 32 and the helper's endpoint markers.
It then retunes one common positive state potential at that center.
Each potential optimization used 28 evaluations. No other mark or exponent
was tried, and no neighboring weight was fitted.

The common fractional exponent is
`3602879701896397/36028797018963968`. The other exact parameters are:

| Center | Mark `z` | State potential `v` | Markers |
|---|---|---|---:|
| 2208 | `7797861747750693/36028797018963968` | `(1,6476607833312521/9223372036854775808)` | 40 |
| 2752 | `668536658701055/4503599627370496` | `(1,8806658302716409/18446744073709551616)` | 49 |

Both banks use grid 128. All marker literals remain in their frozen
proposal artifacts. The FLOAT screens evaluated every integer weight
in each assigned interval, using log-factorial cyclic-run counts. Their
union bounds were -1208.0382048094161 and -2717.263358310497 bits,
respectively. These screen values alone are not certificates.

## Directed replay

The 2208 bank has matching 256- and 384-bit replays, including every
per-weight rational ceiling. Its directed log2 union ceiling is
`-1208038932506483/1000000000000` (-1208.038932506483 bits).
The interval contains 121 weights and 148,517,723 feasible matching-profile
triples.

The 2752 bank also has matching 256- and 384-bit replays, including every
per-weight rational ceiling. Its directed log2 union ceiling is
`-131291896793403/40000000000` (-3282.297419835075 bits). The interval
contains 144 weights and 273,699,888 feasible matching-profile triples.

Each replay reconstructs all operators from the frozen exact witness,
normalizes each histogram before taking atlas minima, and evaluates every
integer weight. These bounds cover all scalar supports and even/odd splits
within their assigned intervals, not the whole construction.

The original [2496 receipt](HD_SHIFTED_2496_BANK.md) remains separate and
unchanged on `2275..2683`. The three intervals are pairwise disjoint and
partition `2154..2827`. The final ledger must sum the three probability
bounds; it must not reuse the middle interval's union bound as the union
of all three. Together, the separate receipts cover all 674 weights and
1,052,774,183 matching-profile triples in this interval.

## Artifacts and reproduction

The frozen proposal hashes are:

- [2208 proposal](theory/hd_shifted_potential_band_2208_2154_2274.json):
  `a8190c89dd626bcefd5013cf1b8364921b469f5d6d6fb27a22c947dd372fca1f`.
- [2752 proposal](theory/hd_shifted_potential_band_2752_2684_2827.json):
  `a21a2e0afbe90dd755cada097e586dafbc400effae4d7b74fad0dc2ffb40bda0`.

Receipt hashes:

- [2208, 256 bits](theory/hd_shifted_potential_band_2208_2154_2274_directed256.json):
  `43393f9a56fb8ba86a47acd6d1eb128370fdbb7493028a5029caecd9b25f2c92`.
- [2208, 384 bits](theory/hd_shifted_potential_band_2208_2154_2274_directed384.json):
  `f21236e742d112f3061365346f73ea509e7c0b293f78730fea7a5527872d1f50`.
- [2752, 256 bits](theory/hd_shifted_potential_band_2752_2684_2827_directed256.json):
  `41e777761b48ef9f65bef87e954b999b169130f1e9a6b5bd7249cd52d8e95f98`.
- [2752, 384 bits](theory/hd_shifted_potential_band_2752_2684_2827_directed384.json):
  `39cebc1883c5e5f9b21e64583d8d21f135808f7896d42ee04771e18bbc14b045`.

From the theory directory, replay each proposal at precision 256 or 384:

```text
python -B hd_shifted_potential_band.py --proposal hd_shifted_potential_band_2208_2154_2274.json --max-weights 121 --precision 384
python -B hd_shifted_potential_band.py --proposal hd_shifted_potential_band_2752_2684_2827.json --max-weights 144 --precision 384
```

`hd_shifted_2208_2752_banks.py` reproduces the bounded proposal workflow.
Its tests verify frozen source parameters, complete interval membership,
receipt hashes, scope, and cross-precision agreement without rebuilding
any bank. All 29 focused and related tests pass. No native source, shared
proof helper, or benchmark was changed. The next step is disjoint ledger
integration, not another fit for these intervals.
