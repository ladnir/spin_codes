# Directed coverage of the alternating-pair support domain

The faster alternating-pair route must handle every outer support, not
only favorable examples. The new coverage ledger records accepted and
unresolved support classes separately. It currently closes all classes
with at least 5,120 active outer groups. Smaller supports remain open.
This is a partial distance proof, not a whole-code certificate.

## Construction and counted objects

This note concerns K=2^20, rate 1/2, the GF256 Toeplitz grouped outer,
the h4 inner, and alternating adjacent-pair routing. It does not concern
the separate XOR outer experiment. The forbidden output weight is 734,003.
Probability is over the setup randomness of this fixed construction.

There are 8,192 outer groups and 4,096 pairs in each matching. Write
`(r,a,b)` for the number of active groups and the numbers of fully active
pairs under the two matchings. Each triple includes every even/odd lane
balance. Feasible triples satisfy

    1 <= r <= 8192,
    max(0,r-4096) <= a,b <= floor(r/2).

For r<8192, exclude the point 2a=2b=r: a nonconstant cycle cannot be
fully paired under both alternating matchings. There are exactly
11,465,830,401 feasible nonzero triples in the full domain.

## Cached bounds and exact accounting

The [backend](theory/hd_shifted_cover_backend.py) caches the directed
regional matrices for fixed numerical witnesses. It separates those
matrices from support-dependent normalizers and support counts.
The scalar exponent is convex on the box intersected with the
nonnegative-histogram domain. At most 24 rational vertices bound it.
Half-integer vertices use directed log-Gamma evaluation; they are not
rounded inward. One witness must bound the entire accepted box.

The [ledger](theory/hd_shifted_cover_ledger.py) maintains disjoint boxes,
counts their feasible integer triples exactly, and sums directed upper
bounds. Empty geometric boxes do not count as missing support classes.
Conversely, neither a search budget nor a failed bound discards a
feasible class. A whole-code flag requires complete global coverage and
a summed bound below 2^-40.

Floating-point search proposes boxes and witnesses only. Directed replay
rebuilds every accepted bound from the witness parameters. The search's
reported probability is never used as a certificate.

## Replayed coverage

Four frozen proposal plans have passed complete replay at both 256-bit
and 384-bit precision. Their rational upper bounds agree across precisions.

| Plan | Requested r range | Accepted boxes | Accepted triples | Log2 bound on accepted classes |
|---|---:|---:|---:|---:|
| 1 | 4096..7648 | 251 | 1,538,866,458 | -210.836753731405 |
| 2 | 6144..7167 | 95 | 627,922,944 | -210.836753731405 |
| 3 | 5120..6143 | 383 | 1,458,395,260 | -104.530220525834 |
| 4 | 5120..6143 | 542 | 1,702,713,344 | -104.530220525834 |

Plans 2 and 4 cover their entire requested domains. Plan 3 retains
244,318,084 unresolved triples, all with r<=6015; plan 4 closes these
using the adaptive-alpha witnesses. The plans overlap through inherited
leaves, so their accepted counts must not be added as a coverage count.

For a contiguous completed range, use the following disjoint weight bands:

- 5120..6143: plan 4;
- 6144..7167: plan 2;
- 7168..7648: plan 1, whose band bound is at most 2^-937.079734692701;
- 7649..8192: the inverse-image ball bound, at most 2^-293.684921198294.

Their union contains 2,420,639,745 feasible triples. Directed addition gives

    log2 Pr[any failing message in these classes]
        <= -52265110262917/500000000000.

No independence between the four events is assumed. Restricting a leaf
to a weight band retains its full probability bound; no contribution is
prorated by the number of retained profiles. The remaining domain below
5120 has 9,045,190,656 triples and still prevents a whole-code claim.

The composite replay imports freshly checked receipts on disjoint
assigned weight bands. It preserves every source gap, rejects overlap,
and updates the destination only after validating the complete import.
The 384-bit joined ledger has 664 accepted boxes, including the dense
tail. Its geometric and feasible-profile accounting both pass. Source
JSON contains proposals only; no saved claimed receipt is trusted.

## Reproduction

The bounded proposal commands are recorded in
[the unsigned analysis](HD_SHIFTED_UNORIENTED.md). Local JSON plans are
ignored research artifacts. The exact plans replayed here have SHA256:

    plan 1: 384071145a74b64e99e0cd81d474b096e441d5937a4a24dd9db9ed3a9a338314
    plan 2: 563aacd1f3d90fea2368edce7a6f5d3b2cee2dac8d03401450cf11352c9c997b
    plan 3: d24864da7c99d094a295d307d34cd79ba73b1d1e10ea35dd6ca175e5b27d6499
    plan 4: c8ce97119f9c1422a95f724bba1ca7000e9e4282caa5b09bf1d2ce52dbbec4a2

From the theory directory, replay each plan separately:

    python -B hd_shifted_cover_ledger.py --plan hd_shifted_bulk_plan.json --precision 384
    python -B hd_shifted_cover_ledger.py --plan hd_shifted_bulk_plan2.json --precision 384
    python -B hd_shifted_cover_ledger.py --plan hd_shifted_bulk_plan3.json --precision 384
    python -B hd_shifted_cover_ledger.py --plan hd_shifted_bulk_plan4.json --precision 384

Replay and combine the current disjoint bands in one command:

    python -B hd_shifted_cover_ledger.py --band-plan hd_shifted_bulk_plan4.json 5120 6143 --band-plan hd_shifted_bulk_plan2.json 6144 7167 --band-plan hd_shifted_bulk_plan.json 7168 7648 --precision 384

The CLI reports unresolved counts, a conservative completed upper weight
range, exact rational aggregate bounds, and the input plan's hash.
Missing optional plans cause their integration tests to skip explicitly;
ordinary unit tests never regenerate a long search implicitly.

No native encoding measurement, production change, or paper change was
made during this proof iteration.
