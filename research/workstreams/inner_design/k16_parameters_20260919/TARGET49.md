# K16 search for a full 49-bit margin

## Follow-up: full bound closed

The full calculation now closes at **49.3275773868 bits** for `(64,12,r=2)`.
All occupancies and their 512-bit replays passed. Matched runtime is
0.342–0.343 ms, essentially equal to the old `(128,19,1)` implementation.
See [FULL49_RESULT.md](FULL49_RESULT.md) for coverage, timings, and reproduction.
The remainder records the earlier screen and its then-incomplete status.

## Outcome of the screening round

The objective is a full setup-failure bound below 2^-49 at K=65536,
N=131072, and bad-word cutoff 13107. The BCH [256,128] outer and routing
distribution remain fixed. A one-active-row bound is only one component
of this objective.

The most promising new lead is the selected (64,12) subspace map with two
independent transvections per update. Its one-row contribution has an outward
49.3295959765-bit bound, checked by 512-bit replay with linear region products.
The ten largest rectangles from the old dense cover also pass after
recomputation under the two-round model. Their smallest individual margin
is 54.8650522520 bits. Those ten checks also passed 512-bit replay.

This is **not a full 49-bit certificate**. Every dense rectangle and the
complete sparse occupancy range still need evaluation and replay. Runtime
for two rounds has not been measured. No default, production kernel, paper
claim, or submitted artifact was changed.

## What the one-round search found

The bounded screen used the existing fixed map recipe, plus the two selected
t64 subspaces. All numerical diagnostics used the injection-shell transfer.

| Inner | Binary64 one-row diagnostic, bits |
|---|---:|
| Existing (128,19) | 44.361697 |
| Selected subspace (64,12) | 45.821646 |
| Selected subspace (64,13) | 46.864714 |
| (16,10) | 43.541086 |
| (32,11) | 45.016119 |
| (32,12) | 46.943630 |
| (32,13) | 48.119232 |
| (32,14) | 48.696599 |
| (32,15) | 49.037743 |
| Diagnostic (64,16) | 45.387107 |
| Diagnostic (64,19) | 45.525780 |

The existing six-tilt outward calculation gives 46.3774183397 bits for the
(64,13) subspace. It therefore does not clear 49 with that calculation.

For (32,15), a sharper outward adapter and different tilt bank give
49.0650099224 bits for one active row. The 512-bit replay passed. This
improves the earlier 48.714-bit one-row result without changing the code.
However, five dense points at occupancies 64,128,218,384,512 and density
coordinate 21/128 still have vacuous bounds. Retuning the actual fixed-input
bound at q=218 did not resolve this failure. This is a limitation of the
tested bounds, not a discovered low-distance codeword.

The historical (32,15) implementation took about 0.352–0.358 ms, slightly
longer than the existing (128,19) control. No new timing was run here.

## Why revisit two rounds

Shrinking t increases update frequency but also changes both fixed maps.
Increasing the number of transvections changes the state mixing without
changing the maps, routing, input partition, or emitted coordinates.

Write P for the marginal action of one transvection on a nonzero state and
J for uniform nonzero refresh. The existing IMT analysis gives
P=(I+J)/2 and J^2=J. Two independent transvections therefore give
P^2=I/4+3J/4=P/2+J/2. Zero remains zero.
Both transvections act before the single feedback addition; there is no
intermediate output or extra permutation.

| Fixed maps | One round: diagnostic | Two rounds: diagnostic |
|---|---:|---:|
| Existing (128,19) | 44.361697 | 51.197261 |
| Selected subspace (64,12) | 45.821646 | 49.332336 |
| Selected subspace (64,13) | 46.864714 | 50.634988 |

The older mixing study already certifies the (128,19,r=2) one-row
contribution at 51.2066374038 bits, with 512-bit replay. It also identifies
a separate dense bottleneck that barely improves with additional mixing.
See `../finite_migration/MIXING_ROUNDS.md`. That earlier result is not new
and does not establish a full 49-bit margin.

The current (64,12) map is attractive because its one-round full certificate
already has sparse and dense component margins of 61.978 and 54.989 bits.
Those numerical bounds cannot simply be reused for two rounds. They supply
candidate witnesses and a complete partition to reevaluate.

## Two-round dense checks

`target49_rounds_model.py` constructs a new fixed-weight transfer. It retains
the zero-source row and combines each nonzero-source row of the one-round
majorant with the corresponding uniform-refresh majorant, each with weight
one half. Both majorants use the same representation of the source measure.
This implements P^2=P/2+J/2 through positive upper bounds.

For the separate Bernoulli alternative, let m=2^s-1 and let c>=1 denote
the unnormalized Fourier cap. The nonzero-row coefficient ratio satisfies

```text
(c/(4(m+1)) + 3/(4m)) / (c/(2(m+1)) + 1/(2m))
    <= (4m+3)/(4m+2).
```

The zero row is unchanged. The adapter conservatively scales every other
row by this factor. It does not assume that an arbitrary one-round bound
remains valid. A first attempt using only the Bernoulli bridge failed on
the tested rectangles; retaining the fixed-weight transfer was essential.

The old cover has 1112 leaves. The probe selects its ten largest retained
contributions, then evaluates new bounds with their saved witnesses.
The weakest two recomputed leaf margins are 54.8650522520 and
62.3974735743 bits. The other eight exceed 72 bits. Selection by old
contribution does not guarantee that these remain the worst new leaves.

Five new adapter tests passed, covering exact baseline regression, an
independent binary64 transfer comparison, the rational Bernoulli ratio,
unchanged zero-source rows, and absence of mutations to historical functions.
The existing four mixing-law tests also passed. High-precision replay checks
arithmetic; completing the full certificate still requires complete coverage
and review of the new transfer adapter.

## Next step proposed by the screening round (now completed)

Reevaluate and replay all sparse occupancies and all 1112 dense leaves for
(64,12,r=2), then sum their bounds with the new one-row contribution.
The one-row bound leaves approximately 0.330 bits above the target; the
combined other contributions must fit in the remaining probability budget.
If the full sum clears 49, benchmark two rounds against both current options
in a matched serial batch. Keep (64,13,r=2) as the higher-headroom alternative.

Sources are the `target49*.py` scripts and `test_target49.py` in this folder.
Receipts are ignored under `measurements/target49_*`. No data files were
committed. The successful partial replay receipts explicitly set
`full_distance_proved` to false.
