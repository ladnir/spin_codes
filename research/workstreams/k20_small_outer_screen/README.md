# K20 small outer: complete 58.75-bit certificate

The small outer now closes at K=2^20: **58.7503767402 bits of setup-failure
margin at 10% relative distance**. The proof covers every one of the 8192
active-group occupancies. It meets the requested 50-bit target without
changing the construction, inner maps, or setup distribution.

The [initial implementation benchmark](../../../spin/experiments/k20_small_outer/README.md)
takes 5.213749 ms. The subsequent
[exact-map optimization](../../../spin/experiments/k20_small_outer/OPTIMIZATION.md)
reduces this to **2.962029 ms**, versus **3.299650 ms** for the larger-outer
library profile in the same holdout. The code and certificate are unchanged.
The small outer is now the research performance winner at the accepted
50-bit target. The subsequent [library promotion](../../../spin/CERTIFIED22.md)
adds this code as `PacketRsT64S22`, alongside the larger-outer profile.

## Complete result

Good setups have minimum distance at least **209716**. Under the fixed
ideal independent setup ensemble, the probability of failing that guarantee
is at most the exact dyadic sum recorded in the whole receipt.

| Component | Margin |
|---|---:|
| q=1, exact expected shells | 58.7503769443 bits |
| q=2, exact expected shell pairs | 122.9604491274 bits |
| All q>=3, exact sum of selected endpoints | 81.5031655382 bits |
| Full union q=1,...,8192 | **58.7503767402 bits** |

The certificate combines 11 authenticated component receipts. Ordinary
small-tilt bounds and two fractional ranges cover the lower occupancies.
Three fractional ranges cover the middle. A marker coefficient bound covers
every q=616,...,8192, with a 6670.925-bit union for that dense interval.
The [method note](METHODS.md) states the conditioning and dense coefficient
argument. Neither method resets the inner state at region boundaries.

The whole receipt is `tmp/k20-small-outer-whole-v1.json`, SHA256
`90492fc15834d2e28a9ddc8289a68dbf7bd17984146524016f9361a9ec54fe65`.
Its component manifest pins every input receipt and inherited mathematical
source. `assemble_whole.py` takes exact per-occupancy minima and sums all
8192 dyadics using integer arithmetic. Decimal displays do not determine
whether the target passes.

`verify_whole.py` is independent of the numerical producers. It checks the
explicit map arrays, source and receipt hashes, full coverage, sparse shell
folds, selected minima, and exact union. It authenticates the local census
and matrix endpoints rather than independently regenerating them.
The saved audit `tmp/k20-small-outer-whole-v1-independent-audit.json` passes:
138 files, 11 component receipts, 10 physical receipts, and all 8192
occupancies authenticate. All **36 workstream tests pass**, including
mutation tests against the actual complete receipt.

```sh
python -B research/workstreams/k20_small_outer_screen/verify_whole.py \
  tmp/k20-small-outer-whole-v1.json
python -B -m unittest discover -s research/workstreams/k20_small_outer_screen -p 'test_*.py' -v
```

## Fixed construction

Keep the exact K18 code distribution and change only the total length:
K=1,048,576, N=2,097,152, and bad-output cutoff D=209,715.
There are 8192 outer groups. Each maps 128 input bits to 256 output bits,
using four GF16 RS[16,8] rows and independent transitive 16-bit randomizers.
Four-bit packets occupy 64 routing regions, each containing 8192 slots.

The inner uses the retained t64/s22 maps, independent uniform invertible
state updates, zero initial state, continuous state, and no final flush.
Its map identity is
`e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581`.
No implementation or setup-distribution change was made.

The [K18 conditioning argument](../k18_margin60/SLACK_REVIEW.md) applies
with 8192 groups instead of 2048. A region now contains 256 two-step macros,
instead of 64. The number of regions remains 64. The outer comparison factor
remains beta=2^256/(65535)^8. The union over active-group sets stays outside
the fractional power; thinning by uniform packet labels occurs before it.

## Initial bounded screen

For q=1, `sparse.py` regenerates the exact s22 census and evaluates the
expected RS shells with outward arithmetic. Four tilts, .0002, .0004,
.0006, and .0008, give a combined bound of **58.7503769443 bits**.
This includes the union over all 8192 groups. The shell fold was independently
checked with exact rationals, and all 93 source pins were authenticated.

For selected larger q, `outward_screen.py` reuses authenticated, length-independent
physical operators. It recomputes regional placement and the global bound
at 256-bit Arb precision for the new geometry. It does not reuse K18 global
endpoints. At tilt .085 the comparison is:

| Active groups q | Ordinary first moment, alpha=1 | Routing-conditioned bound, alpha=.8 |
|---|---:|---:|
| 256 | -2776.26 bits | 27.37 bits |
| 320 | -1324.82 bits | 2035.50 bits |
| 356 | -1057.24 bits | 2617.54 bits |
| 400 | -1196.86 bits | 2863.43 bits |
| 512 | -3541.67 bits | 1505.32 bits |

A negative entry means that particular upper bound is greater than one;
it does not show a bad code. These are separate occupancy bounds, not margins
for the full code. The initial q=256 witness was weak. Changing only the
fractional exponent to .7, at the same tilt .085, gives **666.36 bits** there.

Two additional physical families were exported from one fresh census, at
tilts .02 and .04. They give the following outward selected-occupancy bounds:

| Active groups q | Tilt | Fractional exponent | Margin |
|---|---:|---:|---:|
| 64 | .02 | .8 | 247.75 bits |
| 96 | .02 | .8 | 626.76 bits |
| 128 | .04 | .8 | 299.49 bits |
| 160 | .04 | .8 | 1036.35 bits |
| 192 | .04 | .8 | 1178.49 bits |
| 224 | .04 | .8 | 823.90 bits |
| 256 | .085 | .7 | 666.36 bits |

These witnesses remove the observed obstruction at the sampled occupancies.
They do not cover every integer between them. At q=32 and 48, the tested
.02 tilt still failed; those occupancies needed other witnesses or a sharper
sparse analysis. That initial screen did not evaluate q=2 or the dense tail.

The initial floating placement screen stopped on underflow without producing
a receipt. All reported intermediate-occupancy results use the outward Arb
replacement, which retains arbitrarily small positive terms.

## Records and reproduction

Generated receipts stay in ignored `tmp/`; sources and concise results stay here.

- q=1: `tmp/k20-small-outer-q1-p192-v1.json`, SHA256
  `27f1b3586b7c40b7f9429466f3baf6e3aa5bf751ea86cb5aed1d73f7d400ed29`.
- Intermediate occupancies at .085: `tmp/k20-small-outer-middle-outward-v1.json`, SHA256
  `43e5909d0833ac2e0d8f9ca1b77415622cbb18e5959a8f40931b79fc1b70eb2e`.
- Local .085 operators: the authenticated export recorded in the
  [K18 result](../k18_margin60/README.md#records-and-reproduction).
- Additional physical exports: `tmp/k20-small-outer-physical-p192-v1-t002.json`
  and `tmp/k20-small-outer-physical-p192-v1-t004.json`, respectively SHA256
  `e82e2f4de06973fa66b509e5a21d7d49114c73f7fec2500afea05e27948b4ed9`
  and `1605f1e68643cfb54abdcff61539677c5c6cc62c0051c889473d5d5e39fe88f7`.
- Selected .02 and .04 bounds: `tmp/k20-small-outer-lowmiddle-outward-t002-v1.json`
  and `tmp/k20-small-outer-lowmiddle-outward-t004-v1.json`, respectively SHA256
  `247a614bc8d89a4982c00dded2a2ac87aed4bea3bafdbeebc83d2ed08be45c73`
  and `17f02b65b54e229b81af150e6888af6df219880a2c6fca269d6dc6cd4c9553f0`.
- Refined q=256 witness: `tmp/k20-small-outer-q256-a07-v1.json`.

Run from the repository root, using fresh output names:

```sh
python -B research/workstreams/k20_small_outer_screen/sparse.py \
  --output tmp/k20-small-outer-q1-fresh.json
python -B research/workstreams/k20_small_outer_screen/outward_screen.py \
  --local tmp/k18-margin60-activity-screen-physical-a1_2.json \
  --output tmp/k20-small-outer-middle-fresh.json --alphas .8 .9 1
python -B research/workstreams/k20_small_outer_screen/export_local.py \
  --output-prefix tmp/k20-small-outer-physical-fresh
python -B research/workstreams/k20_small_outer_screen/outward_screen.py \
  --local tmp/k20-small-outer-physical-fresh-t002.json \
  --output tmp/k20-small-outer-t002-fresh.json --alphas .8 1 \
  --qs 32 48 64 96 128 192 256 320
python -B research/workstreams/k20_small_outer_screen/outward_screen.py \
  --local tmp/k20-small-outer-physical-fresh-t004.json \
  --output tmp/k20-small-outer-t004-fresh.json --alphas .8 1 \
  --qs 64 96 128 160 192 224 256 320
python -B research/workstreams/k20_small_outer_screen/outward_screen.py \
  --local tmp/k18-margin60-activity-screen-physical-a1_2.json \
  --output tmp/k20-small-outer-q256-fresh.json --alphas .7 --qs 256
python -B -m unittest discover -s research/workstreams/k20_small_outer_screen -p 'test_*.py' -v
```

The new tests check exact K18 scalar equivalence, K20 cutoff and group-count
scaling, and rejection of invalid geometry before loading the local operators.

## Decision

Keep the larger-outer K20 baseline available: its 68.10-bit whole certificate
is the higher-margin option. Do not treat the old
ordinary-first-moment obstruction as a reason to discard the small outer.

The complete proof above supersedes the initial partial screen. The initial
kernel lost at K20, but streaming stores and exact-map implementation tuning
removed that regression. The optimized inner/routing phase takes about
1.61 ms and its outer takes 1.34 ms in separate phase diagnostics.
These changes preserve the certificate. The later library promotion adds
a distinct profile without changing existing profile identifiers or maps.
The current q=1 endpoint limits this proof to below 58.751 bits. Reaching
60 would require a sharper q=1 bound or a parameter change; it is not needed
for the accepted 50-bit target.

## Replaying the complete components

Use fresh output prefixes. The retained q=1 recipe is given above.
`low_complete.py --output-prefix <prefix>` regenerates q=2, the low ordinary
tilts, and physical families at .0032, .006, and .01. For assembly, keep its
q=2 receipt and ordinary tail receipts at .0008, .0016, and .0032. Replay
q=3,...,63 with `range_replay.py` at alpha=.6 using its .006 and .01 exports.

The middle uses `range_replay.py` with these local tilt / alpha / maximum-q
choices, always starting at q=3: (.02,.8,160), (.04,.8,300), and (.085,.7,768).
The earlier export recipes regenerate .02, .04, and .085 if necessary.

For the dense tail, `dense_complete.py --export-prefix <prefix>` generates
physical families at .15, .25, .5, 1, and 2.1972246 from one fresh census.
Pass those five outputs to the same script with `--local`, a fresh `--output`,
`--min-q 616 --alphas .5 .6 .7 .8 1`. Its default marker grid is j/128 for
j=1,...,128. This checks every integer through 8192.

Finally pass the 11 receipt paths, in the order listed by the saved whole
manifest, to `assemble_whole.py --output <fresh-whole.json> --target-bits 50`.
Run the independent verifier on that output. No K18 global endpoint is used.
Generated receipts remain ignored; the sources and concise results stay here.
