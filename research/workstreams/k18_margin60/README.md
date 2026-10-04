# K18: 60-bit margin without an encoder change

The exact retained K18 construction now has a complete **60.7509985457-bit**
setup-failure bound at **10% relative distance**, improved from 42.0908400604
bits. The outer, routing, inner maps, setup distribution, and implementation
are unchanged. Its previously measured precomputed transpose time remains
about **0.467 ms** for 128-bit elements; no benchmark was run in this workstream.

The code has K=262144 and N=524288. Each of 2048 outer groups maps 128 bits
to 256 bits using four parallel GF16 RS[16,8] rows and independent transitive
16-bit symbol randomizers. Four-bit packets pass through 64 routing regions.
The exact certified t64/s22 maps use independent uniform invertible state
updates, zero initial state, continuous state across all regions, and no flush.
The good setups have minimum distance at least **52429**. The probability
bound concerns the ideal setup ensemble, not the measured distance of each seed.

Map identity:
`e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581`.

## What improved

The old first moment can count many low-weight messages under the same rare
routing arrangement. Counting all those messages is unnecessarily expensive
when the event of interest is whether any such message exists.

Fix a set of active outer groups and condition on their potential slot counts
in each physical step. The conditional failure probability is at most the
smaller of one and the conditional expected number of bad messages. We bound
this clipped quantity by a fractional power before averaging the route.
The [proof review](SLACK_REVIEW.md) gives the conditioning, matrix-path
inequality, and coefficient formula.

The numerical witness uses tilt 17/200 and fractional exponent 4/5. Uniform
packet values, including zero, are averaged before taking fractional powers.
Each pair of chronological physical matrices is multiplied before entrywise
powering. Routing probabilities and the union over active-group subsets stay
outside that power. The matrix calculation carries state across region boundaries.

At exponent one, this expression recovers the old uniform-input bound.
The q=89 regression reproduces 44.841582549 bits. Exponent 4/5 improves that
particular contribution to about 861 bits. This is not the whole-code margin:
the final union is dominated by the unchanged one-active-group contribution.

| Component after combining bounds | Margin |
|---|---:|
| One active group, retained exact-shell bound | 60.751047 bits |
| Two active groups, retained bound | 126.546358 bits |
| All occupancies q>=3, old/new minimum and exact sum | 75.612280 bits |
| Complete union over q=1..2048 | **60.7509985457 bits** |

Fresh outward bounds replace the old endpoints for q=62..143, inclusive.
Every other endpoint is retained and authenticated. The assembler takes
the smaller valid bound for each occupancy and sums all 2048 dyadics exactly.
There is no interpolation over missing occupancies.

## Bounded exploration

Ordinary witness tuning was not sufficient. Six local-row selection activities
at tilt .085 gave no material improvement. Six tilts from .081 through .086
improved the floating combined proposal only to 42.430835903 bits. Those
proposals are not part of the new certificate.

The subsequent two-step fractional screen gave the useful improvement.
Physical comparison operators were freshly generated from the exact s22 maps
using 192-bit Arb arithmetic. Their dyadic upper endpoints were exported and
authenticated. The global replay uses 256-bit Arb, including the fractional
powers and exact-placement matrix polynomial. It does not require another
local census. The successful replay took about 5.4 seconds locally, after
the approximately 66-second census/local-bound baseline run.

The archived 42-bit certificate and its source files remain unchanged.
This is a derived complete certificate using authenticated old endpoints
and fresh improved endpoints, not a fresh recomputation of every old component.

The independent receipt verifier authenticates 112 files and four original
component receipts. It checks all 2048 old/new selections, all 82 improvements,
both exact unions, and the exact comparison with 2^-60. It does not independently
regenerate the local census or the fractional numerical endpoints.
Separate mathematical and code reviews checked the new conditioning argument.
All 21 tests pass, including exact rational checks of fractional-power upper
bounds, exponent-one equivalence, chronological matrix products, preservation
of tiny positive terms, and rejection of altered certificate records.

## Records and reproduction

The selected local receipt is
`tmp/k18-margin60-activity-screen-physical-a1_2.json`, SHA256
`fb5d68b115c8010726a7ff4cbf47da199e5d9cf217a411e1a970327d2a28d478`.
The complete new receipt is
`tmp/k18-margin60-fractional-whole-p256-v2.json`, SHA256
`44fcdb0389d5cc7be870bdfcbeb0ec5500193706b7e690313054971d052a01d0`.
It authenticates 112 source/receipt pins. The v1 attempt failed during
decimal serialization of the exact union and left an incomplete file;
it is not a certificate. V2 stores the large exact union in hexadecimal.
Generated receipts remain ignored and are not added to the artifact.

From the repository root, regenerate the exact local family and one baseline
regression using a fresh output name:

```sh
python -B research/workstreams/k18_margin60/tune.py \
  --output tmp/k18-margin60-reproduction-local.json \
  --tilts .085 --activities 1/2 --min-q 89 --max-q 89 --seconds 300
python -B research/workstreams/k18_margin60/replay_fractional.py \
  --local tmp/k18-margin60-reproduction-local-physical-a1_2.json \
  --output tmp/k18-margin60-reproduction-whole.json
python -B research/workstreams/k18_margin60/verify_fractional.py \
  tmp/k18-margin60-reproduction-whole.json
python -B -m unittest discover -s research/workstreams/k18_margin60 -p 'test_*.py' -v
```

The replay also needs the authenticated original whole receipt identified in
the producer. Its full regeneration recipe remains in the
[original K18 checkpoint](../k16_design/README.md#power-of-two-checkpoint-k--218).
This work does not replace those frozen local models or their assumptions.

## Promotion decision

The certified four-bit lineup now reaches at least 60 bits at all three
selected endpoints: paired-s15 at K16, exact s22 at K18, and larger-outer s20
at K20. Retain the existing K18 kernel and promote that exact map; a larger
outer or extra mixing is unnecessary for this target.
The certificate does not cover intermediate lengths automatically.
