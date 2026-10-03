# Bounded K18/K20 scaling screen

**Status: parked, 2026-10-03.** The user selected the certified four-bit
baseline over further width-eight work. The modest smaller-size gains and
K20 regression do not justify another promoted construction. The results
below remain a research checkpoint; no additional proof run is scheduled.

This study changes the length of the existing width-eight, wider-outer,
24-bit-state construction. It does not change its local maps. The retained
results are **partial outward screens, not whole-code certificates**. K18
clears 40 bits at every sampled occupancy at **10% distance**. K20 clears the sampled occupancies
through eight, but the best tested occupancy-16 witness does not give a useful
probability bound. A negative witness margin does not establish code failure.

The separate [performance campaign](PERFORMANCE.md) measures 0.373 ms at
K18 and 3.389 ms at K20. The older four-bit K20 control is faster at 3.268 ms.
The [native experiment](../../../../spin/experiments/packet8_wider24_scaling/README.md)
contains the store-policy variants and correctness harness.

## Results and coverage

Here `q` is the number of nonzero outer groups. A displayed margin is minus
the base-two logarithm of the corresponding union-bound contribution. All
numbers below are diagnostic displays; exact dyadic endpoints are in
`summary_v1.json`.

| Occupancy q | K18 margin, bits | K20 margin, bits |
|---:|---:|---:|
| 1 | 136.28838 | 106.78409 |
| 2 | 140.05669 | 136.29199 |
| 3 | 61.37227 | 54.88740 |
| 4 | 137.00205 | 127.60221 |
| 8 | 63.10980 | 43.34987 |
| 16 | 54.36513 | -8.44152 |
| 64 | at least 200 | not tested |

The uncapped K18 q64 display is 400.12310 bits; the stored endpoint is
deliberately weakened to `2^-200`. The K20 q16 stored endpoint is the trivial
probability cap one. No result is interpolated to an untested occupancy.
There remain 1,017 untested occupancies at K18 and 4,090 at K20. Even the
positive sampled results therefore do not imply a complete union bound.

The useful witnesses are:

| Scope | theta proposal | alpha |
|---|---:|---:|
| K18 q1 | per-support minimum over the tested local witnesses | not used |
| K20 q1 | 0.000625 | not used |
| K18 / K20 q2 | 0.0025 / 0.000625 | 0.7 |
| K18 / K20 q3, q4 | 0.0075 / 0.001875 | 0.5 |
| K18 / K20 q8 | 0.01 / 0.0025 | 0.4 |
| K18 / K20 q16 | 0.03 / 0.0075 | 0.3 |
| K18 q64 | 0.12 | 0.3 |

`theta` only selects a representable dyadic `z` near `exp(-theta)`. Every
local and macro receipt records and uses the exact selected rational `z`;
the exponential approximation is not a proof assumption.

## Unchanged construction and changed geometry

The outer group remains eight parallel GF16 RS[16,8] rows, with 256 input
bits, 512 output bits, and sixteen aligned 32-bit symbols. Independent
nonzero-transitive symbol randomizers give the same exact expected
64-byte-packet support shells and the same pointwise constant

`beta = 2^512 / (2^32 - 1)^8`.

The local state remains the literal GF256^3 map on eight bytes:
`A(a,b,c)[h] = a + h*b + h^2*c`,
`C(x) = (sum x[h], sum h*x[h], sum h^2*x[h])`, for h=0,...,7 in the
AES polynomial basis. Emission precedes the fresh nonzero-transitive update.
The state starts at zero, persists across all regions, and is not flushed.

| Parameter | K18 | K20 |
|---|---:|---:|
| Input bits K | 262,144 | 1,048,576 |
| Output bits N | 524,288 | 2,097,152 |
| Outer groups G = slots per region | 1,024 | 4,096 |
| Byte-packet regions | 64 | 64 |
| Physical 64-bit steps per region | 128 | 512 |
| Four-step, 32-slot macros per region | 32 | 128 |
| Bad-weight cutoff floor(N/10) | 52,428 | 209,715 |

The authenticated active-label local matrices and G4 matrices are independent
of code length. `global_outward.placement` is called with the changed macro
epoch count; `q1_outward.q1_support_upper` instead receives the changed
physical epoch count. Every call also supplies the changed group count and
cutoff where required. The fixed-width G4 kernel still has 32 potential
slots.

For q>=2, let `R_q` be the exact-placement regional upper matrix assembled
from the authenticated G4 matrices, already fractionally powered inside
each four-step block. The implemented contribution is

`C(G,q) * [beta^q * z^(-cutoff) * a^(64*q)]^alpha * e0 R_q^64 1`,

where `a=((1+z)/2)^8`. The factor `64*q` counts potential packets over all
regions; it does not increase with the number of empty physical steps.
The matrix power carries the same state through region boundaries. There
is no reset to zero between regions and no additional factor `C(G,q)`
inside the fractional power.

For q1, the physical zero/one-active-packet operators are averaged over
128 or 512 steps per region. Their full matrices are composed across 64
regions. The exact outer shell counts are then combined with G times the
per-support witness minimum. This uses the q1 bound, not the pointwise beta
majorant.

## Computation, review, and retained evidence

`scaling_screen.py` imports frozen iteration5 arithmetic. It uses scaled
positive binary64 arithmetic with that arithmetic's outward-error contract,
then makes endpoint cap/floor decisions by exact rational comparison. An
independent reviewer identified an earlier rounded-log cap decision that
could round inward at a boundary. It was fixed before the retained runs;
the regression test uses a value immediately above `2^-200` whose rounded
display margin is exactly 200.

`lowq_v1.json` is preserved only as a **superseded development run**. Its
producer source changed after review, so it is excluded from authentication
and from every reported combined result. `lowq_v2.json` repeats the useful
old-tilt points with the corrected code. These old tilts alone are too large
for the sparse cases at the longer lengths; their negative margins are not
evidence that the construction fails.

`near_one.py` generated four new local witnesses at theta
0.0025, 0.000625, 0.0075, and 0.001875. One exact 2^24-state census was shared
across them; no large census arrays were written to disk. The census took
about nine seconds, and all four outward locals and their G4 macros took
about 93 seconds together. Four additional alphas reused existing local
receipts without a new census. Accumulated numerical wall times, including
the superseded initial run and the corrected replay, were about 4.4 minutes.
The proof worker ran locally without encoder benchmarks or remote commands;
the separate performance campaign was serialized on Peach. Neither part
changed production code.

The retained partial screen receipts are:

- `lowq_v2.json`: existing theta 0.01/0.03 controls;
- `near_one_lowq_v1.json`: the four new near-one witnesses;
- `alpha_k18_v1.json` and `alpha_k20_v1.json`: q8/q16 alpha reuse;
- `matched_k18_q64_v1.json`: one larger, matched-density K18 point;
- `summary_v1.json`: exact endpoint minima and explicit incomplete coverage.

The summary authenticates 91 source/receipt pins. All seven iteration7 tests
pass, including tiny exact noncommuting placement enumeration, explicit
K20 parameter forwarding, the endpoint-boundary regression, and summary
coverage checks. An independent agent reviewed the changed geometry,
endpoint arithmetic, and near-one producer without finding an outstanding
issue. The scope remains partial even though its retained endpoints are
outward under the inherited arithmetic contract.

## Reproduction

Run from the repository root with `C:/Python314/python.exe`. Existing
receipts are never overwritten. Change the output names below when replaying
in a populated checkout; the names shown are the original commands.

```powershell
C:/Python314/python.exe -B -m unittest discover -s research/workstreams/packet8_codesign/iteration7 -p test_*.py -v
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration7/near_one.py
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration7/scaling_screen.py --max-theta .03 --output research/workstreams/packet8_codesign/iteration7/lowq_v2.json
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration7/scaling_screen.py --only-extra --extra-macros research/workstreams/packet8_codesign/iteration7/near_one_v1/macro_theta_1_400_alpha_7_10.json research/workstreams/packet8_codesign/iteration7/near_one_v1/macro_theta_1_1600_alpha_7_10.json research/workstreams/packet8_codesign/iteration7/near_one_v1/macro_theta_3_400_alpha_1_2.json research/workstreams/packet8_codesign/iteration7/near_one_v1/macro_theta_3_1600_alpha_1_2.json --output research/workstreams/packet8_codesign/iteration7/near_one_lowq_v1.json
```

The four reused-local macros are generated by the frozen macro producer:

```powershell
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration5/macro_outward.py --local research/workstreams/packet8_codesign/iteration5/local_theta003_v1.json --alpha 3/10 --output research/workstreams/packet8_codesign/iteration7/macro_theta003_alpha3_10.json
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration5/macro_outward.py --local research/workstreams/packet8_codesign/iteration7/near_one_v1/local_theta_3_400.json --alpha 3/10 --output research/workstreams/packet8_codesign/iteration7/macro_theta0075_alpha3_10.json
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration5/macro_outward.py --local research/workstreams/packet8_codesign/iteration5/cache_v1/local_theta_1_100.json --alpha 2/5 --output research/workstreams/packet8_codesign/iteration7/macro_theta001_alpha2_5.json
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration5/macro_outward.py --local research/workstreams/packet8_codesign/iteration7/near_one_v1/local_theta_1_400.json --alpha 2/5 --output research/workstreams/packet8_codesign/iteration7/macro_theta0025_alpha2_5.json
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration7/scaling_screen.py --powers 18 --q 8 16 --no-q1 --only-extra --extra-macros research/workstreams/packet8_codesign/iteration7/macro_theta003_alpha3_10.json research/workstreams/packet8_codesign/iteration7/macro_theta001_alpha2_5.json --output research/workstreams/packet8_codesign/iteration7/alpha_k18_v1.json
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration7/scaling_screen.py --powers 20 --q 8 16 --no-q1 --only-extra --extra-macros research/workstreams/packet8_codesign/iteration7/macro_theta0075_alpha3_10.json research/workstreams/packet8_codesign/iteration7/macro_theta0025_alpha2_5.json --output research/workstreams/packet8_codesign/iteration7/alpha_k20_v1.json
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration7/scaling_screen.py --powers 18 --q 64 --no-q1 --only-extra --extra-macros research/workstreams/packet8_codesign/iteration5/cache_v1/macro_theta_3_25_alpha_3_10.json --output research/workstreams/packet8_codesign/iteration7/matched_k18_q64_v1.json
C:/Python314/python.exe -B research/workstreams/packet8_codesign/iteration7/summarize_screen.py
```

## Remaining work if this branch is resumed

The K18 evidence supports a targeted coverage study, especially q5..63 and
a few middle/dense occupancies. It does not yet justify launching a complete
certificate. At K20, first optimize the isolated q16 theta/alpha witness:
the present -8.44-bit result is a concrete gap, while q8 has only 3.35 bits
above the nominal target before the full occupancy union is accounted for.
There is no reason from these samples to alter the local construction yet.

A naive full-length replay is substantially more expensive than this screen.
For E G4 epochs, the full-q chronological placement uses
`33 * [32*E*(E-1)/2 + E]` small matrix products per witness. This is 17.6
times the K16 count at K18 and 287.9 times at K20. Scaling the old roughly
180-second, 13-witness complete replay by those counts suggests on the order
of 50 minutes and 14 hours, respectively, before length-dependent overhead.
These are workload estimates, not measured runtimes. Witness-range pruning,
better polynomial assembly, or separate dense-tail bounds should precede a
full K20 run. No such optimization or full replay was attempted here.
