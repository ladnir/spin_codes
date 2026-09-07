# BCH-256 with RM2Sub: current understanding and next decisions

Status: 2026-09-07. This is the current handoff for the BCH-256 workstream.
Historical notes retain the claims and parameter choices made at their dates.

The immediate question is which inner parameters can deliver a 40-bit
distance/setup margin at message lengths K=2^16, 2^18, and 2^20.
The selected t64_s20 map is our proof-backed baseline. The selected
t128_s19 map is the next performance candidate, not another certified choice.

## Construction and meaning of margin

The outer is our fixed binary BCH [256,128] constituent, with minimum distance
at least 38. It is the five-dimensional p37-syndrome restriction described in
[the finite theorem](T64_S20_FULL_CLOSURE.md), not an arbitrary code with these dimensions.
Its exact weight spectrum remains unknown.

For K message bits, encode L=K/128 outer rows, producing N=2K bits.
The setup independently permutes each row and each of the 256 transposed regions.
The inner processes t-bit epochs with an s-bit state, identified with GF(2^s):

    q_0 = 0
    Y_i = X_i + A q_i
    q_(i+1) = alpha_i q_i + B X_i

Here A and B are fixed selected maps with B=A^T and BA=0.
Each alpha_i is a fresh independent uniform nonzero field element.
Multipliers and permutations are mutually independent. State persists across
regions; output precedes update, and there is no final flush.
One sampled setup is shared by all messages.

Let Z_Q count nonzero messages with Q occupied outer rows and output weight
at most H. The proof bounds E_setup[Z_Q] for every Q and sums these bounds.
For an exact upper bound U on the sum, the margin is -log2(U).
Thus a margin of at least 40 bounds the probability that any such message
exists by 2^-40. It is not the total security level of the SPIN protocol.
Replacing fresh randomness with a seeded implementation remains a separate obligation.

New diagnostics use H=floor(N/10). The historical K=2^20 full theorem uses
the slightly stronger cutoff H=209716 rather than floor(N/10)=209715.

## Parameter shortlist at the requested message lengths

The following values are Q1-only nearest-binary64 diagnostics. They use
our certified BCH spectrum constraints through `bridge.bch_bound`, not
the conditioned-binomial spectrum model. They are not outward certificates.

| Selected map | K=2^16 | K=2^18 | K=2^20 | Current role |
|---|---:|---:|---:|---|
| t64_s20 | 54.156357 | 52.422103 | 50.487812 | Proof-backed baseline at K=2^20 |
| t128_s19 | 54.010387 | 52.360418 | 50.448302 | Next performance candidate |
| t64_s16 | 53.698685 | 51.979239 | 50.048776 | Smaller-state exploratory option |
| t128_s15 | 53.186269 | 51.553698 | 49.646891 | Known first-moment obstruction at K=2^20 |
| t256_s14 | 52.343951 | 50.844814 | 48.980155 | Q1-only alternative; not shortlisted for proof |

`parameter_shortlist_q1.py` reproduces this comparison serially. It uses the
refresh transfer, coarse tilts log(a)=-4,-3.9,...,8, and local refinements
around the optima for weights 38,40,...,50. The tilt is lambda=a/L.
The map files are identified by the manifests in `inputs/` and
`generated/larger_state_inputs_v1/`; the script checks those manifests on load.
These dimensions do not certify every map with the same t and s.

The full K=2^20 certificate for t64_s20 is **50.4872982775 bits** after the
[Q1 refresh](DUAL_TRACK_PROGRESS.md). Every one of its 8192 occupancies is covered.
No full certificate has yet been produced here for K=2^16 or K=2^18.
The K=2^20 theorem must not be silently transferred to those shorter encoders.

For all three requested sizes, retain t64_s20 as the first choice to certify.
Then test t128_s19: its Q1 margin is within 0.15 bits of t64_s20 at each size.
It uses half as many state updates, but neither encoder speed nor a full
distance bound has been established for this candidate.

Smaller state should not be selected by subtracting the spare Q1 margin.
Higher occupancies impose a different constraint. In particular, the actual
bad-message expectation for t128_s15 exceeds 2^19600 at Q=2620 and K=2^20.
[The obstruction](ZERO_STATE_FIRST_MOMENT_OBSTRUCTION.md) rules out an
unconditional first-moment closure there; it does not establish large setup-failure probability.
The analogous status of another map or message length requires its own calculation.

## Proved frontier versus the heuristic curve

For the same selected t64_s20 map, the completed all-occupancy bounds are:

| log2(K) | Full certified margin, approximately | Record |
|---:|---:|---|
| 20 | 50.487298 | [Refresh and aggregate](DUAL_TRACK_PROGRESS.md) |
| 24 | 46.508704 | [K24 certificate](K24_DISTANCE_CERTIFICATE.md) |
| 26 | 44.509755 | [K26 certificate](K26_DISTANCE_CERTIFICATE.md) |
| 28 | 42.510028 | [K28 certificate](K28_DISTANCE_CERTIFICATE.md), with retained dense precision |

The K=2^30 work covers only a sparse range; it is not a full 40-bit certificate.
The user has indicated that K=2^28 is already sufficiently large for this exploration.

The certified curve is a conservative guarantee, not a measurement of the
unknown BCH spectrum. The conditioned-binomial reference gives Q1 margins
of about 74.22 bits at K=2^16, 70.55 at K=2^20, and 62.58 at K=2^28.
The roughly 20-bit gap is unresolved spectral uncertainty in these bounds;
it is not 20 bits of established additional security.

[Matched-spectrum backtests](CURVE_SPECTRUM_ASSESSMENT.md) compare exact and
reference spectra for smaller constituents. The largest model overstatement
in the 16 tested cases is below 0.949 bits. This is not a confidence interval
or an error guarantee for BCH-256. Repeated message lengths are not independent
samples of constituent codes. The relevant estimator quantity is the
transfer-weighted spectrum score, not a blind fit through loose BCH-256 caps.

Keep two tracks separate: improve unconditional BCH/inner bounds, and build
evidence for an explicitly labeled reference-spectrum estimate.

## What the bounded low-shell investigation established

The recent search targeted additional BCH facts, not new inner parameters.
[The ROI calculation](LOW_SHELL_ROI.md) found that reoptimizing the unchanged
LP could gain at most about 0.002979 bits with current transfer coefficients.
Material improvement therefore needs stronger constraints.

Halving the current A38 cap would suffice for about 0.969 additional bits
at K=2^28. This is a conditional payoff, not a newly proved cap.
Proving that P dual has no weight-30 words would yield about 1.075 bits.
That stronger code statement remains open.

The [rank investigation](PDUAL_RANK_PROGRESS.md) did prove that the subcode
E={c in P dual : F7(c)=0} has dimension 117 and minimum distance at least 32.
Here F7 is the Fourier coordinate defined in that note. Potential weight-30
words of P dual must lie outside E. The remaining search did not prove absence:
SAT timed out, and randomized searches found weight-36 witnesses but no weight-30 witness.
Search counts do not provide sampling confidence for absence.
[The case-7 note](PDUAL_CASE7_SEARCH.md) records these limitations and the
exact degree-five Boolean description of the remaining coset.

Adding the proved E/coset constraints to the LP gave only 0.004507 bits.
An exact feasible witness limits this augmented relaxation's gain to less
than 0.005635 bits with the current coefficients fixed.
[That route is parked](KERNEL_COSET_PAYOFF.md). We did not replace the
frozen 42.510028-bit aggregate for this negligible gain.

## Integration, reproduction, and next work

GitHub `origin/main` is authoritative. Overleaf is not the integration target.
This workstream owns `workstreams/bch_rm2sub_bridge/`; the separate finite-theory
workstream owns its landscape code. Its nested RM2Sub maps differ from the
selected maps here, so matching t and s does not transfer a certificate.
This handoff does not assert a fresh audit of every incoming landscape result.

The commit preserves compact source code, tests, and notes. Bulk generated
LPs, executables, and search outputs remain local and ignored. Some replays
require those retained local receipts and the inputs described in
[MIGRATION.md](MIGRATION.md). This is not a new self-contained proof bundle.
Existing hash-bound sources and certificate files are left unchanged.

From this directory, the targeted regression suite is:

```text
python -B -m unittest test_kernel_coset_payoff test_pdual_case7_search test_pdual_refinement test_low_shell_roi test_curve_spectrum test_refresh_q1
python -B kernel_coset_payoff.py --verify
```

Before integration, all 21 tests passed, including byte-for-byte verification
of the 804 migrated files. The exact coset replay passed without a solver run.
The new shortlist script reproduced 52.422103 bits for t64_s20 at K=2^18.
The repository hygiene check passed. These checks are not a new full certificate
for that message length or a rerun of every historical numerical calculation.

Recommended next work, in order:

1. Check all occupancies for the selected t64_s20 map at K=2^16 and K=2^18.
2. Screen difficult higher occupancies for t128_s19 at the three requested sizes.
3. Certify passing candidates, then benchmark them serially before ranking speed.
4. Continue weighted low-shell bounds and reference-score calibration when a
   stronger counting fact or independent evidence offers meaningful payoff.

Do not resume the parked coset search merely by increasing its time limit.
