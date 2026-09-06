# BCH-256 frontier at K=2^28

The full computer-assisted distance certificate has margin
42.510012787619644 bits. All 2,097,152 occupancies are covered, all seven
component receipts passed their 512-bit replays, and the exact ledger
`generated/frontier_k28_full_v1.json` passed reconstruction. This leaves
approximately 2.510013 bits above the requested 40-bit setup-failure target.

## Construction and probability space

Keep the fixed BCH [256,128] outer and selected RM2Sub t64_s20 map from
`K26_DISTANCE_CERTIFICATE.md`. This is the same outer restriction and
the same inner map, not the different nested map used by the finite-theory
parameter study. No exact BCH-256 spectrum is assumed.

The selected map has seed 3390111757 and selection-file SHA-256
`23219aa244e3fd6b426f5e5474e3824dff2570459853ed8d3a257a53f408d0f5`.
Here A maps F_2^20 to F_2^64, B=A^T, BA=0, rank(A)=20,
the image of A has distance 16, and ker(B) has distance 8.

Here K=2^28=268,435,456 message bits, L=K/128=2,097,152 outer rows,
N=2K=536,870,912 output bits, and H=floor(N/10)=53,687,091.
An occupancy Q is the number of nonzero outer rows in a message.

The setup theta samples independent row permutations, independent region
permutations, and fresh independent uniform nonzero GF(2^20) multipliers.
The setup is shared by all messages. The state starts at zero, persists
across region boundaries, and is discarded at the end. The epoch rule is

\[
Y_i=X_i+Aq_i,\qquad q_{i+1}=\alpha_iq_i+BX_i.
\]

The output precedes the update. The encoder is injective for every setup
and has rate 1/2. The target setup event is

\[
\mathcal B:=\{\theta:\exists x\ne0,
                  \operatorname{wt}(E_\theta(x))\le H\}.
\]

The calculation bounds Pr_theta[B] by summing first-moment upper bounds
over occupancies. It does not assume that different messages have
independent outputs. Outside B, minimum distance is at least 53,687,092,
so relative distance is greater than 0.1.

The verified rational union bound U satisfies

\[
\Pr_\theta[\mathcal B]\le U<2^{-42}<2^{-40}.
\]

The displayed decimal margin is -log2(U), not an estimate of the actual
failure probability. The strict power-of-two inequalities above are
checked with exact rational arithmetic.

## Results

| Occupancies | Certificate | Dense slope | Dense leaves |
|---|---|---:|---:|
| 1 | Refresh-aware Q1 and the existing exact BCH weighted inequality | — | — |
| 2..7 | Four-state refresh transfer with kernel exclusion | — | — |
| 8..6143 | Recomputed sparse witnesses and frozen-verifier replay | — | — |
| 6144..8191 | Square-root density comparison | 0.75 | 3042 |
| 8192..524287 | Square-root density comparison | 0.75 | 7984 |
| 524288..1572864 | Square-root density comparison | 0.75 | 331 |
| 1572865..2097152 | Square-root density comparison | 0.25 | 890 |

Q1 alone is bounded at 42.510028466018866 bits. The exact sum of the
recorded higher-occupancy bounds has margin 58.999649198045745 bits.
Their union gives the full margin above. The higher-occupancy number
includes deliberate upward rounding in the receipts; it is not a
measurement of the underlying failure probability.

The new curve file, `generated/frontier_curve_v3.json`, contains four
fully certified points:

| Message length | Full certified margin | Modeled-spectrum Q1 only |
|---|---:|---:|
| 2^20 | 50.487298 bits | 70.554413 bits |
| 2^24 | 46.508704 bits | 66.575720 bits |
| 2^26 | 44.509755 bits | 64.576792 bits |
| 2^28 | 42.510013 bits | 62.577055 bits |

From K26 to K28, the certified margin loses 0.9998712875191131 bits per
doubling. This is a finite difference between certified points, not a
uniform theorem over all message lengths. K20 retains its original
slightly stronger cutoff. K30 remains partial, covering Q1..Q32 only;
this checkpoint does not claim a full K30 certificate.

## Computation

The new dense checker uses the pointwise comparison in
`POISSON_DENSITY_REFINEMENT.md`. For a shuffled region of length L,
the density factor is

\[
D_L=\min\{L+1,4\lceil\sqrt L\rceil\}=5796.
\]

For fixed auxiliary band probabilities, each region input is a uniformly
shuffled vector of independent, possibly unequal Bernoulli variables.
The lemma compares its entire vector distribution to iid Bernoulli
inputs of the same mean. Conditional on these auxiliary choices, the
256 region comparisons give D_L^256. The row-density costs and mixture
envelope remain those of the earlier dense argument.

`frontier_dense_density.py` changes only that density cost in the
outward interval bound and uses the existing refined tilt search.
The interval tree covers all real band-mixture densities, including
endpoints, and all integer occupancies in its stated range. A numerical
grid does not fill any gaps. Every accepted leaf bounds each of its
occupancies by 2^-80; leaves partition or overlap the auxiliary parameter
space and are not additional probability events to sum.

The Q1 producer is unchanged. Q2 through Q7 use the parameterized
four-state refresh transfer in `frontier_refresh_sparse.py`. Since the
selected kernel has distance 8, a nonempty epoch of weight at most 7
cannot have zero syndrome. Q2 is outward-bounded at approximately
71.557979 bits. This check avoids relying on the older three-state
bound's small-occupancy witness policy.

The intermediate sparse range uses `frontier_transport_checkpoint.py`.
It transports only auxiliary witnesses from K26; it recomputes all
region polynomials, counting factors, and cutoff corrections at the new
size. It stops each new anchor search once an outward bound satisfies
the existing 70-bit acceptance policy. Write-once discovery checkpoints
allow resumption. Those checkpoints are not certificates and are not
part of the committed evidence. The frozen `frontier_sparse.py` replays
the final sparse receipt independently of the discovery driver.

Three initial attempts to begin dense coverage at Q=4096 did not close:
row-probability slopes 0.75, 0.5, and 1 failed their interval checks.
These are failures of the selected bounds and witnesses, not lower
bounds on code failure. Dense coverage starts at Q=6144 instead; the
sparse calculation must cover every lower occupancy.

## Verification and reproduction

All producers use 256-bit Arb bounds. Replays use 512-bit Arb; sparse
directed binary64 arithmetic remains outward at both precisions.
Q1 and the refresh-aware small occupancies use separate polynomial
powering orders in their replays. Dense replay evaluates every saved
leaf without rerunning witness optimization.

The exact ledger checks source hashes, component replay hashes,
parameters, rational sums, dense tree structure, and an exact partition
of 1..L. It also requires that the 46 shell caps are exactly the set
authenticated by the migration manifest. It does not certify a range
by interpolating between known sizes.

All 29 tests in eleven selected modules passed. These include exact
density-factor checks, comparison against the frozen dense bound,
refresh-receipt rounding checks, ledger regressions, and a curve-upgrade
test that prevents replacing an existing full certificate. The frozen
K26 ledger also passed reconstruction without changes.

Restore the frozen inputs using `MIGRATION.md`. Run numerical jobs
sequentially, from this directory, and use fresh output filenames.

1. Run `frontier_sparse.py --m 28 --q1 --output <q1>` and replay with
   `--q1 --verify --output <q1>`.
2. Run `frontier_refresh_sparse.py --m 28 --output <small>` and replay
   with the same module's `--verify --output <small>`.
3. Run `frontier_transport_checkpoint.py --m 28 --first 8 --last 6143
   --source generated/frontier_k26_q2_q4095_v1.json --output <sparse>`.
   Replay with `frontier_sparse.py --verify --output <sparse>`.
4. For each dense range in the results, run `frontier_dense_density.py
   --m 28 --first <lo> --last <hi> --slope <slope> --output <dense>`,
   then use its `--verify --output <dense>` mode.
5. Aggregate all seven producers with `frontier_ledger.py --m 28
   --certificates <paths> --output <full>`, then reconstruct the ledger
   with `--verify --output <full>`.
6. Run `frontier_curve_upgrade.py --previous generated/frontier_curve_v2.json
   --certificate <full> --output <curve>` to add the newly certified point
   and recompute its separately labeled Q1 model diagnostic.

The regression command is:

```powershell
python -B -m unittest test_frontier_density test_frontier_curve_upgrade test_poisson_density_factor test_frontier_ledger test_frontier test_refresh_q1 test_k30_certificates test_polynomial_regions test_scaled_adaptive test_positive_line_hull test_gap_dyadics
```

`frontier_curve_upgrade.py` upgrades the previously diagnostic-only K28
point in a new curve file and preserves its prior contents. It refuses
to replace an existing full certificate. Its spectrum-model output is
Q1 only: no modeled full margin or statistical confidence interval is
claimed. The independent-multiplier distance result is not a complete
SPIN protocol-security theorem, a PRG-stream theorem, or a performance
measurement.

## The two-track interpretation

The dense refinement helps prove coverage, but it does not remove the
outer-spectrum loss in Q1. The certified frontier and the modeled Q1
frontier therefore remain distinct evidence. A near-unit slope at the
certified points is descriptive evidence about this construction, not
a theorem about every intervening or larger message length.

The current receipts deliberately round most higher-occupancy bounds
up to 2^-80 each. At this L their union is consequently near 2^-59,
even when the underlying inequalities are much stronger. That rounded
union can exceed the modeled Q1 contribution. It cannot establish that
modeled Q1 approximates the entire failure bound.

A useful next step is to retain tighter higher-occupancy bounds while
keeping the same construction. This would test whether the modeled Q1
curve controls the full bound *conditional on the spectrum model*.
It would not validate that spectrum model or turn its approximately
20-bit improvement over the constrained Q1 calculation into a theorem
about the fixed BCH code. Assessing the low-shell spectrum remains a
separate estimator question.
