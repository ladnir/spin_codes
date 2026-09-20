# BCH-256 frontier at K=2^26

The full computer-assisted distance certificate at K=2^26 has margin
44.50975536265787bits. Every occupancy is covered. All five component
certificates passed their512-bit numerical replays, and the exact ledger
`generated/frontier_k26_full_v1.json` passed reconstruction.

## Same construction, larger message length

Keep the fixed BCH [256,128] outer and selected RM2Sub t64_s20 map from
`K24_DISTANCE_CERTIFICATE.md`. The exact outer restriction and map seed
remain unchanged. In particular, the inner is not the different nested
map used by the finite-theory parameter study.

The message length is K=2^26=67,108,864 bits. There are L=524,288 outer
rows and N=2^27=134,217,728 output bits. The bad-output cutoff is
H=floor(N/10)=13,421,772.

The setup samples independent row permutations, independent region
permutations, and fresh independent uniform nonzero GF(2^20) multipliers.
The state starts at zero, persists across regions, and is discarded at
the end. Output precedes each update:

\[
Y_i=X_i+Aq_i,\qquad q_{i+1}=\alpha_iq_i+BX_i.
\]

The setup theta is shared by all messages. The encoder E_theta is
injective for every setup and has rate1/2. The computer-assisted theorem is

\[
\Pr_\theta[\exists x\ne0:\operatorname{wt}(E_\theta(x))\le H]
\le U<2^{-44}<2^{-40}.
\]

Outside this setup event the code has minimum distance at least13,421,773,
and therefore relative distance greater than0.1. The claim concerns code
distance in the fresh-multiplier model, not the full SPIN protocol or a
pseudorandom multiplier stream.

## What is reused, and what is recomputed

`frontier_sparse.py` and the two dense interval checkers remain byte-for-byte
unchanged from the K=2^24 checkpoint. The dense argument still covers the
continuous range of band mixtures, not a sampled approximation.

| Occupancies | New-size calculation |
|---|---|
| 1 | Four-state refresh transfer and the existing exact BCH weighted inequality |
| 2..4095 | Sparse witnesses transported from K24, then recomputed and refined where needed |
| 4096..131071 | Dense interval certificate, row-probability slope0.75 |
| 131072..393216 | Dense interval certificate, row-probability slope0.75 |
| 393217..524288 | Dense interval certificate, row-probability slope0.25 |

Transport changes witness discovery only. For a change from exponent m0
to m, the initial log-tilt witness is shifted by approximately
-(m-m0)*log(2). The auxiliary row probabilities are initially retained.
Every region polynomial, counting factor, cutoff correction, and final
upper bound is recomputed at the new L and H. No old probability bound
is substituted for a new-size probability.

The first transported Q2 witness gave only about68.94bits, below the
sparse checker's70-bit acceptance threshold. A new witness gives about
75.15bits. The driver also refines any other uncovered occupancies before
writing a complete sparse certificate.

Q1 is outward-certified at44.5097710391417bits and passed its512-bit replay.
The sparse certificate covers every Q2..4095. The dense certificates cover
every Q4096..524288, using5501,331, and895 leaves in the three ranges.
No range gap is filled by interpolation or an observed curve.

The exact sum of all higher-occupancy bounds has margin
60.99956804425319bits. Combining it with Q1 gives the full margin
44.50975536265787bits, more than4.5bits above the requested40-bit target.
All25 selected tests in nine modules passed, including the compact-ledger
regressions and the separate density-factor tests.

## Compact full-coverage audit

`frontier_ledger.py` verifies the component receipts and their512-bit
replays. It checks source hashes, exact sums, parameters, and the dense
tree structure. The46 local shell-cap inputs must be exactly the set
authenticated by the migration manifest; extra unauthenticated caps
are not admitted.

The ledger stores ranges rather than a dense dictionary with one object
per occupancy. Adjacent intervals must form an exact partition of1..L.
Its memory use therefore follows the receipt size, not the total number
of dense occupancies. A regression test reproduces both rational sums
from the frozen K24 ledger exactly.

Sparse certificates retain at most80bits per occupancy, and every dense
occupancy receives the recorded bound2^-80. The resulting higher-occupancy
sum is limited by that deliberate receipt rounding. It must not be read
as an estimate of the true higher-occupancy contribution.

## Frontier interpretation

There are now three fully certified points:

| Message length | Full certified margin | Modeled-spectrum Q1 diagnostic |
|---|---:|---:|
| 2^20 | 50.487298 bits | 70.554413 bits |
| 2^24 | 46.508704 bits | 66.575720 bits |
| 2^26 | 44.509755 bits | 64.576792 bits |

The certified margin loses0.9946485964bits per doubling from K20 to K24,
and0.9994742647bits per doubling from K24 to K26. These are descriptive
finite differences, not a theorem about larger sizes. The original K20
certificate uses a slightly stronger distance cutoff than floor(N/10).

`frontier_curve_extend.py` records the completed full point and recomputes
the modeled-spectrum Q1 value at K26 in `generated/frontier_curve_v2.json`.
Full modeled margins remain unset.
In particular, the rounded bound on the higher occupancies can exceed the
modeled Q1 contribution. That is not evidence that dense messages dominate
the actual failure probability; the bound has discarded most of their
certified slack for compactness.

Before treating modeled Q1 as a full-margin estimate, tighten that recorded
higher-occupancy envelope and assess the outer low-shell spectrum model.
The approximately20-bit separation between modeled and constraint-based
Q1 curves remains a heuristic-estimator question, not additional proved
margin.

## Reproduction and next step

Restore inputs with `MIGRATION.md`. Run proof jobs sequentially and use
fresh output filenames; producers and replay receipts are write-once.

For Q1, run `frontier_sparse.py --m 26 --q1 --output <path>`, followed by
the same output path with `--q1 --verify`. For the sparse range, run
`frontier_transport.py --m 26 --inputs <K24-first-range> <K24-second-range>
--output <path>`. The two source receipts cover2..1024 and1025..4095.
Replay the result with `frontier_sparse.py --output <path> --verify`.

For each dense range in the table, run `frontier_dense_refined.py --m 26
--first <lower> --last <upper> --slope <value> --output <path>`, then replay
with the same module's `--output <path> --verify` command.

Finally, call `frontier_ledger.py --m 26 --certificates <five-producer-paths>
--output <ledger-path>`. Its `--output <ledger-path> --verify` mode
reconstructs an existing ledger without writing. The curve extension takes
the completed ledger and the earlier `frontier_curve_v1.json`.

Next, probe where the dense certificate begins at K=2^28 before
committing to another complete sparse calculation. Preserve the K24 and
K26 baselines, and continue the spectrum-model assessment separately.
`POISSON_DENSITY_REFINEMENT.md` derives a smaller comparison factor for
that next step. It is not used in the K26 certificates.
Check Q2 before starting a large sparse search. The inherited70-bit
per-class acceptance rule is stronger than the final40-bit union target;
missing that rule is not itself evidence of a bad code. If necessary,
use the refreshed small-occupancy transfer or a verifier that retains
a sufficient exact bound under an explicitly checked union budget.
