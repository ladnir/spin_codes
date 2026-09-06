# Smaller states and the first Q2 comparison

The first tranche showed little Q1 improvement above s=18. This tranche
extends the same map families to s=12,13,15 and checks Q2 at economical
and larger-state choices. It adds 297 Q1 observations and 36 Q2 observations.
Every new observation uses an exact fixed outer spectrum, except the Q1
random references retained from the first tranche. No random-outer Q2
reference is formed by multiplying expected spectra.

## Smaller states

The new A maps are literal prefixes of the already recorded s=20 generators.
No new quadratic basis is sampled. The audit checks full rank, A^T A=0,
distinct nonzero columns of A^T, exact A spectra, and exact kernel spectra.
Each map receipt records its parent path and hash. Message exponents are
16,18,20, with all eleven outer models from the first tranche.

For exact RM(4,9) [512,256,32], the smallest *tested* state clearing 40 Q1
bits is:

| Message length | t=64 | t=128 | t=256 |
|---|---:|---:|---:|
| k=2^16 | 12 | 12 | 12 |
| k=2^18 | 14 | 13 | 13 |
| k=2^20 | none through 20 | none through 20 | none through 20 |

At k=2^16, the s=12 margins are 41.348, 41.364, and 41.292 bits, respectively.
States below 12 have not been screened. These values therefore identify
useful candidates, not the minimum feasible state or an optimal map.

## Q2 calculation

The Q2 run fixes k=2^16, t in {64,128,256}, and s in {12,18}. It uses
exact BCH [8,4,4], [32,16,8], [64,32,12], [128,64,22], RM(3,7), and RM(4,9).
The two active outer rows use the same fixed constituent. Their coordinate
permutations are independent, even if their local messages are equal.

Conditional on the two outer weights, the row supports are independent
uniform subsets of the B regions. A region containing two active bits has
a uniform two-subset of its L positions. The calculation preserves that
without-replacement law, including same-epoch and different-epoch placements.

`activation_q2.py` constructs three activation-aware epoch transfers for
input weights zero, one, and two. It powers their degree-two matrix
polynomial with binomial support counts. The normalized region transfers
then feed the coefficient polynomial

    (R_0 + x R_1 + y R_1 + xy R_2)^B.

The initial class is zero and the terminal weight sums all three classes.
Dividing each coefficient by choose(B,w1) choose(B,w2) gives its support
average. A separate witness is selected for each pair of outer weights,
before summing with choose(L,2) A_w1(C) A_w2(C). This counts ordered local
messages in a chosen pair of row positions. It does not sample new row codes.

All configurations use the same log-surprisal grid:

    {-8,-7,-6,-5.5,-5,-4.5,-4,-3,-2,-1,0}.

The calculation is nearest binary64 throughout. The finite grid may leave
optimization slack. A failed upper bound does not establish code failure.

The native kernel keeps the three-state contractions explicit, evaluates only
one triangle of the symmetric pair table, and reuses two contiguous buffers.
It performs no allocation inside the coefficient loops. Log-domain evaluation
avoids silently dropping small positive coefficients through underflow.
The kernel is compiled with precise floating-point semantics, without fast
math. Python calls it once per witness; there are no Python callbacks in the
inner recurrence.

All 20 tests pass, including the map-prefix and database witness checks.

The exact tests enumerate region placements and pairs of outer support sets.
They compare the native recurrence to rational arithmetic, check its Q1 axes,
and check row-exchange symmetry. A GF(16) example independently evaluates
actual weight-two activation and following epochs against the envelope.

## Reading the combined results

For RM(4,9) at k=2^16, the completed Q2 margins are:

| State s | t=64 | t=128 | t=256 |
|---|---:|---:|---:|
| 12 | 83.693 | 83.616 | 83.366 |
| 18 | 89.066 | 88.736 | 88.425 |

Each dominant weight pair is (32,32), with log-surprisal witness -5.
For s=12, the Q1+Q2 diagnostic union therefore still has about 41.3 bits
at each epoch size. This is a sparse result, not a full-distance result.

`sparse_comparison.csv` and `sparse_comparison.json` join Q1 and Q2 only when
outer, message length, cutoff, and map identity match exactly. Their combined
margin describes the sum of those two occupation contributions. All
occupations from three onward remain open for these map choices.

The distinction matters already at occupation four. The s=12 maps have
minimum kernel weight four, with 304, 4,960, and 86,720 kernel words of that
weight for t=64,128,256, respectively. The corresponding probabilities use
these counts divided by choose(t,4); raw counts alone should not rank epochs.
Q2 cannot encounter those zero-syndrome inputs. A comfortable Q2 margin
therefore does not settle the higher-occupation trade-off.

The s=18 controls also differ: t=64 has minimum kernel weight six, whereas
t=128 and t=256 still have weight-four kernel words in these unselected
chains. This distinction is a reason to retain multiple epoch choices.

## Reproduction and evidence

The committed first tranche is `88aa11a`. Commit `9d3b7be` corrects its JSON
line-ending rules. A reconstruction from Git blobs under the corrected
attributes matched all 26 dependencies in the refined first-tranche manifest.

From this directory, the new sequential producer commands are:

```text
python run_small_state_pilot.py
powershell -ExecutionPolicy Bypass -File build_activation_q2.ps1
python run_activation_q2_pilot.py
```

The output directories are write-once. To reproduce the numerical experiments,
use a fresh scratch copy of the producers with new output directories and
record its source hashes. The native build script uses the installed MSVC
toolchain; binaries are local build products excluded from Git. The Q2 receipt
records both the C++ source hash and the binary hash used for its run.

Rebuild the index, reports, and tests with:

```text
python build_landscape_db.py
python query_landscape.py export landscape_export.csv
python summarize_activation_pilot.py
python summarize_sparse_comparison.py
python -m unittest -v test_landscape_db.py test_extrapolate_parameters.py test_activation_q1.py test_activation_q2.py
```

The next useful tranche is activation-aware Q3/Q4 and selected middle
occupations for the same s=12 versus s=18 pairs. That will test the kernel
events invisible to Q1/Q2. Lower states can be screened afterwards if the
higher-occupation margin leaves room. BCH-256 remains a secondary
bounded-spectrum comparison, and no implementation benchmark has been run.
