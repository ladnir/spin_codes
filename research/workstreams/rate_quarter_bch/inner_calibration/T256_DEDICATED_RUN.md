# Dedicated t256 run: a deterministic distance ceiling

The dedicated search exposed a structural limitation, not merely loose
numerical witnesses. At K=2^20 with the [128,32] outer, **every t256_s18 setup
has a nonzero output of weight at most 589,888**, or 14.06403% of N=2^22.
It cannot meet either the 16.5% or 19% target. This applies to every inner
map with the stated dimensions and interface, regardless of kernel distance.

The same argument requires s>=22 for 16.5%, and s>=25 for 19%. These are
necessary conditions, not achievable-margin claims. The existing t128_s19
numerical baseline is unaffected.

## A zero-state subspace must exist

Fix any setup. The outer is a binary linear [b,d] code used on L independent
rows. Its output undergoes coordinate permutations before the inner. Write
N=bL and E=N/t, assuming whole epochs. The inner has the interface

    q_0=0,   Y_i=X_i+A_i q_i,   q_(i+1)=alpha_i q_i+B_i X_i,

where each B_i has at most s binary output coordinates. The actual fixed
RM2Sub pair is a special case. The argument needs neither BA=0 nor randomness.

Restrict messages to any Q outer rows. Before the inner, they form a subspace
of dimension dQ, supported on at most bQ coordinates. Coordinate permutations
preserve this support size. Impose B_i X_i=0 in every epoch. There are at
most sE scalar linear equations, so the solution subspace has dimension

    h >= r := dQ-sE.

When r>0 this subspace is nonzero. On it, induction gives q_i=0 and Y_i=X_i
for every epoch. Its outputs therefore remain supported on at most bQ positions.

Every nonzero coordinate functional on a binary h-dimensional subspace is
one on exactly 2^(h-1) words. Averaging over its 2^h-1 nonzero words gives

    d_min <= floor(bQ * 2^(h-1)/(2^h-1))
          <= floor(bQ * 2^(r-1)/(2^r-1)).                       (1)

The second inequality uses h>=r and the decreasing fraction. The low-weight
message can depend on the setup: minimum distance quantifies over all nonzero
messages after fixing the setup. A target exceeding (1) consequently fails
for every setup, not just with a probability inferred from a first moment.

## Quarter-rate parameters

For b=128, d=32, L=32768, and t=256, there are E=16384 epochs. Choose
Q=512s+1, giving r=32. For the states below, (1) simplifies exactly to

    d_min <= 64Q,   d_min/N <= s/128 + 1/65536.

| State s, with t=256 | Q | Exact d_min upper bound | Relative ceiling |
|---:|---:|---:|---:|
| 18 | 9217 | 589888 | 14.06403% |
| 19 | 9729 | 622656 | 14.84528% |
| 20 | 10241 | 655424 | 15.62653% |
| 21 | 10753 | 688192 | 16.40778% |
| 22 | 11265 | 720960 | 17.18903% |
| 24 | 12289 | 786496 | 18.75153% |
| 25 | 12801 | 819264 | 19.53278% |

This ceiling depends on outer dimension and support, not its weight spectrum.
Changing maps or permutations cannot evade it. A pre-inner linear transform
beyond coordinate permutations would require a new support argument.

## Numerical run and diagnosis

The fresh search used log tilts -8,...,1 in steps of 0.25, optimized proposals,
and new partitions over zero, nontrivial, and all-one outer rows. It did not
use the baseline's boxes or witnesses as seeds.

The initial 1023-node run reached its budget at 16.5%, but a replay assertion
rejected a valid improvement: tightening redundant box coordinates decreased
the lattice-count penalty. Replay now checks equality of the represented
simplex vertices and explicitly accounts for that penalty change.

After finding the obstruction, the retained run used 255 nodes per target.
Both complete dense covers give useless upper bounds: logarithmic margin
diagnostics -234172.28 bits at 16.5% and -311333.25 bits at 19%. These numbers
are not failure probabilities. More partition search cannot beat (1).

Singleton checks also failed at intermediate occupancies, without box-width
loss. At the half-active zero/nontrivial-row point, the relaxed transfer was
dominated by its all-zero-state path. Its tilted epoch weights were about
42--49, with negligible weight-four contribution. This supersedes the earlier
suggestion to prioritize eliminating weight-four words as the remedy for t256.

## Independent first-moment certificates

Before deriving (1), we adapted the existing zero-state lower-bound argument
to the exact weight-52 outer shell, with multiplicity 126,946,176. Exact
integer checks prove the following lower bounds on expected bad-message counts:

| Fixed map | At 16.5% | At 19% |
|---|---:|---:|
| t256_s18_nested | >=2^47256 | >=2^101241 |
| t256_s19_nested | >=2^30872 | >=2^84857 |
| t256_s20_nested | >=2^14488 | >=2^68473 |

For even Q and even shell weight, Fourier inversion gives probability at
least 2^-127 that all 128 region counts are even. Exact rational Chernoff
comparisons leave probability at least 2^-128 inside an additional count
interval. For each allowed even regional weight j, an explicit coefficient
term in the kernel polynomial power K(u)^128 gives

    beta_j >= choose(128,h)*k_a^(128-h)*k_(a+2)^h / choose(32768,j),
    a=2*floor(j/256),  h=(j-128a)/2.

Integer comparisons establish beta_j>=2^-R throughout the interval.
Conditional independence of region permutations then gives a zero-state
probability lower bound 2^(-128-128R). Multiplying by choose(L,Q)*A_52^Q
proves the table. Floating-point calculations only propose intervals and
exponents; every bound used in the certificate is checked with integers.

These first-moment results alone would not prove frequent setup failure.
The separate subspace argument (1) does prove failure for every setup at
the ruled-out targets. The two arguments have different logical strength.

## Reproduction and next step

Run from the repository root with the existing NumPy/SciPy environment:

```sh
python workstreams/rate_quarter_bch/inner_calibration/dedicated_dense.py --nodes 255
python workstreams/rate_quarter_bch/inner_calibration/probe_dense_points.py
python workstreams/rate_quarter_bch/inner_calibration/diagnose_dense_paths.py
python workstreams/rate_quarter_bch/inner_calibration/zero_state_obstruction.py
python workstreams/rate_quarter_bch/inner_calibration/zero_state_rank_ceiling.py
python -m unittest discover -s workstreams/rate_quarter_bch/inner_calibration -p 'test_*.py'
```

`ZERO_STATE_RANK_CEILING.json` retains the exact finite-size ceilings;
`ZERO_STATE_OBSTRUCTION.json` retains the first-moment certificates. Numerical
searches and diagnostics have separate receipts. No numerical upper bound
was promoted to an outward certificate. No encoder or benchmark was changed.

Next, retain t=256 but increase state. Screen s=22--24 for 16.5% and s=25--28
for 19%, applying the rank ceiling and zero-state lower tests before expensive
upper-bound searches. These choices still halve the epoch count relative to
t128_s19, but a runtime advantage must be measured.
