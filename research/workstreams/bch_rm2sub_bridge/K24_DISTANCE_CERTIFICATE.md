# BCH-256 size frontier: the K=2^24 certificate

The full computer-assisted distance certificate at K=2^24 has margin
46.50870389202146 bits. Every occupancy is covered. All six component
certificates passed their512-bit numerical replays, and the exact ledger
`generated/frontier_k24_full_v1.json` passed reconstruction.

## Fixed construction and target

Keep the BCH [256,128] outer code and selected RM2Sub t64_s20 map from
`T64_S20_FULL_CLOSURE.md`. The outer is the specified p37-syndrome
restriction over GF(256), not an unspecified code with the same parameters.
The selected inner map has seed3390111757, state dimension20, and epoch
length64. Its A spectrum and B-kernel spectrum are audited exactly.

Encode K=2^24 message bits as L=131,072 outer rows, giving N=2^25 output
bits. Independently permute each encoded row and each transposed region.
The inner state starts at zero and persists across regions. At epoch i,

\[
Y_i=X_i+Aq_i,\qquad q_{i+1}=\alpha_iq_i+BX_i.
\]

Each alpha_i is an independent uniform nonzero GF(2^20) element,
independent of the permutations. Output precedes the state update, and
there is no final flush. The complete setup theta is shared by all messages.
For every setup the resulting encoder E_theta is injective and has rate1/2.

Let H=floor(N/10)=3,355,443. For the fixed construction and setup above,
the computer-assisted distance theorem gives

\[
\Pr_\theta[\exists x\ne0:\operatorname{wt}(E_\theta(x))\le H]
\le U <2^{-46}<2^{-40}.
\]

Outside this bad-setup event, the code has minimum distance at least
3,355,444 and relative distance greater than0.1. This is a distance claim
under fresh independent multipliers, not a full SPIN security theorem or
a claim about pseudorandom multiplier generation.

## Coverage and arithmetic

Let Z_Q count bad messages with exactly Q nonzero outer rows. The argument
bounds E_theta[Z_Q] for every Q and sums these expectations. It does not
assume independence between messages.

| Occupancies | Certificate method |
|---|---|
| 1 | Four-state refresh transfer and the existing exact BCH weighted inequality |
| 2..1024 | Fixed-weight region polynomials and the adaptive outer-band recurrence |
| 1025..4095 | The same recurrence with separately discovered witnesses |
| 4096..32767 | Continuous-density interval certificate |
| 32768..98304 | Continuous-density interval certificate |
| 98305..131072 | Continuous-density interval certificate with different row probabilities |

`frontier_sparse.py` parameterizes the existing region calculation by L.
It retains the selected map and the certified BCH shell caps. Q1 uses
256-bit Arb and a 512-bit replay with a different polynomial-power
representation. The sparse ranges use the existing directed binary64
scaled recurrence between Arb calculations. Their replay uses512-bit Arb
and does not repeat witness optimization.

Sparse results retain at most80 bits per occupancy to keep receipts small.
Every dense occupancy receives the uniform bound2^-80. Receipt precision
therefore limits the recorded higher-occupancy sum; these are not estimates
of the actual contributions.

All producers and numerical replays verify the 804-file migration manifest.
`audit_frontier_k24.py` checks source hashes, replay provenance, parameters,
exact sums, and disjoint coverage of every integer Q=1,...,131072. The
ledger checks the strict rational inequality U<2^-46. Decimal margins
are summaries of the rational bound.

The Q1 bound has margin46.50871955263324bits. The combined bound on all
higher occupancies has margin62.99997798638864bits. Their exact sum gives
the full margin46.50870389202146bits, leaving more than6.5bits above the
requested40-bit threshold. The dense trees have3148 leaves in total.
All19 selected tests in seven modules passed, including the five new
frontier tests and the frozen arithmetic-component tests.

## How the dense screen becomes an interval certificate

The iid comparison and outer-band cost inequality are derived in
`K30_CLOSURE_PROGRESS.md`, under **A scalable route for dense occupancies**.
They apply with the current value L=131072. The encoder is unchanged;
independent Bernoulli inputs are only an auxiliary bounding distribution.

For each band g, fix a reference probability p_g and certified density
cost Gamma_g. The all-one band has p_g=1 and Gamma_g=1. Let h be the upper
concave hull of (p_g,log Gamma_g). For Q active rows, put x=Q/L. The
mean reference input density r lies between x*min_g p_g and x.

The dense bound pays binomial(L,Q) for the active row positions,13^Q
for band assignments, and (L+1)^256 for comparison of shuffled
Poisson-binomial regions with independent Bernoulli(r) bits. The remaining
cost is exp(Q*h(r/x)) times an inner Chernoff moment.

Fix the Chernoff parameter lambda>0 and put z=exp(-lambda). Let T_j(z)
be the existing three-state epoch transfer for a uniform weight-j input.
For independent Bernoulli(r) input bits, its average is

\[
T(r,z)=\sum_{j=0}^{64}{64\choose j}r^j(1-r)^{64-j}T_j(z).
\]

The challenge is to bound every entry of T throughout an interval of r.
Sampling r does not suffice. A global derivative bound proved inefficient,
so the checker uses the following monotone decomposition.

For one matrix entry, let c_j be an exact nonnegative rational upper bound
on the corresponding entry of T_j(z). Define

\[
a_0=c_0,\quad d_0=0,\qquad
a_j=a_{j-1}+\max(c_j-c_{j-1},0),\quad
d_j=d_{j-1}+\max(c_{j-1}-c_j,0).
\]

Then c_j=a_j-d_j, and both sequences a and d are nondecreasing.
If J_r has distribution Binomial(64,r), a coupling with shared uniform
random variables makes J_r nondecreasing in r. Hence, for r_-<=r<=r_+,

\[
\mathbb E[c_{J_r}]
\le \mathbb E[a_{J_{r_+}}]-\mathbb E[d_{J_{r_-}}].
\]

Arb encloses the two endpoint expectations and their difference. The
decomposition itself uses exact rational arithmetic. An outward upper
endpoint, clamped below at zero, supplies the matrix-entry bound.
At a zero-width interval the expression equals the original Bernstein
polynomial, so no fixed derivative penalty remains.

The checker subdivides a rectangle in integer Q and normalized active-row
density v, where r/x=p_min+(1-p_min)*v and 0<=v<=1. For a box
Q in [a,b] and v in [u,w], it uses

\[
r_-=(a/L)(p_{\min}+(1-p_{\min})u),\qquad
r_+=(b/L)(p_{\min}+(1-p_{\min})w).
\]

It bounds log binomial(L,Q) using L times binary entropy, maximized over
the box's Q interval. It bounds the nonnegative band cost by b times
the maximum of h over the box's active-density interval. The entrywise
transfer bound is raised to N/64=2^19 by positive matrix squaring.

Each accepted leaf bounds the dense inequality's right-hand side by2^-80
for every Q and active-density value in that box. The complete tree
therefore supplies U_Q<2^-80 for every Q in its range.
Integer-Q splits are disjoint and exhaustive. Density splits
share their boundary and cover the entire parent interval. The saved
binary tree specifies all splits and the tilt witness at every leaf.
The 512-bit replay checks the complete tree without optimization.

Leaves are not separate bad events: they cover the supremum over possible
band mixtures. Thus the dense range contributes at most its number of
occupancies times2^-80, not its number of leaves times2^-80.

The refined driver changes only witness discovery. A coarser tilt grid
missed passing witnesses in portions of the range. Its failed attempts
produced no completed certificate for those ranges.

## What this adds to the size–margin curve

The comparison now has two fully certified points, K=2^20 and K=2^24.
The K=2^30 result remains partial; this work does not upgrade its coverage.

| Message length | Fully certified margin | Modeled-spectrum Q1 diagnostic |
|---|---:|---:|
| 2^20 | 50.487298 bits | 70.554413 bits |
| 2^24 | 46.508704 bits | 66.575720 bits |

Across these two certified points, the recorded margin loses0.9946485964
bits per doubling. This is a descriptive two-point slope, not a theorem
about intermediate or larger sizes. The older K=2^20 certificate uses a
slightly stronger cutoff209716; its diagnostic uses floor(N/10)=209715.

`frontier_curve.py` exports the full certificates separately from partial
certificates and the existing modeled-spectrum Q1 diagnostics. It leaves
modeled full margins unset. Those modeled Q1 values concern a hypothetical
outer spectrum, not an experimentally measured failure probability.

A difference of about20 bits between the constraint-based Q1 calculation
and the modeled-spectrum calculation remains an estimator question.
It must not be interpreted as20 additional proved bits. Conversely, fitting
only the constraint-based curve as though it were the true margin could
turn proof looseness into an apparent property of BCH-256.

## Reproduction

Restore the local proof inputs as described in `MIGRATION.md`. Each producer
takes `--output <fresh-path>`; append `--verify` to replay its saved
certificate. Replay writes a new adjacent `_replay.json` file, so it also
requires a destination that does not already exist.

Use `frontier_sparse.py --m 24 --q1` for Q1. Use the same module with
`--first 2 --last 1024` and `--first 1025 --last 4095` for the sparse ranges.
For the dense ranges, use:

```text
frontier_dense_refined.py --m 24 --first 4096 --last 32767
frontier_dense.py --m 24 --first 32768 --last 98304
frontier_dense_refined.py --m 24 --first 98305 --last 131072 --slope 0.25
```

Run these Python scripts sequentially, with an output path for each command.
The retained v1/v2 receipt names appear in `audit_frontier_k24.py`.
Run that audit with a fresh `--output` path after numerical replay; its
`--verify` option reconstructs an existing full ledger without writing.

The new `test_frontier.py` checks size parameterization, exact monotone
decomposition, the binomial interval envelope, Arb enclosure of the exact
expression, and maxima of the rational cost hull. Existing tests cover
the frozen polynomial and scaled-recurrence components.

Next candidate: K=2^26, beginning by inspecting where its dense range
closes. Keep the heuristic estimator work separate
and preserve the full K=2^24 result as a fixed baseline.
