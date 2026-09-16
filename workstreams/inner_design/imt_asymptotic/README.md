# IMT asymptotic investigation

The current weight-five IMT inner now has a **complete proof draft at
11% asymptotic relative distance**, with all numerical gates passing.
The [11% update](d11/PROOF_UPDATE.md) tightens the outer bound and supplies
a new dense certificate. The preserved [base draft](ASYMPTOTIC_PROOF.md)
supplies the sparse routing comparison, fixed-occupancy limit, and uniform
union argument; its sparse inequalities already used 11%.
The subsequent [in-session analytic review](d11/PAPER_REVIEW.md) found no
unresolved obstruction, and the manuscript now uses IMT for its 11%
asymptotic theorem. This is not independent external review or formal proof
verification. The encoder default and finite RM2Sub claims are unchanged.

The inner uses the same fixed balanced A and `weight5_seed0` feedback B as
the finite half-rate certificate: t=128, s=19, and one sampled transvection
per epoch. The asymptotic outer is the paper's shared Golay--BA-3 constituent,
not the fixed BCH-256 outer. This investigation retains the randomized
bit-transpose permutation and the certified concave outer-spectrum majorant.

All new JSON records are local and uncommitted. The code needs the existing
map data and outer-majorant files; this directory alone is not a complete
evidence bundle.

## Current results

Let L be the number of outer rows, Q the number of nonzero rows, alpha=Q/L,
and x their mean relative weight. The outer constituent length is b, so N=Lb.

| Regime | Result | Remaining obligation |
|---|---|---|
| alpha in [10^-4,1], x in [0.104,0.896], distance 0.11 | Refined outer bound; all 1,023 boxes pass at 256 bits and replay at 512 bits | Review the imported route reduction and its integration |
| Same domain, distance 0.1099 | Earlier 658-box certificate retained unchanged | Historical baseline; superseded in distance by the refined result |
| Q>=4096, alpha<=10^-4 | Exact rational polynomials prove contraction 1-96 alpha on the whole interval; the marked-input comparison gives a uniform summable bound | Review the analytic conditioning and union argument |
| Every fixed Q>=3 | 62 outward endpoint checks pass, retaining b=(39/4) log2(N)+O(1); the draft derives the actual full-state continuum limit | Review the limit and product-norm argument |
| Q=1,2 | New outward inequalities pass at 256 and 512 bits; the draft uses the same full-state limit | Review the limit and fair-row comparison |

`d11/ASSEMBLY_FINAL.json` assembles the 11% update, including a fresh replay
of all 5,363 outer-support leaves and authentication of the dense replay.
The earlier `d11/ASSEMBLY.json` predates a decimal correction in the update
document and is superseded by this final binding.
`ASSEMBLY_D1099.json` authenticates the base source bindings and the
31 imported outer dependencies. Its producer reconstructs the sparse
polynomials, verifies dense coverage against the retained replay, and
recomputes the likelihood, fixed-Q, and cutoff inequalities at 512 bits.
It binds the proof draft and tests. The aggregate is an evidence check, not
a formal proof verifier. The local receipts retain their original scope.

## Positive occupancy: base certificate

The screen substitutes IMT into the paper's existing exponent. For each
box, fix reference parameters p,y in (0,1) and lambda>0. Set beta=py and
z=exp(-lambda). The exponent per output bit is

```text
alpha * a_hat_BA(x)
  + KL(alpha || p) + alpha * KL(x || y)
  + log(rho(T(beta,z))) / 128 + delta * lambda.
```

Here a_hat_BA is the existing certified concave majorant, in natural units.
It bounds mixtures of row weights; this is not a single-weight-class guess.
The screen chooses among complete IMT transfers, never an entrywise minimum:

1. The binomial mixture of the fixed-weight envelope used in the finite proof.
2. The independent-map Fourier envelope used in that proof.
3. A state-independent scalar moment bound, described below.

`screen.py` implements vectorized binary64 versions of the first two bounds.
Regression tests compare them against the existing implementations. Discovery
uses numerical spectral radii. `certify_imt_dense.py` instead reconstructs
the maps and outward transfers, then checks positive rational Collatz vectors.
The vectors are discovery proposals; no floating-point eigenvalue is trusted
by the certificate.

For a fixed witness and one affine outer support, the exponent is convex in
(alpha,alpha*x). Four vertex checks therefore cover the box. The verifier
checks containment, nonoverlap, and total area using exact rational arithmetic
for every outer-majorant segment. It does not trust the discovery coverage flag.

The retained 10.99% cover uses 311 occupation-transfer boxes, 283 Fourier
boxes, and 64 scalar-bound boxes. Its largest outward exponent endpoint is
approximately -4.3631360128e-7. This is an asymptotic exponent per output bit,
not a finite security margin. The records are `DENSE_OUTWARD_D1099_v3.json`
and `DENSE_REPLAY_D1099_v3.json`. Earlier v1/v2 records are superseded by the
IMT-specific verifier module and final source binding.

### Uniform-input bound and the resolved 11% gap

For a fixed entering state q and iid Bernoulli-beta input X, let w=wt(Aq).
The one-epoch output moment is exactly

```text
g0^(128-w) * g1^w,
g0 = 1-beta+beta*z,    g1 = beta+(1-beta)*z.
```

Taking the maximum over w in {0,48,56,64,72,80} gives a scalar envelope.
At beta=1/2, the expression is independent of q. Consequently the full
inner moment is exactly ((1+z)/2)^N. Equivalently, the invertible inner maps
a uniform input to a uniform output for every realized setup.

At alpha=1,x=1/2, this bound yields h(0.11)-ln(2)/2, approximately
-5.82533603e-5. It removes the apparent failure at the center of the initial
54-point grid. That corrected grid passes, but finer checks reveal nearby
failures. For example, alpha=1,x=0.48 gives about +9.5173e-6 with the current
best witness. The bounded box search retains 32 unresolved boxes; a sampled
point at x=0.5187116976 has exponent about +3.8553e-5.

The old piecewise-affine outer majorant has excess over the random-code
exponent near these central weights. The new `d11/outer_majorant.py` adds
five certified left-half supports and their reflections. Its interval proof
processed 10,721 boxes and closed with no unresolved domain. The resulting
39-segment majorant is concave and no larger than the old one.

With this refinement, `d11/dense_refined.py` closes the complete 11% domain
in 1,023 boxes: 433 occupation, 307 Fourier, and 283 scalar witnesses. All
outward bounds pass 512-bit replay. The maximum exponent upper bound is
approximately -4.1033541376e-7. No IMT map, state size, epoch size, or outer
growth constant changes. The old 32 unresolved boxes describe the coarse
bound only, not the current proof status.

## Sparse mechanism

For a nonzero state, an empty epoch applies the exact transition
P=(I+Pi)/2, where Pi maps every state to the uniform nonzero distribution.
Thus P^r=2^-r I+(1-2^-r)Pi. Every coordinate of the fixed A is nonzero,
so its stationary mean output weight per bit is exactly

```text
p0 = 2^18 / (2^19-1).
```

The proof draft bounds the variance of the accumulated empty-gap output
using this exact geometric mixing law. A Lipschitz estimate then controls
the output-weighted kernel uniformly over all entering and exiting states.
Short spacings and collisions have vanishing probability for fixed Q.

### Fixed Q

The new analytic reduction works with the actual full-state process. Its
two-state impulse matrix has bottom-right entry 1-1/(2^19-1). We then enlarge
that matrix entrywise to the conservative matrix

```text
P_plus = [[0, 1], [1/(2^19-1), 1]],
```

The upper matrix is not stochastic. The draft derives the pathwise
beta-integral bound for nonnegative transition coefficients, so this
enlargement does not assume that the old constants transfer unchanged.

`fixed_limit.py` replaces that entry in the product-norm inequality and
checks all 62 outer-segment endpoints with Arb. At distance 0.11 and the old
39/4 growth constant, the largest upper endpoint remains negative:
-0.000867978620634 per active row bit. This establishes the proposed limit's
numerical inequalities. The written proof now supplies the separate
region-to-continuum reduction; the historical receipt remains conditional.

`continuum.py` compares the finite seven-state envelope with this limit for
empty and one-impulse regions. At tilt theta=8, increasing L from 2^12 to
2^28 reduces the one-impulse relative discrepancy from about 0.117 to
1.90e-6. The intermediate values decrease approximately as 1/L. These checks
start from zero or the stationary live distribution. They are historical
diagnostics, not the proof of convergence. `test_homogenization.py` also
checks the actual small-state kernels for zero, one, and two impulses from
every initial state. The proof, unlike these tests, applies to every fixed Q.

`certify_imt_one_two.py` closes Q=1,2 with the conservative continuum matrix.
Its natural exponent bounds per active-row bit are below -0.10918 and
-0.10914, including the existing outer likelihood and active-position costs.

### Growing sparse occupancy

Set beta=(4/5)alpha and z=1-(8/5)alpha. `sparse.py` constructs a positive
affine Collatz vector with zero-state coordinate 1 and live coordinates
(1/1024)(1+h_i alpha). Shell-dependent corrections account for the different
expansion weights and cancellation probabilities. The arbitrary-live-state
coordinate uses h_D=900.

The certified inequality is

```text
T(beta,z) w(alpha) <= (1-96 alpha) w(alpha),
0 < alpha <= 10^-4.
```

`certify_imt_sparse.py` builds degree-257 rational polynomial upper bounds
for the seven residuals. After removing their exact zero at alpha=0, a
positive-tail coefficient bound proves all seven strictly negative on the
remaining interval. It also certifies the maximum replacements used to
construct the bounds. `SPARSE_EXACT.json` and `SPARSE_EXACT_REPLAY.json`
retain these checks; `SPARSE_G96.json` is only the earlier discovery record.

The gamma=96 inequality leaves a conservative
natural exponent coefficient of about -0.01340 after paying the existing
outer likelihood excess and support cost at c=39/4. This estimate uses
distance 0.11, so it also suffices for the lower target. The draft replaces
the expensive sparse route conditioning by marked Bernoulli inputs. This
costs only ln(8 sqrt(Q))/Q in the exponent per active-row bit. A fixed
cutoff Q=4096 separates the uniform sparse bound from finitely many fixed-Q
arguments, making the final union summable without a hidden growth-rate
assumption on Q.

## Reproduce and inspect

From the repository root, with the existing NumPy, SciPy, mpmath, and
python-flint dependencies and local map data:

```text
python -B -m unittest discover -s workstreams/inner_design/imt_asymptotic -p "test_*.py" -v
python -B -m unittest discover -s workstreams/inner_design/imt_asymptotic/d11 -p "test_*.py" -v
python -B workstreams/inner_design/imt_asymptotic/d11/verify_d11.py --output workstreams/inner_design/imt_asymptotic/d11/ASSEMBLY_NEW.json
python -B workstreams/inner_design/imt_asymptotic/verify_imt_asymptotic.py --output workstreams/inner_design/imt_asymptotic/ASSEMBLY_NEW.json
python -B workstreams/inner_design/imt_asymptotic/certify_imt_sparse.py --verify workstreams/inner_design/imt_asymptotic/SPARSE_EXACT.json --output workstreams/inner_design/imt_asymptotic/SPARSE_EXACT_NEW.json
python -B workstreams/inner_design/imt_asymptotic/certify_imt_one_two.py --output workstreams/inner_design/imt_asymptotic/ONE_TWO_NEW.json
python -B workstreams/inner_design/imt_asymptotic/certify_imt_dense.py --verify workstreams/inner_design/imt_asymptotic/DENSE_OUTWARD_D1099_v3.json --output workstreams/inner_design/imt_asymptotic/DENSE_REPLAY_NEW.json
python -B workstreams/inner_design/imt_asymptotic/fixed_limit.py --output workstreams/inner_design/imt_asymptotic/FIXED_LIMIT_NEW.json
python -B workstreams/inner_design/imt_asymptotic/sparse.py --gamma 96 --output workstreams/inner_design/imt_asymptotic/SPARSE_NEW.json
python -B workstreams/inner_design/imt_asymptotic/continuum.py --output workstreams/inner_design/imt_asymptotic/CONTINUUM_NEW.json
```

Outputs are write-once; use unused paths for another run. Do not use Python
with `-O`. The dense replay checks source hashes and every retained exponent
bound at 512 bits. It reuses the producer's formulas; the regression tests
and derivation review remain separate checks.

To repeat discovery, run `screen.py --mode grid --output NEW_GRID.json`, then
`cover.py --grid NEW_GRID.json --delta .1099 --output NEW_COVER.json`.
Use `certify_imt_dense.py --cover NEW_COVER.json --output NEW_CERTIFICATE.json`
to produce outward bounds. The bounded cover may leave unresolved boxes;
the outward producer rejects such a cover. `OPENBLAS_NUM_THREADS=1` avoids
unnecessary overhead for the small discovery matrices. No encoder benchmark
is run by these scripts.

## Recommended next step

Complete the required finite-size IMT cells and implementation integration
before removing the remaining RM2Sub presentation or switching the global
default. The asymptotic proof has been integrated with the unchanged 39/4
growth constant; independent external review and complete artifact packaging
remain worthwhile before publication.
