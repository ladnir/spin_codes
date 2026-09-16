# BCH-256 overlap audit and feedback diagnosis

The overlap audit did not close the BCH-256 migration. It identifies the
zero-state contribution as the main obstruction at the retained difficult
witness. Two weight-five feedback candidates pass that point. Neither has
a full distance certificate, an optimized implementation, or measured timings.
The supported encoder and the completed quarter-rate certificate are unchanged.

## Retained point results

The point uses K=2^20, t=128, s=19, Q=3482, and normalized density coordinate
v=21/128. The outer bound, reference probability, and output tilt remain fixed
across these comparisons. `DENSE_GAP_POINT.json` in the parent directory
defines the witness. Probability in the underlying distance argument is over
the shared setup, as specified in the parent README.

| Feedback and transfer | Diagnostic margin at this point (bits) |
|---|---:|
| Current weight-three feedback, original cap | -1536.958 |
| Current weight-three feedback, audited overlap cap | -1535.280 |
| Balanced feedback B=A^T, original cap | 862.757 |
| Weight-five seed 0, original cap | 720.623 |
| Weight-five seed 1, original cap | 705.637 |

A positive entry bounds this point only, not the union over all occupancies
and densities. A negative entry means the upper bound is uninformative;
it is not evidence of a bad code realization.

## Exact overlap audit

For each nonzero dual state a, let D be the support of B^T a. The native
audit packs disjoint full-rank bases among A's columns restricted to D and
its complement. Let d and e denote the resulting numbers of bases.
Every nonzero state q has a nonzero evaluation on each basis. Therefore,
with h=wt((Aq)|D), v=wt(Aq), and w=|D|, the overlap satisfies

    max(0, v+w-128, d) <= h <= min(v, w, v-e).

The audit enumerates all 524287 nonzero dual states and retains 161 groups.
The tighter Fourier bound maximizes over the two feasible endpoints.
`PACKING_AUDIT_replay.json` records an exact repeat of the full audit.
Small exhaustive tests independently check the packing and weighted transfer.

## Why this tightening is insufficient

Fix the retained reference input probability theta and output tilt lambda.
Write z=exp(-lambda), g0=1-theta+theta*z, p=theta*z/g0, and rho=1-2p.
Let K_j count weight-j words in ker(B). The exact zero-to-zero contribution is

    T00 = sum_j K_j (theta*z)^j (1-theta)^(128-j)
        = g0^128 * 2^-19 * sum_a rho^wt(B^T a).

There are E=16384 epochs. The all-zero-state trajectory contributes T00^E,
regardless of any overlap cap. The tightened moment exceeds this contribution
by only about 46.94 bits at the retained witness. Substituting the zero-state
contribution alone into the fixed proof expression still gives -1488.34 bits.

For any feedback map with 128 weight-three columns, dual words of input
weight j have a fixed sum of output weights:

    n_j = binom(19,j),
    W_j = 128 * sum_{h in {1,3}} binom(3,h) binom(16,j-h).

For 0<rho<1, convexity bounds the sum of rho raised to these integer weights
from below by distributing W_j as evenly as possible among n_j words.
This gives at least 0.209625 excess bits per epoch over the ideal kernel
term g0^128/2^19. The current map has 0.275778 excess bits per epoch.
Even this optimistic weight-three floor leaves the fixed proof expression
at most about -404.48 bits.

These limits concern the retained outer scalar, reference probability, and
output tilt. They do not rule out other witnesses or a tighter outer bound.
Multiplying a moment lower bound by an outer upper bound does not give a
lower bound on the actual failure probability.

## Files, reproduction, and next step

`POINTS.json` and its replay record the overlap-cap comparisons.
`KERNEL_DIAGNOSIS.json` retains the exact candidate columns and point bounds;
its replay checks the numerical bounds at 512-bit precision. Four tests in
`test_overlap.py` passed, including a rational check of the convexity floor.

The next experiment is a bounded certificate search for weight-five seed 0,
with balanced feedback as a control. Rebuild each candidate's own spectra
and cancellation bounds before checking sparse and dense coverage. Existing
weight-three certificates must not be reused as certificates for new maps.
Performance measurements follow proof coverage; XOR counts alone are not timings.

Run from the repository root with the inner-design Python dependencies:

```text
python -m unittest discover -s workstreams/inner_design/asymmetric/bch256/overlap -p "test_*.py" -v
```

The native audit is built by `build.cmd` from this directory, or by compiling
`packing_audit.cpp` with a C++20 compiler. Build products are deliberately
not tracked. The historical audit authenticates the original executable's
hash as well as its source. A fresh build may have a different binary hash;
it is not a byte-identical replay of that historical run.

For a fresh native audit, run `audit_overlaps.py --output <new-path>` and
then the same command with `--verify`. Compare its `audit` payload with
`PACKING_AUDIT.json`. The historical point tools currently authenticate the
original binary through that receipt; replaying those tools from a clean
checkout needs a new audit binding and regenerated point receipts. This
portability work remains open. Do not overwrite retained receipts silently.
