# Initial landscape at four-bit packets

## Model

Fix a `2^20`-bit message. Partition the message into blocks of `B/2`
bits. Each block uses an independently sampled random linear `[B,B/2]` code.
The encoder globally permutes the resulting four-bit packets and applies one
lap of RandomStepConv with `sigma` retained state bits. The base model has no
additional parity blocks.

For each parameter triple, the evaluator estimates

```text
lambda(B,4,sigma) = -log2 E[X(B,4,sigma)],
```

where `X(B,4,sigma)` counts nonzero outputs below 9% relative weight.
Thus positive values are useful, and the target is 40 bits.

## Initial results

The following table reports representative memory values. A row labeled
"floor" uses `sigma=40`; its dominant inner profile has zero terminations,
so additional memory has negligible effect.

| Outer length `B` | `sigma` | Estimated `lambda` | Dominant support | Episode termination diagnostic |
|---:|---:|---:|---:|---:|
| 256 | 40, floor | -80.78 | 1663 | 0 |
| 512 | 40, floor | -10.51 | 94 | 0 |
| 1024 | 40, floor | 3.86 | 187 | 0 |
| 2048 | 40, floor | 32.62 | 374 | 0 |
| 4096 | 18 | -1.21 | 756 | 34 |
| 4096 | 19 | 27.18 | 753 | 24 |
| 4096 | 20 | 47.25 | 752 | 17 |
| 4096 | 40, floor | 91.53 | 747 | 0 |

The first observed 40-bit crossing occurs between `sigma=19` and
`sigma=20` at `B=4096`. The estimated margins are 27.18 and 47.25 bits,
respectively.

## Interpretation

Memory suppresses paths that contain state terminations. Its benefit vanishes
when a zero-termination path dominates. The floor therefore measures the
late-start obstruction that remains after the inner state becomes persistent.

Increasing `B` improves this floor. One active random outer block then
contains more active packets. A global permutation is less likely to place all
of those packets late enough to produce a low-weight output.

The initial `g=4` slice shows three regimes:

1. At small `sigma`, repeated state terminations dominate.
2. At intermediate `sigma`, each extra state bit gives a large margin gain.
3. At large `sigma`, a zero-termination late-start profile sets the floor.

The floor rises from about -81 bits at `B=256` to about 92 bits at
`B=4096`. This relationship is invisible in a pass/fail table but visible
when `lambda` remains the dependent variable.

## Numerical scope

The random-outer generating polynomial is exact. Sparse outer coefficients
use an occupation-separated inclusion-exclusion sum. Multi-block coefficients
use a lattice saddle approximation. The inner calculation exponentiates the
exact two-state moment matrix and applies a joint Chernoff bound for support
placement and output weight. The support total uses adaptive sampling and
log-linear interpolation.

An exact small-instance check preserves the total message count. Three small
inner instances confirm that the matrix value upper-bounds the exact
low-weight probability. At the multi-block `B=512`, support-785 test point,
an FFT coefficient differs from the lattice saddle by 1.81 bits. Doubling the
support-refinement radius at the former `B=4096,sigma=20` calculation
changed its reported margin by less than `10^-6` bits.

These checks support landscape comparisons. They do not establish a
40-bit certificate. The cells near the target contour require exact central
coefficients and a complete tail cover.

## Next slice

Refine the cells adjacent to each observed 40-bit contour. Then compute enough
intermediate `sigma` values to draw the surface without interpolating across
a change in the dominant profile.
