# Feasibility of one Turbo-style inner in Scalable SPIN

## Construction held fixed

This investigation retains the paper's one-bit regional route, one shared
growing Golay--BAA outer, rate one half, and native schedule
`b >= (39/4) log2 N`. It replaces the `t=128,s=19` IMT inner by Turbo's
actual `t=64,s=22` map: expansion `A`, feedback `C=A^T`, and an independent
uniform `GL(22,2)` state update after each step. The manuscript and the
finite Turbo encoder have not changed.

The map in `spin/src/certified22/productionMap.h` has the nonzero expansion
weights and multiplicities

| weight | 16 | 24 | 28 | 32 | 36 | 40 | 48 | 64 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| count | 2,604 | 291,648 | 888,832 | 1,828,134 | 888,832 | 291,648 | 2,604 | 1 |

`verify_pair_return_map.py --spectrum` independently re-enumerates these
weights from the production map. The MacWilliams transform gives the
weight enumerator of `ker(A^T)`, with minimum nonzero weight 8. Thus each
one-bit impulse has nonzero feedback.

## A three-state moment bound

The useful distinction is between a live state immediately after an
input impulse and one that has passed through a uniform update. Call them
`F` and `U`; `0` denotes the zero state. `F` carries arbitrary live mass.
`U` carries mass whose density at each of the `M=2^22-1` nonzero states
is at most `U/M`. These are upper-envelope coordinates, not a partition
of the actual state space.

For independent `Bernoulli(beta)` step input and `z` in `(0,1]`, put

```
g0 = 1-beta+beta*z,
g1 = beta+(1-beta)*z,
f_w = g0^(64-w) g1^w,
F = max_{w>0 in the expansion spectrum} f_w,
U = (sum_{w>0} A_w f_w)/M.
```

Let `K_j` be the exact weight enumerator of `ker(A^T)`. Define

```
H = sum_j K_j (beta*z)^j (1-beta)^(64-j),
u = 1 - sum_j K_j beta^j (1-beta)^(64-j).
```

Here `H` is the weighted probability of zero feedback from a zero state;
`u` is the unweighted probability of nonzero feedback. With rows and
columns ordered `(0,F,U)`, the following nonnegative matrix dominates
the step's weighted state transition:

```
T(beta,z) = [[H,                   g0^64-H,         0],
             [min(F,u)/M,          0,               F],
             [min(U,u)/M,          0,               U]].
```

To justify the two live rows, fix a nonzero entering state `q` and an
input `x`. The sampled map sends `q` uniformly over all nonzero states,
so each nonzero successor has probability at most `1/M`, while a zero
successor is possible only when `A^T x` is nonzero. Averaging the output
weight over `x` bounds the first live row by `F`; averaging over a
uniform live entering envelope gives `U`. The `min` reset entries use
both the weighted moment and the nonzero-feedback probability. Induction
on the unnormalized state measure proves the matrix domination. It does
not assume the input syndromes are uniform, nor independence between
output weight and the next state.

The existing route comparison turns fixed outer-row weights into the
same Bernoulli input moment used in the current positive-occupation
proof. For a row fraction `alpha`, average row weight `x`, proposal
probability `p`, value probability `y`, and `z=exp(-lambda)`, the
resulting candidate exponent per output bit is

```
alpha * a_hat_BAA(x) + KL(alpha||p) + alpha*KL(x||y)
    + log(rho(T(p*y,z)))/64 + 0.11*lambda.
```

`a_hat_BAA` is the already certified refined outer majorant. Its active
intervals must be used as given; extending an affine support outside
its certified interval is invalid.

## Results and their scope

* `shared_inner_dense_screen.py --summary` gives 54 negative sampled
  positive-occupation exponents. Its largest is about `-1.4871e-5` at
  `alpha=1e-4,x=0.16`.
* `shared_inner_dense_cover.py --max-boxes 100000 --margin 1e-7` gives
  a 30,241-box binary64 discovery cover. The follow-up
  `shared_inner_dense_certificate.py --max-boxes 150000` starts from
  exact rational outer-segment endpoints, recursively partitions exact
  rational boxes, and checks the selected witnesses with 256-bit Arb
  arithmetic. It processed 92,595 boxes and accepted 46,317; every
  accepted box has negative outward upper bound, with maximum below
  `-1.00147e-7`. For a fixed witness the exponent is convex in
  `(alpha,alpha*x)` on one affine outer segment, so four vertex checks
  cover the entire box. This certifies the positive-occupation
  **exponent inequality**, assuming the three-state transfer lemma and
  the imported route comparison. The same exact partition and witness
  search replayed at 512-bit Arb precision with all 46,317 leaves
  passing. It is not by itself a distance theorem.
* `shared_inner_dense_screen.py --sparse-gate` verifies with exact
  rational arithmetic that the positive vector
  `v=(1,(1+24*alpha)/1024,1/1024)` obeys
  `T v <= (1-48*alpha)v` for `0<alpha<=1e-4`. It uses kernel distance 8,
  expansion minimum weight 16, and mean live expansion weight at least
  32. The three normalized residual upper bounds at the right endpoint
  are `-3.020856`, `-1.439532`, and `-2.663212`. Combining this with
  `p=(8/5)alpha`, `y=1/2`, `z=1-(8/5)alpha`, the outer majorant, and
  the existing active-row-position cost gives an outward upper exponent
  below `-0.0134042*alpha`. The exact Collatz residuals and 256-bit
  logarithmic check are replayed by `shared_inner_remaining_certificate.py`.
* For each fixed number of active rows, empty gaps refresh a live state
  to uniform in one step. The limiting zero/live impulse matrix is
  `[[0,1],[1/M,1-1/M]]`, as in the current proof with a smaller reset
  probability. The copied one- and two-row continuum formulas give
  outward exponents below `-0.109827` and `-0.109812` per active-row bit
  (`shared_inner_remaining_certificate.py`). For fixed `Q>=3`,
  `shared_inner_remaining_certificate.py` now recomputes all 62
  inherited fugacity endpoints directly with the Turbo reset
  probability and live density. Their largest outward upper bound is
  below `-0.00086852` per active-row bit. These are continuum
  inequalities; the
  finite-region reduction must still be written for the new inner.

The fixed-row finite-region reduction is favorable here. An empty step
from a live state makes its successor exactly uniform nonzero, regardless
of the entering state. Conditional on that live state, the update is also
independent of the emitted weight. For any fixed number `Q` of active
rows, all but `O_Q(1/L)` of the routed one-bit impulses occur in distinct,
nonadjacent 64-bit steps of a region. On the remaining trajectories, every
long empty gap factors after its first step. With step tilt
`z_L=exp(-theta/L)`, that first step changes the gap moment by
`1+O_Q(1/L)`; the other steps contribute powers of the exact uniform-live
moment. Riemann sums over the uniformly placed impulses then yield the
same two-state simplex kernels as the current proof, with
`p0=2^21/(2^22-1)` and reset probability `1/(2^22-1)`. The entries that
vanish in the limit vanish in the finite transfer for the same reasons:
an empty region cannot change zero/live status, and one isolated bit
cannot take zero to zero. This gives the route to a componentwise
comparison across `b=O(log L)` regions, but its uniform constants and
boundary cases still need lemma-level verification.

These computations make a unified inner **strongly feasible at the same 11%
target and growth schedule**. The numerical inequalities now close with
outward arithmetic, but they do not alone establish the asymptotic theorem.
The three-state domination and imported route comparison need a written
lemma-level audit, the fixed-occupation finite-region limit needs its
uniform constants and boundary cases checked, and the first-moment
assembly must be redone for this setup law. The existing IMT Lean theorem
does not transfer.

## Recommended proof order

1. Audit the three-state domination and Bernoulli route comparison with
   explicit quantifiers over outer realization and routed supports.
2. Review the exact rational positive-occupation cover and its spectral
   Collatz vectors independently; replay it at higher precision.
3. Prove the fixed-row finite-region comparison with uniform constants
   and verify its zero entries and shared-state boundary cases.
4. Assemble the shared-outer first moment and reverify the theorem in
   Lean before changing the paper's construction claims.

No implementation timing was inferred from these proof screens. The
finite Turbo encoder already uses this inner at K18/K20; the one-bit
scalable variant has different routing work and requires measurement if
presented as a practical operating point.
