# Keeping active rows and input weight jointly

This experiment asks whether two separate relaxations cause the rate-three-quarters proof gap.
The middle bound fixes total input weight but sums over active outer rows.
The sparse bound fixes active rows but sums over total input weight.
Both quantities can be retained without changing the construction or its randomness.

## Conditional source factor

Fix a GF(256) MDS `[b,k]` outer, with parity `p=b-k`, repeated in `L` rows.
Let `A` be the number of nonzero rows and `W` their total support weight.
Feasibility requires `(p+1)A<=W<=bA`.
Each of the `b` regions contains `L` positions and receives one independent uniform permutation.

For a fixed vector of regional support sizes, introduce independent Bernoulli coordinates with probabilities `alpha_j` equal to those sizes divided by `L`.
This is a comparison distribution, not the distribution of the actual outer codeword.
Conditioning on the regional counts recovers independent uniform supports within regions.
The existing modal-conditioning bound contributes at most

\[
\left(2^L/\binom{L}{L/2}\right)^b.
\]

For positive marks `tau,lambda`, put `beta=tau*lambda`.
Let `H` be the sum of the `b` Bernoulli variables in one row.
Its empty and occupied source factors are

\[
E=\Pr[H=0],\qquad
O=\tau^{-p}\mathbb E[\beta^H\mathbf1_{H\ge p+1}].
\]

The exact-`A` source expectation is `binom(L,A) E^(L-A) O^A`.
An activation mark `v` produces the equivalent generating polynomial `(E+vO)^L`.
Using the exact coefficient avoids the small loss from a Chernoff extraction in `v`.

When `beta>=2`, both row factors increase under pairwise averaging of the regional probabilities.
For `E`, the coefficient sequence is `(1,0,...,0)`.
For `O`, it is zero through degree `p`, then `tau^-p beta^h`.
Both sequences have nonnegative second differences.
At the occupied jump, the only additional condition is `beta-2>=0`.
Thus both factors are bounded by their values at `alpha=W/(bL)`:

\[
E(\alpha)=(1-\alpha)^b,\qquad
O(\alpha)=\tau^{-p}\sum_{h=p+1}^b
\binom bh(\beta\alpha)^h(1-\alpha)^{b-h}.
\]

For `beta<2`, the script uses the established relaxation
`O(alpha)=tau^-p(1-alpha+beta*alpha)^b`.
Its geometric coefficient sequence is discretely convex for every positive `beta`.
The relaxation includes empty and infeasibly light occupied rows, so it remains an upper bound.

Under the unchanged regional rank-cover argument, the previous middle source factor `(E+O)^L` can therefore be replaced by

\[
\binom LA E(\alpha)^{L-A}O(\alpha)^A.
\]

The inner matrix power, modal penalty, output mark, and projective rank factor are unchanged.
The total-weight extraction retains `lambda^-W` and the inner input mark `x^-W`.
No independence between actual input forms is asserted or used.

## Fixed-weight sparse extraction

The sparse regional polynomial also accepts a total-weight mark.
Write `C_h` for its positive matrix coefficient at support size `h` before inserting the regional rank mark.
Replace its regional matrix by

\[
R_A(\lambda)=\sum_{h=0}^A C_h
\frac{\binom Ah}{\binom Lh}(\tau\lambda)^h.
\]

The positive coefficient bound at total weight `W` contributes `lambda^-W` to `e_0^T R_A(lambda)^b 1`.
Here `lambda` may be any positive number.
Unlike the earlier minimum-weight-only mark, this step does not require `lambda>=1`.

## Results: the same obstruction survives conditioning

The tested construction is full-Cauchy width eight, state 32, with outer `[128,96]`, requested `K=2^20`, and output distance 10%.
It uses the conservative return cover from `gf8_mds_geometry_screen.py`.
These are floating point screens, not interval certificates or statements about realized code distance.

| Class and bound | Unconditioned margin | Jointly conditioned margin |
|---|---:|---:|
| `W=70042`, middle bound; dominant `A=1524` | -17,188.5148 bits | -17,182.0103 bits |
| `A=256`, sparse bound; tilted `W=13146` | -1,131.10564 bits | -1,131.10560 bits |

At the first class, reoptimization leaves all four marks unchanged to roughly seven decimal places.
The 6.50-bit improvement is the binomial modal-mass factor from fixing `A`.
At the second class, the optimized weight mark is approximately `1.000114`, effectively the original neutral mark.
The joint middle cover for `(A,W)=(256,13146)` remains weaker: its 65-iteration witness gives about -5,842 bits and does not report convergence.

Retaining both quantities is valid, but it does not remove these tested saddle-point classes.
The experiment therefore argues against spending another broad sweep on activation and weight marks alone.
A changed local construction or a stronger coupling between input weight, state evolution, and accumulated output weight is the next useful direction.

## Reproduction

`theory/gf8_joint_source_screen.py` implements both conditional bounds.
Its five tests check exact pairwise homogenization, the activation polynomial, exact toy column conditioning, feasible classes, and recovery of the old sparse bound at the neutral mark.
The homogenization and the unrestricted fixed-weight mark were also independently reviewed.
No existing certificate or implementation was changed.

```text
python research/workstreams/char2_field_spin/theory/gf8_joint_source_tests.py
python research/workstreams/char2_field_spin/theory/gf8_joint_source_screen.py --mode middle --active 1524 --weight 70042 --iterations 80
python research/workstreams/char2_field_spin/theory/gf8_joint_source_screen.py --mode sparse --active 256 --weight 13146 --iterations 65
```
