# Goal 25: replace the remaining cases by one local bound

## Result

The remaining boundary calculation has a single probabilistic interface. It
is a point probability for a finite-state Markov bridge. The packet
histograms do not require separate formulas.

A factor-four Gaussian lattice bound passes every remaining exact comparison.
The audit contains 44 coefficients from dimensions three and four. Its worst
case retains 0.886 bits of slack.

The factor-four inequality is a candidate theorem, not a proved bound. The
remaining task is a finite Fourier certificate on a compact saddle region.

## Exact Markov-bridge representation

Fix an active packet-weight face $S$ and set absent packet variables to zero.
Use $x_0=1$. For log-fugacities $\boldsymbol\ell$, define

\[
B_S(\boldsymbol\ell):=L_S(e^{\boldsymbol\ell}).
\]

Restrict this matrix to state weights that are reachable from zero and can
return to zero. Let $\rho>0$ be its Perron eigenvalue. Fix a positive right
eigenvector $r$. Define

\[
P_{a,b}
:=
\frac{B_S(\boldsymbol\ell)_{a,b}r_b}
{\rho r_a}.
\tag{1}
\]

Equation (1) defines a stochastic matrix on the active states. For every
path from state zero back to state zero,

\[
\prod_{t=1}^N B_S(\boldsymbol\ell)_{q_{t-1},q_t}
=
\rho^N\prod_{t=1}^N P_{q_{t-1},q_t}.
\tag{2}
\]

The eigenvector factors cancel because both endpoints are zero. Therefore
the tilted path law from Goal 22 is exactly the Markov chain (1), started at
zero and conditioned to return to zero at time $N$.

Let $\mathbf K$ record its nonzero packet counts. Let

\[
F_S(\boldsymbol\ell)
:=e_0^{\mathsf T}B_S(\boldsymbol\ell)^N e_0.
\]

For every attainable histogram $\mathbf h$ on the face,

\[
A_{N,4}(\mathbf h,D)
=
F_S(\boldsymbol\ell)e^{-\langle\boldsymbol\ell,\mathbf h\rangle}
\Pr_{\boldsymbol\ell}
[\mathbf K=\mathbf h\mid q_0=q_N=0].
\tag{3}
\]

At a saddle $\boldsymbol\ell_*$, the conditional mean in (3) equals
$\mathbf h$. The exact bridge covariance is

\[
\Sigma_N
:=
\nabla^2\log F_S(\boldsymbol\ell_*).
\tag{4}
\]

## Candidate local bound

Let $d=|S|-1$. Let $\nu_S$ be the index of the reachable histogram lattice.
The candidate inequality is

\[
\boxed{
\Pr_{\boldsymbol\ell_*}
[\mathbf K=\mathbf h\mid q_0=q_N=0]
\le
4\,
\frac{\nu_S}
{(2\pi)^{d/2}\sqrt{\det\Sigma_N}}.
}
\tag{5}

If the right side exceeds one, the trivial probability bound should replace
it. The active face is selected before applying (5), so absent packet weights
do not create singular covariance coordinates.

Equation (5) would recover the polynomial factor omitted by Goal 18. It would
also cover every histogram in an active face at once.

## Exact structural checks

After Goals 23 and 24, the proof-gym shell contains three active face graphs.

| Active packet weights | Dimension | Active state weights | Primitivity exponent | Lattice index |
|---|---:|---|---:|---:|
| $\{0,1,2,3\}$ | 3 | $\{0,1,2,3\}$ | 2 | 2 |
| $\{0,1,3,4\}$ | 3 | $\{0,1,2,3,4\}$ | 3 | 2 |
| $\{0,1,2,3,4\}$ | 4 | $\{0,1,2,3,4\}$ | 2 | 2 |

The graph and lattice calculations use exact integer arithmetic. The lattice
index equals two at every checked length from six through twenty.

The same index holds for all larger lengths. Every boundary path has even
total input weight, which gives an index-two parity constraint. A zero-packet
self-loop embeds every length-six histogram difference at every larger
length. The length-six differences already generate an index-two lattice.

The reachable histogram differences have full rank on each active face.
Every supported path has positive probability at finite fugacities. Hence the
finite bridge covariance is positive definite. On a compact fugacity set,
continuity gives a positive uniform lower eigenvalue.

## Exact-data audit

The audit uses the nine open types at each of four exact shell sizes. It also
uses the three- and four-dimensional proportional profiles through scale
four. This gives 44 comparisons.

The largest observed deficit of the unscaled Gaussian approximation is
1.114 bits. The factor four in (5) supplies two bits. Thus the smallest
observed slack is 0.886 bits.

After inserting the exact faces from Goals 23 and 24, candidate (5) gives:

| Binary length | Candidate loss over exact shell |
|---:|---:|
| 48 | 1.714 bits |
| 80 | 1.483 bits |
| 112 | 1.362 bits |
| 136 | 1.300 bits |

These values are diagnostics. They do not certify (5).

## Why a uniform theorem is plausible

The increments are bounded and the active chains are primitive. The
characteristic function of the bridge is a ratio of matrix powers:

\[
\phi_N(\boldsymbol\theta)
=
\frac{F_S(\boldsymbol\ell_*+i\boldsymbol\theta)}
{F_S(\boldsymbol\ell_*)}.
\tag{6}
\]

Standard spectral proofs split the Fourier domain near and away from the
dual lattice. Near a dual-lattice point, analytic perturbation gives the
Gaussian term and an error controlled by higher derivatives. Away from those
points, primitivity gives a strict spectral-radius gap.

This proof pattern is established for finite Markov-additive processes.
Hervé and Ledoux also give a version uniform over compact sets of transition
matrices. Their theorem treats densities rather than this lattice bridge, so
it supports the method but does not prove (5). Lattice local limit theorems
for additive functionals of mixing Markov chains provide additional context.

References:

- Loïc Hervé and James Ledoux, [Additional material on local limit theorem
  for finite Additive Markov Processes](https://arxiv.org/abs/1305.5644),
  especially the compact-family theorem in Section 6.
- Florence Merlevède, Magda Peligrad, and Costel Peligrad,
  [On the local limit theorems for lower psi-mixing Markov
  chains](https://arxiv.org/abs/2110.10193), which includes lattice cases.

## Remaining certificate

A proof of (5) needs explicit finite constants. It should use Fourier
inversion over one fundamental domain of the index-two lattice.

1. Choose a compact box of saddle parameters required by the outer sum.
2. Bound third and fourth derivatives of $\log F_S$ on a central Fourier
   neighborhood.
3. Bound the complex spectral radius away from the two dual-lattice points.
4. Integrate both bounds and verify that their sum is at most four times the
   Gaussian lattice factor.

The exact data use saddle coordinates between approximately -3.0 and -0.68.
This range is only a calibration box. The outer proof must determine the
box required at target length.

The next goal should compute that target saddle box. Without it, a Fourier
certificate would either cover too little or waste margin on an unnecessarily
large parameter range.

## Reproduction

Run

```powershell
python scripts/audit_riffle_boundary_uniform_local_candidate.py --factor 4
```

The receipt is `receipts/goal25_boundary_uniform_local_candidate.json`.
