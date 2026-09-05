# Two-stage sparse-EA route

## Objective

The closed degree-33 constituent costs 16,895 XORs. This route asks whether
a two-stage constituent can satisfy the same realized-spectrum caps at lower
cost. Reusing the same caps also reuses the complete RandomStepConv-M22 SPIN
transfer without changing its probability space.

## Construction

Let \(E_0:\mathbb F_2^{256}\to\mathbb F_2^{512}\). Each row of \(E_0\)
independently samples a uniform support of size \(r_0\). Let
\(E_1:\mathbb F_2^{512}\to\mathbb F_2^{512}\) have the analogous row law
with degree \(r_1\). Setup samples \(E_0\) and \(E_1\) independently and
defines

\[
  C_2:=A E_1 A E_0,
\]

where \(A\) is the zero-initialized length-512 accumulator. A bounded setup
samples at most 16 independent pairs and accepts the first pair for which
\(C_2\) has rank 256. Every SPIN outer position reuses the accepted pair.

The XOR count is

\[
  512(r_0-1)+511+512(r_1-1)+511
  =512(r_0+r_1)-2.
\]

## Reusable theorem interface

Let \(T_w\) be the integer caps from the closed degree-33 theorem. The new
outer task is to prove

\[
  \Pr[
    \operatorname{rank}(C_2)<256
    \text{ or }
    (\exists w)\ A_w(C_2)>T_w
  ]
\]

with enough margin after bounded rank-tested setup. Conditional on this cap
event, the existing SPIN transfer applies verbatim. No new occupation-\(Q\)
calculation is required.

The complete binary64 one-word calculation gives the following initial gate.

| Degrees | XORs | Zero-cap tail margin | Minimum cap/mean headroom | Largest uniform variance factor compatible with 40 end-to-end bits |
|---:|---:|---:|---:|---:|
| \((17,3)\) | 10,238 | 43.9504 bits | 0.04213 bits | 53.81 |
| \((19,3)\) | 11,262 | 49.4771 bits | 0.04213 bits | 57.46 |
| \((17,5)\) | 11,262 | 45.9211 bits | 0.04213 bits | 57.04 |

Every expected spectrum lies below every frozen positive-shell cap. The last
column assumes
\(\operatorname{Var}(A_w)\le F\mathbb E[A_w]\) with the same \(F\) on every
shell. It is a sensitivity calculation, not a proved variance bound.

The two candidates with 11,262 XORs expose a proof tradeoff. The
\((19,3)\) candidate has 3.56 more bits of zero-cap-tail margin. The
\((17,5)\) candidate has a stronger second mixer. Krawtchouk orthogonality
gives

\[
 \mathbb E_{W\gets\operatorname{Bin}(512,1/2)}
   [\beta_r(W)^2]
 =\binom{512}{r}^{-1}.
\]

Increasing the second degree from 3 to 5 reduces this dense squared-character
average by 13.6583 bits. This identity motivates testing \((17,5)\), but it
does not by itself bound the composed variance.

The next equal-raw-cost line is

| Degrees | XORs | Zero-cap tail margin |
|---:|---:|---:|
| \((17,7)\) | 12,286 | 47.4939 bits |
| \((19,5)\) | 12,286 | 51.1684 bits |
| \((21,3)\) | 12,286 | 52.5371 bits |

These values are also binary64 diagnostics. Degree 7 gives the strongest
second square mixer in this line.

## Exact-circuit cost check

Raw row-by-row XOR counts substantially overstate the cost of a sampled
sparse matrix. A deterministic Paar common-subexpression heuristic was run
on exact sampled matrices. The resulting straight-line circuits were checked
exactly in both the forward and transposed directions.

| Map | Raw forward XORs | Optimized forward XORs | Optimized transposed XORs | Forward saving |
|---|---:|---:|---:|---:|
| one-stage \(E_0\), \(256\to512\), degree 33 | 16,384 | 8,966 | 9,222 | 45.28% |
| first-stage \(E_0\), \(256\to512\), degree 17 | 8,192 | 5,478 | 5,734 | 33.13% |
| square \(E_1\), \(512\to512\), degree 7 | 3,072 | 2,800 | 2,800 | 8.85% |

Including accumulators, the one-stage degree-33 circuit costs 9,477 forward
XORs. The depth-two \((17,7)\) circuit costs 9,300. Thus depth two saves only
177 XORs, or 1.87%, after optimizing both exact sampled circuits. The
corresponding transposed totals are 9,733 and 9,556, again a difference of
177. The one-stage optimized circuit has an upper bound of 2,025 live
intermediates under the generic schedule. A greedy transposed schedule
reduces the measured surrogate to 1,305 packed values.

A matched Peach benchmark maps \(2^{21}\) input blocks to \(2^{20}\) output
blocks. The one-stage path takes 8.703 ms, compared with 4.896 ms for two
optimized BCH constituents. These values are the means of two run medians.
Thus the one-stage outer is 1.777 times slower in this isolated benchmark.

The sampled circuit does not certify the realized spectrum of that fixed
matrix. It is an exact implementation-cost experiment within the same random
matrix ensemble. The closed theorem remains a probability statement over
the setup sample.

## Exact moment decomposition

For a fixed message \(x\ne0\), the first stage induces a distribution on
\(U_x:=AE_0x\). Conditional on \(U_x=u\), the rows of \(E_1u\) are
independent Bernoulli variables whose parameter depends only on
\(\operatorname{wt}(u)\). This gives an exact composition of the one-word
weight transition and produces the diagnostic table above.

For distinct nonzero messages \(x,y\), put

\[
  U:=AE_0x,
  \qquad
  V:=AE_0y.
\]

Conditional on the four-symbol type of \((U,V)\), every row of
\((E_1U,E_1V)\) has the exact hypergeometric pair law. Consequently, the
second factorial moment of every final shell is determined by a composition
of two four-symbol type transitions. The existing small-model verifier
evaluates this composition in exact rational arithmetic.

The law of total variance exposes the new term:

\[
 \operatorname{Var}(A_w(C_2))
 =\mathbb E_{E_0}[
    \operatorname{Var}_{E_1}(A_w(C_2)\mid E_0)]
  +\operatorname{Var}_{E_0}(
    \mathbb E_{E_1}[A_w(C_2)\mid E_0]).
\]

The one-stage proof controls only the analogue of the first term when its
input messages are deterministic. The second term records the common
randomness in \(E_0\). Thus replacing degree 33 by a product of degrees 17
and 3 inside the old formula is invalid.

There is also a direct dual formulation. For a final Fourier character
\(z\in\mathbb F_2^{512}\), define

\[
  S_z(E_1):=A^{\mathsf T}E_1^{\mathsf T}A^{\mathsf T}z.
\]

If \(R_0(s)\) is the indicator of \(E_0^{\mathsf T}s=0\), Fourier inversion
gives

\[
 A_w(C_2)+\mathbf1\{w=0\}
 =2^{-256}\sum_z K_w^{(512)}(\operatorname{wt}(z))
   R_0(S_z(E_1)).
\]

The additive term at \(w=0\) removes the zero message and does not affect the
variance. Therefore the shell variance is the quadratic form induced by

\[
 \mathbb E_{E_1}[
   \Pr_{E_0}[E_0^{\mathsf T}S_z=0,
             E_0^{\mathsf T}S_u=0]]
 -\mathbb E_{E_1}[
   \Pr_{E_0}[E_0^{\mathsf T}S_z=0]]
  \mathbb E_{E_1}[
   \Pr_{E_0}[E_0^{\mathsf T}S_u=0]].
\]

The exact verifier verify_two_stage_dual_covariance_small.py enumerates the
\((K,B,r_0,r_1)=(2,3,1,1)\) ensemble. For every shell, it checks this dual
quadratic form and the law of total variance against direct enumeration.
All comparisons use exact rational arithmetic.

## Proposed proof reduction

The exact pair law suggests the following certificate.

1. Express both terms in the variance decomposition as quadratic forms of
   the first-stage ordered-pair type measure.
2. Apply the second-stage likelihood comparison conditionally on each
   intermediate pair type.
3. Contract the resulting functions through the first-stage four-state
   accumulator transfer without materializing all
   \(\binom{515}{3}=22{,}632{,}705\) intermediate types.
4. Bound the final shell variances against the frozen integer caps.
5. Add rank conditioning, 16-attempt abort, and the unchanged conditional
   SPIN failure.

The principal technical choice is the compression in step 3. The preferred
route is a weighted spectral or generating-function contraction of the
four-symbol transfer. Enumerating the complete intermediate type table for
every input pair is not acceptable.

The direct dual formula offers a second implementation route. It moves the
second stage into a random transformation of the final Fourier characters
and then applies the known return kernel of \(E_0\). This representation may
permit a block or tensor contraction before any 22.6-million-entry table is
formed.

## Square-mixer obstruction and invertible conditioning

The one-stage outward certifier was also applied to the square second mixer
alone at shell 80. This is a diagnostic of the reference covariance, not a
certificate for the two-stage constituent.

| Square degree | Normalized positive pair excess | Variance/mean upper bound |
|---:|---:|---:|
| 5 | \(9.83582\mathbin{\cdot}10^{-4}\) | \(1.11852\mathbin{\cdot}10^{92}\) |
| 7 | \(1.15379\mathbin{\cdot}10^{-5}\) | \(1.31208\mathbin{\cdot}10^{90}\) |

Degree 7 improves both quantities by a factor of 85.25, but remains far from
the factor-512 target. Fixed-degree escalation is therefore not the missing
argument under this positive-only relaxation.

There is a structural way to remove the square mixer's uniform-reference
defect. If \(E_1\) is invertible, then \(AE_1A\) is invertible. For a dense
uniform \(512\times256\) reference map \(G\), \(AE_1AG\) has exactly the
same law as \(G\), for every fixed invertible \(E_1\). In particular, two
distinct nonzero messages still produce independent uniform words. This is
an exact fact, not a heuristic.

A rank experiment found 39 full-rank matrices among 256 independently
sampled degree-7 square matrices. This 15.23% frequency is a measurement,
not a lower bound. Conditioning on full rank couples the rows of \(E_1\), so
the present row-product proof cannot simply be reused. Moreover, the
manageable one-stage likelihood comparison uses a biased
dominant-character reference rather than the dense uniform reference.
Invertibility therefore removes the uniform component of the obstruction,
but a new conditional contraction or signed-defect lemma is still required
for the sparse first stage.

## Current status

Proved facts:

- the one-word and ordered-pair transition formulas compose exactly;
- the small exact verifier confirms multistage means and variances;
- the small dual verifier confirms the composed return-kernel formula;
- the closed caps imply the existing downstream SPIN theorem; and
- the optimized sampled-circuit identities and ranks are checked exactly.

Certified realized-spectrum facts for the target two-stage parameters: none
yet. The present mean and margin values use nearest-binary64 arithmetic.

Open obligations:

- outward certification of the complete one-word spectrum;
- a compressed target-size evaluation of both variance terms;
- simultaneous cap accounting with the frozen integers;
- rank-conditioned bounded setup; and
- a manifest binding the new receipts to the existing transfer.

The cost experiment changes the priority. The optimized \((17,7)\) circuit
saves only 1.87% relative to the already-certified optimized degree-33
one-stage circuit, while requiring a new conditional covariance theorem.
The degree-33 benchmark is 77.7% slower than optimized BCH and shows material
register pressure. The next choice is therefore either a pressure-aware
one-stage circuit search or the new depth-two theorem. If depth two resumes,
the first mathematical target should be the invertible-conditioned degree-7
square ensemble and a bound on the transported sparse-reference defect.
