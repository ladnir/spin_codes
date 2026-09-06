# Handoff: concrete distance certificates for Structured SPIN

> **Superseded for the current finite-parameter lane.** Read
> `PARAMETER_LANDSCAPE_HANDOFF.md` first. It records the RM2Sub landscape
> database, parameter extrapolation, activation-state audit, the proved
> \((t,s)=(128,15)\) first-moment obstruction, and the active
> \((t,s)=(64,20)\) BCH closure route. The material below remains historical.

Date: 2026-09-03.

Superseding status: the one-stage sparse-mixer--accumulator certificate that
was open when most of this handoff was written is now closed.  The definitive
statement is `ONE_STAGE_SPARSE_EA_RANDOMSTEPCONV_CERTIFICATE.md`; its audit
entry point is `audit_one_stage_sparse_ea_certificate.py`.  It proves, for
the close RandomStepConv-M22 SPIN variant at \((K,N)=(2^{20},2^{21})\),
minimum distance at least 228,590 except with probability below
\(2^{-41.3621359295}\).  The historical sections below retain the derivation
and rejected approaches, but any sentence saying that this one-stage route is
still incomplete is superseded by the closure section below.

Second superseding status: the fixed-RM(4,9), RM2Sub-t64,s14 small-k
multiband gap is now closed as a binary64 diagnostic. The refined groups
1--95, 96--416, and 417--512 cover every composition for Q=30 through 256.
Their union has 1228.638 bits. The audit checks 2,857,249 compositions. See
small_k_replay/RM2SUB_REFINED_BAND_BRIDGE.md. Directed outward rounding and
combination with the Q=1--29 receipts remained open at that stage. The later
combined Q=1--256 diagnostic has 42.577982 bits, with Q=1 weakest. Directed
outward rounding is the remaining certificate step.

This document is the starting point for a new chat. It records the global
research objective, the frozen implementation, the theorem interfaces, the
current sparse-mixer concentration task, and the next useful proof steps.

## The real objective

The primary objective is a concrete, independently checkable distance
certificate for Structured SPIN, or for a close variant that preserves its
important operational properties.

The preferred finite target is:

- message dimension \(K=2^{20}\);
- output length \(N=2^{21}\);
- minimum distance at least
  \[
    D=\lfloor 0.11N\rfloor+1=230687;
  \]
- setup-failure probability at most \(2^{-40}\);
- one local outer constituent sampled once and reused at every outer-block
  position;
- practical setup, without enumerating \(2^{128}\) or \(2^{256}\) messages;
- linear ordinary and transposed encoding work; and
- approximately 11 ms or less on the Peach Ryzen 9 7950X reference machine.

A literal 10.9-percent result remains useful when it materially improves
parameters or proof closure. At \(N=2^{21}\), that threshold is \(D=228590\).
Every result must state its actual distance and failure bound.

The phrase “one constituent” is important. Sampling an independent outer code
at every block position gives an easier ensemble, but it is not the desired
construction. A bounded setup may sample several complete candidate
constituents, test an efficient property such as rank, and retain one accepted
constituent. The accepted constituent is then reused everywhere.

The frozen implementation describes lengths in 128-bit implementation
blocks, while the mathematical code documents use binary-coordinate notation.
Any final theorem must give the exact conversion and must not silently identify
these units.

## Required startup reading and scope

The new chat should read these files before changing anything:

1. collaboration/README.md
2. collaboration/INTEGRATION_CONTRACT.md
3. collaboration/TASK_ASSIGNMENTS.md
4. SPIN_NAMING.md
5. RIFFLE_NEXT_WORK_ROADMAP.md
6. workstreams/finite_asymptotic_theory/BRIEF.md
7. this handoff
8. workstreams/finite_asymptotic_theory/FINITE_K20_ATTEMPT_LOG.md
9. workstreams/finite_asymptotic_theory/pure_expander_accumulate/README.md
10. workstreams/finite_asymptotic_theory/pure_expander_accumulate/ONE_STAGE_VARIANCE_ROUTE.md

Use the Controlled Writing for Cryptography skill for mathematical prose.
Write only within workstreams/finite_asymptotic_theory/. Read-only inspection
outside that directory is permitted.

The frozen baseline is not at repository root. It is under:

    constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/

Do not modify its frozen_source directory, MAIN_CODE_FREEZE.md, proof receipts,
or SOURCE_MANIFEST.json. The manifest is at:

    constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/frozen_source/SOURCE_MANIFEST.json

Do not reinterpret a proposed construction without first obtaining user
approval. In particular, do not replace a single reused constituent by
independently resampled row constituents.

Performance is a first-class constraint. Never run two benchmarks or large
numerical certificate jobs simultaneously. At handoff time no Python process
is running.

## Frozen Structured SPIN

The frozen candidate is Structured SPIN \((B=256,t=128,s=19)\). Its
implementation applies:

1. a fixed 256-to-128 extended-BCH-based outer map;
2. an independently sampled disjoint ParityFanout-31x33 transform per outer
   block;
3. per-block coordinate routing;
4. a bit transpose;
5. independent region permutations; and
6. the fixed RM2Sub-S19 inner with step size 128.

The frozen correctness evidence passes both dense-inner and staged
end-to-end checks. Its checksum is 0x95c9d722a9539fef. The frozen 21-trial
Peach median is 10.823872 ms.

The displayed end-to-end distance ledger is not a proof. At relative distance
0.11 it reports a nearest-binary64 first-moment margin of 55.864642 bits.
The occupation partition is

\[
  q=1,\qquad 2\le q\le100,\qquad 101\le q\le8192,
\]

with displayed margins 55.8646, 108.8244, and 121.0192 bits. The one-active
class limits the result.

Two gates prevent this diagnostic from becoming the desired certificate:

1. the local outer spectrum uses a modeled even-floor spectrum rather than a
   proved spectrum or sufficient envelope for the actual constituent; and
2. the final ledger uses nearest binary64 rather than independently
   checkable outward arithmetic.

The frozen proof status also describes outer maps sampled independently at
different outer-block positions. That probability space does not meet the
desired single-reused-constituent interface. A final theorem must resolve
this mismatch explicitly.

## General theorem interfaces already completed

The workstream deliverables are present:

- FINITE_LENGTH_FRAMEWORK.md;
- ARBITRARY_LENGTH_WRAPPER.md;
- ASYMPTOTIC_SCALING.md;
- STRUCTURED_SPIN_THEOREM_TARGET.md; and
- MERGE_SUMMARY.md.

The exact finite first-moment theorem uses a setup experiment that samples an
outer map \(O\), interleaver \(\Pi\), and inner map \(I\), and sets

\[
  E=I\circ\Pi\circ O.
\]

For

\[
  Z_D
  =
  |\{x\ne0:\operatorname{wt}(E(x))<D\}|,
\]

the theorem proves

\[
  \Pr[E\text{ is noninjective or }d_{\min}(E)<D]
  \le \mathbb E[Z_D].
\]

With a deterministic outer spectrum and a uniform global permutation,

\[
  \mathbb E[Z_D]
  =
  \sum_h A_h^{\mathrm{out}}p_h^{\mathrm{in}}(D).
\]

For a factored structured interleaver, Hamming weight need not determine an
orbit. The exact replacement is the type formula

\[
  \mathbb E[Z_D]
  =
  \mathbb E_O\!\left[
    \sum_{q,\tau}A_{O,q,\tau}P_{O,q,\tau}(D)
  \right].
\]

Separate average bounds on \(A_{O,q,\tau}\) and
\(P_{O,q,\tau}\) do not bound their expected product. A valid modular proof
needs a pointwise transfer bound, independence, or a directly proved
conditional product bound.

The arbitrary-length wrapper starts from a certified admissible
\([N^+,K^+,D^+]\) code. Input restriction preserves distance. Puncturing
\(q=N^+-N\) fixed output coordinates gives distance at least \(D^+-q\).
If \(q<D^+\), puncturing preserves injectivity. The wrapper adds no failure
probability.

The asymptotic Structured SPIN target uses outer blocks

\[
  B_N=\Theta(\log N)
\]

when the weighted one-block contribution is exponentially small in \(B_N\).
The one-block class then generally has only \(\Theta(\log N)\) failure
margin. Linear-occupation classes can have \(\Theta(N)\) margin. If sparse
suppression is only polynomial in \(B_N\), a sublinear schedule
\(B_N=N^{\beta+o(1)}\) may be necessary.

## Existing reference results

These results should prevent the next chat from repeating old searches.

### Closed theorems or certificates

- Random outer plus fixed RM2Sub-S19 has a proved asymptotic rate-one-half,
  relative-distance-0.11 theorem with
  \(B_N=9\log_2N+O(1)\). This proves that the frozen inner can support the
  desired distance when paired with an ideal outer.
- One sampled Golay--BA-3 constituent, reused everywhere and paired with
  RM2Sub-S19, has a proved asymptotic rate-one-half, relative-distance-0.11
  theorem with
  \[
    B_N=(39/4)\log_2N+O(1).
  \]
  Ordinary and transposed work are linear. This is the strongest proved
  scalable structured reference, but its attempted finite \(B=240\) dense
  proof did not close at 11 or 10 percent with the tested bounds.
- A global two-sided regular expander followed by memory-15 wrapped random
  convolution has a complete finite certificate at
  \((K,N)=(2^{20},2^{21})\), distance 228590, and failure below
  \(2^{-50.2054}\). It has a linear operation count, but it is not the
  Structured SPIN architecture, has no decoder theorem, and has not been
  benchmarked.
- A repeated exact EBCH128 outer with RandomStepConv-M30 has a complete
  finite 11-percent comparator with about 26.192 failure bits.
- The unpadded repeated EBCH128 plus RandomStepConv-M20 comparator has a
  complete 10.9-percent certificate with about 7.584 failure bits.

### Important failures or limitations

- Pure bit transpose followed by one accumulator has a proved adjacency
  obstruction and cannot yield positive relative distance when the local
  block is sublinear.
- Region permutations remove that exact obstruction, but their structured
  transfer still needs a complete finite certificate.
- Repeated Golay--BA-3 plus RM2Sub-S19 is promising and asymptotically proved,
  but the attempted finite dense relaxations did not close.
- Independent row-local Fanout-56 can close a finite proof, but it violates
  the one-repeated-code interface and approximately doubles the measured
  implementation cost.
- Random Toeplitz and RandomStepConv inners are useful proof comparators.
  Their representations or online costs are not the preferred final
  implementation.
- Expected spectrum alone is insufficient when one sampled constituent is
  reused everywhere. The proof needs a high-probability realized-spectrum
  event or a direct transfer-weighted concentration theorem.
- Minimum distance and dimension alone are insufficient. The inner proof is
  sensitive to the complete shell shape and to mixed block compositions.

FINITE_K20_ATTEMPT_LOG.md contains the complete route-by-route history,
classification, numerical values, and supporting files.

## Current local task: one sparse mixer followed by an accumulator

The active local task is not yet the frozen Structured SPIN proof. It is a
proof-oriented candidate for the repeated local outer constituent.

Fix

\[
  k=256,\qquad n=512,\qquad r=33.
\]

Setup independently samples rows

\[
  R_1,\ldots,R_n
  \gets
  \{R\in\mathbb F_2^k:\operatorname{wt}(R)=r\}.
\]

Define the sparse linear map

\[
  (Ex)_j=\langle R_j,x\rangle.
\]

Let \(A\) be the zero-initialized prefix accumulator on \(n\) bits. The
candidate constituent is

\[
  C=A E:\mathbb F_2^{256}\longrightarrow\mathbb F_2^{512}.
\]

Setup samples \(E\) once. Every SPIN outer-block position reuses the same
realized \(C\). The rows of one \(E\) are independent; different SPIN
outer-block positions do not resample \(E\).

The current proposal may sample at most 16 independent candidate matrices and
accept the first full-rank matrix. Gaussian elimination performs the rank
test. This rank-conditioned setup is proposed, not yet fully certified.

Because \(r\) is odd,

\[
  E\mathbf 1_k=\mathbf 1_n,
  \qquad
  AE\mathbf 1_k=(1,0,1,0,\ldots,1,0).
\]

Every realized constituent therefore contains a deterministic weight-256
word. Any proof that replaces the row law by a uniform law on all binary
vectors is invalid.

For the realized constituent, define

\[
  N_w(E)
  =
  |\{x\in\mathbb F_2^k:\operatorname{wt}(AEx)=w\}|.
\]

When \(E\) is injective and \(w>0\), \(N_w(E)\) is the ordinary shell
multiplicity.

The desired output is a simultaneous integer cap event

\[
  \mathcal E_{\mathrm{out}}
  =
  \{\operatorname{rank}(E)=256
    \text{ and }N_w(E)\le T_w\text{ for all }w>0\},
\]

with a certified failure probability. The existing 10.9-percent SPIN cap
interface places zero caps outside weights 42 through 470. The current
variance target is

\[
  \operatorname{Var}(N_w)
  \le 512\,\mathbb E[N_w]
\]

for the defect shells 42 through 79 and, by an explicitly proved symmetry,
their complements 433 through 470. This factor 512 is sufficient for the
current three-stage cap budget. It is not claimed to be necessary.

The one-stage construction costs 16,895 XORs per local constituent in the
existing operation model. Lower-cost multi-stage diagnostics include:

| Stages | Sparse degrees | XORs | Extreme-shell diagnostic bits |
|---:|---:|---:|---:|
| 1 | \(33\) | 16,895 | 44.6379 |
| 2 | \((17,3)\) | 10,238 | 43.9504 |
| 3 | \((13,3,3)\) | 9,725 | 42.4540 |
| 3 | \((15,2,3)\) | 10,237 | 46.5549 |

These values are nearest-binary64 first-moment diagnostics for rank and
extreme shells. They do not prove simultaneous central caps. The one-stage
route is being solved first because it avoids interstage covariance
composition. If its mechanism closes, transfer it to \((15,2,3)\).

## Exact results for the current task

The row-character bias is

\[
  \beta(x)
  =
  \frac{K_r^{(k)}(\operatorname{wt}(x))}{\binom kr}.
\]

For a row-index subset \(t\subseteq[n]\), define

\[
  I_t(E)
  =
  \mathbf 1\left\{\bigoplus_{j\in t}R_j=0\right\},
  \qquad
  c_w(t)=K_w^{(n)}(\operatorname{wt}(A^{-\mathsf T}t)).
\]

Fourier inversion proves

\[
  N_w(E)
  =
  2^{k-n}\sum_t c_w(t)I_t(E)
\]

and

\[
  \operatorname{Var}(N_w)
  =
  2^{2(k-n)}c_w^{\mathsf T}C c_w,
\]

where \(C\) is the covariance matrix of the indicators \(I_t\).

The covariance kernel depends only on
\((|s|,|t|,|s\cap t|)\). A 257-state radial walk evaluates it exactly. The
kernel is invariant under permutations of the 512 row indices and therefore
decomposes into the standard \(S_{512}\) sectors.

Walsh conjugation moves the low-distance calculation to a sparse primal
level profile. Define

\[
  f_w(z)=\mathbf 1\{\operatorname{wt}(Az)=w\},
  \qquad
  D=2^{-n}WCW.
\]

At rate one half,

\[
  \operatorname{Var}(N_w)=f_w^{\mathsf T}Df_w.
\]

Let

\[
  g_{w,p}
  =
  |\{z:\operatorname{wt}(z)=p,\operatorname{wt}(Az)=w\}|.
\]

The accumulator run count computes every \(g_{w,p}\) exactly, and

\[
  p\le\min\{2w,2(n-w)+1\}.
\]

For an ordered message pair \((x,y)\), let \(P_{x,y}\) be its one-row
two-output probability matrix. The exact sector block is

\[
  (P_{x,y}^{\otimes n})^{(j)}
  =
  \det(P_{x,y})^j\operatorname{Sym}^{n-2j}(P_{x,y}).
\]

An exact diagonal insertion inequality proves that every sector diagonal is
bounded by the maximum of sectors zero, one, and two. Positive
semidefiniteness then gives

\[
  f_w^{\mathsf T}Df_w
  \le
  \left(\sum_p\sqrt{g_{w,p}d_p}\right)^2,
  \qquad
  d_p=\max_{j\in\{0,1,2\}}D^{(j)}_{p,p}.
\]

Here \(w\) is the final local codeword weight
\(\operatorname{wt}(AEx)\). The index \(p\) is the Hamming weight of the
primal variable \(z\) before applying \(A\). A bound at \(p=100\) is not a
bound on the final shell \(w=100\).

Sector zero initially suffered catastrophic cancellation. Define

\[
  q_x=\frac{1-\beta(x)}2,
  \qquad
  \delta_{x,y}
  =
  \frac{\beta(x+y)-\beta(x)\beta(y)}4,
\]

and

\[
  h_{j,p}(q)
  =
  [z^p](1-z)^j(1-q+qz)^{n-j}.
\]

The exact Hoeffding identity is

\[
  \operatorname{Cov}\!\left(
    \mathbf1\{\operatorname{wt}(Ex)=p\},
    \mathbf1\{\operatorname{wt}(Ey)=p\}
  \right)
  =
  \sum_{j=1}^n
  \binom nj\delta_{x,y}^{\,j}
  h_{j,p}(q_x)h_{j,p}(q_y).
\]

The determinant matrix \((\delta_{x,y})\) is a covariance matrix. The Schur
product theorem makes every entrywise power
\((\delta_{x,y}^{\,j})\) positive semidefinite. Consequently, every
fixed-order aggregate after summing over \(x,y\) is nonnegative.

The radial compression defines

\[
  F_{j,\ell,p}(t)
  =
  \sum_{a=0}^k
  K_a^{(k)}(t)h_{j,p}(q_a)\beta(a)^{j-\ell}.
\]

If \(v_t(\ell)\) is the per-vector probability that the XOR of \(\ell\)
sampled sparse rows has weight \(t\), the order-\(j\) aggregate equals

\[
  \binom nj4^{-j}
  \sum_{\ell=0}^j(-1)^{j-\ell}\binom j\ell
  \sum_{t=0}^k
  \binom kt\,v_t(\ell)F_{j,\ell,p}(t)^2.
\]

This identity replaces \(2^{2k}\) message pairs by 257 radial states. The
alternating sum over \(\ell\) must still be evaluated with enough precision.
Fixed-order nonnegativity does not permit dropping its negative terms.

## What is fully certified locally

Small exact verifiers establish:

- four exhaustive one-word and pair-kernel checks;
- 165 exact radial triple identities at \((k,r)=(4,3)\);
- exact covariance-block reconstruction at two small sizes;
- exact sector projection energies at two small sizes;
- 2,304 pair-shell Hoeffding identities at \((k,n,r)=(4,8,3)\);
- 72 nonnegative fixed-order aggregates; and
- 72 exact radial fixed-order identities; and
- 1,980 exact polynomial radial-transform identities.

At the target parameters, sector zero is completely certified at primal
level \(p=100\):

\[
  D^{(0)}_{100,100}
  \in
  1.0000035547087893253531579783751399012502409269071873\ldots
  \mathbin{\pm}2.14\cdot10^{-111}.
\]

The complete receipt checks one finite Arb enclosure for every Hoeffding
order \(j=1,\ldots,512\):

    workstreams/finite_asymptotic_theory/pure_expander_accumulate/sector_zero_hoeffding_radial_K256_B512_r33_p100_complete_arb.json

The first full run used 1,024-bit precision and produced valid rows through
order 511. At order 512, Python-FLINT's generic Arb exponentiation returned
nan when squaring a zero-centered interval. This was a library edge case, not
a failed inequality. The evaluator now computes a square as \(x\cdot x\).
A separate 2,048-bit order-512 run gives a contribution below
\(4.27\cdot10^{-46}\). The merger combines the finite orders and verifies
complete coverage.

The corrected evaluator reproduces the exact small-model value

\[
  0.61776065826416015625
\]

for \((k,n,r,p)=(4,8,3,2)\).

The same transform has the exact polynomial form

\[
 F_{j,\ell,p}(t)=2^{k-(n-j)}\sum_{s=0}^{n-j}
 \binom{n-j}{s}K_p^{(n)}(j+s)v_t(s+j-\ell).
\]

The exact small verifier checks all 1,980 instances. The polynomial Arb mode
matches the scalar mode exactly on the complete small model, but is slower at
the target size: 81.4 seconds versus 44.2 seconds through order 128 at level
100. Keep `--transform-mode scalar` as the production default.

The complement-orbit route does not require orbit positivity. Enclose each
signed four-term orbit and sum only its positive upper endpoint. The directed
prototype returns 0.8136835098268285 in the small model, versus nondirected
positive mass 0.8136835098266602 and exact diagonal 0.61776065826416015625.
The target implementation now quotients both message complements and message
swap. It covers all 2,862,209 pair types with 361,985 representatives.
Construct the determinant directly from its exact integer numerator: the
earlier subtraction of rounded biases inflated the upper bound to
\(7.51\cdot10^{18}\). The corrected outward receipt proves

\[
 D^{(0)}_{100,100}\le2.45522037042628
\]

in 479.6 seconds. This is looser than the independent Arb value near
1.00000355 but has ample factor-512 slack. A compiled implementation should
use a recurrence for the symmetric-power coefficients and amortize the type
traversal across multiple levels.

Do not use the full sector-zero trace as a shortcut. Low-dual-level row
collisions make it enormous, and direct long-double evaluation has
catastrophic cancellation.

## Diagnostics that must not be promoted

The following numbers guide the proof but are not complete certificates:

- The length-512 shell-\(w=42\) bound from the sector-diagonal route is
  31.0225, compared with the sufficient factor 512.
- Substituting outward endpoints only at primal levels 80 and 84 gives a
  mixed-arithmetic shell-\(42\) value 33.8177.
- The old pairwise-positive sector-zero relaxation gives bounds
  \(2.72\cdot10^6\) at \(p=90\) and \(6.97\cdot10^{12}\) at \(p=100\).
  The exact complete \(p=100\) value near one proves that this growth was
  relaxation loss.
- Four-message complement grouping gives about 1.0000032 at \(p=100\) in
  nondirected arithmetic. Small exact examples have negative complement
  orbits, so orbit positivity is not a theorem.
- Moderate-size exact variances are encouraging, but they do not extrapolate
  to \(k=256,n=512\).

Do not use the rejected weight-level maximum bound. It gives factors in the
thousands at moderate sizes because it maximizes sectors independently and
destroys necessary cancellation.

## Historical open obligations, now discharged

The following list records the obligations that led to the final proof.  All
items required for the one-stage RandomStepConv-M22 variant have now been
discharged.  They remain useful as a map of the proof dependencies.

1. Obtain complete outward sector-zero bounds for every required primal level
   through \(p=159\). Only \(p=100\) is complete.
2. Obtain outward sector-one and sector-two bounds through \(p=159\).
3. Use the exact accumulator counts and the diagonal merger to certify every
   defect shell \(w=42,\ldots,79\). State and verify the symmetry used for
   shells \(433,\ldots,470\).
4. If the diagonal bound exceeds the available cap slack, use the exact
   accumulator correlation table \(Q_w(p,\ell)\) or signed block quadratic
   forms. Do not return to the rejected weight-level maximum.
5. Convert shell means and variances into simultaneous integer caps with an
   explicit allocation of the outer failure budget.
6. Certify one-draw rank failure. For rank-tested resampling, divide the
   unconditioned spectrum failure by a proved lower bound on rank success and
   add the bounded-attempt abort probability.
7. Attach the resulting cap event to a compatible SPIN transfer proof.
   A weight-only uniform-permutation transfer cannot be applied silently to
   the frozen factored interleaver.
8. Only after the one-stage mechanism closes, transfer it to the lower-XOR
   \((15,2,3)\) candidate and compare the optimized online cost with the
   10.823872-ms frozen implementation.

Even a complete local spectrum event does not by itself prove the final SPIN
distance theorem. The final proof must bound the conditional inner failure
for every constituent satisfying the event and add all setup-failure terms.

## Final closure and recommended next action

The proof uses exact sector-zero coverage through level 159, outward sectors
one and two through level 159, direct dominant-character defect-shell bounds,
central-shell likelihood bounds, explicit simultaneous integer caps, rank
conditioning with at most 16 setup attempts, and a complete finite transfer.
The final setup margin is 41.3632966090 bits; the conditional transfer margin
is 51.6422972013 bits; their union has 41.3621359295 bits of margin.

The production central-shell command accepts `--workers`.  The eight-worker
Peach run used one FLINT thread per process and finished in 84.7095 seconds,
versus 771.2017 seconds for the serial local run.  Every one of the 353 shell
rows and the final claim matched exactly.  Never run two certificate or
benchmark jobs concurrently.

The recommended next task is to preserve this cap interface and replace
RandomStepConv-M22 with RM2Sub-S19.  If that transfer does not close, the next
choices are to certify the frozen structured routing or to port the
concentration method to the lower-XOR \((15,2,3)\) sparse-EA candidate.  These
are new theorem tasks, not unfinished steps in the certificate reported here.

## Reproduction commands

Run these commands from:

    C:\Users\peter\.codex\worktrees\3061\permute_conv

The target sector-zero evaluator requires Python-FLINT.

To evaluate one complete target level after the square fix:

    python workstreams/finite_asymptotic_theory/pure_expander_accumulate/evaluate_sector_zero_hoeffding_radial_arb.py --message-bits 256 --output-bits 512 --right-degree 33 --level P --maximum-order 512 --precision 1024 --progress-every 32 --suppress-json --output RECEIPT.json

If only a failed order needs more precision, use a strict order range:

    python workstreams/finite_asymptotic_theory/pure_expander_accumulate/evaluate_sector_zero_hoeffding_radial_arb.py --message-bits 256 --output-bits 512 --right-degree 33 --level P --minimum-order J --maximum-order J --precision 2048 --progress-every 1 --suppress-json --output ORDER_RECEIPT.json

Merge enough receipts to cover every order exactly:

    python workstreams/finite_asymptotic_theory/pure_expander_accumulate/merge_sector_zero_hoeffding_order_receipts.py RECEIPT_A.json RECEIPT_B.json --precision 2048 --output COMPLETE.json

The merger rejects missing orders and disjoint duplicate enclosures. It may
ignore a non-finite row only when another source supplies a finite enclosure
for the same order.

The existing positive diagonal verifier accepts entries as SECTOR:LEVEL:

    python workstreams/finite_asymptotic_theory/pure_expander_accumulate/certify_primal_schur_diagonal_outward.py --message-bits 256 --output-bits 512 --right-degree 33 --entry 1:P --entry 2:P --output DIAGONAL_RECEIPT.json

It can also evaluate a band:

    python workstreams/finite_asymptotic_theory/pure_expander_accumulate/certify_primal_schur_diagonal_outward.py --message-bits 256 --output-bits 512 --right-degree 33 --level-min LO --all-through HI --output BAND_RECEIPT.json

The current sector-zero implementation inside that binary64 verifier uses the
old pairwise-positive relaxation. Do not use it beyond the range where it is
actually smaller than the desired bound. Replace sector-zero entries with
the complete Arb receipts before final shell assembly.

Merge outward bands and request final local shells with:

    python workstreams/finite_asymptotic_theory/pure_expander_accumulate/merge_primal_schur_diagonal_outward.py --receipt BAND_A.json --receipt BAND_B.json --shell-weight W --output SHELL_RECEIPT.json

The current shell merger expects all controlling sector entries. It does not
yet accept a separate Arb sector-zero receipt, so that integration is an
immediate code task.

## Important local files

The main proof narrative is:

    workstreams/finite_asymptotic_theory/pure_expander_accumulate/ONE_STAGE_VARIANCE_ROUTE.md

The construction and decision gates are:

    workstreams/finite_asymptotic_theory/pure_expander_accumulate/README.md
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/PAIR_KERNEL_AND_FIRST_GATE.md

The complete history of finite alternatives is:

    workstreams/finite_asymptotic_theory/FINITE_K20_ATTEMPT_LOG.md

The overall workstream handoff is:

    workstreams/finite_asymptotic_theory/MERGE_SUMMARY.md

The central scripts are:

    workstreams/finite_asymptotic_theory/pure_expander_accumulate/evaluate_sector_zero_hoeffding_radial_arb.py
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/merge_sector_zero_hoeffding_order_receipts.py
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/certify_primal_schur_diagonal_outward.py
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/merge_primal_schur_diagonal_outward.py
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/verify_sector_zero_hoeffding_small.py
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/verify_primal_schur_power_factorization_small.py
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/verify_subset_covariance_blocks_small.py

The most important receipts are:

    workstreams/finite_asymptotic_theory/pure_expander_accumulate/sector_zero_hoeffding_K4_B8_r3_exact.json
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/sector_zero_hoeffding_radial_K4_B8_r3_p2_complete_arb_v2.json
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/sector_zero_hoeffding_radial_K256_B512_r33_p100_complete_arb.json
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/primal_schur_diagonal_K256_B512_r33_p80_84_outward.json
    workstreams/finite_asymptotic_theory/pure_expander_accumulate/primal_schur_shell_bound_K256_B512_r33_w42_probe.json

## Current worktree state

The finite-theory log and merge summary are modified. The complete
pure_expander_accumulate directory is currently untracked in this worktree.
Preserve all of it. Do not reset or discard these files.

Before handing results back, run:

    python -c "import ast,json,pathlib; root=pathlib.Path(r'workstreams/finite_asymptotic_theory/pure_expander_accumulate'); [ast.parse(p.read_text(encoding='utf-8'),filename=str(p)) for p in root.glob('*.py')]; [json.loads(p.read_text(encoding='utf-8')) for p in root.glob('*.json')]; print('parse ok')"
    git diff --check -- workstreams/finite_asymptotic_theory
    git status --short -- workstreams/finite_asymptotic_theory

## Suggested prompt for the new chat

Continue the finite/asymptotic SPIN theory work from
workstreams/finite_asymptotic_theory/NEXT_CHAT_HANDOFF.md. Read every required
startup file listed there and use Controlled Writing for Cryptography. Work
only inside workstreams/finite_asymptotic_theory/. Preserve the frozen
Structured SPIN baseline and the single-reused-constituent semantics. The
one-stage RandomStepConv-M22 certificate is complete, with distance 228,590
and 41.3621 bits of total margin. Audit it first from
ONE_STAGE_SPARSE_EA_CERTIFICATE_MANIFEST.json. Then attempt the smallest next
extension: preserve the certified outer caps and replace RandomStepConv-M22
by RM2Sub-S19. Keep the real objective visible: a practical finite distance
certificate for frozen Structured SPIN, or a close linear-time variant, at
\(K=2^{20}\), \(N=2^{21}\), approximately 11-percent distance, at least 40
bits of total margin, and competitive implementation cost.
