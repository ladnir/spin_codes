# One-stage sparse-EA SPIN finite certificate

## Result

This certificate proves a finite distance statement for a new SPIN variant.
It does not certify the frozen Structured SPIN construction.

Set

\[
 k=2^{20},\qquad N=2^{21},\qquad
 D=\lceil0.109N\rceil=228{,}590.
\]

The setup below returns a binary linear \([N,k]\) code or aborts. Counting an
abort as failure, its output satisfies

\[
 \Pr[d_{\min}<D\text{ or setup aborts}]
 <2^{-41.3621359294}<2^{-40}.
 \tag{1}
\]

The probability in (1) includes the outer constituent, bounded rank-tested
setup, routing permutations, and RandomStepConv maps.

## Construction

Let \(K=256\), \(B=512\), and \(L=4096\). View a message as a matrix
\(U\in\mathbb F_2^{L\times K}\).

One outer-setup attempt samples independent sets

\[
 S_j\gets\binom{[K]}{33}
 \qquad(0\le j<B)
\]

and defines \(E:\mathbb F_2^K\to\mathbb F_2^B\) by

\[
 (Eu)_j:=\bigoplus_{i\in S_j}u_i.
\]

Let \(A\) be the zero-initialized prefix accumulator,

\[
 (Av)_j:=\bigoplus_{i=0}^{j}v_i.
\]

The attempted constituent is \(C:=AE\). Setup makes at most 16 independent
attempts and accepts the first constituent of rank 256. It aborts if all
attempts have smaller rank. Every outer row uses the same accepted map \(C\).

For each outer row, setup samples one uniform permutation of its 512 output
coordinates. The encoder applies these permutations and transposes the
\(4096\)-by-\(512\) array. Setup then samples one uniform permutation of the
4096 positions in each transposed region.

The inner encoder is RandomStepConv with memory 22. At each of the \(N\)
positions, setup samples an independent uniform binary \(23\)-by-\(23\)
linear map. The state starts at zero, and the final state is discarded. All
sampled maps are fixed and shared by every encoded message.

## Realized-spectrum event

For a constituent \(C\), define

\[
 A_w(C):=
 \left|\{u\in\mathbb F_2^{256}\setminus\{0\}:
              \operatorname{wt}(Cu)=w\}\right|.
\]

The cap verifier constructs integers \(T_0,\ldots,T_{512}\). Their positive
support is exactly \(42\le w\le470\). The good event is

\[
 \mathcal G:=
 \{\operatorname{rank}(C)=256\text{ and }
   A_w(C)\le T_w\text{ for every }w\}.
\]

For one unconditioned attempt, Markov's inequality bounds the kernel and the
zero-cap tails. Cantelli's inequality bounds every positive-cap shell. After
conditioning on full rank and adding the 16-attempt abort event, the outward
receipt proves

\[
 \Pr[\neg\mathcal G\text{ or setup aborts}]
 <2^{-41.3632966090}.
 \tag{2}
\]

The separate contributions have margins 93.9393 bits for rank failure,
44.6380 bits for zero-cap tails, and 41.5206 bits before rank conditioning
for the positive-cap shells.

## Moment bound

For a nonzero message \(x\), define the one-row character bias

\[
 \beta(x):=
 \frac{K_{33}^{(256)}(\operatorname{wt}(x))}{\binom{256}{33}}.
\]

For distinct nonzero \(x,y\), put \(z=x+y\). The exact output law of one
sparse row is

\[
 \Pr[(Ex,Ey)=(a,b)]
 =\frac14\left(
  1+(-1)^a\beta(x)+(-1)^b\beta(y)
   +(-1)^{a+b}\beta(z)
 \right).
 \tag{3}
\]

Retain one of the three biases with maximum absolute value. The corresponding
reference law keeps that marginal or correlation and sets the other two
biases to zero. If the retained magnitude is \(b<1\) and the other magnitudes
sum to \(s\), every atom of (3) is at most

\[
 1+\delta,\qquad \delta:=\frac{s}{1-b},
\]

times its reference atom. Independence of the 512 sparse rows gives a
likelihood factor \((1+\delta)^{512}\). The verifier also handles the
all-ones message and complement pairs separately.

The accumulator is invertible. Fourier inversion therefore reduces each
reference shell probability to a one-variable polynomial in the retained
bias. Arb encloses these polynomial values at 1,536-bit precision.

On the defect shells, the verifier subtracts certified lower endpoints for
the independent marginals. It proves

\[
 \operatorname{Var}(A_w)
 \le F_w\,\mathbb E[A_w],
 \qquad
 \begin{cases}
  42\le w\le79,\\
  433\le w\le470,
 \end{cases}
 \tag{4}
\]

where \(\max_w F_w\le245.145100\). The high shells were computed directly;
the proof assumes no complement symmetry of a realized constituent.

Near the central shells, direct subtraction of two ratios close to one loses
the covariance scale in binary64. The central verifier instead asks Arb for
the absolute deviations of the reference ratios from one. It rewrites each
covariance excess as a sum of nonnegative deviation and likelihood terms.
This calculation covers every \(80\le w\le432\). It does not claim the
factor-512 bound there. The larger central variance bounds still yield caps
within the available transfer envelope.

An exact rational small model checks every ordered message pair and shell at
\((K,B,r)=(4,8,3)\). It verifies 1,638 likelihood inequalities and 1,638
covariance-excess inequalities.

## Cap allocation and transfer

The cap for each positive shell is the smaller of two integers:

1. the cap with Cantelli failure at most \(2^{-51}\); and
2. the largest cap allowed by the frozen low, central, or high band envelope.

The second constraint preserves the existing transfer hypotheses for every
occupation \(Q\ge3\). The verifier recomputes occupations \(Q=1\) and
\(Q=2\) with the exact integer caps. It obtains

\[
\begin{aligned}
 \log_2\Pr[\text{bad},Q=1\mid\mathcal G]&\le-61.3895125311,\\
 \log_2\Pr[\text{bad},Q=2\mid\mathcal G]&\le-105.2734740108,\\
 \log_2\Pr[\text{bad},3\le Q\le159\mid\mathcal G]
   &\le-51.6796100317,\\
 \log_2\Pr[\text{bad},160\le Q\le4096\mid\mathcal G]
   &\le-57.0011587683.
\end{aligned}
\]

Their outward sum is below \(2^{-51.6422972013}\). Combining this value with
(2) proves (1).

## Independent sector receipts

The earlier Schur-sector route is also complete at the requested primal
levels. Separate outward receipts cover:

- sector zero at every level 1 through 159;
- sector one at every level 1 through 159; and
- sector two at every level 2 through 159.

The diagonal shell relaxation built from these entries closes weights 42
through 53 but fails thereafter. A lower-bound receipt proves that this
relaxation first exceeds 512 at weight 65. The final variance theorem uses
the dominant-character likelihood bound instead. The sector receipts remain
valid independent bounds and document why the diagonal relaxation was not
used to claim closure.

## Parallel validation

The central verifier supports independent shell workers. Each worker uses one
FLINT thread, and the launch sets the BLAS and OpenMP thread counts to one.
On Peach, eight workers reduced wall time from 771.202 seconds to 84.710
seconds. The 353 shell rows and the claim object matched the serial receipt
exactly.

## Reproduction

Run the following commands from the repository root:

```text
python workstreams/finite_asymptotic_theory/pure_expander_accumulate/verify_dominant_character_deviation_small.py --message-bits 4 --output-bits 8 --right-degree 3 --output workstreams/finite_asymptotic_theory/pure_expander_accumulate/dominant_character_deviation_K4_B8_r3_all_shells_exact.json
python workstreams/finite_asymptotic_theory/pure_expander_accumulate/certify_dominant_character_deviation_outward.py --message-bits 256 --output-bits 512 --right-degree 33 --shell-min 80 --shell-max 432 --precision 1536 --workers 8 --progress-every 0 --output workstreams/finite_asymptotic_theory/pure_expander_accumulate/dominant_character_deviation_K256_B512_r33_w80_432_outward.json
python workstreams/finite_asymptotic_theory/certify_one_stage_sparse_ea_caps_outward.py
python workstreams/finite_asymptotic_theory/certify_one_stage_sparse_ea_transfer_outward.py
python workstreams/finite_asymptotic_theory/audit_one_stage_sparse_ea_certificate.py
```

The low and high defect receipts are already listed in the manifest. Their
generation commands are recorded by the companion scripts.

## Scope

This theorem concerns the one-stage sparse-EA outer, uniform routing, and
RandomStepConv-M22. It is a close SPIN proof variant, not the frozen
Structured SPIN construction. It does not prove the same result for RM2Sub,
the frozen factored interleaver, or the frozen BCH-based outer constituent.

The RandomStepConv setup is information-theoretic and impractically large.
The next mathematical task is to transfer the realized-spectrum mechanism to
a practical inner or to the lower-XOR multistage sparse-EA outer.
