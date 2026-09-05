# Finite random bounded-memory convolution baseline

## Result

At message dimension

\[
k=2^{20}=1{,}048{,}576,
\]

output length

\[
N=240\cdot 8832=2{,}119{,}680,
\]

and cutoff

\[
D=\lfloor 0.11N\rfloor=233{,}164,
\]

the random-outer/RandomStepConv ensemble has a large first-moment distance
margin. Nearest-binary64 evaluation gives 1,422.7138 bits at memory 10 and
11,421.6932 bits at memory 19. Memory 9 does not close.

This result isolates the inner question. A random bounded-memory linear
state machine is strong enough at this length. The result does not yet give
the desired repeated-BA outer, and the numerical calculation is not an
outward-rounded certificate.

## Encoder and probability space

Fix a memory size \(M\geq 1\). Sample a uniform binary matrix

\[
G\gets\mathbb F_2^{k\times N}.
\]

Independently, for every bit position \(t\in\{0,\ldots,N-1\}\), sample

\[
A_t\gets\mathbb F_2^{(M+1)\times(M+1)}.
\]

All matrix entries are independent fair bits. The setup samples \(G\) and
the \(A_t\)'s once. The resulting matrices are shared by all messages.

For a message \(u\in\mathbb F_2^k\), first compute \(x=uG\). Set
\(s_0=0\in\mathbb F_2^M\). Then compute

\[
  (y_t,s_{t+1})=A_t(x_t,s_t)
  \qquad (0\leq t<N).
\]

The encoder returns \(y=(y_0,\ldots,y_{N-1})\) and discards \(s_N\). This is
the bitwise, or \(g=1\), RandomStepConv model. It is a time-varying linear
finite-state encoder. The word *random function* means a random linear map
here; a uniformly random arbitrary function would define a nonlinear code
and a different distance problem.

The maps \(A_t\) are independent across positions. Reusing one map at every
position would also define a different ensemble.

## Exact two-state transfer

For \(0<z<1\), put

\[
q=2^{-M},\qquad b(z)=\frac{1+z}{2}.
\]

Classify the state as zero or live. If \((x_t,s_t)\ne0\), then
\(A_t(x_t,s_t)\) is uniform in \(\mathbb F_2^{M+1}\). Consequently, its output
bit is fair, its next state is uniform, and those two values are independent.
The output-weight-tilted transfers for a fixed zero or one input bit are

\[
T_0(z)=
\begin{pmatrix}
1&0\\
qb&(1-q)b
\end{pmatrix},
\qquad
T_1(z)=
\begin{pmatrix}
qb&(1-q)b\\
qb&(1-q)b
\end{pmatrix}.
\]

Rows give the current zero/live class and columns give the next zero/live
class. In particular, zero input from zero state remains deterministically
zero. Every other input-state pair invokes a uniform random image.

The existing RandomStepConv verifier exhaustively enumerates small maps and
reports exact agreement for \((g,M,N,a)=(1,1,3,1)\), \((2,1,2,1)\), and
\((2,1,2,3)\). The command is

```text
python scripts/analyze_riffle_randomstepconv_g4_sigma20_goal01.py --mode check --no-write
```

## All-message random-outer bound

Fix a nonzero message \(u\). Because \(G\) is uniform, the coordinates of
\(x=uG\) are independent fair bits. Averaging \(T_0\) and \(T_1\) gives

\[
P_M(z)=\frac{T_0(z)+T_1(z)}2
=
\begin{pmatrix}
\frac12+\frac{qb}{2}&\frac{(1-q)b}{2}\\
qb&(1-q)b
\end{pmatrix}.
\]

Let \(e_0=(1,0)\) and \({\bf 1}=(1,1)^\mathsf T\). For the number \(Z_D\) of
nonzero messages whose encoded word has weight at most \(D\), the Chernoff
bound and linearity of expectation give

\[
\mathbb E[Z_D]
\leq
(2^k-1)z^{-D}e_0P_M(z)^N{\bf 1}.
\tag{1}
\]

No independence between different messages is used. The setup matrices are
shared exactly as they would be in one sampled code. Markov's inequality
then gives

\[
\Pr[d_{\min}\leq D]\leq\mathbb E[Z_D].
\tag{2}
\]

This proof also counts a nonzero message that the random outer maps to zero.
Thus a successful bound implies injectivity as well as the distance claim.

## Finite evaluation

The evaluator minimizes the right side of (1) over \(z\). It powers the
positive two-state matrix with scaling and uses nearest binary64 arithmetic.

| Memory \(M\) | \(\log_2\mathbb E[Z_D]\) upper diagnostic | Margin (bits) |
|---:|---:|---:|
| 8 | 23,944.8599 | -23,944.8599 |
| 9 | 7,680.1952 | -7,680.1952 |
| 10 | -1,422.7138 | 1,422.7138 |
| 12 | -8,834.2152 | 8,834.2152 |
| 16 | -11,277.0897 | 11,277.0897 |
| 19 | -11,421.6932 | 11,421.6932 |
| 23 | -11,441.0779 | 11,441.0779 |

At \(M=19\), the minimizing parameters are approximately

\[
-\log z=2.0906618952,
\qquad z=0.1236052950.
\]

Optimization is not needed to see closure. The fixed rational tilt
\(z=1/8\) gives 1,332.6908 bits at \(M=10\) and 11,402.7900 bits at
\(M=19\). This fixed-tilt form is the natural starting point for an outward
or exact-rational certificate.

The ideal random-code comparator at this length has 11,454.52 bits. Thus
RandomStepConv with \(M=19\) loses only about 32.83 bits against that
comparator.

For comparison, before shortening, the parent dimension is
\(k_+=1{,}059{,}840=N/2\). At that dimension, the same calculation gives
157.6932 bits at \(M=19\) and 177.0779 bits at \(M=23\). The large target-code
margin therefore comes mainly from the 11,264 shortened message coordinates,
as it does for the ideal random comparator.

## One-active repeated-BA diagnostic

The first architecture-faithful calculation keeps the repeated
Golay--BA-3 outer with \(B=240\). A coordinate permutation is sampled for
the active outer row. After transpose, every one of the 240 regions receives
an independent uniform permutation of its \(L=8832\) positions. The
RandomStepConv maps are sampled independently of those objects.

For a constituent word of weight \(w\), exactly \(w\) regions contain one
active input bit. In each such region, that bit has a uniform position. The
calculation averages

\[
\frac1L\sum_{j=0}^{L-1}T_0(z)^jT_1(z)T_0(z)^{L-1-j}
\]

for an active region and extracts the exact coefficient for \(w\) active
regions. It then sums over the expected BA spectrum.

Let \(\mathcal G_{240}\) be the event that the sampled BA constituent has no
nonzero word outside weights 23 through 217. The expected BA tail count and
Markov's inequality give

\[
\Pr[\mathcal G_{240}]\geq0.9936608994938328.
\]

Dividing the allowed-shell expectation by this lower bound gives the
following valid conditional first-moment diagnostics for occupation one.

| Memory \(M\) | Unconditional BA margin | Margin conditional on \(\mathcal G_{240}\) |
|---:|---:|---:|
| 19 | 4.3351 | 21.5360 |
| 20 | 5.6836 | 29.5748 |
| 21 | 6.7111 | 35.0931 |
| 22 | 7.4452 | 38.9445 |
| 23 | 7.9273 | 41.6332 |

At \(M=23\), weight 23 is the dominant allowed constituent shell. Rare BA
codes with words of weights 2 through 5 dominate the unconditional ensemble
average. This is the same heavy-tail phenomenon found in the Toeplitz study;
it is not a defect of a typical selected BA constituent.

This section covers only messages with one active outer row. It is not a
distance certificate for the repeated-BA construction.

## Interpretation

The full random-outer result answers the baseline question: a fresh random
bounded-memory step is already almost as strong as a fully random inner at
the finite target. The RM2Sub difficulty is therefore not caused by bounded
memory alone. It comes from the structured state interface, especially the
large outer subcode that is silent to every RM2Sub syndrome observation.

RandomStepConv avoids a fixed silent syndrome kernel. At every nonzero
input-state pair, a fresh random map produces a uniform next output-state
pair. A fixed outer word can still experience random returns to zero state;
the two-state transfer accounts for those returns exactly.

The model is asymptotically linear time for fixed \(M\), but it is only a
proof baseline. Storing unrestricted step maps takes \(N(M+1)^2\) random
bits: about 101.1 MiB at \(M=19\) and 145.5 MiB at \(M=23\). A direct dense
implementation also performs \(O(NM^2)\) bit operations. Replacing the maps
by pseudorandomly generated maps changes the information-theoretic ensemble
unless the proof states an additional computational assumption.

## Status and remaining obligations

### Proved symbolically

- the exact two-state transfer for the stated random linear map;
- the all-message first-moment inequality (1);
- the conversion from expected bad messages to setup-failure probability;
- the exact one-active region average and coefficient formula.

### Supported numerically

- the finite margins in both tables;
- the threshold between \(M=9\) and \(M=10\) for the global random outer;
- the 40-bit occupation-one crossing between \(M=22\) and \(M=23\) for the
  conditioned repeated-BA outer.

### Open

1. Replace binary64 optimization and matrix powering by an outward-rounded
   interval cover.
2. Cover occupations \(Q=2,\ldots,L\) for the one-repeated-BA construction.
3. Freeze or sample one BA constituent under a setup rule that proves the
   required spectrum event without treating independent rows as independent
   BA draws.
4. Decide whether the random-step model is only a proof comparator or a
   candidate construction, and account for its setup representation and
   online cost accordingly.
5. Use the successful transfer as the reference law for a modified RM2Sub
   design, then prove a comparison theorem for the modified structured step.

## Reproduction

The all-message evaluator and receipts are

- `evaluate_random_outer_randomstepconv_g1.py`;
- `random_outer_randomstepconv_g1_k20_d11.json`;
- `random_outer_randomstepconv_g1_kparent_d11.json`.

The occupation-one evaluator and receipts are

- `evaluate_ba_randomstepconv_g1_one_active.py`;
- `ba240_randomstepconv_g1_s19_q1_d11.json` through
  `ba240_randomstepconv_g1_s23_q1_d11.json`.

## Follow-up closure with one repeated random constituent

`FINITE_K20_REPEATED_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md` closes the
all-occupation finite comparison. It uses one sampled and repeated
\([512,256]\) constituent, 4,096 active information rows, 44 zero-shortened
parent rows, and RandomStepConv-M30. The outward verifier proves

\[
 \Pr[d_{\min}<233{,}165]<2^{-108}
\]

conditional on the stated constituent-spectrum event. One random constituent
satisfies that event with certified probability greater than
0.6365388449160215. A bounded 28-attempt selection procedure gives a combined
setup-abort-or-distance-failure probability below $2^{-40}$.

This follow-up resolves the first outward-arithmetic obligation for the
random-outer comparison model. It does not resolve the repeated-BA
occupations or supply an efficient outer-spectrum test. The next design task
is to compare an implementable structured inner to the successful transfer
law, while retaining the single-repeated-constituent interface.
