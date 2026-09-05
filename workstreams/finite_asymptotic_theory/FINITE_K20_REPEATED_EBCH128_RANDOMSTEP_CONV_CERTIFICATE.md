# Historical padded repeated-EBCH128 RandomStepConv certificate

## Result

This certificate records the earlier padded instance with \(L=16560\). It
does not certify the current power-of-two instance with \(L=16384\).

One fixed \([128,64,22]\) extended BCH constituent suffices for the
RandomStepConv comparison model. The same constituent is used in every outer
row. No outer code is sampled, and setup does not enumerate its codewords.

At message dimension \(k=2^{20}\), the outward verifier proves

\[
 \Pr_\omega[d_{\min}(C_\omega)<233165]
 <2^{-26.1921836158}<2^{-20}.
 \tag{1}
\]

Consequently, outside an event of probability less than \(2^{-20}\),

\[
 \frac{d_{\min}(C_\omega)}{2119680}
 \ge \frac{233165}{2119680}
 =0.11000009435386474\ldots.
\]

The theorem concerns RandomStepConv-M30. It does not establish the same
claim for RM2Sub.

## Construction and probability space

Set

\[
 B=128,\quad K=64,\quad L=16560,\quad L_0=16384,
\]

and hence

\[
 k=KL_0=2^{20},\qquad N=BL=2119680.
\]

Use the same concrete constituent and coordinate convention as
`FINITE_K20_EBCH128_OUTER_DEFINITION.md`, but retain the historical padding
specified below. Its cyclic generator polynomial is
`0xf4845518b9582a1f`; coordinate 127 is the extension-parity coordinate. The
complete spectrum consumed by the proof is stored in
`scripts/EBCH128_64.wd`. Its nonzero minimum weight is 22. The spectrum is
complement symmetric, has total mass \(2^{64}\), and contains one all-one
word. The theorem uses the explicit premise that this imported enumerator is
the enumerator of the stated concrete generator. The current local verifier
authenticates the table but does not derive that binding.

Encode each of the first \(L_0\) information rows with this fixed code. Set
the remaining \(L-L_0=176\) information rows to zero. For every outer row,
sample an independent uniform permutation of its 128 code coordinates. For
every coordinate region, sample an independent uniform permutation of its
\(L\) row positions. Concatenate the 128 regions. Equations (1)--(3) of the
outer-definition document fix the input blocking, permutation convention,
and every output coordinate.

RandomStepConv has state \(s_t\in\mathbb F_2^{30}\) and initial state
\(s_0=0\). At each position \(t\), sample an independent uniform binary
\(31\)-by-\(31\) matrix \(M_t\) and set

\[
 (y_t,s_{t+1}) := M_t(x_t,s_t).
\]

The encoder emits \(y_t\) and discards the terminal state. The random object
\(\omega\) consists of all row permutations, region permutations, and
matrices \(M_t\). The encoder samples \(\omega\) once and reuses it for every
message.

## RandomStepConv transfer

Fix a Chernoff parameter \(z\in(0,1)\). Put

\[
 q_0:=2^{-30},\qquad b:=\frac{1+z}{2}.
\]

The two transfer states record whether the two convolution states agree.
For input difference zero or one, respectively, the tilted transfers are

\[
 T_0=
 \begin{pmatrix}
 1&0\\ q_0b&(1-q_0)b
 \end{pmatrix},
 \qquad
 T_1=
 \begin{pmatrix}
 q_0b&(1-q_0)b\\ q_0b&(1-q_0)b
 \end{pmatrix}.
 \tag{2}
\]

For a fair candidate input, define \(C:=(T_0+T_1)/2\).

The transfer is monotone under input deletion. To see this, process an
arbitrary suffix backwards. Its continuation vector satisfies
\(v_{\mathrm{zero}}\ge v_{\mathrm{live}}\). Equation (2) then gives

\[
 T_0v\ge T_1v
\]

coordinatewise. Replacing an input difference by zero can therefore only
increase the expected tilted moment.

## Spectrum comparison

Let \(A_w\) denote the exact BCH spectrum. Excluding the zero and all-one
words, define

\[
 \Gamma:=
 \max_{0<w<128}
 \frac{A_w2^{127}}{\binom{128}{w}}.
 \tag{3}
\]

The exact spectrum gives

\[
 \Gamma=
 \frac{82008050427946169694673280391056138960896}
      {4242647163381783985125},
 \qquad
 \log_2\Gamma=64.0674346406\ldots.
\]

Weights 24 and 104 maximize (3). Thus the permuted counting measure of the
ordinary BCH words is pointwise bounded by \(\Gamma\) times a uniform even
row. The all-one word is handled separately. Monotonicity permits its input
to be thinned to the same reference row. Hence every active row contributes
at most

\[
 \Lambda:=\Gamma+1
 \tag{4}
\]

times the reference distribution.

## Dispersed parity pivots

A uniform even row admits the following equivalent generation procedure.
Sample a pivot \(J\) uniformly from the 128 coordinates. Sample the other
127 bits independently and uniformly. Set the pivot bit to their parity.
The resulting row is uniform over the even-weight vectors.

The proof deletes the pivot bit. This deletion is valid by monotonicity. The
pivot is an auxiliary proof variable; the encoder does not sample it.

Consider a message with \(Q\) active outer rows. Sample the proof pivots
independently. Let \(M_j\) count the rows whose pivot equals coordinate
region \(j\). Then

\[
 (M_1,\ldots,M_{128})
 \sim\operatorname{Multinomial}(Q;1/128,\ldots,1/128).
\]

Region \(j\) contains \(Q-M_j\) fair candidate positions. After the region
permutation, those positions form a uniform subset of the \(L\) positions.
For \(0\le a\le L\), define

\[
 R_a:=\binom La^{-1}[u^a](T_0+uC)^L.
 \tag{5}
\]

Conditional on the pivot loads, the reference moment is

\[
 e_0^{\mathsf T}
 \left(\prod_{j=1}^{128}R_{Q-M_j}\right)\mathbf1.
 \tag{6}
\]

This formula keeps the parity-dependent coordinate dispersed across the
regions. It avoids the artificial zero-input region in the earlier proof.

## Positive coefficient bound

Fix positive witnesses \(x\) and \(r\). Positivity in (5) gives

\[
 R_a\le
 \frac{x^{-a}}{\binom La}(T_0+xC)^L
 \tag{7}
\]

entrywise. Put

\[
 A_x:=(T_0+xC)^L,
 \qquad
 G_Q(r):=
 \sum_{m=0}^Q
 \frac{r^m}{m!\binom L{Q-m}}.
\]

Applying (7) in (6), averaging the multinomial pivot loads, and bounding a
positive degree-\(Q\) coefficient at \(r\) gives

\[
 \mathbb E[z^W\mid Q]
 \le
 \frac{Q!}{128^Q}
 x^{-127Q}r^{-Q}G_Q(r)^{128}
 e_0^{\mathsf T}A_x^{128}\mathbf1.
 \tag{8}
\]

The verifier further uses

\[
 G_Q(r)le
 (Q+1)\max_{0\le m\le Q}
 \frac{r^m}{m!\binom L{Q-m}}.
 \tag{9}
\]

It locates the maximum in (9) from exact dyadic ratio comparisons. Thus the
dense proof does not rely on floating-point summation of \(G_Q(r)\).

## Occupation sum

For each \(Q\), the proof multiplies the conditional bound by

\[
 \binom{L_0}{Q}\Lambda^Q.
\]

It then applies the Chernoff factor \(z^{-D}\), where

\[
 D=\lceil0.11N\rceil=233165.
\]

Occupation one uses the exact BCH spectrum rather than (3). Occupations 2
through 99 use the exact region coefficient recurrence and a fixed aligned
parity deletion. Occupations 100 through 16384 use (8) and (9). Summing all
occupation bounds gives (1). Occupation one is dominant; occupation two is
the next-largest contribution.

## Verification

The diagnostic optimizer stores every Chernoff and coefficient witness as an
exact binary64 hexadecimal value. The outward verifier treats each stored
value as an exact dyadic rational. It uses Arb for exponentials, logarithms,
and factorial expressions. Every addition and multiplication in a positive
binary64 matrix recurrence advances by one ULP toward positive infinity.

The verifier reports

\[
 \log_2\Pr_\omega[d_{\min}<D]
 <-26.1921836158846.
\]

The receipt is
`ebch128_randomstepconv_g1_s30_parity_pivot_outward_d11.json`. The manifest
`FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_MANIFEST.json` binds the exact
concrete forward circuit, imported spectrum, evaluator, verifier, and
receipt. It labels the spectrum-generator binding as an imported premise
rather than a locally derived fact.

## Scope

This result is an information-theoretic certificate for the stated random
linear inner. A literal setup stores \(N\) independent 31-by-31 matrices and
performs dense state updates. The theorem does not claim that this inner is
an efficient replacement for RM2Sub. Its role is to prove that the fixed
repeated BCH outer has sufficient spectrum once the inner handles the state
interface uniformly.
