# Repeated EBCH128x4--BA constituent: concentration status

## Result

There are now two distinct targets. Three accumulator stages remain the
performance-oriented target, but their concentration proof is open. A
universal 128-stage fallback now closes the complete finite transfer,
including every occupation (1\le Q\le4096), conditional on one explicit
accumulator singular-value lemma. Neither route is yet unconditional.

The three-stage implication is as follows.

Let \(C_0\) be the direct sum of four copies of the fixed extended BCH
\([128,64,22]\) code. Sample three independent uniform permutations of 512
coordinates. Before each of three accumulator stages, apply one sampled
permutation. Let \(C_3\) be the resulting \([512,256]\) code, and define

\[
 A_w(C_3)=|\{c\in C_3\setminus\{0\}:\operatorname{wt}(c)=w\}|.
\]

Suppose that, for every \(1\le w\le512\),

\[
 \operatorname{Var}(A_w(C_3))
 \le 512\,\mathbb E[A_w(C_3)]. \tag{1}
\]

Then a shellwise Markov--Cantelli argument produces a spectrum event whose
failure probability is at most \(2^{-42.255553}\). Conditional on that event, the
nearest-binary64 RandomStepConv-M22 transfer has at least 46.11499 bits of
margin for every occupation \(1\le Q<160\). Sampled dense compositions have
at least 559.78 bits of margin. These results make (1) a plausible sufficient
target, but two obligations remain:

1. prove (1) for the stated three-stage ensemble; and
2. replace the dense samples and binary64 calculations by a complete outward
   all-occupation verifier.

If the dense cover retains the observed slack, the union of the spectrum-event
failure and a \(2^{-46.11499}\) conditional distance bound would be below
\(2^{-42.1594}\). Thus this route has enough room for the requested 40-bit
failure claim.

## Construction and probability space

The fixed local BCH code has length 128, dimension 64, and minimum distance
22. Its imported full weight enumerator is `scripts/EBCH128_64.wd`. The base
code

\[
 C_0=C_{\rm BCH}^{\oplus4}
\]

has length 512 and dimension 256.

For \(j\in\{1,2,3\}\), setup independently samples
\(P_j\gets S_{512}\). Let \(\operatorname{Acc}\) be the invertible binary
prefix accumulator with zero initial state. Recursively define

\[
 C_j=\operatorname{Acc}(P_j C_{j-1}).
\]

The phrase *three accumulator stages* always means the three maps in this
display. Some coding papers would call the same construction BA-4 because
they count the base block code as one constituent.

Setup samples \(P_1,P_2,P_3\) once. The same realized code \(C_3\) is repeated
in all \(L=4096\) SPIN outer rows. There is no rowwise resampling.

After outer encoding, each row receives an independent uniform
512-coordinate permutation. The array is bit-transposed, and each of the 512
length-4096 regions receives an independent uniform permutation. The inner is
RandomStepConv with memory \(M=22\). All permutations and inner step maps are
sampled once and are shared by every message. The finite target is

\[
 k=2^{20},\qquad N=2^{21},\qquad
 D=\lceil0.109N\rceil=228{,}590.
\]

## Proved algebraic facts

The fourth power of the imported BCH weight-enumerator polynomial gives the
base spectrum exactly. The accumulator input-output enumerator gives the
ensemble expected spectrum after each uniform interleaver and accumulator.
The implemented recurrence preserves total expected mass \(2^{256}\).

The imported BCH enumerator is formally self-dual. The dual expected-spectrum
recurrence uses

\[
 \Pr[\operatorname{wt}(A^{-\mathsf T}Px)=w
       \mid\operatorname{wt}(x)=h]
 =\frac{N_{512}(w,h)}{\binom{512}{h}},
\]

where \(N_{512}(w,h)\) is the ordinary accumulator input-output enumerator.
This identity follows by bijection counting. The numerical spectra in the
receipts use binary64 log arithmetic and are therefore diagnostics rather
than outward certificates.

The exact MacWilliams transform also proves that the local dual distance is
22. The direct sum \(C_0\) therefore has dual distance 22. If \(U,V\) are
independent uniform messages for \(C_0\), then for any set of at most 21
coordinates the pair symbols

\[
 ((UG_0)_i,(VG_0)_i)\in\mathbb F_2^2
\]

are independent and uniform. Equivalently, the base pair-type distribution
matches every multinomial falling-factorial moment of total degree at most
21. This is a proved restriction on the initial pair modes, not a heuristic
random-code approximation.

The variance sum uses distinct nonzero message pairs. Removing the events
\(U=0\), \(V=0\), and \(U=V\) from the independent-message law changes its
mass by exactly

\[
 3\cdot2^{-256}-2\cdot2^{-512}.
\]

Thus the rank-two pair distribution inherits the same low-degree moment
identities up to an explicit correction of this size. The diagonal pair term
in (2) is already \(\mathbb E[A_w]\).

The ordinary uniform-random \([512,256]\) comparator does not require a
biased sampler. Its already certified one-shot spectrum event holds with
probability greater than \(1-2^{-42.2614729204}\). That fact concerns the
uniform random-code ensemble; it does not itself imply concentration for
the BA ensemble.

A compact comparator now gives the same variance premise exactly. Over
\(\mathbb F_{2^{256}}\), sample \(a,b\) and encode
\(u\mapsto(u,au+bu^2)\). Distinct nonzero messages have independent uniform
parity halves because the relevant determinant is \(uv(u+v)\ne0\). The
resulting 512-bit sampler closes the same finite theorem with
42.2593129905 bits overall. This algebraic result also does not imply the BA
concentration lemma. It isolates pairwise independence as the exact property
that the BA proof must replace.

There is also a certified fixed-base fallback. Apply
\(x\mapsto ax+bx^2\) over \(\mathbb F_{2^{512}}\) after the four EBCH blocks.
This map gives pairwise-uniform images for every fixed pair of distinct
nonzero base words. Its restriction loses rank with probability below
\(2^{-256}\). The prior cap-conditional transfer then closes with the same
42.2593129905-bit margin and needs no accumulator stage. This fallback adds
two field multiplications per row and is not the BA-only construction studied
below.

## Conditional spectrum event

Write \(\mu_w=\mathbb E[A_w(C_3)]\), and set
\(\delta=2^{-51}\). Under (1), construct an integer cap \(T_w\) as follows.
If \(\mu_w\le\delta\), set \(T_w=0\) and apply Markov's inequality. Otherwise,
choose the smallest cap used by the receipt so that Cantelli's inequality,
with variance upper bound \(512\mu_w\), bounds
\(\Pr[A_w>T_w]\) by at most \(\delta\).

For every shell, the cap is the smaller valid integer cap supplied by
Markov's inequality or Cantelli's inequality. This correction matters in the
rare shells: at weight 42, Markov gives cap 2 while Cantelli alone gives cap
33. The union of the 512 shell-failure bounds is at most
\(2^{-42.255553}\). The nonzero caps have support from weight 42 through
weight 470. The caps imply the two retained categories used by the transfer:

| Category | Weight interval | Bernoulli parameter | Log2 majorant |
|---|---:|---:|---:|
| defect | symmetric tails | 0.0438836791 | 64.697919 |
| central | 80 through 432 | 0.50039779997 | 257.274022 |

The defect category merges the low and high tails using the proved
monotonicity of the RandomStepConv tilted moment.

## Conditional finite transfer

The following table reports nearest-binary64 calculations under (1). Exact
means exact enumeration of shell weights or category compositions within the
stated interface; it does not mean outward-rounded arithmetic.

| Occupations | Method | Aggregate margin |
|---|---|---:|
| \(Q=1\) | exact shell transfer | 71.741239 bits |
| \(Q=2\) | exact shell-pair transfer | 152.726349 bits |
| \(Q=3,\ldots,16\) | exact shared two-category recurrence | 46.114992 bits |
| \(Q=17,\ldots,32\) | exact shared two-category recurrence | 340.286124 bits |
| \(Q=33,\ldots,64\) | exact shared two-category recurrence | 619.816215 bits |
| \(Q=65,\ldots,96\) | exact shared two-category recurrence | 1294.098372 bits |
| \(Q=97,\ldots,128\) | exact shared two-category recurrence | 2205.619520 bits |
| \(Q=129,\ldots,159\) | exact shared two-category recurrence | 2902.566478 bits |
| sampled \(Q\ge160\) | sampled dense compositions | 559.784886 bits |

The all-defect composition at \(Q=3\) is the current bottleneck. The dense
number is evidence only. A convex cover must certify every integer pair of
defect and central row counts for every \(160\le Q\le4096\).

## Why the variance lemma is the right missing object

The expected weight enumerator alone cannot bound the spectrum of one sampled
constituent tightly enough. For the exact second moment,

\[
 \mathbb E[A_w^2]
 =\mathbb E[A_w]
  +\sum_{c\ne d\in C_0\setminus\{0\}}
    \Pr[\operatorname{wt}(T c)=w,
        \operatorname{wt}(T d)=w], \tag{2}
\]

where \(T=\operatorname{Acc}P_3\operatorname{Acc}P_2
\operatorname{Acc}P_1\). The pair probability in (2) depends on the four
coordinate counts of \((c_i,d_i)\), equivalently on

\[
 (\operatorname{wt}(c),\operatorname{wt}(d),
   \operatorname{wt}(c+d)).
\]

Consequently, the ordinary BCH weight enumerator does not determine (2). A
direct proof needs a joint weight enumerator, or an operator argument that
bounds the nonconstant pair modes without materializing that enumerator.

The one-stage pair kernel is now exact. Index the four pair symbols by
\(\mathbb F_2^2\). For input type \(n=(n_a)_{a\in\mathbb F_2^2}\) and
output type \(m=(m_s)_{s\in\mathbb F_2^2}\), define the four-by-four matrix

\[
 K(\boldsymbol x,\boldsymbol y)_{s,s'}
 =x_{s+s'}y_{s'}.
\]

Starting from state zero, every length-512 path through this matrix records
the input difference symbol \(s+s'\) and the new accumulated output symbol
\(s'\). Therefore the exact transition probability is

\[
 \Pr[n\longmapsto m]
 =\frac{[\boldsymbol x^n\boldsymbol y^m]
   e_0^{\mathsf T}K(\boldsymbol x,\boldsymbol y)^{512}\mathbf1}
  {512!/\prod_a n_a!}. \tag{3}
\]

An independent integer dynamic program and exhaustive labelled-permutation
enumeration agree for four representative pair types at length 8. This check
does not bound the length-512 operator, but it removes ambiguity from the
kernel that the proof must analyze.

The most promising formulation is an association-scheme contraction bound.
A uniform interleaver averages each shell in the Johnson scheme. The
accumulator maps that shell through its run enumerator. Composing the three
stages defines a positive operator on pair types. It is enough to bound the
nonconstant part of this operator so that every diagonal shell variance is at
most 512 times its mean. This is a finite, verifier-ready target; it does not
require a claim that the complete BA spectrum is pointwise random.

The dual-distance fact sharpens this route. The operator bound only needs to
control input modes of degree at least 22, plus the explicit rank-deficient
correction above. A generic second-singular-value bound would ignore this
annihilation and is likely too weak.

The simple high-degree norm is in fact too weak. For the length-512 one-word
weight chain, projecting away polynomial degrees zero through 21 leaves
binary64 induced norms of \(2^{-3.9971}\), \(2^{-8.4864}\), and
\(2^{-15.2928}\) after one, two, and three stages. This is nowhere near the
scale required by a worst-input L2 argument. The proof must also use the
specific base distribution or the target output shell; dual distance alone
does not close (1).

### Exact scalar-character decomposition

The pair chain has an exact reducing sector that explains both the repeated
singular values and the failure of marginal-only bounds. For a rank-two pair
\((x,y)\), let the three nonzero binary characters be \(x\), \(y\), and
\(x+y\). If \(L_i\varphi\) is a function of the Hamming weight of character
\(i\), then

\[
 P_B L_i\varphi=L_i Q_B\varphi, \tag{4}
\]

where \(P_B\) is the common-interleaver pair operator and \(Q_B\) is the
one-word weight operator. This follows because the same linear accumulator
map acts on both words and therefore commutes with every binary character.
The adjoint channel is induced by the inverse linear maps and has the same
property. Hence the span of the three lifted one-word sectors is reducing,
not merely invariant. Its orthogonal complement is the interaction sector
that retains shared-overlap information.

Binary64 diagnostics verify the decomposition at lengths 8 and 16. The
scalar-sector ranks are 22 and 46, the invariance residuals are below
\(5.8\mathbin\cdot10^{-16}\) and \(4.8\mathbin\cdot10^{-15}\), and the
interaction-block norms are 0.279680 and 0.185857. The reverse leakage is at
roundoff scale. These values are useful representation evidence, but they do
not close (1). A generic L2 bound still pays for the enormous density of a
fixed rate-half code before mixing. The successful proof must use the target
shell pointwise, as the factorization below does.

### Triangle-convolution reduction

A more targeted reduction avoids the unknown joint BCH enumerator. Let
\(M=2^{256}\), let \(H(u)=\operatorname{wt}(uG_0)\), and let
\(\Phi_{t,w}(h_1,h_2,h_3)\) be the probability that both words have output
weight \(w\) after \(t\) common interleaver--accumulator stages, starting from
a pair type with triple weights \((h_1,h_2,h_3)\). The exact kernel (3)
defines \(\Phi\).

Suppose nonnegative functions \(f_w,g_w\) satisfy

\[
 \Phi_{3,w}(h_1,h_2,h_3)
 \le f_w(h_1)f_w(h_2)g_w(h_3) \tag{5}
\]

on the support of the base spectrum. For independent uniform messages
\(U,V\), the three random variables \(U,V,U+V\) are pairwise independent and
uniform. Fourier expansion on the message group, Parseval's identity, and
\(\lVert\widehat g\rVert_\infty\le \mathbb E[g]\) give

\[
 \mathbb E[f_w(H(U))f_w(H(V))g_w(H(U+V))]
 \le \mathbb E[f_w(H)^2]\mathbb E[g_w(H)]. \tag{6}
\]

Set (f_w(0)=g_w(0)=0), so the factorization covers only the rank-two
pairs (U\ne0), (V\ne0), and (U\ne V). The omitted diagonal contribution
to the second moment is exactly \(\mathbb E[A_w]\). Consequently,

\[
 \mathbb E[A_w(C_3)^2]
 \le \mathbb E[A_w(C_3)]
       +M^2\mathbb E[f_w(H)^2]\mathbb E[g_w(H)]. \tag{7}
\]

Every expectation on the right side of (7) is determined by the ordinary
base spectrum. Thus (5)--(7) can prove (1) without a genus-two BCH
enumerator. For fixed values of \(\Phi\), finding the best factorization is a
convex problem in \(\log f_w\) and \(\log g_w\).

The small RM(1,3) models pass this proof-of-concept. After separating the
degenerate pairs exactly, the length-8 model gives worst variance-to-mean
upper bounds 3.0755, 1.0330, and 0.9416 after one, two, and three stages.
For direct sums of two and three constituents, the bounds are respectively
19.2869, 2.7092, 1.2642 and 107.1594, 7.7268, 1.5681. Thus two and three
stages are below the analogue target \(F=B\) through length 24, while one
stage is not at lengths 16 and 24. The exact length-24 ratios are 1.6389,
2.1517, and 1.2607. These binary64 results are evidence for the reduction,
not for the length-512 inequality.

### Genus-two MacWilliams LP

The ordinary spectrum need not be the only admissible input to a scalable
proof. Let (J_C(n)) count ordered pairs ((x,y)in C^2) of pair type
(n=(n_{00},n_{01},n_{10},n_{11})). Let (mathcal T_B) be the genus-two
MacWilliams transform induced by the four-by-four binary character matrix.
Then

\[
 J_{C^\perp}=|C|^{-2}\mathcal T_B J_C. \tag{8}
\]

The known ordinary spectra of (C) and (C^\perp) fix all three character
marginals of these two nonnegative tables. Rank-one pair counts are also
fixed by the ordinary spectrum. For the current base code, the primal and
dual ordinary spectra agree because the imported enumerator is formally
self-dual. This fact does not assert that the concrete code is self-dual.

For each stage count (t) and shell (w), define

\[
 c_{t,w}(n):=
 \Pr[\operatorname{wt}(T_t x)=\operatorname{wt}(T_t y)=w
      \mid \operatorname{type}(x,y)=n], \tag{9}
\]

where the probability is over the (t) independent interleavers and (T_t)
includes their accumulators. The exact second moment is the linear functional

\[
 \mathbb E[A_w(C_t)^2]=\sum_n J_C(n)c_{t,w}(n). \tag{10}
\]

Consequently, maximizing (10) over nonnegative primal and dual tables that
satisfy (8) and the exact marginal constraints gives a valid upper bound for
the fixed BCH code. An outward dual LP solution would be a certificate. This
route requires neither message-pair enumeration nor an assumed genus-two BCH
enumerator.

The small tests justify pursuing this interface. At length 16, ordinary
marginals alone permit a worst variance factor of 22.9738. Adding only
nonnegativity of the MacWilliams-transformed table reduces every tested
one-, two-, and three-stage optimum to the exact value within binary64. The
worst exact value is 1.4622. Exact self-duality adds no visible improvement in
this test.

The existing character-projection kernel bound does not combine successfully
with the LP. On the actual length-16 code, that pointwise envelope already
gives a three-stage variance factor of 115779.95. The MacWilliams LP cannot
remove loss that occurs on the true pair table. A successful LP certificate
therefore needs the exact (c_{t,w}), or a materially tighter upper bound.

### Universal singular-mixing fallback

The long-chain fallback avoids the unknown BCH pair enumerator entirely. Let

\[
 \Omega_1=\mathbb F_2^{512}\setminus\{0\},\qquad
 \Omega_2=\{(x,y)\in\Omega_1^2:x\ne y\}.
\]

Both spaces carry the uniform measure. One stage samples
(P\gets S_{512}) and applies (x\mapsto\operatorname{Acc}(Px)), with the
same (P) on both coordinates of (Omega_2). This defines one-word and
ordered-distinct-pair Markov operators. The uniform measures are stationary
because both the permutation and accumulator are bijections.

The open lemma is the following exact finite statement.

> **Accumulator singular-value lemma.** On the orthogonal complement of the
> constants in stationary (L_2), the second singular values of both stage
> operators are at most (63/1000).

Assume this lemma. Let (C\subseteq\mathbb F_2^{512}) be any fixed binary
([512,256]) code. For (t) independently sampled stages,
stationary-(L_2) contraction starts from chi-square divergences

\[
 \chi_1^2=\frac{2^{512}-1}{2^{256}-1}-1,
 \qquad
 \chi_2^2=
 \frac{(2^{512}-1)(2^{512}-2)}
 {(2^{256}-1)(2^{256}-2)}-1. \tag{11}
\]

For the weight-(w) shell, set

\[
 \beta_w=\frac{\binom{512}{w}}{2^{512}-1},\qquad
 \pi_w=\frac{\binom{512}{w}(\binom{512}{w}-1)}
 {(2^{512}-1)(2^{512}-2)}.
\]

Writing (lambda=63/1000), Cauchy--Schwarz gives errors

\[
 \epsilon_{1,w}=\lambda^t
   \sqrt{\chi_1^2\beta_w(1-\beta_w)},\qquad
 \epsilon_{2,w}=\lambda^t
   \sqrt{\chi_2^2\pi_w(1-\pi_w)}. \tag{12}
\]

The shell-count mean lies between
((2^{256}-1)(\beta_w-\epsilon_{1,w})) and
((2^{256}-1)(\beta_w+\epsilon_{1,w})). Its second factorial moment is at
most

\[
 (2^{256}-1)(2^{256}-2)(\pi_w+\epsilon_{2,w}). \tag{13}
\]

Equations (12)--(13) give an explicit variance upper bound. Markov's
inequality for zero caps and Cantelli's inequality for positive caps then
test the exact integer cap vector used by the certified random
([512,256]) transfer.

At (t=127), this cap-event union bound has only 36.7702 bits of margin. At
(t=128), the outward Arb verifier gives 40.3474252482 bits. Therefore 128
is the first passing stage count under this universal (63/1000) argument.
Conditional on the same cap vector, the existing outward sparse and dense
receipts cover every occupation and give 51.6439589890 bits. Their union with
the BA-128 cap-event failure has 40.3468518020 bits of margin and proves

\[
 \Pr[d_{\min}<228{,}590]<2^{-40}. \tag{14}
\]

The probability in (14) is over 128 constituent-stage permutations, the
row and transposed-region routing permutations, and the RandomStepConv-M22
maps. The 128 constituent-stage permutations are sampled once; the resulting
single constituent is repeated in all 4096 rows. This is not rejection
sampling and does not bias the base BCH code.

Statement (14) remains conditional solely because the displayed
singular-value lemma is open. All later cap and distance arithmetic is
outward certified. The 128-stage construction is a proof fallback, not a
performance-competitive proposal. In coding conventions that count the base
block code as a constituent, this may be called BA-129.

The numerical evidence isolates a plausible proof decomposition. Functions
of the three character weights form a reducing scalar sector. The one-word
second singular value at length 512 is 0.0626241431 in binary64. On the
interaction complement, exact-kernel/binary64 calculations for even lengths
4 through 20 give norms below (3/B); the length-20 value is
0.141808598 versus (3/20=0.15). These observations motivate the sufficient
pair of sublemmas

\[
 \lVert Q_{512}\rVert_{1^\perp}\le63/1000,
 \qquad
 \lVert P_{512}|_{\rm interaction}\rVert\le3/512. \tag{15}
\]

The data do not prove either inequality in (15), nor do they by themselves
prove that the scalar-sector norm equals the one-word norm at length 512.

The one-word side now has a weaker but fully outward certificate. After
conjugating the weight chain into stationary \(L_2\), subtract the constant
singular direction and sum the squares of every remaining matrix entry. Exact
integer path counts and 256-bit Arb arithmetic give

\[
 \lVert Q_{512}-\Pi\rVert_{\mathrm{HS}}^2
 <0.011892115495235 < (11/100)^2.
\]

Hence the one-word second singular value is strictly below \(0.11\). This
does not establish the \(0.063\) hypothesis used by BA-128, and it says
nothing about the pair interaction sector. If the complete ordered-pair
operator could also be bounded by \(0.11\), the same diagnostic cap
calculation first passes at 161 accumulator stages with 41.5690 event bits.
That BA-161 implication is still conditional and is not a performance route.
The one-word receipt and verifier are
`accumulator_oneword_singular_B512_outward.json` and
`certify_accumulator_oneword_singular_outward.py`.

## Two-stage optimization diagnostic

The same cap interface was tested with two accumulator stages. With the
ordinary \(F=512\) shell caps, RandomStepConv-M24 gives 43.3264 bits for
\(Q=3,\ldots,8\), but the \(Q=1\) term has only 16.3383 bits because the
cap event permits very low constituent weights.

A stronger proof event sets every shell cap through a cutoff \(c\) to zero
and charges \(\sum_{w\le c}\mathbb E[A_w]\) directly by Markov's inequality.
This changes only the analysis; it does not reject, resample, or bias the BA
constituent. For \(c=23,24,25,26\), the outer-event margins are respectively
41.5402, 41.2674, 40.9189, and 40.5002 bits. Even with inner memory increased
to 40 bits, the corresponding \(Q=1\) margins are only 32.4548, 34.2823,
36.1182, and 38.1036 bits. The best union is therefore at most 37.8528 bits.
This rejects two stages under the present \(F=512\) shell-cap interface; it
does not prove that the two-stage construction itself fails.

An inverse sweep identifies the constant that BA-2 would actually need.
With RandomStepConv-M64 and a forced-zero cutoff at weight 26, the
event-plus-\(Q=1\) gate has 40.2090 bits under \(F=1\). Among the tested
factors, the largest passing value is \(F=2^{1.5}=2.8285\), with only
40.0260 bits before adding \(Q\ge2\). The next value, \(F=2^{1.75}\), fails.
Thus a complete BA-2 proof through this interface needs a uniform variance
constant smaller than 2.8285, and probably smaller still. The exact
length-24 ratio 2.1517 is below this necessary gate, but the length-24
triangle-convolution upper bound 7.7268 and the earlier length-32 sample maximum
10.8456 are not. BA-2 remains unproved and appears fragile.

## Rejected proof interfaces

The following failures reject bounds, not constructions.

- A minimum-distance-only wrapper fails for dense occupations because it
  discards the central spectrum shape.
- A primal/dual-distance wrapper remains far too weak. The terminated
  accumulator creates a weight-one dual word after the first stage, and even
  after enough stages to recover moderate dual distance, limited-independence
  moment tails allow too many low-weight codewords.
- Replacing the central-moment tail bound by the optimal Christoffel bound
  does not rescue that wrapper. At twelve stages, a 43.7852-bit event gives
  primal support from weight 44 through 468 and dual distance 38. The exact
  Christoffel consequence still permits about (2^{145}) words of weight 44.
  Its RandomStepConv-M64 (Q=1) bound is (2^{65.34}), rather than below
  (2^{-40}). At thirteen stages, the event margin rises to 45.1895 bits and
  the dual distance reaches 44, but the (Q=1) bound remains (2^{53.25}).
  Thus primal and dual distances alone fail even after optimizing the moment
  polynomial and increasing the inner memory.
- A spectrum-free column-Holder bound for arbitrary full-support codes loses
  hundreds of thousands of bits.
- A generic L2 contraction argument remains too weak even after deleting the
  first 21 polynomial modes. The one-word proxy contracts by only 15.29 bits
  after three stages.
- Projecting the pair kernel onto its three binary character weights gives
  an exact one-stage inequality. Exact rational checks at lengths 8 and 16
  verify all three projection inequalities and find equality witnesses.
  Factoring and iterating the inequality is nevertheless far too loose. With
  the sharper convolution average, greedy propagation at length 8 gives
  variance ratios 64.04, 474.28, and 8232.77 after one, two, and three
  stages. Jointly optimizing all intermediate factors only lowers the
  three-stage value to 8086.06. This rejects the projection majorant, not the
  greedy optimizer or the four-symbol pair kernel.
- A direct marginal-intersection bound also fails. If two output words have
  weight (w), their difference has weight at most
  (2\min(w,512-w)). Bounding the joint event by an optimized weighted
  geometric mean of these three one-word events gives an exact
  Triangle--Holder reduction. At length 512, its best worst-shell variance
  factor is (2^{256}) for every tested depth through sixteen stages. The
  bound loses all shared-interleaver overlap information.
- Splitting the defect tail into severity categories worsens the \(Q=3\)
  aggregate margin from 46.11 bits to at most 44.03 bits in the tested cuts.
  Moving the single tail boundary inward can improve sparse margin, but it
  inflates the central envelope by up to 32 bits per row. Splitting that
  inflated shoulder into a third category fails already for sparse pure-
  shoulder compositions. These are rejected change-of-measure interfaces.
- The universal hope \(\operatorname{Var}(A_w)\le\mathbb E[A_w]\) is false
  in a small exact analogue. In 500 sampled \([32,16]\) codes built from four
  RM(1,3) blocks, the maximum observed variance-to-mean ratios after one
  through four accumulator stages were 2.5517, 10.8456, 2.2244, and 1.1548.
  This experiment supports a moderate-inflation target but proves nothing at
  length 512.

## Evidence and reproducibility

Expected primal and dual spectra are generated by
`evaluate_ebch128x4_ba_expected_spectrum.py`. The principal spectrum receipt
is `ebch128x4_ba0_16_B512_expected_spectra.json`.

Conditional shell caps are generated by
`evaluate_ebch128x4_ba_variance_caps.py` and recorded in
`ebch128x4_ba_variance_cap_requirements.json`. Exact sparse transfer receipts
for variance inflation \(2^9=512\) are:

- `ebch128x4_ba3_variance9_q1q2_m22.json`;
- `ebch128x4_ba3_variance9_sparse_q3_16_m22.json`;
- `ebch128x4_ba3_variance9_sparse_q17_32_m22.json`;
- `ebch128x4_ba3_variance9_sparse_q33_64_m22.json`;
- `ebch128x4_ba3_variance9_sparse_q65_96_m22.json`;
- `ebch128x4_ba3_variance9_sparse_q97_128_m22.json`; and
- `ebch128x4_ba3_variance9_sparse_q129_159_m22.json`.

`ebch128x4_ba3_variance9_dense_probe_m22.json` contains the sampled dense
diagnostic. `small_ba32_spectrum_variance_probe.json` records the small-code
experiment. The limited-independence and column-Holder failures are recorded
in `ebch128x4_ba_limited_independence_probe_s11_s13.json` and
`full_support_B512_randomstepconv_m22_column_holder_probe.json`.

`verify_accumulator_pair_type_kernel.py` implements the exact one-stage path
count and compares it with exhaustive permutations. Its receipt is
`accumulator_pair_type_kernel_small_exact.json`.

`analyze_accumulator_pair_chain_small.py` materializes the exact rank-two
pair-type chains at lengths 8 and 16. The multinomial orbit measure is
stationary to binary64 precision. At both lengths, the leading nontrivial
pair-chain singular values equal the one-word values and occur with
multiplicity three. The second values are 0.5174040542 and 0.4179626228.
For the length-16 direct-sum base, the exact worst shell variance-to-mean
ratios after one through four stages are 1.4622, 1.1628, 1.0415, and 1.0034.
These numbers suggest a representation reduction, but they do not prove the
length-512 constant. The receipts are
`accumulator_pair_chain_B8_exact.json` and
`accumulator_pair_chain_B16_exact.json`.

`analyze_accumulator_pair_interaction_small.py` constructs the exact
scalar-character sector and measures its interaction complement. Its receipts
cover even lengths 4 through 20 in the files
`accumulator_pair_interaction_B{4,6,8,10,12,14,16,18,20}_probe.json`.

`evaluate_ba_pair_singular_mixing_sweep.py` and
`ba_pair_singular_mixing_sweep_B512.json` record the universal variance-factor
sweep. `evaluate_ba_pair_singular_random_cap_event.py` and
`ba_pair_singular_random_cap_event_B512.json` test the exact random-code cap
vector. The outward cap-event verifier and combined conditional certificate
are:

- `certify_ba128_pair_singular_random_cap_event_outward.py` and
  `ba128_pair_singular_random_cap_event_outward_B512.json`; and
- `certify_ebch128x4_ba128_randomstepconv_combined.py` and
  `ebch128x4_ba128_randomstepconv_combined_outward_s22.json`.

`analyze_accumulator_weight_chain_singular.py` records the failed generic
high-degree norm test in
`accumulator_weight_chain_B512_singular_probe.json`. Its optional multi-stage
diagnostic is `accumulator_weight_chain_B512_singular_long_probe.json`; values
after eight stages reach a binary64 floor and are not certificate evidence.
The defect-split,
boundary, and shoulder receipts are produced by
`evaluate_ebch128x4_ba_variance_proxy_split_sparse.py`, the `--low-upper`
option of `evaluate_ebch128x4_ba_variance_proxy_sparse.py`, and
`evaluate_ebch128x4_ba_variance_proxy_shoulder_sparse.py`.

`probe_triangle_holder_pair_bound_small.py` implements both the original
Holder average and the sharper convolution average in (5)--(7) through
length 24. The current receipts are
`triangle_young_pair_bound_B8_probe.json`,
`triangle_young_pair_bound_B16_probe.json`, and
`triangle_young_pair_bound_B24_probe.json`. The two-stage forced-zero
diagnostics are recorded in the
`ebch128x4_ba2_variance9_zero{23,24,25,26}_{caps,q1_m24,q1_m40}.json`
receipts. The ordinary-cap two-stage memory sweep is recorded in
`ebch128x4_ba2_variance9_q1_m{24,28,32,40,64}.json`,
`ebch128x4_ba2_variance9_q1q2_m24.json`, and
`ebch128x4_ba2_variance9_sparse_q3_*.json`.
`sweep_ebch128x4_ba2_variance_threshold.py` and its receipt give the inverse
BA-2 constant. `probe_pair_projection_holder_bound_small.py`,
`pair_projection_young_B8_probe.json`, and
`probe_pair_projection_global_small.py` with
`pair_projection_global_young_B8_probe.json` record the rejected character-
projection reduction. `verify_pair_projection_kernel_inequality.py` verifies
its one-stage premise exactly; the receipts are
`pair_projection_kernel_B8_exact.json` and
`pair_projection_kernel_B16_exact.json`.

`build_ebch128x4_ba_christoffel_caps.py` constructs the optimal shell caps
implied by the high-probability primal and dual distance events. The stage-12
and stage-13 cap receipts and the corresponding RandomStepConv-M64 (Q=1)
receipts record the rejection of this interface.

`evaluate_ebch128x4_ba_marginal_geometric_variance.py` records the rejected
three-marginal intersection bound through sixteen stages.

`probe_genus2_macwilliams_lp_small.py` verifies the genus-two transform and
solves the small marginal, formal-dual, and self-dual LPs. Its receipts are
`genus2_macwilliams_lp_B8_probe.json` and
`genus2_macwilliams_lp_B16_probe.json`.
`probe_projection_macwilliams_lp_small.py` and
`projection_macwilliams_lp_B16_probe.json` record the rejected combination
with the loose character-projection envelope.

## Remaining proof obligations

1. Prove the accumulator singular-value lemma. The outward Hilbert--Schmidt
   argument proves the one-word value below \(0.11\), but BA-128 needs the
   sharper \(0.063\) bound. It also needs an analytic interaction bound and a
   proof that the scalar/interaction decomposition controls the complete pair
   operator. Alternatively, a complete pair bound below \(0.11\) would give
   the slower conditional BA-161 fallback.
2. After item 1, rerun the two outward verifiers and bind their hashes in a
   final manifest. No further all-\(Q\) transfer work is needed for BA-128.
3. Separately optimize toward BA-3 by constructing a scalable exact
   second-order input-output representation. Use it for an outward
   genus-two LP dual or a triangle-convolution witness for every shell.
4. A BA-3 witness must imply variance factor at most 512. A BA-2 witness must
   meet the stricter factor below 2.8285 and is unlikely to retain enough
   failure budget under the present interface.
5. Only after a short-stage mathematical closure, measure its implementation
   cost against the available performance budget.

The next proof task is item 1. It is narrower than the earlier genus-two
problem and already suffices for a complete, though impractical, finite
theorem. A direct length-512 pair-type matrix remains too large; the verifier
must exploit the scalar/interaction decomposition.
