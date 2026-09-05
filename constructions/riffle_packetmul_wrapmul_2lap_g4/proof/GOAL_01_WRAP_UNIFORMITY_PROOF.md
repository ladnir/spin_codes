# Goal 01: wrap-state uniformity

## Result

The WrapMul interface removes the conditioning obstruction at the first-lap
terminal state. After fixing the complete inherited setup, every nonzero
terminal state becomes uniform over the \(2^{64}-1\) nonzero states.

This result does not make the second-lap output uniform. It reduces the
nonzero-terminal distance event to a finite list-size question.

For the 26 authenticated support-33 outer words, a uniform bad-state list cap
of 316,505 fits the exact remaining ledger budget. A cap of 316,506 does not.

## Setup experiment

Fix an outer word \(x\). The inherited setup samples the packet permutation
\(\Pi\) and packet multipliers \(A\). Define

\[
y:=\Pi(D_Ax),
\qquad
\ell:=L(y).
\]

The new setup independently samples

\[
B\gets\mathbb F_{2^{64}}^\times.
\]

The encoder uses the wrapped terminal state

\[
Z:=B\ell.
\]

The same realized value \(B\) is public and is reused for every outer word.
No runtime randomness is sampled between laps.

## Conditional-uniformity lemma

Let \(\mathcal F\) contain the fixed outer word and all inherited setup
randomness. On the event \(\ell\ne0\), for every
\(z\in\mathbb F_{2^{64}}^\times\),

\[
\Pr[Z=z\mid\mathcal F]
=\frac1{2^{64}-1}.
\]

### Proof

Fix a realization of \(\mathcal F\) for which \(\ell\ne0\). Multiplication by
\(\ell\) is a bijection of \(\mathbb F_{2^{64}}^\times\). For each nonzero
\(z\), exactly one scalar,

\[
B=z\ell^{-1},
\]

satisfies \(B\ell=z\). The scalar \(B\) is uniform and independent of
\(\mathcal F\). The displayed conditional probability follows.

The statement remains valid after conditioning on any event determined by
\(\mathcal F\), including a first-node stratum. Such conditioning does not
inspect \(B\).

## Zero-terminal and autonomous transfers

If \(\ell=0\), then \(Z=0\) for every \(B\). The new encoder output equals the
parent output on this event. Therefore, the authenticated PacketMul
support-33 terminal-zero bound transfers without a probability loss.

If \(\ell\ne0\), then \(Z\ne0\). The three-node and four-node autonomous
certificates hold for every nonzero entering state. They therefore transfer
pointwise. In particular, the sufficient zero-prefix cutoff remains node
20,976.

Packet support and its placement law also transfer. WrapMul acts after the
first lap and does not change the driven packet word.

## Exact list reduction

Fix \(x\) and a realization of \(\mathcal F\). Define the fixed driven offset

\[
c_y:=F(F(y)).
\]

For \(z\ne0\), define the corresponding two-lap output

\[
Y_y(z):=c_y\mathbin\oplus J(z).
\]

The bad-state list is

\[
\mathcal B_y
:=\left\{
z\in\mathbb F_{2^{64}}^\times:
\operatorname{wt}(Y_y(z))\le188{,}765
\right\}.
\]

The conditional-uniformity lemma gives the exact identity

\[
\Pr_B[
\operatorname{wt}(G_B(y))\le188{,}765
\mid\mathcal F,\ell\ne0
]
=\frac{|\mathcal B_y|}{2^{64}-1}.
\]

This identity is pointwise in the driven offset. It does not average over the
packet permutation or packet multipliers.

## Support-33 list budget

Let \(R\) be the exact remaining numerical budget from the parent Goal 04
ledger. Its log-base-two interval is

\[
[-41.027690825931,-41.027690825930].
\]

Suppose every realization associated with any of the 26 authenticated
support-33 words satisfies

\[
|\mathcal B_y|\le K.
\]

A union bound over the 26 words gives

\[
\Pr[\text{a bad support-33 output with }\ell\ne0]
\le\frac{26K}{2^{64}-1}.
\]

The same public scalar \(B\) is reused across the words. The union bound does
not require their bad events to be independent.

Exact rational arithmetic gives

\[
\max\left\{
K\in\mathbb Z_{\ge0}:
\frac{26K}{2^{64}-1}\le R
\right\}
=316{,}505.
\]

For this value,

\[
\log_2\left(\frac{26K}{2^{64}-1}\right)
\approx-41.0276915166.
\]

The cap 316,506 exceeds \(R\). Thus 316,505 is the largest uniform cap that
could use the entire remaining budget.

Spending the entire budget on support 33 would leave essentially no margin
for other open rows. A proof should seek a smaller list or combine list sizes
with placement probabilities.

## First boundary list

For first occupied node

\[
r\in\{20{,}972,20{,}973,20{,}974,20{,}975\},
\]

the complete four-node autonomous blocks contribute at least 188,748. A bad
output must therefore have weight at most 17 outside those blocks.

For a fixed boundary drive, let \(\mathcal T_{y,r}(17)\) be the nonzero states
whose remaining output has weight at most 17. Then

\[
\mathcal B_y\subseteq\mathcal T_{y,r}(17).
\]

The next bounded goal should determine or bound
\(|\mathcal T_{y,r}(17)|\) over the four boundary strata. This count directly
tests whether WrapMul provides enough margin where the parent proof first
stops.

## Limitations

Goal 01 does not bound \(|\mathcal B_y|\) or
\(|\mathcal T_{y,r}(17)|\). It does not close the support-33 row.

The uniformity lemma applies only when the terminal state is nonzero. The
transferred PacketMul argument remains necessary for the zero-terminal event.

The report makes no claim about integrated performance or the complete
construction.
