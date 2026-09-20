# Distance-push obstruction

The 9% theorem has 80.1685 bits of margin. This note asks what happens when
the target distance increases while the outer block remains the modeled
binary `[256,128]` code.

## Zero-state class

Write the transposed input as 16,384 epochs

\[
X_0,\ldots,X_{16383}\in\mathbb F_2^{128}.
\]

Let \(C_\pi\) denote the outer code after the sampled coordinate and region
permutations. Define

\[
\mathcal Z_\pi
:=
\left\{
x\in C_\pi:
B(X_i)=0\text{ for }0\leq i<16383
\right\}.
\]

The final epoch has no constraint because the encoder does not emit its
outgoing state.

Fix any \(x\in\mathcal Z_\pi\). Since \(Q_0=0\), induction gives
\(Q_i=0\) at the start of every epoch. Therefore

\[
Y_i=X_i
\]

for every epoch. The inner encoder is the identity on
\(\mathcal Z_\pi\). The sampled field scalars do not affect these words.

This class is large for every fixed setup. The map

\[
x\longmapsto
\bigl(B(X_0),\ldots,B(X_{16382})\bigr)
\]

has rank at most \(16\cdot16383\). Hence

\[
\dim(\mathcal Z_\pi)
\geq
2^{20}-16\cdot16383
=786448.
\]

The dimension bound does not by itself prove that \(\mathcal Z_\pi\)
contains a word below the target weight. It shows that zero-state words are
an exact structural subcode, rather than an artifact of the Chernoff bound.

## Dominant occupation at 9.35%

Set

\[
D=196130,
\qquad
\delta=D/2^{21}=0.0935220718.
\]

The regular calculation is dominated by approximately 2,650 active outer
blocks. Before the value bits are exposed, these blocks contribute
\(256\cdot2650=678400\) candidate coordinates. A typical unconditioned
assignment has 339,200 one bits.

The optimized Chernoff parameter has log-surprisal \(-0.205\). Along the
zero-state path, one transposed region has the following tilted profile.

- It contains 2,650 candidate coordinates across 64 epochs.
- One epoch contains 41.406 candidate coordinates on average.
- One epoch emits 11.999 bits on average.
- Approximately 28.98% of the candidate values are one.

Thus the bad profile does not rely on a live state canceling the input.
Instead, the outer word has low value weight and every nonfinal epoch
satisfies the 16 checks imposed by \(B\).

A Bernoulli-shell approximation makes the epoch class more concrete. After
conditioning on \(B(X_i)=0\) and tuning the mean to 12, the most likely
epoch weights are 12, 14, 10, and 16, with respective probabilities 21.5%,
19.9%, 16.3%, and 13.5%. Weight eight contributes 8.2%, and the zero word
contributes 7.4%. The 64 exceptional weight-four kernel words contribute
only 0.08%. Thus the obstruction uses ordinary medium-weight kernel words,
not a rare minimum-weight pattern.

For one ordinary region, the base-two logarithm of the exact zero-to-zero
moment is

\[
-2265.043938.
\]

The final region may leave the state live after its last epoch. Its moment is

\[
-2249.148836.
\]

Omitting the unused final syndrome therefore saves 15.895102 exponent bits.
The resulting pointwise margin at occupation 2,650 is 43.495984 bits.

The full two-state transfer gives 43.493341 bits at the same occupation and
tilt. The difference is less than 0.003 bits. Therefore the explicit
zero-state class accounts for essentially the complete dominant transfer.

After summing the neighboring occupations, the regular margin is 40.1680
bits at distance 196,130. It falls to 38.9885 bits at distance 196,131. The
one-active class still has 76.6013 bits at distance 196,130.

At 11% distance, the same class moves to approximately 3,094 active outer
blocks. Its explicit zero-state term has base-two exponent
\(42544.7812\), so the first-moment calculation fails by more than 42,000
bits. The full transfer has the same exponent to two decimal places.

## Interpretation

The current certificate stops near 9.3522% because the triangular inner map
has a large identity subcode and the regular outer bound loses one parity
bit per active block. The obstruction is not live-state cancellation, state
termination, or insufficient randomness in the field scalars.

## Outer parity audit

The modeled outer spectrum contains only even-weight regular words. The
regular analyzer replaces each such word by fair independent bits and pays
a density factor of two. At occupation 2,650, this envelope loses 2,650
exponent bits.

Conditioning a 256-bit value column on even parity should recover almost one
bit under the tilted zero-state law. Subtracting this ideal parity loss moves
the projected 40-bit crossing to

\[
D=198379,
\qquad
D/2^{21}=0.0945944786.
\]

This projection is not yet a certificate because the current recurrence
does not track each outer block's parity across the 256 regions. It shows
that the 9.3522% threshold is partly an envelope artifact.

The parity correction does not remove the zero-state obstruction. At 11%,
the explicit class uses approximately 3,094 active outer blocks. Granting a
full 3,094-bit parity recovery still leaves a failed exponent of about
39,451 bits. Thus exact parity accounting may move the threshold by roughly
one tenth of a percentage point, but it cannot plausibly move this inner to
the rate-one-half GV distance.

The calculation identifies the class that controls the first-moment bound.
It does not yet prove that a typical sampled setup contains a word from this
class at the predicted weight. A matching lower bound or a second-moment
calculation would turn the obstruction into refutation evidence about the
sampled construction itself.

The next design comparison should vary the number of independent checks per
128-bit epoch. Increasing the state width, decreasing the epoch width, or
adding an output transformation that is nontrivial on \(\ker B\) directly
targets this class.
