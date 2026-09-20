# Goal 02: three-block zero-gap transfer

## Question

The field outer code has block distance three.  Its minimum nonzero profile
therefore contains three BCH words.  Does the profile with binary weights
\((22,22,22)\) contract at output weight \(D=188{,}766\)?

The calculation below answers this question affirmatively.  It does not yet
cover larger BCH weights or the field-outer weight spectrum.

## Probability spaces

Set

\[
N=524{,}352,
\qquad
D=188{,}766.
\]

The real three-block experiment samples three independent binary vectors
\(X_1,X_2,X_3\).  Each \(X_j\) is uniform on the weight-22 slice of
\(\{0,1\}^{128}\).  This is the distribution produced by independently
permuting the bits of three fixed BCH words of weight 22.

For the comparison experiment, sample \(X\) uniformly from the weight-66
slice of \(\{0,1\}^{384}\).  Partition its coordinates into three labeled
parts of 128 bits.  Define the balanced event

\[
\mathcal B:=\{\operatorname{wt}(X_1)=
\operatorname{wt}(X_2)=\operatorname{wt}(X_3)=22\}.
\]

The real law equals the comparison law conditioned on \(\mathcal B\).  The
balance probability is

\[
p_{\mathrm{bal}}
=\Pr[\mathcal B]
=\frac{\binom{128}{22}^3}{\binom{384}{66}}
=2^{-6.0614583497}.
\]

The global packet permutation is independent of \(X\) and \(\mathcal B\).
Thus, for every event \(E\) determined after that permutation,

\[
\Pr_{\mathrm{real}}[E]
=\Pr[E\mid\mathcal B]
\le \frac{\Pr[E]}{p_{\mathrm{bal}}}.
\tag{1}
\]

Equation (1) is the superblock transfer.  It replaces three balanced slices
with one exchangeable 384-bit slice at a cost of 6.061 bits in the exponent.

## Active paths and gaps

Divide \(X\) into 96 four-bit packets.  Delete its zero packets and retain the
order of the \(H\) nonzero packets.  The global packet permutation makes their
positions a uniform \(H\)-subset of the \(N\) packet positions.  Exchangeability
of the comparison slice makes the retained packet order uniform as well.

Let \(q_i\in\{0,1\}^4\) be the accumulator state after the \(i\)-th nonzero
packet.  Define

\[
a_i:=\operatorname{wt}(q_i),
\qquad
m_0:=|\{i:a_i=0\}|.
\]

Let \(y_i\ge1\) be the number of packet positions for which state \(q_i\)
persists.  A leading zero-state gap absorbs all positions before the first
nonzero packet.  The binary output weight is

\[
W=\sum_{i=1}^{H}a_i y_i.
\tag{2}
\]

There are \(\binom NH\) possible placements of the active packets.

## Exact treatment of zero-state gaps

Fix an active path with parameters \((H,m_0)\).  There are \(H-m_0\) positive
state gaps.  Replacing every positive state weight by one gives

\[
\#\{\text{placements with }W\le D\}
\le
\binom{D}{H-m_0}
\binom{N-H+m_0}{m_0}.
\tag{3}
\]

The first factor counts positive gap lengths whose sum is at most \(D\).  The
second factor lets the leading gap and all \(m_0\) zero-state gaps absorb every
remaining position.  Consequently,

\[
\Pr[W\le D\mid H,m_0]
\le
\min\!\left\{
1,
\frac{
\binom{D}{H-m_0}\binom{N-H+m_0}{m_0}
}{\binom NH}
\right\}.
\tag{4}
\]

An exact 16-state dynamic program computes the joint law of \((H,m_0)\) over
the weight-66 comparison slice.  Averaging (4) gives

\[
\Pr[W\le D]\le 2^{-30.8994148543}.
\]

The transfer (1) therefore gives

\[
\Pr_{\mathrm{real}}[W\le D]
\le 2^{-24.8379565046}.
\tag{5}
\]

This bound already passes Goal 02.

## Positive-state-weight refinement

Bound (3) discards the difference between accumulator weights one through
four.  Fix \(u>0\).  For a positive state of weight \(a_i\),

\[
\sum_{y_i\ge1}e^{-u a_i y_i}
=\frac{e^{-u a_i}}{1-e^{-u a_i}}.
\]

The exponential counting bound and the same exact treatment of zero-state
gaps give

\[
\#\{\text{placements with }W\le D\}
\le
e^{uD}
\binom{N-H+m_0}{m_0}
\prod_{i:a_i>0}
\frac{e^{-u a_i}}{1-e^{-u a_i}}.
\tag{6}
\]

The dynamic program can average the product in (6) without enumerating active
paths.  It tracks the used input weight, packet support, accumulator state,
and number of zero prefix states.  For every positive prefix state, it
multiplies the path weight by the corresponding factor in (6).

Using

\[
u=\frac{33.7}{188{,}766}
=0.0001785279129
\]

gives the comparison-slice bound

\[
\Pr[W\le D]\le2^{-42.5086746381}.
\]

After the transfer penalty,

\[
\boxed{
\Pr_{\mathrm{real}}[W\le188{,}766
\mid(22,22,22)]
\le2^{-36.4472162884}
}.
\tag{7}
\]

The parameter 33.7 was selected by a small one-dimensional numerical search.
The validity of (7) does not depend on proving that this parameter is optimal.

## Refutation search

A reproducible search sampled 20,000 paths from the real balanced law.  The
mean packet support was 51.2875, and the mean number of zero prefix states was
3.2534.  At least one zero return occurred in 96.305% of samples.  The largest
observed count was 13; the 99.9th percentile was 10.

Thus zero returns remain common, but the search found no distinct three-block
obstruction.  The proof succeeds because it sums zero-state gaps separately
and retains the weights of positive states.  The search is numerical evidence
only and is not part of bound (7).

## Validation and scope

The implementation checks three independent gates:

1. The dynamic program sums to \(\binom{384}{66}\).
2. Its packet-support marginal matches the inclusion-exclusion formula.
3. Exhaustive enumeration verifies both gap inequalities on 128,670 small
   instances.

The result proves contraction only for the fixed profile \((22,22,22)\).  A
complete distance proof must cover every attainable BCH weight triple and
combine those bounds with the field-outer spectrum.
