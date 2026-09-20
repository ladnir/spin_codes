# Goal 01: one active outer block

Build an exact first-moment evaluator for messages that activate one outer
block. The evaluator accepts a weight spectrum \(A_w\) and performs these
steps:

1. Condition on the outer weight \(w\).
2. Average over the uniform \(w\)-subset of active transposed rows.
3. Average over the bit permutation within every row.
4. Apply a Chernoff bound to the RandomStepConv output weight.
5. Multiply by \(L A_w\) and sum over \(w\).

For one active block, an inactive row has transfer \(R_0(z)\). An active row
contains one nonzero input bit at a uniform inner position and has transfer
\(R_1(z)\). The complete transfer conditioned on weight \(w\) is

\[
 \frac{[u^w](R_0(z)+uR_1(z))^B}{{B\choose w}}.
\]

The coefficient is a matrix coefficient. Matrix multiplication retains the
order of the transposed rows.

The initial spectrum is the expected spectrum of a uniform random
\([B,B/2]\) linear injection:

\[
 \overline A_w=
 (2^{B/2}-1)\frac{{B\choose w}}{2^B-1},
 \qquad w>0.
\]

This fractional spectrum is an ideal benchmark, not a constituent code.

Success criteria:

- exhaustive small cases validate the row coefficient calculation;
- the ideal spectrum gives a plausible baseline against the random-outer
  calculation;
- the receipt exposes every weight contribution, so a BCH spectrum or an
  upper envelope can replace the benchmark without changing the inner model.

Recommended next goal: evaluate two active outer blocks from their weight
pair \((w_1,w_2)\), then identify a compression that extends to higher
occupations.
