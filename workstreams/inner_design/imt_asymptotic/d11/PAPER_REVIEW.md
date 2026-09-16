# IMT proof review and manuscript integration

Reviewed on 2026-09-15. The review found no unresolved analytic obstruction
to the assembled 11% claim for the native-length, rate-half family. The paper
now instantiates the asymptotic theorem with IMT. This was an in-session
derivation review, not independent external review or machine-checked proof.

The reviewed argument is in `paper/structured_appendix.tex` (outer and route)
and `paper/structured_imt_appendix.tex` (IMT transfers and occupancy analysis).
The construction and adjoint are in `paper/structured_spin.tex`. The exact
expansion and feedback generators are transcribed in the appendix and checked
against the certified maps by `paper/check_imt_integration.py`.

## Analytic checks

- The transvection marginal is derived for each fixed nonzero entering state;
  epoch independence is what permits its conditional use. Shared setup does
  not supply independent coins to different messages, nor does the first
  moment require them.
- The occupation transfer preserves a dominating weighted measure on zero,
  arbitrary live states, and uniform image shells. Cancellation bounds use
  the feedback spectrum separately from the expansion spectrum. The Fourier
  transfer bounds each target of the tilted syndrome law. No entrywise
  minimum between those representations is used.
- Dense routing domination is applied to a moment restricted to the fixed
  total input weight before changing the iid reference parameter. Concavity
  combines unequal outer-row weights only after the region comparison.
- Sparse counting-measure domination by fair rows is pointwise after summing
  messages. Marked-input conditioning costs at most `(8 sqrt(Q))^b` times
  the KL factor. This avoids an unjustified `b log(L)` sparse loss.
- The empty-gap covariance identity is valid from every initial live state.
  The resulting weighted-kernel bound controls both endpoint states, not just
  the stationary row sum. Short spacings and epoch collisions have vanishing
  probability for fixed Q. The omitted impulse durations cost O_Q(1/L).
- Componentwise convergence is used only after checking structural zeros.
  It averages within-region placements and is uniform over active-row bit
  patterns and boundary states. It does not claim uniformity over prescribed
  microscopic impulse locations. This is sufficient for coefficient extraction.
- The product-norm argument works for the nonstochastic upper impulse matrix
  because its beta-integral comparison is pathwise with nonnegative coefficients.
- The final union uses fixed Q<4096, a uniform geometric bound for Q>=4096
  in the sparse range, and a uniform dense exponent. No interchange of a
  pointwise limit and an unbounded sum is needed.

## Corrections and scope

The old manuscript defined the expected BA enumerator to count nonzero words
but included the zero-input term in its exact formula. Its sum now starts at
input weight one. The imported outer proof already treats nonzero weights;
its retained sources and receipts are unchanged.

The earlier field-multiplier transfers and constants were replaced, not
renamed. The manuscript explicitly gives the conservative impulse matrix,
the new sparse contraction `1-96 alpha`, and the refined 39-segment outer
majorant. The new theorem keeps t=128, s=19, and c=39/4.

The finite BCH ladder, parameter plots, and timing tables remain RM2Sub
results. This review does not migrate those configurations or switch any
encoder default. Smaller finite IMT cells and the desired performance scope
remain separate work.

Historical proof drafts, numerical producers, and receipts remain byte-for-byte
unchanged so their bindings still replay. Their original `review_pending`
fields describe their creation-time status; this review record and the paper
record the subsequent review. None of these fields asserts formal verification.

## Checks and next step

The new manuscript checker authenticates the retained 11% assembly and dense
replay and verifies all 38 generator words. The existing finite checker still
checks the five RM2Sub margins, selected map, plots, and historical timings.
The numerical evidence already includes a 512-bit dense replay and exact
sparse polynomial reconstruction; no benchmark is part of this update.

Next: migrate the required finite IMT cells and measurements before removing
the remaining RM2Sub presentation. External artifact packaging must include
the new local evidence; no data files were committed by this integration.
