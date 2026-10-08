# High-distance proof closure

The proof goal is complete for the explicitly documented s=16 refinement:
K20, rate 1/2, 40% GF65536-symbol distance, and more than 59.87 bits of
setup-failure margin. The native-performance target remains separate.
The randomized RS[256,128] outer, one regional route, and shared beta64 inner
are the baseline. Compatible refinements must be stated explicitly and
proved before implementation optimization or promotion.

## Certified coverage at the start

- All low-weight words with 1 through 5 active outer rows: union below 2^-56.5182.
- All low-weight words with 7548 through 8192 active rows: below 2^-326.6848.
- All nonzero quiet words with 1 through 962 active rows: covered by the
  sparse quiet certificate and its exact-coefficient extension.
- Quiet low-weight words with 5865 through 8192 active rows: below 2^-355.6458.

Remaining classes: quiet 963 through 5864 and nonquiet 6 through 7547.
The goal is not complete until their union is covered with the required margin.

## Iteration log

1. Started the joint input/output-weight calculation. For a Cauchy tile with
   a potential inputs, the systematic map [I_a;F_U] is MDS[a+8,a,9]. Its split
   weight enumerator is exact. The first feedback-aware upper calculation
   distinguishes zero and nonzero input tiles before imposing the support
   mark. This avoids charging all potential positions whenever a packet is
   nonzero. Numerical outcomes and checked coverage will be appended after
   verification; this entry is not a certificate.

All accepted endpoints require directed arithmetic. Floating searches select
witnesses only. Preserve unsuccessful tests when they identify a specific
loss, and use those findings to choose the next proof step.

2. Closed quiet active-row counts 1508 through 8192 by tile conditioning,
   the distinct-hyperplane rank refinement, and verified regional coefficient
   caps. Witnesses and formulas are in
   [the joint-input note](HD_RANDOMIZED_OUTER_JOINT_INPUT.md). The unchanged
   construction now leaves quiet counts 963--1507 and nonquiet 6--7547 open.
   Larger input marks did not close the selected low-middle witnesses.
3. Derived a compatible but distinct [sum-rank inner](HD_SUMRANK_INNER.md).
   Its local quiet and flat-state mask ranks are exact. The finite-field
   regression checks pass. Initial packetwise fractional screens still
   leave intermediate counts open; they are not certificates. Current
   work retains input weights and reduces state-history splitting loss.
4. Closed every quiet active-row class for the s=4 sum-rank candidate.
   The directed combined failure bound is below 2^-93.9900. This is a
   complete quiet-class result, not a whole-code certificate. A new
   pointwise live-output bound uses all occupied tiles, not just untouched
   ones. Its exact combinatorial formula and regressions are recorded.
5. Derived and checked a two-component density operator for independent
   [feedback refresh](HD_FEEDBACK_REFRESH.md). It removes arbitrary
   residual mass at birth. State sizes 4, 8, and 12 are proof-side screens;
   none is implemented or promoted. Regional short-block bounds are being
   tested to close nonquiet classes. Initial FLOAT underflow outputs were
   rejected; the short-block products now use log arithmetic.
6. Retaining profiles across pairs and enforcing each region's exact
   occupancy gives a new [s=16 proof candidate](HD_SUMRANK_FEEDBACK_CERTIFICATE.md).
   Directed all-word intervals now cover every count from 1 through 8192.
   The clean pinned replay reports 53 receipts, no coverage gaps, and
   log2(failure) <= -59.87186558535774. This candidate uses a larger internal-state multiplier;
   it is not the old native implementation and carries no timing claim.
   The s=4 quiet-class certificate remains independently closed. Numerical
   positive bounds for s=4,8,12 have not been mislabeled as counterexamples.
7. Verified the complete replay at 256-bit Arb precision and selected weak
   witnesses at 384 bits. All 375 high-distance regression tests pass.
   The four frozen native-source hashes are unchanged. No native benchmark
   or implementation work was performed. The next recommendation is to
   reduce state/refresh cost while retaining a complete certificate, then
   return to native performance engineering.

No native implementation, timing, public API, or prior certificate was changed.
