# Goal 01: initial cyclic live-coverage report

## Result

The exact gap identity passes the initial checks, but the first full-size
diagnostic is negative for the proposed 9% target. Two laps remove the fixed
late-start boundary. They do not remove the corresponding clustered-support
event on the cycle.

At packet support \(h=16384\), the floating two-tilt calculation gives

\[
\log_2 p_h(188766)\le -41480.35.
\]

The one-lap diagnostic at the same support was \(-41470.64\). The measured
improvement is only about 9.7 bits.

This comparison is diagnostic. Neither value is an outward-rounded
end-to-end certificate.

## Explicit cyclic-cluster witness

Fix a cyclic interval of \(R\) packet positions that crosses the boundary
between the end and start of the retained lap. Require all \(h\) active
packets to lie in this interval. Also require the two interval endpoints next
to the complementary gap to be active.

Suppose the state resets after \(k\) live zero packets in the complementary
gap. The retained lap then has at most \(R+k\) live packets. Conditional on
this event, its weight is stochastically dominated by

\[
\operatorname{Bin}(g(R+k),1/2).
\]

Consequently,

\[
p_h(D)\ge
\frac{\binom{R-2}{h-2}}{\binom Nh}
\sum_{k=0}^{L-1}a^kq\,
\Pr[\operatorname{Bin}(g(R+k),1/2)\le D],
\tag{3}
\]

where \(L=N-R\). Truncating the positive sum preserves a lower bound.

For \(h=16384\), the initial evaluator uses \(R\approx99000\) and reports a
lower-bound exponent near \(-41755\). Thus the explicit witness lies within
about 275 bits of the upper bound on a 41,000-bit scale.

## Interpretation

The burn-in lap changes the bad event rather than eliminating it:

1. A rare global permutation clusters all active packets into about 19% of
   the packet cycle.
2. The complementary zero-input gap occupies about 81% of the cycle.
3. One random-map reset turns the state off inside that gap.
4. The next active packet does not arrive until the other end of the gap.

For \(\sigma=20\), a gap of this length has a substantial reset probability.
The reset therefore adds little exponent beyond the support-clustering cost.

Typical supports have short gaps, but the end-to-end first moment must also
pay for exponentially many outer words. The rare clustered supports are the
relevant configurations.

## Status and next decision

Goal 01 remains open because the target calculations are not outward rounded.
The structural question is already narrow: increasing the number of laps
alone does not appear sufficient while zero remains absorbing between active
packets.

The next useful comparison should change the state-transition law. Candidate
options include a nonzero-preserving state transition or a small independent
state injection at each step. Either change requires a distinct construction
name.
